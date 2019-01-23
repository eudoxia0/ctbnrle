/* Continuous Time Bayesian Network Reasoning and Learning Engine
 * Copyright (C) 2009 The Regents of the University of California
 *
 * see docs/AUTHORS for contributor list
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "importancesampler.h"
#include "ctbn.h"
#include "utils.h"
#include "samplequeue.h"



namespace ctbn {

//#define DEBUG
using namespace std;

ImportanceSampler::ImportanceSampler(const Process *pr,
		const Trajectory *traj,
		const VarSample *m) {
	SetProcess(pr);
	SetTrajectory(traj);
	SetMethod(m);
	discard_bn_weight = false;
}

ImportanceSampler::~ImportanceSampler() {}

ImportanceSampler *ImportanceSampler::Clone() const {
	return new ImportanceSampler(*this);
}

void ImportanceSampler::SampleTrajectories(vector<Trajectory> &tr, 
		vector<double> &w, int numsamples, Random &rand) {
	MakeIndex();
	const Markov* pr = dynamic_cast<const Markov*>(p);
	Instantiation x = evid->Values(pr->GetDynamics()->Domain() +
			pr->GetDynamics()->CondDomain(), begintime);
	for (int i=0; i<numsamples; i++) {
		Trajectory traj;
		traj.SetBeginTime(begintime);
		traj.SetEndTime(endtime);
		double weight = SampleSingleTrajectory(traj);
		tr.push_back(traj);
		w.push_back(weight);
	}
}

void ImportanceSampler::MakeCache(const Instantiation &inst, 
		const vector<Trajectory::Index> &nextevid) {
	nodeindex.clear();
	nodeval.clear();
	nextevidval.clear();

	for (unsigned int i=0; i<nodes.size(); i++) {
		nodeindex.push_back(nodes[i]->CondDomain().Index(inst));
		nodeval.push_back(inst.Value(str.Node2Var(i)));
		nextevidval.push_back(nextevid[i].Values().Value(str.Node2Var(i)));
	}
}

double ImportanceSampler::SampleDyn(Trajectory &traj, Random &rand) {
	unsigned int numofnodes = nodes.size();
	vector<Trajectory::Index> nextindex = evidindex; //index for each node
	changeindex = vector<bool>(numofnodes, 0); 
	int * nextevidtype = new int [numofnodes];
	double weight = 0;

	double currt = begintime;
	const Dynamics* dyn = dynamic_cast<const Markov*>(p)->GetDynamics();
	//Get initial values
	Instantiation curri = traj.Values(dyn->Domain()+dyn->CondDomain(), currt);

	for (unsigned int i=0; i<numofnodes; i++)
		nextevidtype[i] = nextindex[i].TestInc(nodes[i]->Domain());
	//check conddomain
	if (nextindex.size()>numofnodes) ++nextindex[numofnodes];
	MakeCache(curri, nextindex);
	int count = 1;
	int movecount = 0;
	vector<int> undefinelist; //save the nodes that need to be resampled
	for (unsigned int i=0; i<numofnodes; i++) 
		undefinelist.push_back(i);
	//SampleQueue here is only used to save the transition time 
	//of each node (variable), the transition value is not used
	//in each event (they are all set to be 0 in our implemetation).
	SampleQueue events(numofnodes);
	while (currt<endtime) {
		for (unsigned int k=0; k<undefinelist.size(); k++) {
			unsigned int i = undefinelist[k];
			nodeindex[i] = nodes[i]->CondDomain().Index(curri);
			double newt = SampleTime(curri, i, nodeindex[i], nodeval[i], 
				nextevidtype[i], nextindex[i], nextevidval[i], currt, rand);
			events.Add(SampleQueue::Event(i, 0, newt));
		}	
		SampleQueue::Event e;
		events.Head(e); //Get the next transition variable
		int nextid = e.var;
		double nextt = e.time;
		//check transitions in conddomain of the process
		if (nextindex.size()>numofnodes && nextt>nextindex[numofnodes].Time()) {
			nextt = nextindex[numofnodes].Time();
			nextid = -1;
		}
		//##############################################
		//This part checks if the sampling process falls into
		//some state that it will never be able to jump out
		if (fabs(nextt-currt)<1e-12) movecount++;
		else movecount = 0;
		if (movecount>100) exit(1);
		//##############################################
		if (nextt>=endtime) {
			nextt = endtime;
			for (unsigned int i=0; i<numofnodes; i++)
				weight += SampleTimeWeight(curri, i, nextid, nodeindex[i], 
							nodeval[i], nextevidtype[i], nextindex[i], 
							nextevidval[i], currt, nextt);
			return weight;
		}
		//Calculate the weight contribution of sampling time
		//for each node
		for (unsigned int i=0; i<numofnodes; i++)
			weight += SampleTimeWeight(curri, i, nextid, nodeindex[i], 
						nodeval[i], nextevidtype[i], nextindex[i], 
						nextevidval[i], currt, nextt);
//cout << "weight: " << weight << endl;
//cout << "--------------------------------------------" << endl;
		//int oldval = nodeval[nextid];
		if (nextid!=-1) {
			int oldval = nodeval[nextid]; 
			int newval = SampleTransition(curri, nextid, nodeindex[nextid], 
					nodeval[nextid], nextevidtype[nextid], nextindex[nextid], 
					nextevidval[nextid], nextt, weight, rand);
			if (oldval!=newval) {
				traj.AddTransition(str.Node2Var(nextid), nextt, newval);
				nodeval[nextid] = newval;
				/*
				weight += SampleTransitionWeight(curri, nextid, oldval, nodeindex[nextid], 
						nodeval[nextid], nextevidtype[nextid], nextindex[nextid], 
						nextevidval[nextid], currt, nextt);
		
					*/
				//assert(!changeindex[nextid] || nextevidtype[nextid]==2);
			} //else assert(changeindex[nextid] && nextevidtype[nextid]==1);

			currt = nextt;
#ifdef DEBUG
			cout << "nextid: " << nextid << ", nextt: " << nextt << endl;
			cout << "nextindex time: " << nextindex[nextid].Time() << endl;
			cout << "nexti: ";
			curri.PrintVal(cout);

#endif
			if (changeindex[nextid]) {
				nextevidtype[nextid] = nextindex[nextid].TestInc(nodes[nextid]->Domain());
				// switch from below to next line reported as bug by Asela Gunawardana
				//nextevidval[nextid] = nextindex[nextid].Values().Value(nextid);
				nextevidval[nextid] = nextindex[nextid].Values().Value(str.Node2Var(nextid));
#ifdef DEBUG
				cout << "Change index: nextid: " << nextid << endl;
				cout << "next evid type: " << nextevidtype[nextid] << endl;
				cout << "nextevid val: " << nextevidval[nextid] << endl;
#endif
			}
		} else {
			//Add the transition of the conddomain to the trajectory
			//we have to assume that the evidence of conddomain is fully observed
			//otherwise, we can't sample any trajectory
			traj.AddTransition(nextindex[numofnodes].Values(), nextindex[numofnodes].Time());
			++nextindex[numofnodes];
			currt = nextt;
		}
		UndefineTime(nextid, undefinelist);
		for (unsigned int k=0; k<undefinelist.size(); k++) {
			events.Remove(undefinelist[k]);
		}
	}
	delete [] nextevidtype;
	return weight;
}

double ImportanceSampler::SampleSingleTrajectory(Trajectory &traj, bool logscale, Random &rand) {
	double currt = begintime;
	Instantiation curri = evid->Values(
		dynamic_cast<const Markov*>(p)->GetDynamics()->Domain(), currt);

	//Sample initial values
	// notes by Jing: for network application: only use the dynamic's weight
	double weight = 0.0;
	if (discard_bn_weight) SampleInitial(curri);
	else weight = SampleInitial(curri);
	traj.AddTransition(curri, begintime);
	weight += SampleDyn(traj);
	//cout << weight << endl;
	return logscale ? weight : exp(weight);
}

} // end of ctbn namespace
