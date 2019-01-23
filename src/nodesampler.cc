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
#include "nodesampler.h"

#define DEBUG

namespace ctbn {

using namespace std;

NodeSampler::NodeSampler(const Process *pr, const Trajectory *traj, const VarSample *m) {
	SetProcess(pr);
	SetTrajectory(traj);
	SetMethod(m);
	GetChildrenContext();
	MakeIndex();
}

NodeSampler *NodeSampler::Clone() const {
	return new NodeSampler(*this);
}

double NodeSampler::SampleNodeTrajectory(const Trajectory &tr, Trajectory &nodetraj, int nodeid, Random &rand) {
	Trajectory::Index nextindex = evidindex[nodeid];
	changeindex = vector<bool>(nodes.size(), 0);
	double currt = tr.TimeBegin();
	nodetraj.SetBeginTime(begintime);
	nodetraj.SetEndTime(endtime);
	double weight = 0.0;
	//remove the content of nodeid in tr
	//get the index of it's parent for sampling
	//cout << "begin remove" << endl;
	//tr.SetUnknown(nodeid);
	//cout << "remove node" << endl;
	//tr.Draw(cout);

	Context c = nodes[nodeid]->CondDomain();
	Trajectory::Index parenttraj = tr.Begin(c);

	Instantiation curri = tr.Values(dynamic_cast<const Markov*>(p)->GetDynamics()->Domain(), currt);
	//set evidence of initial value of nodeid
	curri.SetVal(nodeid, nextindex.Values().Value(nodeid));
	//sample initial value;
	const BN *bn=dynamic_cast<const BN*>(dynamic_cast<const Markov*>(p)->GetStartDist());
	bn->Sample(curri, rand);

	//At this monent, Let's just use the network that
	//the initial distribution of each node is independent
	//we only consider the weight contribution of the sampling node
	Instantiation tmpi(curri, -1);
	tmpi.SetVal(nodeid, nextindex.Values().Value(nodeid));
	//weight += bn->ImportanceSampleWeight(curri, tmpi, true); //Feb-12 comment only for now

	int nextevidtype = nextindex.TestInc(nodes[nodeid]->Domain());
	int nextevidval = nextindex.Values().Value(nodeid);
	int nodeval = curri.Value(nodeid);

	nodetraj.AddTransition(nodeid, currt, nodeval);
	//++parenttraj;
	double nextt = currt;
	int nodeparindex = parenttraj.Values().Index();
	while (currt<endtime) {
		double newt = SampleTime(curri, nodeid, nodeparindex, nodeval, nextevidtype, nextindex, nextevidval, currt, rand);
		double trajtime = parenttraj.Time() + parenttraj.DeltaT();
		if (newt>trajtime) nextt = trajtime;
		else nextt = newt;
		int nextvar = parenttraj.NextVar();
		if (nextt>=endtime) {
			currt = nextt;
			continue;
		}
		if (newt>trajtime) {
			++parenttraj;
			nodeparindex = parenttraj.Values().Index();
			curri.SetVal(nextvar, parenttraj.Values().Value(nextvar));
			continue;
		}
		int oldval = nodeval;
		int newval = SampleTransition(curri, nodeid, nodeparindex, nodeval, nextevidtype, nextindex, nextevidval, newt, weight, rand);

		if (oldval!=newval) {
			nodetraj.AddTransition(nodeid, newt, newval);
			nodeval = newval;
			//comment only for now
			//weight += SampleTransitionWeight(curri, nodeid, oldval, nodeparindex, nodeval, nextevidtype, nextindex, nextevidval, currt, nextt);
		}
		currt = nextt;
		if (currt>=nextindex.Time()) {
			nextevidtype = nextindex.TestInc(nodes[nodeid]->Domain());
			nextevidval = nextindex.Values().Value(nodeid);
		}
	}
	//tr.Draw(cout);
	//calculate weight contribution of nodeid's children
	//weight += ChildrenSampleWeight(tr, nodeid); //Feb-12, comment only for now

	weight = 0.0; //For now, it only works for Jing's special case.
	return weight;
}

double  NodeSampler::ChildrenSampleWeight(const Trajectory &tr, int nodeid) const {
	double retweight = 0.0;
	vector<int> childrenid = str.GetChildren(nodeid);
	int numchildren = childrenid.size();
	if (numchildren==0) return 0.0;

	Trajectory::Index trajindex = tr.Begin(childrencontext[nodeid]);

	//llh of initial values
	const BN *bn=dynamic_cast<const BN*>(dynamic_cast<const Markov*>(p)->GetStartDist());
	Instantiation inst = trajindex.Values();
	for (int i=0; i<numchildren; i++)
		retweight += bn->Node(childrenid[i])->Prob(inst, true);

	map<int, SS *> nodess;
	for (int i=0; i<numchildren; i++)
		nodess.insert(make_pair(childrenid[i], nodes[childrenid[i]]->BlankSS()));
	double currt = trajindex.Time();
	while (!trajindex.Done()) {
		Instantiation curri = trajindex.Values();
		double deltat = trajindex.DeltaT();
		int nextvar = trajindex.NextVar();
		for (int i=0; i<numchildren; i++)
			nodes[childrenid[i]]->AddSS(curri, currt, deltat, nodess[childrenid[i]]);
		++trajindex;
		map<int, SS*>::const_iterator vit = nodess.find(nextvar);
		currt = trajindex.Time();
		if (vit!=nodess.end()) {
			Instantiation nexti = trajindex.Values();
			nodes[nextvar]->AddTransSS(curri, nexti, currt, vit->second);
		}
	}

	for (map<int, SS*>::const_iterator i=nodess.begin(); i!=nodess.end(); ++i) {
		retweight += nodes[i->first]->LLH(i->second);
		delete i->second;
	}
	return retweight;
}


// void NodeSampler::SampleNodeTrajectory(Trajectory &tr, int nodeid, Random &rand)
// {
//   Trajectory::Index nextindex = evidindex[nodeid];
//   int nextevidtype = nextindex.TestInc(nodes[nodeid]->Domain());
//   int nextevidval = nextindex.Values().Value(nodeid);
//   double currt = tr.TimeBegin();

//   //remove the content of nodeid in tr
//   //get the index of it's parent for sampling
//   tr.SetUnknown(nodeid);
//   Context c = nodes[nodeid]->CondDomain();
//   Trajectory::Index parenttraj = tr.Begin(c);

//   Instantiatoin curri = tr.Values(dynamic_cast<const CTBN*>(p)->Domain(), currt);
//   //sample initial value;
//   SampleInitial(curri, rand);

//   int nodeval = curri.Value(nodeid);
//   tr.AddTransition(nodeid, currt, nodeval);
//   parenttraj.TestInc();
//   while (currt<tr.TimeEnd())
//     {
//       int nodeindex = parenttraj.Value().Index();
//       double newt = SampleTime(curri, nodeid, nodeindex, nodeval, nextevidtype, nextindex, nextevidval, currt, rand);

//       if (newt>endtime) return;
//       if (newt>parenttraj.Time()) continue;

//       int oldval = nodeval;
//       int newval = SampleTransition(curri, nodeid, nodeindex, nodeval, nextevidtype, nextindex, nextevidval, newt, rand);
//       if (oldval!=newval)
//         {
//           tr.AddTransition(nodeid, newt, newval);
//           nodeval = newval;
//         }
//       currt = newt;
//       if (currt>=nextindex.Time())
//         {
//           nextevidtpe = nextindex.TestInc(nodes[nodeid]->Domain());
//           nextevidval = nextindex.Values().Value(nodeid);
//         }
//     }
// }

// double  NodeSampler::SampleWeight(const Trajectory &tr, int nodeid, bool logscale)
// {
//   double retwegiht = 0.0;
//   Context condcontext = nodes[nodeid]->CondDomain();
//   Context nodecontext = nodes[nodeid]->Domain();

//   vector<int> childrenid = str.GetChildren(nodeid);
// //   map<int, int> childrenval;
// //   map<int, int> childrenindex;
//   vector<int, int> childrenval;
//   vector<int, int> childrenindex;
//   int numchild = childrenid.size();

//   Trajectory::Index trajindex = tr.Begin(mbcontext);
//   Trajectory::Index nodeindex = evidindex[nodeid];

//   Instantiation curri = tr.Values(globalcontext, tr.TimeBegin());
//   int nodeval = curri.Value(nodeid);
//   int nextevidtype = nodeindex.TestInc(nodecontext);
//   int nextevidval = nextindex.Values().Value(nodeid);

//   for (int i=0; i<numchild; i++)
//     {
//       chidrenval.push_back(curri.Value(childrenid[i]));
//       childrenindex.push_back(nodes[childrendi[i]]->CondDomain().Index(curri));
//     }
// //   for (int i=0; i<numchild; i++)
// //     {

// //     }

//   int parindex = condcontext.Index(curri);

//   while (!trajindex.Done())
//     {
//       int nextvar = trajindex.NextVar();
//       double nextt = trajindex.Time()+trajindex.DeltaT();
//       bool inc = false;
//       if (nextt>nodeindex.Time()+nodeindex.DeltaT())
//         //evidence begin before next transition
//         {
//           nextt = nodeindex.Time()+nodeindex.DeltaT();
//           inc = true;
//         }
//       //calculate weight contribution of nodeid
//       retweight += SampleTimeWeight(curri, nodeid, nextvar, parindex, nodeval, nextevidtype, nodeindex, nextevidval, currt, nextt);

//       if (inc)
//         nodeindex.TestInc(nodes[nodeid]->Domain());
//       else
//         {
//           //contribution of nodeid's children
//           for (int i=0; i<numchild; i++)
//             retweight += SampleTimeWeight(curri, childrenid[i], nextvar, childrenindex, childrenval[i], trajindex, childrenval[i], currt, nextt); //fill the parameters

//           ++trajindex;
//           if (nextvar==nodeid)
//             {
//               retweight += SampleTransitionWeight(curri, nodeid, oldval, nodeindex, nodeval[nextid], nextevidtype[nextid], nextindex[nextid], nextevidval[nextid], currt, nextt);
//             }
//           if (IsChild(nextvar, nodeid)) //transition var is child of nodeid
//             retweight += SampleTransitionWeight(curri, nextvar, childrenval[nextvar], childrenindex[nextvar], childrenval[nextvar], 2, trajindex, nextevidval[nextid], currt, nextt);


//         }
//     }
//   if (logscale) return retweight;
//   else return exp(retweight);
// }

void NodeSampler::GetChildrenContext() {
	childrencontext.clear();
	int n = nodes.size();
	for (int i=0; i<n; i++) {
		Context c;
		vector<int> child = str.GetChildren(i);
		for (unsigned int j=0; j<child.size(); j++) {
			c = c + nodes[child[j]]->Domain(); //add children
			c = c + nodes[child[j]]->CondDomain(); //add children's parents
		}
		childrencontext.push_back(c);
	}
}

} // end of ctbn namespace
