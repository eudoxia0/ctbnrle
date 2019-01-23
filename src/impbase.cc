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
#include "impbase.h"
#include "nullptr03.h"


namespace ctbn {

using namespace std;

void ImpBase::SetProcess(const Process *pr) {
	p = pr;
	const CTBNDyn *ctbndyn = 
		dynamic_cast<const CTBNDyn*>(dynamic_cast<const Markov *>(pr)->GetDynamics());
	ctbndyn->GetStructure(str);
	int n = ctbndyn->NumofNodes();
	nodes.clear();
	for (int i=0; i<n; i++)
		nodes.push_back(dynamic_cast<const MarkovDyn*>(ctbndyn->Node(i)));
}

double ImpBase::SampleTime(const Instantiation &curri, int id,
		int parindex, int nodeval,
		int evidencetype, Trajectory::Index &nextindex, int e,
		double t, Random &rand) {
	double ret = 0.0;
	const matrix &Q = (*nodes[id])(parindex)->Intensity();
	if (evidencetype==0) {
		//no upcoming evidence
		ret = t + method->SampleTime(Q, id, nodeval, curri, t, evid, -1, -1, rand);
	} else if (evidencetype==1 && e!=-1) {
		//there is an upcoming evidence
		double te = nextindex.Time();
		ret = t + method->SampleTime(Q, id, nodeval, curri, t, evid, te, e, rand);
		if (ret>te&&nodeval==e) {
			ret = te;
			changeindex[id] = 1;
		}
	} else {
		//evidence change value
		ret = nextindex.Time();
		changeindex[id] = 1;
	}
	return ret;
}

double ImpBase::SampleTimeWeight(const Instantiation &curri, int nodeid, int transitionid,
		int parindex, int nodeval,
		int evidencetype, Trajectory::Index &nextindex, int e,
		double currt, double nextt) {
	const matrix &Q = (*nodes[nodeid])(parindex)->Intensity();
	double te = nextindex.Time();
	if (evidencetype==2 || (evidencetype==1 && e==-1)) {
		//node is forced to agree the evidence,
		//calculate the likelihood
		return Q[nodeval][nodeval] * (nextt-currt);
	} else if (evidencetype==1 && e!=-1) {
		//there is an upcoming evidence
		return method->TimeWeight(Q, curri, nodeid, transitionid, 
				nodeval, currt, nextt-currt, evid, te, e);
	} else {
		return method->TimeWeight(Q, curri, nodeid, transitionid, 
				nodeval, currt, nextt-currt, evid, -1, -1);
	}
}

double ImpBase::SampleInitial(Instantiation &initval, Random &rand) {
	const BN* bn = dynamic_cast<const BN*>(dynamic_cast<const Markov*>(p)->GetStartDist());
	return bn->ImportanceSample(initval, true,rand);
}

int ImpBase::SampleTransition(Instantiation &inst, int id,
		int parindex, int nodeval,
		int evidencetype, Trajectory::Index &nextindex, int e,
		double t, double &weight, Random &rand) {
	forceflip = 0;
	int retval;
	const matrix &Q = (*nodes[id])(parindex)->Intensity();
	if (evidencetype==2) {
		//force a transition according to the evidence
		retval = e;
		//matrix Q = (*nodes[id])(parindex)->Intensity();
		//weight += log(Q[nodeval][retval]) - log(-(Q[nodeval][nodeval]));
		// cshelton 8/20/12: removed second part of term above
		//  (fixes weight calculation) -- also removed similar below in 
		//     corner case of too close to evidence
		weight += log(Q[nodeval][retval]);
	} else if (evidencetype==1) {
		if (e!=-1) {
			//there is an upcoming evidence
			double te = nextindex.Time();
			if (te>t) {
				//matrix Q = (*nodes[id])(parindex)->Intensity();
				if (te-t>1e-6 || nodeval==e) {
					//sample next transition according to the evidence
					retval = method->SampleTransition(Q, id, nodeval, inst, t, evid, te, e, weight, rand);
				} else {
					//too close to the upcoming evidence,
					//force the variable transition to the 
					//value of the upcoming evidence and 
					//update the correponding weight
					retval = e;
					forceflip = 1;
					//matrix Q = (*nodes[transitionid])(parindex)->Intensity();
					weight += log(Q[nodeval][e]); // - log(-(Q[nodeval][nodeval]));
					forceflip = 0;
				}
			} else {
				//we are at the time where evidence is a transition that
				//the value changes from unknown to known. At this time,
				//the current value of the variable should agree with the
				//evidence. So we do nothing.
				retval = nodeval;
			}

		} else {
			//the evidence changes from known value to unknow value
			retval = nodeval;
		}
	} else {
		//there is no upcoming evidence
		//matrix Q = (*nodes[id])(parindex)->Intensity();
		retval = method->SampleTransition(Q, id, nodeval, inst, t, evid, -1, -1, weight, rand);
	}
	//update the instantiation
	if (nodeval != retval) // below should probably be left to caller, as std.Node2Var(id)
			// might already have been computed there
		inst.SetVal(str.Node2Var(id), retval);
	return retval;
}

/*
double ImpBase::SampleTransitionWeight(const Instantiation &nextinst, int transitionid, int oldval,
		int parindex, int nodeval,
		int evidencetype, Trajectory::Index &nextindex, int e,
		double currt, double nextt) {
	double retweight = 0.0;
	if (evidencetype==2) {//force a transition according to the evidence
		int val = oldval;
		int nextval = nodeval;
		matrix Q = (*nodes[transitionid])(parindex)->Intensity();
		retweight = log(Q[val][nextval]) - log(-(Q[val][val]));
	} else if (evidencetype==1&&e==-1)
		//evidence  change to known value to unknown
		retweight = 0;
	else if (evidencetype==1&&e!=-1) {
		//retweight = method->TransitionWeight();
		if (forceflip==1) {
			int val = oldval;
			int nextval = nodeval;
			matrix Q = (*nodes[transitionid])(parindex)->Intensity();
			retweight = log(Q[oldval][e]) - log(-(Q[val][val]));
			forceflip = 0;
			//cout << "force flip weight: " << retweight << endl;
		} else {
			double te = nextindex.Time();
			matrix Q = (*nodes[transitionid])(parindex)->Intensity();
			retweight = method->TransitionWeight(Q, transitionid, oldval, nextinst, nextinst, currt, nextt-currt, evid, te, e);
		}
	} else {
		retweight = 0;
	}
	return retweight;
}
*/

void ImpBase::UndefineTime(int id, vector<int> &undefinelist) {
	undefinelist.clear();
	if (id!=-1) {
		changeindex[id] = 0;
		undefinelist = str.GetChildren(id);
		undefinelist.push_back(id);
		for (unsigned int i=0; i<undefinelist.size(); i++) {
			changeindex[undefinelist[i]] = 0;
		}
	} else {
		//conddomain changes
		Context conddomain = dynamic_cast<const Markov*>(p)->GetDynamics()->CondDomain();
		vector<int> varlist = conddomain.VarList();
		vector<bool> clist(nodes.size(), false);
		for (unsigned int i=0; i<varlist.size(); i++) {
			vector<int> c = str.GetChildren(str.Var2Node(varlist[i]));
			for (unsigned int j=0; j<c.size(); j++) {
				clist[c[j]] = true;
			}
		}
		for (unsigned int i=0; i<clist.size(); i++)
			if (clist[i]) 
				undefinelist.push_back(i);
	}
}

void ImpBase::MakeIndex() {
	if (p == nullptr03 || evid == nullptr03) return;
	evidindex.clear();
	for (unsigned int i=0; i<nodes.size(); i++) {
		Context nodecontext = nodes[i]->Domain();
		Trajectory::Index it = evid->Begin(nodecontext);
		evidindex.push_back(it);
	}
	//Add index of conddomain of the whole process if has one
	Context conddomain = dynamic_cast<const Markov *>(p)->GetDynamics()->CondDomain();
	if (conddomain.NumVars()!=0) {
		Trajectory::Index it = evid->Begin(conddomain);
		evidindex.push_back(it);
	}
}

} // end of ctbn namespace
