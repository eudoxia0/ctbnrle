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
#include <algorithm>
#include "gibbsbase.h"
#include "ctbndyn.h"
#include "markov.h"
#include "markovdyn.h"
#include "multisimple.h"
#include "bn.h"
#include "rk.h"
#include "utils.h"
#include "params.h"



namespace ctbn {

using namespace std;

// "process" is a pointer to a (subclass of ) Process object.
// "evidence" is a pointer to Trajectory object that holds the context.
// Variables with value -1 means no evidence.
// Interval evidence X(t:t+dt)=x is represented as X(t)=x, X(t+dt)=-1.
// Point evidence X(t)=x is represented as X(t)=x, X(t+epsilon)=-1.
GibbsBase::GibbsBase(const Process *process, const Trajectory *evidence,
				 int burnin) {
	p = process;
	init_traj = NULL;
	evid = evidence;
	numBurninIter = burnin;
	Initialize();
}

void GibbsBase::Initialize() {
	const Markov *ctbn = dynamic_cast<const Markov *>(p);
	const CTBNDyn *ctbndyn = dynamic_cast<const CTBNDyn *>
		(ctbn->GetDynamics());
	own_var_list = ctbndyn->Domain().VarList();
	context = ctbndyn->Domain() + ctbndyn->CondDomain();
	p0 = ctbn->GetStartDist();

	vector<int> var_list = context.VarList();
	for (size_t j=0; j!=var_list.size(); ++j) {
		// markovdyn for each variable
		const MarkovDyn *markov_dyn = dynamic_cast<const MarkovDyn *>
			(ctbndyn->NodeByVar(var_list[j]));
		Context mb = markov_dyn->Domain(); // *will* become Markov Blanket
		vector<int> domain_var_list = mb.VarList();
		mb = mb + markov_dyn->CondDomain();
		// domain_var_list should have exactly one variable
		mds[domain_var_list[0]] = markov_dyn;

		// parents and children set for each variable
		Context parent_context = Context();
		Context children_context = Context();
		for (int i=0; i<ctbndyn->NumofNodes(); ++i) {
			Context c = ctbndyn->Node(i)->Domain();
			Context cc = ctbndyn->Node(i)->CondDomain();
			if (c.HasId(var_list[j]))
				parent_context = parent_context + cc;
			if (cc.HasId(var_list[j])) {
				children_context = children_context + c;
				mb = mb + c + cc;
			}
		}
		parents[var_list[j]] = parent_context;
		children[var_list[j]] = children_context;
		markovblanket[var_list[j]] = mb;
	}
	burntin = false;
}

void GibbsBase::SampleTrajectories(vector<Trajectory> &traj,
					vector<double> &w,
					int numsamples, Random &rand) {
	BurnIn(rand);
	for (int i=0; i<numsamples; ++i) {
		traj.push_back(Get());
		w.push_back(0.0);  // log weight
		Next();
	}
}

void GibbsBase::BurnIn(Random &rand) const {
	if (burntin) return;
	if (!init_traj) {
		SampleInitialTrajectory(rand);
	} else {
		tr = *init_traj;
	}
	Next(numBurninIter, rand);
	burntin = true;
}

void GibbsBase::Next(int num_iter, Random &rand) const {
	for (int i=0; i<num_iter; ++i) {
		for (size_t var=0; var!=own_var_list.size(); ++var)
			SampleVariable(own_var_list[var], rand);
	}
}

void GibbsBase::SetInitTraj(const Trajectory &t) {
	if (init_traj)
		delete init_traj;
	init_traj = new Trajectory(t);
	tr = *init_traj;  // also reset the current sample to init_traj
}

void GibbsBase::ClearInitTraj() {
	if (init_traj) {
		delete init_traj;
		init_traj = NULL;
	}
}

void GibbsBase::SampleInitialTrajectory(Random &rand) const {
	tr = Trajectory();
	tr.SetBeginTime(evid->TimeBegin());
	tr.SetEndTime(evid->TimeEnd());

	// Sample a start distribution from BN. The start distribution should
	// be consistent with the initial values in evidence.
	// NOTE: start with the next event after tr.Timebegin(), 
	// since the first one is already handled in BN's sample; 
	Trajectory::Index i = evid->Begin(context);
	Instantiation inst(i.Values());
	Instantiation myinst(inst);
	p0->Sample(myinst, rand);
	tr.AddTransition(myinst, tr.TimeBegin());
	double dt = i.DeltaT();

	// loop over evidence
	for (++i; !i.Done(); ++i) {
		Instantiation newinst(i.Values());
		Instantiation myoldinst(myinst);
		myinst.SetVal(Instantiation(newinst.KnownVars(),newinst),true);

		Context diff(myinst, myoldinst);
		vector<int> diff_vars = diff.VarList();
		for (vector<int>::iterator varid=diff_vars.begin();
							varid!=diff_vars.end();++varid) {
			double event_time;
			if (inst.Value(*varid)==-1)
				event_time = i.Time() - randomizer.RandReal() * dt;
			else event_time = i.Time(); // observed transition!
			tr.AddTransition(*varid, event_time, 
						newinst.Value(*varid));
		}
		inst = newinst;
		dt = i.DeltaT();
	}
}

void GibbsBase::SampleVariable(int v, Random &rand) const {
	// Remove all v's transitions from tr.
	oldtr = tr;
	tr.SetUnknown(v,true);

	vector<int> x0, xT;
	vector<double> t0, tT;
	SplitFreeIntervals(v, x0, xT, t0, tT, 
				tr.Values(context, tr.TimeBegin()), rand);
	if (x0.empty())  // no free interval for v
		return;

	for (size_t i=0; i!=x0.size(); ++i)
		SampleVariableInterval(v, x0[i], xT[i], t0[i], tT[i], rand);
}

void GibbsBase::SplitFreeIntervals(int v, 
					vector<int> &x0, 
					vector<int> &xT,
				 	vector<double> &t0, 
					vector<double> &tT, 
					const Instantiation &inst0, 
					Random &rand) const {
	// NOTE: if no evidence for v specifies the value of v at 
	// tr.TimeBegin() and tr.TimeEnd(), then we must sample 
	// x0 and xT for them.
	Context vcontext;
	vcontext.AddVar(v,context.Cardinality(v));
	Trajectory::Index i = evid->Begin(vcontext);
	Instantiation inst(i.Values());
	int oldval = inst.Value(v);
	if (oldval == -1) {
		// open free interval
		t0.push_back(i.Time());
		x0.push_back(-1);
	} else {
		tr.AddTransition(v, tr.TimeBegin(), oldval);
	}

	for (++i; !i.Done(); ++i) {
		inst = i.Values();
		int newval = inst.Value(v);

		// "&& oldval!=-1" added (cshelton 9/18/12) b/c I think that
		// will be fixed when sampling
		if (newval!=oldval && newval!=-1 && oldval!=-1)
			tr.AddTransition(v, i.Time(), newval);

		if (oldval==-1 && newval!=-1) {
			// close free interval
			tT.push_back(i.Time());
			xT.push_back(newval);
		} else if (oldval!=-1 && newval==-1) {
			// open free interval
			t0.push_back(i.Time());
			x0.push_back(oldval);
		}
		oldval = newval;
	}
	if (oldval == -1) {
		// close free interval
		tT.push_back(evid->TimeEnd());
		xT.push_back(-1);
	}
}

void GibbsBase::GetTau(const Trajectory &t, 
				double begin_time, 
				double end_time, 
				int v, const Context &varset,
				vector<double> &tau, 
				vector<Instantiation> &inst) {
	Trajectory::Index index = t.Begin(varset, begin_time, end_time);
	tau.push_back(begin_time);
	inst.push_back(index.Values());

	for (++index; !index.Done(); ++index) {
		Instantiation newinst(index.Values());
		inst.push_back(newinst);
		tau.push_back(index.Time());
	}

	tau.push_back(end_time);
}


// assumes input vectr T is right sized
bool GibbsBase::GetTMatrix(int v, const Instantiation &oldinst, 
					const Instantiation &newinst, vectr &T) const {
	int vdim = context.Cardinality(v);
	Context diff(oldinst, newinst);

	// NOTE: assume there is exactly one variable changed
	int changeid = diff.VarList()[0];

	// if the changed variable is not a child of v; 
	// simply return an identity matrix;
	if (!children[v].HasId(changeid)) return false;
	// o.w. return a diagonal matrix whose T[i][i] = q(a->a')|i
	int a = oldinst.Value(changeid);
	int b = newinst.Value(changeid);
	const MarkovDyn *md = mds[changeid];
	Instantiation inst = oldinst;

	for (int i=0; i<vdim; ++i) {
		inst.SetVal(v, i);
		T[i] = (*md)(inst)->Intensity()[a][b]; //get the CIM
	}
	return true;
}

} // end of ctbn namespace
