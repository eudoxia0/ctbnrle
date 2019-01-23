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
#include "gibbssampler.h"
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
GibbsSampler::GibbsSampler(const Process *process, const Trajectory *evidence,
				 int burnin) : GibbsBase(process,evidence,burnin) {
}

GibbsSampler::~GibbsSampler() {
}


void GibbsSampler::SampleVariableInterval(int v, int x0, int xT, 
						double t0, double tT, 
						Random &rand) const {
	static double kEpsilon = ParamDouble("GibbsEpsilon",1e-5);

	if (x0==-1) { // must be initial interval and var not known at beginning
		Instantiation init_val = SampleConditionalRv(tr.Values(context,t0),
								t0,v,rand);
		x0 = init_val.Value(v);
		tr.AddTransition(v,t0,x0);
	}

	// interval too small
	if (fabs(tT-t0) < 2.0*kEpsilon) {
		if (x0==xT || xT==-1)  // already consistent
			return;
		tr.AddTransition(v, (t0+tT)/2.0, xT);
		return;
	}
	
	vector<double> tau;
	vector<Instantiation> inst;
	GetTau(tr, t0, tT, v, markovblanket[v], tau, inst);

	// Sample the time to stay at state x0.
	// Compute P(future) vectors for {tau} from K back to 0.

	vector<MultiSimple> p_future;  // normalized distribution vectors
	vector<double> log_z;  // log of the normalization factors of p_future
	double logfactor;  // required by GetDist but not used
	int tau_size = static_cast<int>(tau.size());
	int vdim = context.Cardinality(v);
	for (int i=0; i!=tau_size; ++i) {
		p_future.push_back(MultiSimple(vdim));
		log_z.push_back(0.0);
	}
	vectr dist(vdim, 0.0);
	if (xT == -1)
		dist = 1.0;
	else
		dist[xT] = 1.0;
	p_future[tau_size-1].SetDist(dist);
	vectr Tvec(vdim);
	for (int i=tau_size - 2; i>=0; --i) {
		// compute p_future[i] from p_future[i+1]
		matrix R = GetRMatrix(v, inst[i]);
		p_future[i+1].GetDist(dist, logfactor);
		logfactor = expmtv(dist, R, tau[i+1] - tau[i]);
		if (i > 0 && GetTMatrix(v,inst[i-1],inst[i],Tvec))
			dist = dist.dotstar(Tvec);
		p_future[i].SetDist(dist);
		logfactor += log(p_future[i].Normalize());
		log_z[i] = logfactor + log_z[i+1];
	}
	double log_p_future0 = p_future[0].Prob(x0, true) + log_z[0];
	double ksi = rand.RandReal();
	double log_p_past = 0.0;
	int which_interval;
	for (which_interval=0; which_interval<tau_size-1; ++which_interval) {
		double new_log_p_past = log_p_past;
		matrix R = GetRMatrix(v, inst[which_interval]);
		new_log_p_past += (tau[which_interval+1] - 
					tau[which_interval]) * R[x0][x0];
		double log_p_future = 
			p_future[which_interval+1].Prob(x0, true) +
			log_z[which_interval+1];
		double Ft = 1.0 - 
			exp(new_log_p_past + log_p_future - log_p_future0);
		if (ksi < Ft)
			break;
		if (which_interval < tau_size-2) {
			if (GetTMatrix(v, inst[which_interval], 
						inst[which_interval+1],Tvec))
				new_log_p_past += log(Tvec[x0]);
		}
		log_p_past = new_log_p_past;
	}
	// This happens only when xT == x0, in which case F(tT) < 1.
	if (which_interval >= tau_size - 1)
		return;
	matrix R = GetRMatrix(v, inst[which_interval]);
	double t = SearchPoint(ksi, x0, R, log_p_future0, log_p_past, 
				p_future[which_interval+1], 
				log_z[which_interval+1],
				tau[which_interval], 
				tau[which_interval+1], 
				kEpsilon);

	// stay till the end of trajectory
	if (fabs(tT-t) < 2.0*kEpsilon) {
		if (x0==xT || xT==-1)  // already consistent
			return;
		tr.AddTransition(v, (t+tT)/2.0, xT);
		return;
	}

	// Sample the next state to transition to.
	p_future[which_interval+1].GetDist(dist, logfactor);
	expmtv(dist, R, tau[which_interval+1] - t);
	for (int i=0; i<vdim; ++i)
		dist[i] *= R[x0][i];
	dist[x0] = 0.0;
	MultiSimple pft(vdim);
	pft.SetDist(dist);
	pft.Normalize();
	int next_state = pft.Sample(rand);
	// add the event to tr
	tr.AddTransition(v, t, next_state);
	SampleVariableInterval(v, next_state, xT, t, tT, rand);
}

matrix GibbsSampler::GetRMatrix(int v, const Instantiation &inst) const {
	const MarkovDyn *md = mds[v];
	matrix R = (*md)(inst)->Intensity();
	vector<int> children_id = children[v].VarList();
	for (size_t i=0; i!=children_id.size(); ++i) {
		int val = inst.Value(children_id[i]);
		Instantiation pinst = inst;
		const MarkovDyn *cmd = mds[children_id[i]];
		for (int j=0; j<R.getm(); ++j) {
			pinst.SetVal(v, j);
			R[j][j] += (*cmd)(pinst)->Intensity()[val][val];
		}
	}
	return R;
}

double GibbsSampler::SearchPoint(double ksi, 
					int x0, 
					const matrix &R, 
					double log_p_future0, 
					double log_p_past, 
					MultiSimple p_future, 
					double log_z, 
					double t0, 
					double t1, 
					double epsilon) const {
	// This loop is guaranteed to terminate since each iteration we will
	// reduce (t1 - t0) by a half.
	for (;;) {
		double t = (t0+t1)/2.0;
		double dt = (t1-t0)/2.0;
		if (dt < epsilon)
			return t;
		// Compute F(t).
		vectr dist;
		double logfactor;
		p_future.GetDist(dist, logfactor);
		logfactor = expmtv(dist, R, dt);
		logfactor += log(dist.normalize());
		logfactor += log_z;

		double new_log_p_past = log_p_past;
		new_log_p_past += dt*R[x0][x0];
		MultiSimple new_p_future = p_future;
		new_p_future.SetDist(dist);
		double Ft = 1.0 - exp(new_log_p_past + log(dist[x0]) 
							  + logfactor - log_p_future0);
		if (fabs(Ft-ksi) < epsilon)
			return t;
		if (Ft < ksi) {  // search (t, t1)
			t0 = t;
			log_p_past = new_log_p_past;
		} else {
			t1 = t;  // search (t0, t)
			p_future = new_p_future;
			log_z = logfactor;
		}
	}
}

Instantiation GibbsSampler::SampleConditionalRv(const Instantiation &partial,
						double t, int v, Random &rand) const {
	double log_prob_e = p0->Prob(partial, true);
	int card_v = context.Cardinality(v);
	vector<double> posterior_prob(card_v);
	double max_log_prob = -1e20;
	for (int val=0; val<card_v; ++val) {
		Instantiation inst(partial);
		inst.SetVal(v, val);
		posterior_prob[val] = p0->Prob(inst, true);
		if (posterior_prob[val] > max_log_prob)
			max_log_prob = posterior_prob[val];
	}
	// now do "beta pass" for evidence after this time (on this var
	// 	or others)
	Context vc;
	vc.AddVar(v,context.Cardinality(v));
	Trajectory::Index vi = evid->Begin(vc);
	double t0 = vi.Time();
	assert(vi.Values().Value(v)==-1); // otherwise why are we here?
	++vi; 
	vectr beta(card_v);
	if (vi.Done() || vi.Values().Value(v)==-1) {
		beta = 1.0;
	} else {
		beta = 0.0;
		beta[vi.Values().Value(v)] = 1.0;
	}
	vector<double> tau;
	vector<Instantiation> tauinst;
	GetTau(tr,t0,vi.Time(),v,markovblanket[v],tau,tauinst);
	vectr Tvec(card_v);
	for(int i=tau.size()-2;i>=0;--i) {
		expmtv(beta,GetRMatrix(v,tauinst[i]),tau[i+1]-tau[i]);
		if (i>0 && GetTMatrix(v,tauinst[i-1],tauinst[i],Tvec))
			beta = beta.dotstar(Tvec);
	}

	// put the two together:
	double summation = 0.0;
	for (int val=0; val<card_v; ++val) {
		posterior_prob[val]= exp(posterior_prob[val]-max_log_prob)*beta[val];
		summation += posterior_prob[val];
	}
	int sample_val = rand.SampleMultinomial(posterior_prob, summation);
	Instantiation inst = partial;
	inst.SetVal(v, sample_val);
	return inst;
	
}

} // end of ctbn namespace
