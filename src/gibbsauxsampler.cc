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
#include "gibbsauxsampler.h"
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
GibbsAuxSampler::GibbsAuxSampler(const Process *process, const Trajectory *evidence,
				 int burnin) : GibbsBase(process,evidence,burnin) {
}

GibbsAuxSampler::~GibbsAuxSampler() {
}

void GibbsAuxSampler::SampleVariableInterval(int v, int x0, int xT, 
						double t0, double tT, Random &rand) const {
	const MarkovDyn *md = mds[v];
	int nv = context.Cardinality(v);
	vector<double> tau;
	vector<Instantiation> inst;
	GetTau(oldtr, t0, tT, v, markovblanket[v], tau, inst);

	vector<double> newtau;
	vector<matrix> prop;
	vector<vectr> lhood;
	matrix eye(nv,nv,0.0);
	for(int i=0;i<nv;i++) eye[i][i] = 1.0;

	vector<int> chvar = children[v].VarList();
	vector<const MarkovDyn *> chmd;
	for(size_t j=0;j<chvar.size();j++) chmd.push_back(mds[chvar[j]]);
	
	vectr Tvec(nv);
	for(size_t i=0;i<tau.size()-1;i++) {
		int val = inst[i].Value(v);
		const matrix &Q = (*md)(inst[i])->Intensity();
		double alpha = 2*Q.diag().absmax();
		double auxrate = alpha+Q[val][val];
		matrix M = eye+Q/alpha;
		vectr qs(nv,0.0);
		for(int k=0;k<nv;k++) {
			inst[i].SetVal(v,k);
			for(size_t j=0;j<chvar.size();j++) {
				int chval = inst[i].Value(chvar[j]);	
				qs[k] += (*chmd[j])(inst[i])->Intensity()[chval][chval];
			}
		}
		inst[i].SetVal(v,val);
		newtau.push_back(tau[i]);
		if (i==0) prop.push_back(eye);
		else if (inst[i-1].Value(v)!=val)
			prop.push_back(M);
		else if (GetTMatrix(v,inst[i-1],inst[i],Tvec))
			prop.push_back(matrix(nv,nv,Tvec));
		else prop.push_back(eye);
		double dt,currt=tau[i],lastt=tau[i];
		while((currt+=(dt=rand.SampleExp(auxrate)))<tau[i+1]) {
			newtau.push_back(currt);
			prop.push_back(M);
			vectr lh(nv); for(int k=0;k<nv;k++) lh[k] = exp(qs[k]*dt);
			lhood.push_back(lh);
			lastt = currt;
		}
		dt = tau[i+1]-lastt;
		vectr lh(nv); for(int k=0;k<nv;k++) lh[k] = exp(qs[k]*dt);
		lhood.push_back(lh);
	}
	newtau.push_back(tau[tau.size()-1]);

	vectr curra(nv,0.0);
	vector<vectr> alpha;
	if (x0!=-1) curra[x0] = 1.0;
	else {
		double mlp = -INFINITY;
		Instantiation i0 = tr.Values(context,t0);
		for(int k=0;k<nv;k++) {
			i0.SetVal(v,k);
			curra[k] = p0->Prob(i0,true);
			if (curra[k]>mlp) mlp = curra[k];
		}
		for(int k=0;k<nv;k++) curra[k] = exp(curra[k]-mlp);
		curra.normalize();
	}
	for(size_t i=0;i<newtau.size()-1;i++) {
		curra = curra*prop[i];
		curra = curra.dotstar(lhood[i]);
		curra.normalize();
		alpha.push_back(curra);
	}
	int currv = xT!=-1 ? xT : rand.SampleMultinomial(curra);
	for(int i=newtau.size()-2;i>0;i--) {
		vectr w(alpha[i-1]);
		for(int k=0;k<nv;k++) w[k] *= prop[i][k][currv];
		w.normalize();
		int prevv = rand.SampleMultinomial(w);
		if (prevv != currv) {
			tr.AddTransition(v,newtau[i],currv);
			currv = prevv;
		}
	}
	if (x0==-1) tr.AddTransition(v,newtau[0],currv);
}

} // end of ctbn namespace
