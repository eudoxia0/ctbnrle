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
#include "samplinginf.h"
#include "markov.h"
#include "markovdyn.h"
#include "multirv.h"
#include "ctbndyn.h"
#include "utils.h"



namespace ctbn {

using namespace std;

SamplingInfAbs::SamplingInfAbs() {
	p = NULL;
	numofsamples = 100;
	timelimit = 0.0;
}

SamplingInf::SamplingInf() {
}

SamplingPreInf::SamplingPreInf(const vector<QueryCalculator *> &queries,
		int workingsetsize)
			: qs(queries), ans(queries.size(),0) {
	maxlnw = -INFINITY;
	ttlw = 0;
	sampled = false;
	setsize = workingsetsize;
}

SamplingInfAbs::~SamplingInfAbs() {
	if(sampler!=NULL) delete sampler;
}

SamplingInf::~SamplingInf() {
}

SamplingPreInf::~SamplingPreInf() {
}

SamplingInf* SamplingInf::Clone() const {
	SamplingInf* siClone = new SamplingInf(*this);
	siClone->SetSampler(this->sampler->Clone());
	return siClone;
}

SamplingPreInf *SamplingPreInf::Clone() const {
	SamplingPreInf* siClone = new SamplingPreInf(*this);
	siClone->SetSampler(this->sampler->Clone());
	return siClone;
}

void SamplingInfAbs::SetProcess(const Process *pr) {
	if(sampler==NULL) {
		cout << "Sampler not set!" << endl;
		exit(1);
	}
	sampler->SetProcess(pr);
	p = dynamic_cast<const Markov*>(pr);
	Reset();
}

void SamplingInfAbs::SetSampler(const Sampler *s) {
	sampler = s->Clone();
	if(p!=NULL) sampler->SetProcess(p);
	Reset();
}

void SamplingInfAbs::SetTrajectory(const Trajectory *tr) {
	if(sampler==NULL) {
		cout << "Sampler not set!" << endl;
		exit(1);
	}
	if(p==NULL) {
		cout << "Process not set!" << endl;
		exit(1);
	}
	sampler->SetTrajectory(tr);
	Reset();
}

void SamplingInfAbs::SetNumSample(int n) {
	numofsamples = n;
}

void SamplingInfAbs::SetTimeLimit(double t, int batchsize) {
	timelimit = t;
	bsize = batchsize;
}

// Sample completions of given evidence trajectory (stored in sampler), and
// associate them with their weights.
void SamplingInfAbs::Sample() {
	if(sampler==NULL || p==NULL) {
		cout << "Sampler or process not set!" << endl;
		exit(1);
	}
	if (timelimit == 0.0) {
		AddSamples(numofsamples,true);
	} else {
		double begin_time = getcputime();
		double end_time = begin_time + timelimit;
		do {
			AddSamples(bsize,false);
		} while (getcputime() < end_time);
	}
}

double SamplingInfAbs::NormalizeSet(vector<double> &w) {
	double maxweight = w[0];
	for(unsigned int i=1; i<w.size(); i++)
		if(w[i]>maxweight) maxweight = w[i];
	for(unsigned int i=0; i<w.size(); i++)
		w[i] = exp(w[i] - maxweight);
	return maxweight;
}

void SamplingInf::Sample() {
	SamplingInfAbs::Sample();
	NormalizeSet(weights);
}

void SamplingInf::AddSamples(int nsamp, bool all) {
	if (all)
		sampler->SampleTrajectories(traj, weights, numofsamples);
	else {
		vector<Trajectory> t;
		vector<double> w;
		sampler->SampleTrajectories(t, w, nsamp);
		traj.insert(traj.end(),t.begin(),t.end());
		weights.insert(weights.end(),w.begin(),w.end());
	}
}

void SamplingPreInf::AddSamples(int nsamp, bool all) {
	while(nsamp) {
		vector<Trajectory> t;
		vector<double> w;
		int n = nsamp>setsize ? setsize : nsamp;
		sampler->SampleTrajectories(t,w,n);
		nsamp -= n;

		double lnw = NormalizeSet(w);
		double f;
		if (lnw > maxlnw) {
			f = 1.0;
			double oldf = exp(maxlnw-lnw);
			ttlw *= oldf;
			maxlnw = lnw;
			for(size_t i=0;i<qs.size();i++)
				ans[i] *= oldf;
		} else
			f = exp(lnw-maxlnw);
		for(int j=0;j<n;j++) {
			double wf = w[j]*f;
			ttlw += wf;
			for(size_t i=0;i<qs.size();i++)
				ans[i] += qs[i]->Calculate(t[j])*wf;
		}
	}
}

void SamplingInf::Reset() {
	traj.clear();
	weights.clear();
}

void SamplingPreInf::Reset() {
	ans = vector<double>(0,qs.size());
	maxlnw = -INFINITY;
	ttlw = 0.0;
	sampled = false;
}

double SamplingInf::Smooth(const Instantiation &x, double t, bool log) {
	QueryProb query(x,t);
	double ret = CalcQuery(query);
	return log ? ::log(ret) : ret;
}

double SamplingPreInf::Smooth(const Instantiation &,double,bool) {
	cout << "not implemented" << endl;
	return 0;
}

double SamplingInf::Filter(const Instantiation &x, double t, bool log) {
	cout << "not implemented" << endl;
	return 0;
}

double SamplingPreInf::Filter(const Instantiation &x, double t, bool log) {
	cout << "not implemented" << endl;
	return 0;
}


double SamplingInf::Prob(double t, bool log) {
	cout << "not implemented" << endl;
	return 0;
}

double SamplingPreInf::Prob(double t, bool log) {
	cout << "not implemented" << endl;
	return 0;
}

void SamplingInf::AddExpSuffStats(SS *ss, double w) {
	if(sampler==NULL) {
		cout << "Sampler not set!" << endl;
		exit(1);
	}

	if(traj.size()==0) Sample();
	double sumweights = 0.0;
	for(unsigned int i=0; i<traj.size(); i++) {
		if(isnan(weights[i])) {
			cerr << "Suff Stats NaN:" << endl;
			//traj[i].Draw(cout);
			cerr << endl;
			cerr << weights << endl;
			exit(0);
		}
		sumweights += weights[i];
	}

	SS *newss = p->SuffStats(traj, weights);
	newss->Scale(sumweights); //by Yu
	ss->AddSS(newss);
	delete newss;
}

void SamplingPreInf::AddExpSuffStats(SS *ss, double w) {
	assert(0); // not implementable
}

void SamplingInf::AddExpSuffStats(const Dynamics *dyn, SS *ss, double w) {
	if(sampler==NULL) {
		cout << "Sampler not set!" << endl;
		exit(1);
	}
	if(traj.size()==0) Sample();
	double sumweights = 0.0;
	for(unsigned int i=0; i<traj.size(); i++) {
		if(isnan(weights[i])) {
			cerr << "Suff Stats NaN:" << endl;
			//traj[i].Draw(cout);
			cerr << endl;
			cerr << weights << endl;
			exit(0);
		}
		sumweights += weights[i];
	}
	SS *newdynss = dyn->SuffStats(traj, weights);
	newdynss->Scale(sumweights);
	ss->AddSS(newdynss);
	delete newdynss;
}

void SamplingPreInf::AddExpSuffStats(const Dynamics *dyn, SS *ss, double w) {
	assert(0); // not implementable
}

void SamplingInf::AddExpSuffStats(const RV *p0, SS *ss, double w) {
	if(sampler==NULL) {
		cout << "Sampler not set!" << endl;
		exit(1);
	}
	if(traj.size()==0) Sample();
	double sumweights = 0.0;
	for(unsigned int i=0; i<traj.size(); i++) {
		if(isnan(weights[i])) {
			cerr << "Suff Stats NaN:" << endl;
			//traj[i].Draw(cout);
			cerr << endl;
			cerr << weights << endl;
			exit(0);
		}
		sumweights += weights[i];
	}
	SS *newp0ss = p0->SuffStats(traj, weights);
	newp0ss->Scale(sumweights);
	ss->AddSS(newp0ss);
	delete newp0ss;
}

void SamplingPreInf::AddExpSuffStats(const RV *p0, SS *ss, double w) {
	assert(0);
}

double SamplingInf::CalcQuery(QueryCalculator &calc) {
	if(sampler==NULL) {
		cout << "Sampler not set!" << endl;
		exit(1);
	}

	if(traj.size()==0) Sample();

	double sumweights = 0;
	for(unsigned int i=0; i<traj.size(); i++)
		sumweights += weights[i];
	//   for(int i=0; i<traj.size(); i++)
	//     weights[i] /= sumweights;

	double retval = 0.0;

	int n = traj.size();
	for(int i=0; i<n; i++)
		retval += weights[i] * calc.Calculate(traj[i]);
	return retval/sumweights;
}

double SamplingPreInf::CalcQuery(QueryCalculator &calc) {
	if (!sampled) {
		Sample();
		sampled = true;
	}
	for(size_t i=0;i<qs.size();i++)
		if (qs[i] == &calc) return ans[i]/ttlw;
	cerr << "SamplePreInf::CalcQuery called with calculator not shown at construction" << endl;
	assert(0);
}

} // end of ctbn namespace
