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
#include "multisimple.h"
#include "extramath.h"
#include <numeric>
#include <algorithm>



namespace ctbn {

using namespace std;

SOBJCLASSDEF(MultiSimple)
SOBJCLASSDEF(MultiLogSimple)
SOBJCLASSDEF(MultiZSimple)
SOBJCLASSDEF(SparseMultiZSimple)
SOBJCLASSDEF(MultiSimpleSS)

MultiSimple::MultiSimple(int n) : RVSimple(), theta(n,0.0) {
}

MultiLogSimple::MultiLogSimple(int n) : RVSimple(), logtheta(n,0.0) {
}

MultiZSimple::MultiZSimple(int n) : RVSimple(), theta(n,0.0) {
	logz = 0.0;
}

SparseMultiZSimple::SparseMultiZSimple(int n) : MultiZSimple(n), ix(n) {
	//iota(ix.begin(),ix.end(),0); // set to 0, 1, 2, ...
	for(int i=0;i<n;i++) ix[i] = i;
}

//---

MultiZSimple::MultiZSimple(const vectr &dist, double logfactor) : RVSimple(),
			theta(dist), logz(logfactor) {
}

SparseMultiZSimple::SparseMultiZSimple(const vectr &dist,
		const vector<int> &indexes, double logfactor) :
			MultiZSimple(dist,logfactor), ix(indexes) {
}

//---

MultiSimple::MultiSimple(istream &is) : RVSimple(is) {
	is >> theta;
}

MultiLogSimple::MultiLogSimple(istream &is) : RVSimple(is) {
	is >> logtheta;
}

MultiZSimple::MultiZSimple(istream &is) : RVSimple(is) {
	is >> theta >> logz;
}

SparseMultiZSimple::SparseMultiZSimple(istream &is) : MultiZSimple(is) {
	is >> ix;
}

//---

MultiSimple::~MultiSimple() {
}

MultiLogSimple::~MultiLogSimple() {
}

MultiZSimple::~MultiZSimple() {
}

SparseMultiZSimple::~SparseMultiZSimple() {
}
//---

void MultiSimple::LoadOld(istream &is) {
	RVSimple::LoadOld(is);
	is >> theta;
}

void MultiLogSimple::LoadOld(istream &is) {
	RVSimple::LoadOld(is);
	is >> logtheta;
}

void MultiZSimple::LoadOld(istream &is) {
	RVSimple::LoadOld(is);
	is >> theta >> logz;
}

void SparseMultiZSimple::LoadOld(istream &is) {
	MultiZSimple::LoadOld(is);
	is >> ix;
}

//---

void MultiSimple::SaveOld(ostream &os) const {
	RVSimple::SaveOld(os);
	os << os.fill() << theta;
}

void MultiLogSimple::SaveOld(ostream &os) const {
	RVSimple::SaveOld(os);
	os << os.fill() << logtheta;
}

void MultiZSimple::SaveOld(ostream &os) const {
	RVSimple::SaveOld(os);
	os << os.fill() << theta << os.fill() << logz;
}

void SparseMultiZSimple::SaveOld(ostream &os) const {
	MultiZSimple::SaveOld(os);
	os << os.fill() << ix;
}

//---

MultiSimple *MultiSimple::Clone() const {
	return new MultiSimple(*this);
}

MultiLogSimple *MultiLogSimple::Clone() const {
	return new MultiLogSimple(*this);
}

MultiZSimple *MultiZSimple::Clone() const {
	return new MultiZSimple(*this);
}

SparseMultiZSimple *SparseMultiZSimple::Clone() const {
	return new SparseMultiZSimple(*this);
}
//---

double MultiSimple::Prob(int ind, bool log) const {
	if (log) return ::log(theta[ind]);
	else return theta[ind];
}

double MultiLogSimple::Prob(int ind, bool log) const {
	if (log) return logtheta[ind];
	else return exp(logtheta[ind]);
}

double MultiZSimple::Prob(int ind, bool log) const {
	if (log) return logz+::log(theta[ind]);
	else return exp(logz)*theta[ind];
}

double SparseMultiZSimple::Prob(int ind, bool log) const {
	vector<int>::const_iterator spot = lower_bound(ix.begin(),ix.end(),ind);
	if (spot==ix.end() || *spot!=ind) return log ? -INFINITY : 0;
	return MultiZSimple::Prob(spot-ix.begin(),log);
}

//---

double MultiSimple::Normalize() {
	return theta.normalize();
}

double MultiLogSimple::Normalize() {
	double lmax = logtheta.max();
	double sum = 0.0;
	int n = logtheta.getm();
	for(int i=0;i<n;i++)
		if (isfinite(logtheta[i])) sum += exp(logtheta[i]-lmax);
	lmax -= log(sum);
	for(int i=0;i<n;i++)
		if (isfinite(logtheta[i])) logtheta[i] += lmax;
	return exp(-lmax);
}

double MultiZSimple::Normalize() {
	/*
	double sum = theta.sum();
	double ret = exp(logz)*sum;
	logz = -log(sum);
	*/
	double ret = theta.normalize();
	logz = 0.0;
	return ret;
}

//---

double MultiSimple::Sum(bool log) const {
	if (log) return ::log(theta.sum());
	else return theta.sum();
}

double MultiLogSimple::Sum(bool log) const {
	double lmax = logtheta.max();
	double sum = 0.0;
	int n = logtheta.getm();
	for(int i=0;i<n;i++)
		if (isfinite(logtheta[i])) sum += exp(logtheta[i]-lmax);
	lmax -= ::log(sum);
	if (log) return lmax;
	else return exp(lmax);
}

double MultiZSimple::Sum(bool log) const {
	if (log) return logz+::log(theta.sum());
	else return exp(logz)*theta.sum();
}

//---

void MultiSimple::Restrict(const vector<int> &ind) {
	int n = theta.getm();
	int indn = ind.size();
    unsigned int j=0;
	for(int i=0;i<n;i++) {
		if (j<ind.size() && i==ind[j]) j++;
		else theta[i] = 0.0;
	}
}

void MultiLogSimple::Restrict(const vector<int> &ind) {
	int n = logtheta.getm();
	int indn = ind.size();
    unsigned int j=0;
	for(int i=0;i<n;i++) {
		if (j<ind.size() && i==ind[j]) j++;
		else logtheta[i] = -INFINITY; 
	}
}

void MultiZSimple::Restrict(const vector<int> &ind) {
	int n = theta.getm();
	int indn = ind.size();
    unsigned int j=0;
	for(int i=0;i<n;i++) {
		if (j<ind.size() && i==ind[j]) j++;
		else theta[i] = 0;
	}
	SetZ();
}

void SparseMultiZSimple::Restrict(const vector<int> &ind) {
	int n = ix.size();
	int indn = ind.size();
	vectr oldt(theta);
	theta.resize(indn);
	theta = 0.0;
	int i=0,j=0;
	while(i<n && j<indn) {
		if (ind[j]==ix[i]) {
			theta[j] = oldt[i];
			++j; ++i;
		} else if (ind[j]<ix[i]) {
			theta[j] = 0.0;
			++j;
		} else {
			++i;
		}
	}
	ix = ind;
	SetZ();
}

//---

void MultiSimple::Reindex(const vector<vector<int> > &ind) {
	vectr oldt(theta);
	int n = ind.size();
	theta.resize(n);
	for(int i=0;i<n;i++) {
		double sum = 0.0;
		int nn = ind[i].size();
		for(int j=0;j<nn;j++)
			sum += oldt[ind[i][j]];
		theta[i] = sum;
	}
}

void MultiLogSimple::Reindex(const vector<vector<int> > &ind) {
	double lmax = logtheta.max();
	vectr oldt(logtheta);
	int n = ind.size();
	logtheta.resize(n);
	for(int i=0;i<n;i++) {
		double sum = 0.0;
		int nn = ind[i].size();
		for(int j=0;j<nn;j++)
			sum += exp(oldt[ind[i][j]]-lmax);
		logtheta[i] = lmax + log(sum);
	}
}

void MultiZSimple::Reindex(const vector<vector<int> > &ind) {
	vectr oldt(theta);
	int n = ind.size();
	theta.resize(n);
	for(int i=0;i<n;i++) {
		double sum = 0.0;
		int nn = ind[i].size();
		for(int j=0;j<nn;j++)
			sum += oldt[ind[i][j]];
		theta[i] = sum;
	}
	SetZ();
}

void SparseMultiZSimple::Reindex(const vector<vector<int> > &ind) {
	vectr oldt(theta);
	int n = ind.size();
	int oldn = ix.size();
	vector<int> oldix(ix);

	int newn = 0;
	for(int i=0;i<n;i++)
		if (ind[i].size()>0) ++newn;
	theta.resize(newn);
	ix.resize(newn,0);

	int kk = 0;
	for(int k=0;k<n;++k) {
		int indn = ind[k].size();
		if (indn==0) continue;
		double sum = 0.0;
		int i=0,j=0;
		while(i<oldn && j<indn) {
			if (ind[k][j]==oldix[i]) {
				sum += oldt[i];
				++j; ++i;
			} else if (ind[k][j]<oldix[i]) {
				++j;
			} else {
				++i;
			}
		}
		theta[kk] = sum;
		ix[kk] = k;
		++kk;
	}
}


//---

void MultiSimple::MultBy(const RVSimple *x) {
	const MultiSimple *ms = dynamic_cast<const MultiSimple *>(x);
	theta.multby(ms->theta);
}

void MultiLogSimple::MultBy(const RVSimple *x) {
	const MultiLogSimple *ms = dynamic_cast<const MultiLogSimple *>(x);
	logtheta += ms->logtheta;
}

void MultiZSimple::MultBy(const RVSimple *x) {
	const MultiZSimple *ms = dynamic_cast<const MultiZSimple *>(x);
	theta.multby(ms->theta);
	logz += ms->logz;
	SetZ();
}

//---

void MultiSimple::Add(const RVSimple *x) {
	const MultiSimple *ms = dynamic_cast<const MultiSimple *>(x);
	theta += ms->theta;
}

void MultiLogSimple::Add(const RVSimple *x) {
	const MultiLogSimple *ms = dynamic_cast<const MultiLogSimple *>(x);
	int n = logtheta.getm();
	for(int i=0;i<n;i++) {
		if (logtheta[i]<ms->logtheta[i])
			logtheta[i] = ms->logtheta[i]
				+ ::log(exp(logtheta[i]-ms->logtheta[i]) + 1.0);
		else 
			logtheta[i] = logtheta[i]
				+ ::log(exp(ms->logtheta[i]-logtheta[i]) + 1.0);
	}
}

void MultiZSimple::Add(const RVSimple *x) {
	const MultiZSimple *ms = dynamic_cast<const MultiZSimple *>(x);
	if (ms->logz>logz) {
		theta = theta*exp(logz-ms->logz) + ms->theta;
		logz = ms->logz;
	} else {
		theta += ms->theta*exp(ms->logz-logz);
	}
	SetZ();
}

//--- 

void MultiSimple::Add(const RVSimple *x, double w) {
	const MultiSimple *ms = dynamic_cast<const MultiSimple *>(x);
	theta.add(ms->theta,w);
}

void MultiLogSimple::Add(const RVSimple *x, double w) {
	const MultiLogSimple *ms = dynamic_cast<const MultiLogSimple *>(x);
	int n = logtheta.getm();
	double logw = ::log(w);
	for(int i=0;i<n;i++) {
		double xlt = ms->logtheta[i]+logw;
		if (logtheta[i]<xlt)
			logtheta[i] = xlt + ::log(exp(logtheta[i]-xlt) + 1.0);
		else 
			logtheta[i] = logtheta[i] + 
				::log(exp(xlt-logtheta[i]) + 1.0);
	}
}

void MultiZSimple::Add(const RVSimple *x, double w) {
	const MultiZSimple *ms = dynamic_cast<const MultiZSimple *>(x);
	double xlz = ms->logz+::log(w);
	if (xlz>logz) {
		theta = theta*exp(logz-xlz) + ms->theta;
		logz = xlz;
	} else {
		theta += ms->theta*exp(xlz-logz);
	}
	SetZ();
}

//---

void MultiSimple::Mult(double x) {
	theta *= x;
}

void MultiLogSimple::Mult(double x) {
	logtheta += log(x);
}

void MultiZSimple::Mult(double x) {
	logz += log(x);
}

//---

MultiSimpleSS *MultiSimple::BlankSS() const {
	return new MultiSimpleSS(theta.getm());
}

MultiSimpleSS *MultiLogSimple::BlankSS() const {
	return new MultiSimpleSS(logtheta.getm());
}

MultiSimpleSS *MultiZSimple::BlankSS() const {
	return new MultiSimpleSS(theta.getm());
}

//---

void MultiSimple::AddSS(int x, SS *ss, double w) const {
	MultiSimpleSS * msss = dynamic_cast<MultiSimpleSS *>(ss);
	msss->c[x] += w;
}

void MultiLogSimple::AddSS(int x, SS *ss, double w) const {
	MultiSimpleSS * msss = dynamic_cast<MultiSimpleSS *>(ss);
	msss->c[x] += w;
}

void MultiZSimple::AddSS(int x, SS *ss, double w) const {
	MultiSimpleSS * msss = dynamic_cast<MultiSimpleSS *>(ss);
	msss->c[x] += w;
}

void SparseMultiZSimple::AddSS(int x, SS *ss, double w) const {
	MultiSimpleSS * msss = dynamic_cast<MultiSimpleSS *>(ss);
	msss->c[x] += w;
}

//---

void MultiSimple::AddExpSS(SS *ss, double w) const {
	MultiSimpleSS * msss = dynamic_cast<MultiSimpleSS *>(ss);
	msss->c += w*theta;
}

void MultiLogSimple::AddExpSS(SS *ss, double w) const {
	MultiSimpleSS * msss = dynamic_cast<MultiSimpleSS *>(ss);
	int n = logtheta.getm();
	for(int i=0;i<n;i++)
		msss->c[i] += w*exp(logtheta[i]);
}

void MultiZSimple::AddExpSS(SS *ss, double w) const {
	MultiSimpleSS * msss = dynamic_cast<MultiSimpleSS *>(ss);
	int n = theta.getm();
	w *= exp(logz);
	for(int i=0;i<n;i++)
		msss->c[i] += w*theta[i];
}

void SparseMultiZSimple::AddExpSS(SS *ss, double w) const {
	MultiSimpleSS * msss = dynamic_cast<MultiSimpleSS *>(ss);
	int n = theta.getm();
	w *= exp(logz);
	for(int i=0;i<n;i++)
		msss->c[ix[i]] += w*theta[i];
}

//---
void MultiSimple::AddSS(const SS *toadd, const RVSimple *rvs,
		const std::vector<std::vector<int> > &mapping,
		SS *ss, double w) const {
	const MultiSimple *ms = dynamic_cast<const MultiSimple *>(rvs);
	assert(ms!=NULL); // otherwise, need double-dispatch and more code
	MultiSimpleSS *msss = dynamic_cast<MultiSimpleSS *>(ss);
	const MultiSimpleSS *addmsss = 
				dynamic_cast<const MultiSimpleSS *>(toadd);

	int n = theta.getm();
	for(int i=0;i<n;i++) {
		int li = mapping[i].size();
		for(int ki=0;ki<li;ki++) {
			msss->c[i] +=
				addmsss->c[mapping[i][ki]];
		}
	}
}

void MultiLogSimple::AddSS(const SS *toadd, const RVSimple *rvs,
		const std::vector<std::vector<int> > &mapping,
		SS *ss, double w) const {
	const MultiLogSimple *ms = dynamic_cast<const MultiLogSimple *>(rvs);
	assert(ms!=NULL); // otherwise, need double-dispatch and more code
	MultiSimpleSS *msss = dynamic_cast<MultiSimpleSS *>(ss);
	const MultiSimpleSS *addmsss = 
				dynamic_cast<const MultiSimpleSS *>(toadd);

	int n = logtheta.getm();
	for(int i=0;i<n;i++) {
		int li = mapping[i].size();
		for(int ki=0;ki<li;ki++) {
			msss->c[i] +=
				addmsss->c[mapping[i][ki]];
		}
	}
}

void MultiZSimple::AddSS(const SS *toadd, const RVSimple *rvs,
		const std::vector<std::vector<int> > &mapping,
		SS *ss, double w) const {
	const MultiZSimple *ms = dynamic_cast<const MultiZSimple *>(rvs);
	assert(ms!=NULL); // otherwise, need double-dispatch and more code
	MultiSimpleSS *msss = dynamic_cast<MultiSimpleSS *>(ss);
	const MultiSimpleSS *addmsss = 
				dynamic_cast<const MultiSimpleSS *>(toadd);

	int n = theta.getm();
	for(int i=0;i<n;i++) {
		int li = mapping[i].size();
		for(int ki=0;ki<li;ki++) {
			msss->c[i] +=
				addmsss->c[mapping[i][ki]];
		}
	}
}

void SparseMultiZSimple::AddSS(const SS *toadd, const RVSimple *rvs,
		const std::vector<std::vector<int> > &mapping,
		SS *ss, double w) const {
	const SparseMultiZSimple *ms = 
			dynamic_cast<const SparseMultiZSimple *>(rvs);
	assert(ms!=NULL); // otherwise, need double-dispatch and more code
	MultiSimpleSS *msss = dynamic_cast<MultiSimpleSS *>(ss);
	const MultiSimpleSS *addmsss = 
				dynamic_cast<const MultiSimpleSS *>(toadd);

	int n = theta.getm();
	for(int i=0;i<n;i++) {
		int li = mapping[ix[i]].size();
		for(int ki=0;ki<li;ki++) {
			msss->c[ix[i]] +=
				addmsss->c[mapping[ix[i]][ki]];
		}
	}
}





//---

int MultiSimple::Sample(Random &rand) const {
//	return rand.SampleMultinomial(theta,theta.getm(),theta.sum());
	return rand.SampleMultinomial(theta,theta.sum());
}

int MultiLogSimple::Sample(Random &rand) const {
	int n = logtheta.getm();
	vector<double> addw(n,0.0);
	double lmax = logtheta.max();
	double sum = 0.0;
	for(int i=0;i<n;i++)
		if (isfinite(logtheta[i])) addw[i] = sum;
		else addw[i] = (sum += exp(logtheta[i]-lmax));
	return rand.SampleAddWeights(addw);
}

int MultiZSimple::Sample(Random &rand) const {
//	return rand.SampleMultinomial(theta,theta.getm(),theta.sum());
	return rand.SampleMultinomial(theta,theta.sum());
}

int SparseMultiZSimple::Sample(Random &rand) const {
	return ix[MultiZSimple::Sample(rand)];
}

//---

void MultiSimple::Maximize(const SS *ss) {
	const MultiSimpleSS * msss = dynamic_cast<const MultiSimpleSS *>(ss);
	theta = msss->c;
	double sum = theta.sum();
	if (sum>0) theta /= sum;
	else theta = 1.0/theta.getm();
}

void MultiLogSimple::Maximize(const SS *ss) {
	const MultiSimpleSS *msss = dynamic_cast<const MultiSimpleSS *>(ss);
	logtheta = msss->c;
	double sum = logtheta.sum();
	if (sum>0) sum = log(sum);
	else logtheta = 1.0/logtheta.getm();
	int n = logtheta.getm();
	for(int i=0;i<n;i++)
		logtheta[i] = log(logtheta[i])-sum;
}

void MultiZSimple::Maximize(const SS *ss) {
	const MultiSimpleSS *msss = dynamic_cast<const MultiSimpleSS *>(ss);
	theta = msss->c;
	double sum = theta.sum();
	if (sum>0) logz = -log(sum);
	else {
		theta = 1.0;
		logz = -log((double) theta.getm());
	}
}

void SparseMultiZSimple::Maximize(const SS *ss) {
	MultiZSimple::Maximize(ss);
	int n = theta.getm();
	ix.resize(n);
	// iota(ix.begin(),ix.end(),0); // set to 0, 1, 2, ...
	for(int i=0;i<n;i++) ix[i] = i;
}

//---

void MultiSimple::GetDist(vectr &d, double &logfactor) const {
	d = theta;
	logfactor = 1.0;
}

void MultiLogSimple::GetDist(vectr &d, double &logfactor) const {
	d = logtheta;
	int n = d.getm();
	
	logfactor = d[0];
	for(int i=1;i<n;i++)
		if (d[i]>logfactor) logfactor = d[i];
	for(int i=0;i<d.getm();i++)
		d[i] = exp(d[i]-logfactor);
}

void MultiZSimple::GetDist(vectr &d, double &logfactor) const {
	d = theta;
	logfactor = logz;
}

void SparseMultiZSimple::GetDist(vectr &d, double &logfactor) const {
	logfactor = logz;
	int nn = 0;
	int n = ix.size();
	for(int i=0;i<n;i++) if (ix[i]>nn) nn = ix[i];
	if (d.getm()<n) d.resize(n);
	d = 0;
	for(int i=0;i<n;i++)
		d[ix[i]] = theta[i];
}

void SparseMultiZSimple::GetSparseDist(vectr &d, vectr &ind, double &logfactor) const {
	logfactor = logz;
	int n = ix.size();
	if (d.getm()!=n) d.resize(n);
	if (ind.getm()!=n) ind.resize(n);
	for (int i=0; i<n; i++){
		ind[i] = ix[i];
		d[i] = theta[i];
	}
}

//---

void MultiSimple::SetDist(const vectr &d, double logfactor) {
	theta = d*exp(logfactor);
}

void MultiLogSimple::SetDist(const vectr &d, double logfactor) {
	int n = d.getm();
	for(int i=0;i<d.getm();i++)
		logtheta[i] = log(d[i])+logfactor;
}

void MultiZSimple::SetDist(const vectr &d, double logfactor) {
	theta = d;
	logz = logfactor;
}

void SparseMultiZSimple::SetDist(const vectr &d, double logfactor) {
	logz = logfactor;
	int n = d.getm();
	unsigned int nn = 0;
	for(int i=0;i<n;i++) if (d[i]!=0.0) nn = i;
	if (ix.size()!=nn) ix.resize(nn);
	for(int i=0,j=0;i<n;i++)
		if (d[i]!=0.0) {
			ix[j] = i;
			theta[j] = d[i];
			j++;
		}
}

void SparseMultiZSimple::SetSparseDist(const vectr &d, const vectr &ind, double logfactor) {
//	logfactor = logz;
	unsigned int n = d.getm();
	if (ix.size()!=n) {
		ix.resize(n);
		theta.resize(n);
	}
	for (unsigned int i=0; i<n; i++){
		theta[i] = d[i];
		ix[i] = ind[i];
	}
}

//---

void MultiSimple::MakeUniform() {
	int n = theta.getm();
	for(int i=0;i<n;i++) theta[i] = 1.0;
		//if (theta[i]>0) theta[i] = 1.0;
}
void MultiLogSimple::MakeUniform() {
	int n = logtheta.getm();
	for(int i=0;i<n;i++) logtheta[i] = 0.0;
	//if (isfinite(logtheta[i])) logtheta[i] = 0.0;
}
void MultiZSimple::MakeUniform() {
	int n = theta.getm();
	for(int i=0;i<n;i++) theta[i] = 1.0;
		//if (theta[i]>0) theta[i] = 1.0;
	logz = 0.0;
}

//---
void MultiSimple::Scramble(double alpha, double degree, Random &rand) {
	int n = theta.getm();
	double *in = new double[n];
	double *newtheta = new double[n];
	for(int i=0;i<n;i++) in[i] = alpha;
	rand.SampleDirichlet(in,n,newtheta);
	delete []in;
	for (int j = 0; j < n; ++j)
		theta[j] = newtheta[j] * degree + theta[j] * (1.0 - degree);
	delete []newtheta;
}

void MultiZSimple::Scramble(double alpha, double degree, Random &rand) {
	int n = theta.getm();
	double *in = new double[n];
	for(int i=0;i<n;i++) in[i] = alpha;
	double *newtheta = new double[n];
	rand.SampleDirichlet(in,n,newtheta);
	delete []in;

	for (int j = 0; j < n; ++j)
		theta[j] = newtheta[j] * degree + theta[j] * exp(logz)* (1.0 - degree);
	logz = 0.0;
	delete []newtheta;
}

void MultiLogSimple::Scramble(double alpha, double degree, Random &rand) {
	int n = logtheta.getm();
	double *in = new double[n];
	vectr original(logtheta);
	for(int i=0;i<n;i++) in[i] = alpha;
	for(int i=0;i<n;i++) 
		if (isfinite(logtheta[i])) 
			original[i] = exp(original[i]);
		else
			original[i] = 0.0;
	double *newtheta = new double[n];
	rand.SampleDirichlet(in,n,newtheta);
	delete []in;
	for (int j = 0; j < n; ++j)
		logtheta[j] = log(newtheta[j] * degree + original[j] * (1.0 - degree));
	delete []newtheta;
}

//---
double MultiSimple::LLH(const SS *ss) const {
	const MultiSimpleSS *rvss = dynamic_cast<const MultiSimpleSS *>(ss);
	double llh = 0.0;
	int n = theta.length();
	for(int i=0; i<n; i++)
		if(rvss->c[i]!=0.0) llh += rvss->c[i] * log(theta[i]);
	return llh;
}

double MultiZSimple::LLH(const SS *ss) const {
	const MultiSimpleSS *rvss = dynamic_cast<const MultiSimpleSS *>(ss);
	double llh = 0.0;
	int n = theta.length();
	for(int i=0; i<n; i++)
		if(rvss->c[i]>0.0) llh += rvss->c[i] * (log(theta[i]) + logz);
	return llh;
}

double MultiLogSimple::LLH(const SS *ss) const {
	const MultiSimpleSS *rvss = dynamic_cast<const MultiSimpleSS *>(ss);
	double llh = 0.0;
	int n = logtheta.length();
	for(int i=0; i<n; i++)
		if(rvss->c[i]!=0.0) llh += rvss->c[i] * logtheta[i];
	return llh;
}

double SparseMultiZSimple::LLH(const SS *ss) const {
	const MultiSimpleSS *rvss = dynamic_cast<const MultiSimpleSS *>(ss);
	double llh = 0.0;
	int n = rvss->c.length();
	//int m = ix.length();
	int m = ix.size();
	for(int i=0,j=0; i<n; i++) {
		while (j<m && ix[j]<i) ++j;
		if(rvss->c[i]>0.0) {
			if (ix[j]==i) llh += rvss->c[i] * (log(theta[j]) + logz);
			else return log(0.0);
		}
	}
	return llh;
}

//double MultiLogSimple::LLH(const SS *ss) const
//---

void SparseMultiZSimple::SetSparseDist(const vectr &d, double logfactor) {
	MultiZSimple::SetDist(d,logfactor);
}

void SparseMultiZSimple::SetSparseIndexes(const vector<int> &ind) {
	ix = ind;
}

void SparseMultiZSimple::GetSparseDist(vectr &d, double &logfactor) const {
	MultiZSimple::GetDist(d,logfactor);
}

//---

void MultiZSimple::SetZ() {
	// set so that maximum is 1 (could set so that sum=1 instead,
	//     but this seems to work better)
	double s = theta.max();
	if (s>0) {
		theta /= s;
		logz += log(s);
	}
}

//---
double RVSimple::GetScore(double numTrans, const SS* ss) const {
	const MultiSimpleSS *rvss = dynamic_cast<const MultiSimpleSS *>(ss);
	const vectr& c = rvss->c;

	unsigned int size = c.length();

/*
	double alpha_xx_u(numTrans / ((size-1)*size));
	double alpha_x_u(alpha_xx_u*size);
	double localScore(0);
	localScore += lngamma(alpha_x_u) - lngamma(alpha_x_u + c.sum());
*/
    double alpha_xx_u(numTrans / size);
    double alpha_x_u(numTrans);
    double localScore = lngamma(alpha_x_u) - lngamma(alpha_x_u + c.sum());
	
	for(unsigned int i = 0; i < size; i++)
		localScore += lngamma(alpha_xx_u + c[i]) - lngamma(alpha_xx_u);
	return localScore;
}
		
		
//---------

MultiSimpleSS::MultiSimpleSS(int n) : SS(), c(n,0.0) {
}

MultiSimpleSS::MultiSimpleSS(istream &is) : SS(is) {
	is >> c;
}

MultiSimpleSS::MultiSimpleSS(const MultiSimpleSS &msss) : SS(*this), c(msss.c) {
}

MultiSimpleSS::~MultiSimpleSS() {
}

MultiSimpleSS *MultiSimpleSS::Clone() const {
	return new MultiSimpleSS(*this);
}

void MultiSimpleSS::LoadOld(istream &is) {
	SS::LoadOld(is);
	is >> c;
}

void MultiSimpleSS::SaveOld(ostream &os) const {
	SS::SaveOld(os);
	os << os.fill() << c;
}

void MultiSimpleSS::Scale(double w) {
	int n = c.length();
	for(int i=0; i<n; i++)
		c[i] /= w;
}


void MultiSimpleSS::AddSS(const SS* nss, double w) { 
	const MultiSimpleSS *mss = dynamic_cast<const MultiSimpleSS*>(nss);
	int n = c.length();
	for(int i=0; i<n; i++)
		c[i] += mss->c[i] * w;

}

} // end of ctbn namespace
