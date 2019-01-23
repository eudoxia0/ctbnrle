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
#include "markov.h"
#include "ctbndyn.h"
#include "markovdyn.h"
#include "params.h"
#include <vector>



namespace ctbn {

using namespace std;

SOBJCLASSDEF(Markov)

Markov::Markov(RV *startdist, Dynamics *dyn) : Process() {
	p0 = startdist;
	d = dyn;
}

Markov::Markov(istream &is) : Process() {
	p0 = NULL;
	d = NULL;
	LoadOld(is);
}

Markov::Markov(const Markov &m) : Process() {
	p0 = m.p0->Clone();
	d = m.d->Clone();
}

Markov::~Markov() {
	if (p0!=NULL) delete p0;
	if (d!=NULL) delete d;
}

Markov &Markov::operator=(const Markov &m) {
	if (&m==this) return *this;
	if (p0!=NULL) delete p0;
	if (d!=NULL) delete d;
	p0 = m.p0->Clone();
	d = m.d->Clone();
	return *this;
}

Markov *Markov::Clone() const {
	return new Markov(*this);
}

void Markov::Mult(const Process *p) {
	const Markov *m = dynamic_cast<const Markov *>(p);
	p0->MultBy(m->p0);
	d->Mult(m->d);
}

SS *Markov::BlankSS() const {
	return new MarkovSS(p0->BlankSS(),d->BlankSS());
}

void Markov::Maximize(const SS *ss) {
	const MarkovSS *mss = dynamic_cast<const MarkovSS *>(ss);
	p0->Maximize(mss->p0ss);
	d->Maximize(mss->dss);
}

void Markov::LoadOld(istream &is) {
	if (p0!=NULL) delete p0;
	if (d!=NULL) delete d;
	p0 = dynamic_cast<RV *>(StreamObj::LoadOldPtr(is));
	d = dynamic_cast<Dynamics *>(StreamObj::LoadOldPtr(is));
}

void Markov::serial_preload() {
	if (p0!=NULL) delete p0;
	if (d!=NULL) delete d;
	p0 = NULL; d = NULL;
}

void Markov::SaveOld(ostream &os) const {
	p0->SaveOldPtr(os);
	os << os.fill();
	d->SaveOldPtr(os);
}

void Markov::Sample(Trajectory &tr, Random &rand) const {
	Instantiation x0;
	p0->Sample(x0,rand);
	tr.AddTransition(x0,tr.TimeBegin());
	d->SampleTrajectory(tr,tr.TimeBegin(),rand);
	//dynamic_cast<CTBNDyn *> (d)->SampleTrajectory(tr,tr.TimeBegin(),rand);
}

// return value is the total llh on a single trajectory
// also pass out the p0llh and dynllh seperately
double Markov::LLH(const Trajectory &tr, double w, 
			double &p0llh, double &dynllh) const {
	MarkovSS *mss = dynamic_cast<MarkovSS *>(BlankSS());
	const Context &c = d->Domain() + d->CondDomain();
	const Context &owndomain = d->Domain(); 
	Instantiation curri = tr.Values(c, tr.TimeBegin());
	p0->AddSS(curri, mss->p0ss, w);
	Trajectory::Index index = tr.Begin(c);
	while(!index.Done()) { 
		double t = index.Time();
		double deltat = index.DeltaT();
		d->AddSS(curri, t, deltat, mss->dss, w);
		int sign = index.TestInc(owndomain);
		if (!index.Done()) {
			Instantiation nexti = index.Values();
			t = index.Time();
			if (sign==2)
				d->AddTransSS(curri, nexti, t, mss->dss, w);
			curri = nexti;
		}
	}
	mss->Scale(w);
	p0llh = p0->LLH(mss->p0ss);
	dynllh = d->LLH(mss->dss);
	delete mss;
	return p0llh + dynllh;
}

double Markov::LLH(const vector<Trajectory> &tr, 
			const vector<double> &w) const {
	MarkovSS *mss = dynamic_cast<MarkovSS *>(BlankSS());
	const Context &c = d->Domain() + d->CondDomain();
	const Context &owndomain = d->Domain();
	double sumw = 0.0;
	for(unsigned int i=0; i<tr.size(); i++) {
		Instantiation curri = tr[i].Values(c, tr[i].TimeBegin());
		p0->AddSS(curri, mss->p0ss, w[i]);
		Trajectory::Index index = tr[i].Begin(c);
		while(!index.Done()) { 
			double t = index.Time();
			double deltat = index.DeltaT();
			d->AddSS(curri, t, deltat, mss->dss, w[i]);

			int sign = index.TestInc(owndomain);
			if(!index.Done()) {
				Instantiation nexti = index.Values();
				t = index.Time();
				if (sign==2)
					d->AddTransSS(curri, nexti, 
							t, mss->dss, w[i]);
				curri = nexti;
			}
		}
		sumw += w[i];
	}
	mss->Scale(sumw);
	double llh = LLH(mss);
	delete mss;
	return llh;
}


double Markov::LLH(const SS *ss) const {
	const MarkovSS *mss = dynamic_cast<const MarkovSS*>(ss);
	return p0->LLH(mss->p0ss) + d->LLH(mss->dss);
}

SS* Markov::SuffStats(const vector<Trajectory> &tr, 
			const vector<double> &w) const {
	MarkovSS *mss = dynamic_cast<MarkovSS *>(BlankSS());

	mss->p0ss = p0->SuffStats(tr,w);	
	mss->dss = d->SuffStats(tr,w);

	for(unsigned int i=0; i<tr.size(); i++) {
		if(isnan(w[i])!=0) {
			cout << "Suff Stats NaN:" << endl;
			tr[i].Draw(cout);
			cout << endl;
		}
	}
	return mss;
}

void Markov::Scramble(double a, double b, double alpha, double degree, 
			Random &rand) {
	d->Scramble(a, b, alpha, degree, rand);
	p0->Scramble(alpha, degree, rand);
}

double Markov::GetScore(double numTrans, double amtTime, 
			const SS* ss) const {
	double p0Score = p0->GetScore(numTrans, 
		dynamic_cast<const MarkovSS* > (ss)->p0ss);
	double dynScore = d->GetScore(numTrans, amtTime, 
		dynamic_cast<const MarkovSS* > (ss)->dss);
	return p0Score + dynScore;
}

double MarkovSS::NodeSS(int id, int val1, int val2, 
				const Instantiation &cond) const {
	return dynamic_cast<const CTBNDynSS *>(dss)->
						NodeSS(id, val1, val2, cond);
}
SOBJCLASSDEF(MarkovSS)

MarkovSS::MarkovSS(std::istream &is) {
	p0ss = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
	dss = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
}

MarkovSS::MarkovSS(SS *rvss, SS *dynss) {
	p0ss = rvss;
	dss = dynss;
}

MarkovSS::MarkovSS(const MarkovSS &mss) : SS() {
	p0ss = mss.p0ss->Clone();
	dss = mss.dss->Clone();
}

MarkovSS &MarkovSS::operator=(const MarkovSS &mss) {
	if (&mss!=this) {
		if (p0ss != NULL) delete p0ss;
		if (dss != NULL) delete dss;
		p0ss = mss.p0ss->Clone();
		dss = mss.dss->Clone();
	}
	return *this;
}

MarkovSS::~MarkovSS() {
	if (p0ss != NULL) delete p0ss;
	if (dss != NULL) delete dss;
}

MarkovSS *MarkovSS::Clone() const {
	return new MarkovSS(*this);
}

void MarkovSS::LoadOld(istream &is) {
	if (p0ss != NULL) delete p0ss;
	if (dss != NULL) delete dss;
	p0ss = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
	dss = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
}

void MarkovSS::serial_preload() {
	if (p0ss != NULL) delete p0ss;
	if (dss != NULL) delete dss;
	p0ss = NULL;
	dss = NULL;
}

void MarkovSS::SaveOld(ostream &os) const {
	p0ss->SaveOldPtr(os);
	os << os.fill();
	dss->SaveOldPtr(os);
}

void MarkovSS::Scale(double w) {
	p0ss->Scale(w);
	dss->Scale(w);
}

void MarkovSS::AddSS(const SS* nss, double w) {
	const MarkovSS* mss = dynamic_cast<const MarkovSS*>(nss);
	p0ss->AddSS(mss->p0ss, w);
	dss->AddSS(mss->dss, w);
}

} // end of ctbn namespace
