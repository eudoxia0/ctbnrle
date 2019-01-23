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
#include "dyncomp.h"
#include "markovsimple.h"



namespace ctbn {

using namespace std;

SOBJCLASSDEF(DynCompSS)

DynCompSS::DynCompSS(int n) : SS(), ss(n) {
}

DynCompSS::DynCompSS(istream &is) : SS(is) {
	LoadOld(is);
}

DynCompSS::DynCompSS(const DynCompSS &dcss) : SS(dcss) {
	int n = dcss.ss.size();
	for(int i=0;i<n;i++)
		ss.push_back(dcss.ss[i]->Clone());
}

DynCompSS &DynCompSS::operator=(const DynCompSS &dcss) {
	if (&dcss==this) return *this;
	SS::operator=(dcss);
	int n = ss.size();
	for(int i=0;i<n;i++) if (ss[i]!=NULL) delete ss[i];
	if (dcss.ss.size()!=ss.size()) ss.resize(dcss.ss.size());
	n = dcss.ss.size();
	for(int i=0;i<n;i++) ss[i] = dcss.ss[i]->Clone();
	return *this;
}

DynCompSS::~DynCompSS() {
	int n = ss.size();
	for(int i=0;i<n;i++)
		if (ss[i]!=NULL) delete ss[i];
}

DynCompSS *DynCompSS::Clone() const {
	return new DynCompSS(*this);
}

void DynCompSS::LoadOld(istream &is) {
	ss.resize(0);
	int n;
	is >> n;
	for(int i=0;i<n;i++) {
		SS *ssptr = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
		ss.push_back(ssptr);
	}
}

void DynCompSS::SaveOld(ostream &os) const {
	int n = ss.size();
	os << n;
	for(int i=0;i<n;i++) {
		os << os.fill();
		ss[i]->SaveOldPtr(os);
	}
}

void DynCompSS::Scale(double w) {
	int n = ss.size();
	for(int i=0; i<n; i++)
		ss[i]->Scale(w);
}

void DynCompSS::AddSS(const SS *nss, double w) {
	int n = ss.size();
	const DynCompSS *dynss = dynamic_cast<const DynCompSS*>(nss);
	for(int i=0; i<n; i++)
		ss[i]->AddSS(dynss->ss[i], w);
}

double DynCompSS::Element(std::vector<int> index, int s1, int s2) const {
	double ret = 0.0;
	int n = ss.size();
	for(int i=0; i<n; i++) {
		ret += dynamic_cast<const MarkovSimpleSS*>(ss[i])->
							Element(s1, s2);
	}
	return ret;
}

SOBJCLASSDEF(DynComp)

DynComp::DynComp(const Context &var, const Context &cvar, 
			const DynSimple *base) :
	Dynamics(var,cvar), impl(cvar.Size()) {
	for(int i=0;i<cvar.Size();i++)
		impl[i] = base->MakeNew(var.Size());
}

DynComp::DynComp(istream &is) :
	Dynamics(is), impl(CondDomain().Size()) {
	int n = impl.size();
	for(int i=0;i<n;i++)
		impl[i] = dynamic_cast<DynSimple *>(StreamObj::LoadOldPtr(is));
}

DynComp::DynComp(const DynComp &dc) : Dynamics(dc), impl(dc.impl.size()) {
	for(unsigned int i=0;i<dc.impl.size();i++)
		impl[i] = dc.impl[i]->Clone();
}

DynComp::~DynComp() {
	for(unsigned int i=0;i<impl.size();i++)
		if (impl[i]!=NULL) delete impl[i];
}

void DynComp::LoadOld(istream &is) {
	Dynamics::LoadOld(is);
	int n = CondDomain().Size();
	for(unsigned int i=0;i<impl.size();i++)
		if (impl[i]!=NULL) delete impl[i];
	impl.resize(n);
	for(int i=0;i<n;i++)
		impl[i] = dynamic_cast<DynSimple *>(StreamObj::LoadOldPtr(is));
}

void DynComp::SaveOld(ostream &os) const {
	Dynamics::SaveOld(os);
	int n = impl.size();
	for(int i=0;i<n;i++) {
		os << os.fill();
		impl[i]->SaveOldPtr(os);
	}
}

DynComp &DynComp::operator=(const DynComp &dc) {
	if (&dc != this) {
		for(unsigned int i=0;i<impl.size();i++)
			if (impl[i]!=NULL) delete impl[i];
		impl.resize(dc.impl.size());
		for(unsigned int i=0;i<impl.size();i++)
			impl[i] = dc.impl[i]->Clone();
	}
	return *this;
}

DynComp *DynComp::Clone() const {
	return new DynComp(*this);
}

void DynComp::Normalize() {
	int n = impl.size();
	for(int i=0;i<n;i++)
		impl[i]->Normalize();
}

void DynComp::Restrict(const Instantiation &x) {
	if (impl.size()==1) {
		vector<int> ids;
		Domain().ConsistentIndexes(ids,x);
		impl[0]->Restrict(ids);
	} else {
		vector<int> condi;
		CondDomain().ConsistentIndexes(condi,x);
		int ci=0;
		int n = impl.size();
		for(int i=0;i<n;i++) {
			if (condi[ci]==i) {
				vector<int> ids;
				Domain().ConsistentIndexes(ids,x);
				impl[i]->Restrict(ids);
				++ci;
			} else {
				vector<int> ids;
				impl[i]->Restrict(ids);
			}
		}
	}
}

void DynComp::Mult(const Dynamics *x) {
	const DynComp *xx = dynamic_cast<const DynComp *>(x);
	assert(xx != NULL); // otherwise, I'm not sure what to do

	// update domains (and save old ones):
	Context oldv(Domain());
	Context oldcv(CondDomain());
	Dynamics::Mult(x);

	vector<DynSimple *> newimpl(CondDomain().Size());
	for(unsigned int i=0;i<newimpl.size();i++)
		newimpl[i] = impl[0]->MakeNew(Domain().Size());
	
	Instantiation v(Domain()+CondDomain(),0);
	Context anycond(oldcv,xx->CondDomain(),Context::UNION);
	Context otherc = CondDomain() - oldcv;
	Context notc = v - oldcv;
	do {
		int myci = oldcv.Index(v);
		vector<vector<int> > myind(oldv.Size());
		do {
			int myi = oldv.Index(v);
			int i = Domain().Index(v);
			myind[myi].push_back(i);
		} while(v.Inc(notc));
		impl[myci]->Expand(myind,Domain().Size());
		
		do {
			int ci = CondDomain().Index(v);
			newimpl[ci]->Mult(impl[myci]);
		} while(v.Inc(otherc));
	} while(v.Inc(oldcv));
	
	otherc = CondDomain() - xx->CondDomain();
	notc = v - xx->CondDomain();
	do {
		int xci = xx->CondDomain().Index(v);
		vector<vector<int> > xind(xx->Domain().Size());
		do {
			int xi = xx->Domain().Index(v);
			int i = Domain().Index(v);
			xind[xi].push_back(i);
		} while(v.Inc(notc));
		DynSimple *ds = xx->impl[xci]->Clone();
		ds->Expand(xind,Domain().Size());
		do {
			int ci = CondDomain().Index(v);
			newimpl[ci]->Mult(ds);
		} while(v.Inc(otherc));
		delete ds;
	} while(v.Inc(xx->CondDomain()));
	
	for(unsigned int i=0;i<impl.size();i++)
		delete impl[i];
	impl = newimpl;
}

RVCondSimple *DynComp::Cond(double t0, double t1,
		const Instantiation &x) const {
	int condi = CondIndex(x);
	vector<int> i;
    	Domain().ConsistentIndexes(i,x);
	return impl[condi]->CondRestrict(t0,t1,i);
}



RVCondSimple *DynComp::Cond(double t, const Instantiation &from,
		const Instantiation &to, bool transition) const {
	int condi = CondIndex(from);
	vector<int> fromi,toi;
	Domain().ConsistentIndexes(fromi,from);
	Domain().ConsistentIndexes(toi,to);
	return impl[condi]->CondRestrict(t,fromi,toi,transition);
}

DynCompSS *DynComp::BlankSS() const {
	int nc = impl.size();
	DynCompSS *ret = new DynCompSS(nc);
	for(int i=0;i<nc;i++)
		ret->ss[i] = impl[i]->BlankSS();
	return ret;
}

void DynComp::AddSS(const Instantiation &x, double t0, double deltat,
		SS *ss, double w) const {
	DynCompSS *dcss = dynamic_cast<DynCompSS *>(ss);
	int condi = CondIndex(x);
	int i = Domain().Index(x);
	impl[condi]->AddSS(i,t0,deltat,dcss->ss[condi],w);
	
}

void DynComp::AddTransSS(const Instantiation &x1, const Instantiation &x2,
		double t, SS *ss, double w) const {
	DynCompSS *dcss = dynamic_cast<DynCompSS *>(ss);
	// assume conditioning variables are the same...
	int condi = CondIndex(x1);
	int i1 = Domain().Index(x1), i2 = Domain().Index(x2);
//	cout << "condi: " << condi << endl;
//	cout << "i1: " << i1 << endl;
//	cout << "i2: " << i2 << endl;
//	x1.SaveOld(cout); cout << endl;
//	x2.SaveOld(cout); cout << endl;
//	cout << endl;
	if(i1!=i2) //add by Yu
		impl[condi]->AddTransSS(i1,i2,t,dcss->ss[condi],w);
}

void DynComp::AddExpSS(const RV *alpha, const RV *beta,
		double t0, double deltat, SS *ss, double w) const {
	DynCompSS *dcss = dynamic_cast<DynCompSS *>(ss);

	RV *aa = alpha->Clone();
	aa->Project(Domain(),CondDomain());
	RV *bb = beta->Clone();
	bb->Project(Domain(),CondDomain());

	Instantiation ci(CondDomain(),0);
	do {
		// assuming that a and b have the same distribution over
		// the conditioning variables...
		RVSimple *a = aa->MakeSimple(ci,false);
		double ww = w*a->Sum();

		if (ww<=0.0) continue;
		a->Normalize();

		RVSimple *b = bb->MakeSimple(ci,true);

		int cind = ci.Index();
		impl[cind]->AddExpSS(a,b,t0,deltat,dcss->ss[cind],ww);
		delete a;
		delete b;
	} while(ci.Inc());
	delete aa;
	delete bb;
}

void DynComp::AddExpTransSS(const RV *x1, const RV *x2,
		const Context &changevar,
		double t, SS *ss, double w) const {

	if (!changevar.IsOverlap(Domain())) return;
	DynCompSS *dcss = dynamic_cast<DynCompSS *>(ss);

	RV *aa = x1->Clone();
	aa->Project(Domain(),CondDomain());
	RV *bb = x2->Clone();
	bb->Project(Domain(),CondDomain());

	RV *acond = x1->Clone();
	acond->Project(CondDomain(),Context());
	RV *bcond = x2->Clone();
	bcond->Project(CondDomain(),Context());
	acond->MultBy(bcond);
	acond->Normalize();

	Instantiation ci(CondDomain(),0);
	do {
		RVSimple *a = aa->MakeSimple(ci,false);
		double ww = w*acond->Prob(ci);
		if (ww<=0.0) continue;

		a->Normalize();
		RVSimple *b = bb->MakeSimple(ci,true);

		int cind = ci.Index();
		impl[cind]->AddExpTransSS(a,b,t,dcss->ss[cind],ww);
		delete a;
		delete b;
	} while(ci.Inc());
	delete acond;
	delete bcond;
	delete aa;
	delete bb;
}


void DynComp::AddExpSS(const Instantiation &condi,
		const RVSimple *alpha, const RVSimple *beta,
		double t0, double deltat, SS *ss, double w) const {
	AddExpSS(CondDomain().Index(condi),alpha,beta,t0,deltat,ss,w);
}

void DynComp::AddExpTransSS(const Instantiation &condi,
		const RVSimple *x1, const RVSimple *x2,
		const Context &changevar,
		double t, SS *ss, double w) const {
	AddExpTransSS(CondDomain().Index(condi),x1,x2,changevar,t,ss,w);
}

void DynComp::AddExpSS(int condi,
		const RVSimple *alpha, const RVSimple *beta,
		double t0, double deltat, SS *ss, double w) const {
	DynCompSS *dcss = dynamic_cast<DynCompSS *>(ss);
	impl[condi]->AddExpSS(alpha,beta,t0,deltat,dcss->ss[condi],w);
}

void DynComp::AddExpTransSS(int condi,
		const RVSimple *x1, const RVSimple *x2,
		const Context &changevar,
		double t, SS *ss, double w) const {
	if (!changevar.IsOverlap(Domain())) return;
	DynCompSS *dcss = dynamic_cast<DynCompSS *>(ss);
	impl[condi]->AddExpTransSS(x1,x2,t,dcss->ss[condi],w);
}

void DynComp::AddSS(const SS *toadd, const Dynamics *dyn,
		SS *ss, double w) const {
	const DynComp *d = dynamic_cast<const DynComp *>(dyn);
	assert(d!=NULL); 
	assert(CondDomain().IsSubset(d->CondDomain()));

	DynCompSS *dcss = dynamic_cast<DynCompSS *>(ss);
	const DynCompSS *toadddcss = dynamic_cast<const DynCompSS *>(toadd);

	Instantiation i(Domain()+CondDomain(),0);
	do {
		int mycondi = CondDomain().Index(i);
		int dcondi = d->CondDomain().Index(i);
		vector<vector<int> > mapping(Domain().Size());
		do {
			d->Domain().ConsistentIndexes(mapping[Domain().Index(i)],i);
		} while(i.Inc(Domain()));
		impl[mycondi]->AddSS(toadddcss->ss[dcondi],d->impl[dcondi],
				mapping,dcss->ss[mycondi],w);
	} while(i.Inc(CondDomain()));
}

void DynComp::SampleNextEvent(const Instantiation &i, double t,
		double &nextt, Instantiation &nexti, Random &rand) const {
	int condi = CondIndex(i);
	int ii = Domain().Index(i);
	int nextii;
	impl[condi]->SampleNextEvent(ii,t,nextt,nextii,rand);
	Instantiation newi(Domain());
	newi.SetIndex(nextii);
	nexti.SetVal(newi,true);
}

void DynComp::Maximize(const SS *ss) {
	const DynCompSS *dcss = dynamic_cast<const DynCompSS *>(ss);
	int nc = impl.size();
	for(int i=0;i<nc;i++)
		impl[i]->Maximize(dcss->ss[i]);
}

void DynComp::Scramble(double a, double b, double alpha, double degree, 
			Random &rand) {
	//by Yu
	int nc = impl.size();
	for(int i=0;i<nc;i++)
		impl[i]->Scramble(a, b, alpha, degree, rand);
}

double DynComp::LLH(const SS *ss) const {
	const DynCompSS *dynss = dynamic_cast<const DynCompSS*>(ss);
	double llh = 0.0;
	int nc = impl.size();
	for(int i=0; i<nc; i++)
		llh += impl[i]->LLH(dynss->ss[i]);
	return llh;
}

double DynComp::GetScore(double numTrans, double amtTime, 
				const SS* ss) const {
	const DynCompSS *dynss = dynamic_cast<const DynCompSS* > (ss);
	double localScore(0);
	int nc = impl.size();
	// Sum over all possible parent instantiations
	for(int i = 0; i < nc; i++)
		localScore += 
			impl[i]->GetScore(nc, numTrans, 
						amtTime, dynss->ss[i]);
	return localScore;
}

} // end of ctbn namespace
