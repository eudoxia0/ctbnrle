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
#include "rvcomp.h"

namespace ctbn {

using namespace std;

SOBJCLASSDEFTEMPL1(DynComp,DS)

template<class DS>
DynComp<DS>::DynComp(const Context &var, const Context &cvar) :
	Dynamics(var,cvar), impl(cvar.Size(),DS(var.Size())) {
}

template<class DS>
DynComp<DS>::DynComp(std::istream &is) :
	Dynamics(is), impl(CondDomain().Size(),DS(Domain().Size())) {
	int n = impl.size();
	for(int i=0;i<n;i++)
		impl[i].LoadOld(is);
}

template<class DS>
DynComp<DS>::~DynComp() {
}

template<class DS>
void DynComp<DS>::LoadOld(std::istream &is) {
	Dynamics::LoadOld(is);
	int n = CondDomain().Size();
	impl.resize(n);
	for(int i=0;i<n;i++)
		impl[i].LoadOld(is);
}

template<class DS>
void DynComp<DS>::SaveOld(std::ostream &os) const {
	Dynamics::SaveOld(os);
	int n = impl.size();
	for(int i=0;i<n;i++) {
		os << os.fill();
		impl[i].SaveOld(os);
	}
}

template<class DS>
DynComp<DS> *DynComp<DS>::Clone() const {
	return new DynComp<DS>(*this);
}

template<class DS>
void DynComp<DS>::Normalize() {
	int n = impl.size();
	for(int i=0;i<n;i++)
		impl[i].Normalize();
}

template<class DS>
void DynComp<DS>::Restrict(const Instantiation &x) {
	if (impl.size()==1) {
		vector<int> ids;
		Domain().ConsistentIndexes(ids,x);
		impl[0].Restrict(ids);
	} else {
		vector<int> condi;
		CondDomain().ConsistentIndexes(condi,x);
		int ci=0;
		int n = impl.size();
		for(int i=0;i<n;i++) {
			if (condi[ci]==i) {
				vector<int> ids;
				Domain().ConsistentIndexes(ids,x);
				impl[i].Restrict(ids);
				++ci;
			} else {
				vector<int> ids;
				impl[i].Restrict(ids);
			}
		}
	}
}

template<class DS>
void DynComp<DS>::Mult(const Dynamics *x) {
	const DynComp<DS> *xx = dynamic_cast<const DynComp<DS> *>(x);
	assert(xx!=NULL); // otherwise, I'm not sure what to do

	// update domains (and save old ones):
	Context oldv(Domain());
	Context oldcv(CondDomain());
	Dynamics::Mult(x);

	vector<DS> newimpl(CondDomain().Size());

	Instantiation v(Domain()+CondDomain(),0);
	Context notc(v,oldcv,Context::DIFFERENCE);
	do {
		int ci = CondDomain().Index(v);
		int myci = oldcv.Index(v);
		vector<vector<int> > myind(oldv.Size());
		do {
			int myi = oldv.Index(v);
			int i = Domain().Index(v);
			myind[myi].push_back(i);
		} while(v.Inc(notc));
		impl[myci].Expand(myind,Domain().Size());
		newimpl[ci].Mult(&impl[myci]);
	} while(v.Inc(oldcv));

	v.SetIndex(0);
	notc = v-xx->CondDomain();
	do {
		int ci = CondDomain().Index(v);
		int xci = xx->CondDomain().Index(v);
		vector<vector<int> > xind(xx->Domain().Size());
		do {
			int xi = xx->Domain().Index(v);
			int i = Domain().Index(v);
			xind[xi].push_back(i);
		} while(v.Inc(notc));
		DS ds(xx->impl[xci]);
		ds.Expand(xind,Domain().Size());
		newimpl[ci].Mult(&ds);
	} while(v.Inc(xx->CondDomain()));

	impl = newimpl;
}

template<class DS>
RV *DynComp<DS>::Cond(double t0, double t1) const {
	// difficult due to same variable names on "both sides"
	// need a method for remapping variable names...
	assert(0);
}

template<class DS>
RV *DynComp<DS>::CondRestrict(double t0, double t1, const Context &c) const {
	// difficult due to same variable names on "both sides"
	// need a method for remapping variable names...
	assert(0);
}

template<class DS>
DynCompSS *DynComp<DS>::BlankSS() const {
	int nc = impl.size();
	DynCompSS *ret = new DynCompSS(nc);
	for(int i=0;i<nc;i++)
		ret->ss[i] = impl[i].BlankSS();
	return ret;
}

template<class DS>
void DynComp<DS>::AddSS(const Instantiation &x, double t0, double t1,
		SS *ss, double w) const {
	DynCompSS *dcss = dynamic_cast<DynCompSS *>(ss);
	int condi = CondIndex(x);
	int i = Domain().Index(x);
	impl[condi].AddSS(i,t0,t1,dcss->ss[condi],w);
}

template<class DS>
void DynComp<DS>::AddTransSS(const Instantiation &x1, const Instantiation &x2,
		double t, SS *ss, double w) const {
	DynCompSS *dcss = dynamic_cast<DynCompSS *>(ss);
	// assume conditioning variables are the same...
	int condi = CondIndex(x1);
	int i1 = Domain().Index(x1), i2 = Domain().Index(x2);
	impl[condi].AddTransSS(i1,i2,t,dcss->ss[condi],w);
}

template<class DS>
void DynComp<DS>::AddExpSS(const RV *alpha, const RV *beta,
		double t0, double t1, SS *ss, double w) const {
	// need to know how to convert alpha into a RVSimple object
	// probably best way is to make DynComp have two template
	// arguments (the second being the class C for RVComp<C>)
	// to do later...
	assert(0);
}

template<class DS>
void DynComp<DS>::AddExpTransSS(const RV *x1, const RV *x2,
		double t0, SS *ss, double w) const {
	// see above
	assert(0);
}

template<class DS>
void DynComp<DS>::SampleNextEvent(const Instantiation &i, double t,
		double &nextt, Instantiation &nexti, Random rand) const {
	int condi = CondIndex(i);
	int ii = Domain().Index(i);
	int nextii;
	impl[condi].SampleNextEvent(ii,t,nextt,nextii,rand);
	nexti.SetIndex(nextii);
}

template<class DS>
void DynComp<DS>::Maximize(const SS *ss) {
	const DynCompSS *dcss = dynamic_cast<const DynCompSS *>(ss);
	int nc = impl.size();
	for(int i=0;i<nc;i++)
		impl[i].Maximize(dcss->ss[i]);
}

}
