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
#include "rvcomp.h"
#include "multisimple.h"



namespace ctbn {

using namespace std;

SOBJCLASSDEF(RVComp)

RVComp::RVComp(const Context &var, const Context &cvar,
		RVCondSimple *touse) : RV(var,cvar) {
	impl = touse;
}

RVComp::RVComp(istream &is) : RV(is) {
	impl = dynamic_cast<RVCondSimple *>(StreamObj::LoadOldPtr(is));
}

RVComp::RVComp(const RVComp &rvc) : RV(rvc) {
	impl = rvc.impl->Clone();
}

RVComp &RVComp::operator=(const RVComp &rvc) {
	if (&rvc != this) {
		if (impl!=NULL) delete impl;
		impl = rvc.impl->Clone();
	}
	return *this;
}

RVComp::~RVComp() {
	if (impl!=NULL) delete impl;
}

void RVComp::LoadOld(std::istream &is) {
	RV::LoadOld(is);
	if (impl!=NULL) delete impl;
	impl = dynamic_cast<RVCondSimple *>(StreamObj::LoadOldPtr(is));
}

void RVComp::serial_preload() {
	if (impl!=NULL) delete impl;
	impl = NULL;
}

void RVComp::SaveOld(std::ostream &os) const {
	RV::SaveOld(os);
	os << os.fill();
	impl->SaveOldPtr(os);
}

RVComp *RVComp::Clone() const {
	return new RVComp(*this);
}

void RVComp::Normalize() {
	impl->Normalize();
}

double RVComp::Prob(const Instantiation &x, bool log) const {
	int condi = CondDomain().Index(x);
	if (condi==-1) return log ? -INFINITY : 0.0;
	RVSimple *condrv = impl->Condition(condi);
	vector<int> ids;
	Domain().ConsistentIndexes(ids,x);
	condrv->Restrict(ids);
	double ret = condrv->Sum(log);
	delete condrv;
	return ret;
}

void RVComp::Restrict(const Instantiation &x) {
	vector<int> i,ci;
	Domain().ConsistentIndexes(i,x);
	CondDomain().ConsistentIndexes(ci,x);
	impl->Restrict(i,ci);
}

void RVComp::Project(const Context &c, const Context &cc) {
	Context all = c+cc+Domain()+CondDomain();
	Context condc = cc+CondDomain();
	Context nonc = all-condc;

	RVCondSimple::RemapT ind;

	forallassign(Instantiation &vcond,condc) {
		int oldci = CondDomain().Index(vcond);
		int newci = cc.Index(vcond);
		vector<vector<int> > remap(c.Size());
		forallassign_subset(Instantiation &v,nonc,vcond) {
			int oldi = Domain().Index(v);
			int newi = c.Index(v);
			remap[newi].push_back(oldi);
		}
		ind[make_pair(newci,oldci)] = remap;
	}
	impl->Reindex(c.Size(),cc.Size(),ind);
	RV::Project(c,cc);
}

void RVComp::MultBy(const RV *x) {
	const RVComp *xx = dynamic_cast<const RVComp *>(x);
	Context newc(Domain(),xx->Domain());
	Context newcc(CondDomain(),xx->CondDomain());
	newcc = newcc - newc;
	if (newc.Size()!=Domain().Size() || newcc.Size()!=CondDomain().Size())
		Project(newc,newcc);
	if (newc.Size()!=xx->Domain().Size()
			|| newcc.Size()!=xx->CondDomain().Size()) {
		RVComp *cpy = xx->Clone();
		cpy->Project(newc,newcc);
		impl->MultBy(cpy->impl);
		delete cpy;
	} else {
		impl->MultBy(xx->impl);
	}
}

SS *RVComp::BlankSS() const {
	return impl->BlankSS();
}

void RVComp::AddSS(const Instantiation &x, SS *ss, double w) const {
	int condi = CondDomain().Index(x);
	int i = Domain().Index(x);
	impl->AddSS(condi,i,ss,w);
}

void RVComp::AddExpSS(const Instantiation &x, SS *ss, double w) const {
	int condi = CondDomain().Index(x);
	vector<int> ids,cids;
	Domain().ConsistentIndexes(ids,x);
	cids.push_back(condi);
	if (ids.size()==static_cast<unsigned int>(Domain().Size())) {
		impl->AddExpSS(condi,ss,w);
	} else {
		RVCondSimple *condrv = impl->Clone();
		condrv->Restrict(ids,cids);
		condrv->Normalize();
		condrv->AddExpSS(condi,ss,w);
		delete condrv;
	}
}

void RVComp::AddExpSS(const Instantiation &x, const RVSimple *rvs, 
		SS *ss, double w) const {
//	cerr << "here" << endl;
//	cerr << rvs << endl;
	int condi = CondDomain().Index(x);
//	cerr << condi << endl;
//	cerr << "impl of this ";
//	cerr << impl << endl;
	impl->AddExpSS(condi,rvs,ss,w);
}

void RVComp::AddSS(const SS *toadd, const RV* rv, SS *ss, double w) const {
	const RVComp *rvc = 
		dynamic_cast<const RVComp *>(rv);
	assert(rvc!=NULL); 
	assert(CondDomain().IsSubset(rvc->CondDomain()));

	Instantiation i(Domain()+CondDomain(),0);
	forallassign(Instantiation &condi,CondDomain()) {
		vector<vector<int> > mapping(Domain().Size());
		forallassign_subset(Instantiation &i,Domain(),condi)
			rvc->Domain().ConsistentIndexes(mapping[Domain().Index(i)],i);
		impl->AddSS(toadd,rvc->impl,mapping,condi.Index(),
				rvc->CondDomain().Index(condi),ss);
	}
}

void RVComp::Sample(Instantiation &x, Random &rand) const {
	int condi = CondDomain().Index(x);
	vector<int> ids,cids;
	Domain().ConsistentIndexes(ids,x);
	cids.push_back(condi);
	int i;
	if (ids.size()==static_cast<unsigned int>(Domain().Size())) {
        // sample all...
		i = impl->Sample(condi,rand);
	} else {
		RVCondSimple *condrv = impl->Clone();
		condrv->Restrict(ids,cids);
		condrv->Normalize();
		i = condrv->Sample(condi,rand);
		delete condrv;
	}
	Instantiation newi(Domain());
	newi.SetIndex(i);
	x.SetVal(newi,true);
}

void RVComp::Maximize(const SS *ss) {
	impl->Maximize(ss);
}

RVSimple * RVComp::MakeSimple(const Instantiation &x, bool normalize) const {
	RVSimple *ret = impl->Condition(CondDomain().Index(x));
	vector<int> ids;
	Domain().ConsistentIndexes(ids,x);
	ret->Restrict(ids);
	if (normalize) ret->Normalize();
	return ret;
}

void RVComp::Scramble(double alpha, double degree, Random &rand) {
	impl->Scramble(alpha, degree, rand);
}

double RVComp::LLH(const SS *ss) const {
	return impl->LLH(ss);
}

double RVComp::GetScore(double numTrans, const SS *ss) const {
	return impl->GetScore(numTrans,ss);
}

} // end of ctbn namespace
