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
#include "rvsimple.h"

#if 0
#ifndef CTBNRLE_RVSIMPLE_TCC
#define CTBNRLE_RVSIMPLE_TCC

// implementation of RVCondSimpleComp template
#include "rvsimple.h"

#include <iostream>


SOBJCLASSDEFTEMPL1(RVCondSimpleComp,RVS)


template<class RVS>
RVCondSimpleComp<RVS>::RVCondSimpleComp(int n, int nc) :
	RVCondSimple(), impl(nc,RVS(n)){
}

template<class RVS>
RVCondSimpleComp<RVS>::RVCondSimpleComp(std::istream &is) : RVCondSimple(is) {
	int nc;
	is >> nc;
	impl.resize(nc);
	for(int i=0;i<nc;i++)
		impl[i].LoadOld(is);
}

template<class RVS>
RVCondSimpleComp<RVS>::~RVCondSimpleComp() {}


template<class RVS>
void RVCondSimpleComp<RVS>::LoadOld(std::istream &is) {
	RVCondSimple::LoadOld(is);
	int nc;
	is >> nc;
	impl.resize(nc,RVS(0));
	for(int i=0;i<nc;i++)
		impl[i].LoadOld(is);
}

template<class RVS>
void RVCondSimpleComp<RVS>::SaveOld(std::ostream &os) const {
	RVCondSimple::SaveOld(os);
	int nc = impl.size();
	os << os.fill() << nc;
	for(int i=0;i<nc;i++) {
		os << os.fill();
		impl[i].SaveOld(os);
	}
}

template<class RVS>
RVCondSimpleComp<RVS> *RVCondSimpleComp<RVS>::Clone() const {
	return new RVCondSimpleComp<RVS>(*this);
}

template<class RVS>
double RVCondSimpleComp<RVS>::Prob(int ind, int cind, bool log) const {
	return impl[cind].Prob(ind,log);
}

template<class RVS>
void RVCondSimpleComp<RVS>::Normalize() {
	int nc = impl.size();
	for(int i=0;i<nc;i++)
		impl[i].Normalize();
}

template<class RVS>
RVS *RVCondSimpleComp<RVS>::Condition(int ind) const {
	return impl[ind].Clone();
}

template<class RVS>
void RVCondSimpleComp<RVS>::Restrict(const std::vector<int> &ind, const std::vector<int> &cind) {
	int nc = impl.size();
	std::vector<int> blank;
	unsigned int k=0;
	for(int i=0;i<nc;i++) {
		if (k<cind.size() && i==cind[k]) {
			impl[i].Restrict(ind);
			k++;
		} else {
			impl[i].Restrict(blank);
		}
	}
}

template<class RVS>
void RVCondSimpleComp<RVS>::Reindex(int n, int nc,
	const RVCondSimple::RemapT &ind) {

	std::vector<RVS> newimpl(nc,RVS(n));

	for (RVCondSimple::RemapT::const_iterator i=ind.begin();i!=ind.end();++i) {
		int newci = i->first.first;
		int oldci = i->first.second;
		RVS cpy(impl[oldci]);
		cpy.Reindex(i->second);
		newimpl[newci].Add(&cpy);
	}
	impl = newimpl;
}

template<class RVS>
void RVCondSimpleComp<RVS>::MultBy(const RVCondSimple *x) {
	const RVCondSimpleComp<RVS> *xx
		= dynamic_cast<const RVCondSimpleComp<RVS> *>(x);
	int nc = impl.size();
	for(int i=0;i<nc;i++)
		impl[i].MultBy(&(xx->impl[i]));
}

template<class RVS>
SS *RVCondSimpleComp<RVS>::BlankSS() const {
	int nc = impl.size();
	RVCSCompSS *ret = new RVCSCompSS(nc);
	for(int i=0;i<nc;i++)
		ret->ss[i] = impl[i].BlankSS();
	return ret;
}

template<class RVS>
void RVCondSimpleComp<RVS>::AddSS(int cind, int ind, SS *ss, double w) const {
	RVCSCompSS *rss = dynamic_cast<RVCSCompSS *>(ss);
	impl[cind].AddSS(ind,rss->ss[cind],w);
}

template<class RVS>
void RVCondSimpleComp<RVS>::AddExpSS(int cind, SS *ss, double w) const {
	RVCSCompSS *rss = dynamic_cast<RVCSCompSS *>(ss);
	impl[cind].AddExpSS(rss->ss[cind],w);
}

template<class RVS>
void RVCondSimpleComp<RVS>::AddExpSS(int cind, const RVSimple *rvs, SS *ss, double w) const {
	RVCSCompSS *rss = dynamic_cast<RVCSCompSS *>(ss);
	rvs->AddExpSS(rss->ss[cind],w);
}

template<class RVS>
void RVCondSimpleComp<RVS>::AddSS(const SS *toadd, 
		const RVCondSimple *rv,
		const std::vector<std::vector<int> > &mapping,
		int mycondi, int rvccondi,
		SS *ss, double w) const {
	const RVCondSimpleComp<RVS> *rvcsc =
		dynamic_cast<const RVCondSimpleComp<RVS> *>(rv);
	const RVCSCompSS* toaddrvcss = 
		dynamic_cast<const RVCSCompSS *>(toadd);
	RVCSCompSS* rvcss = 
		dynamic_cast<RVCSCompSS *>(ss);

	impl[mycondi].AddSS(toaddrvcss->ss[rvccondi],
			&(rvcsc->impl[rvccondi]),
			mapping,rvcss->ss[mycondi]);
}

template<class RVS>
int RVCondSimpleComp<RVS>::Sample(int cx, Random &rand) const {
	return impl[cx].Sample(rand);
}

template<class RVS>
void RVCondSimpleComp<RVS>::Maximize(const SS *ss) {
	const RVCSCompSS *sss = dynamic_cast<const RVCSCompSS *>(ss);
	int nc = impl.size();
	for(int i=0;i<nc;i++)
		impl[i].Maximize(sss->ss[i]);
}

template<class RVS>
void RVCondSimpleComp<RVS>::Mult(RVSimple *&x) const {
	int nc = impl.size();
	RVSimple *ret = impl[0].Clone();
	ret->Mult(x->Prob(0));
	for(int i=1;i<nc;i++)
		ret->Add(&(impl[i]),x->Prob(i));
	delete x;
	x = ret;
}

template<class RVS>
void RVCondSimpleComp<RVS>::RMult(RVSimple *&x) const {
	assert(0); // not easily doable (and probably not important)
}

template<class RVS>
RVSimple *RVCondSimpleComp<RVS>::Convert(const RVSimple *rvs) const {
	// can't really multiply anyway... so we'll just call Clone
	return rvs->Clone();
}

template<class RVS>
void RVCondSimpleComp<RVS>::Scramble(double alpha, double degree, Random &rand)
{
        int nc = impl.size();
        for(int i=0; i<nc; i++)
                impl[i].Scramble(alpha, degree, rand);
}

template<class RVS>
double RVCondSimpleComp<RVS>::LLH(const SS *ss) const
{
        double llh = 0.0;
	const RVCSCompSS *sss = dynamic_cast<const RVCSCompSS *>(ss);
        int nc = impl.size();
        for(int i=0; i<nc; i++)
                llh += impl[i].LLH(sss->ss[i]);
        return llh;
}

template<class RVS>
double RVCondSimpleComp<RVS>::GetScore(double numTrans, const SS* ss) 
					const {
	double localScore(0);
	const RVCSCompSS *rvss = dynamic_cast<const RVCSCompSS *> (ss);
	int nc = impl.size();
		
	for(int i = 0; i < nc; i++)		
		localScore += impl[i].GetScore(numTrans/nc, rvss->ss[i]);
	return localScore;
}


#endif
#endif
