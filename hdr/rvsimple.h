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
#ifndef CTBNRLE_RVSIMPLE_H
#define CTBNRLE_RVSIMPLE_H

#include "serial/serial.h"

#include "matrix.h"
#include "random.h"
#include "rv.h"
#include "ss.h"

#include <iosfwd>
#include <utility>
#include <vector>



namespace ctbn {

// RVSimple is the ABC (Abstract Base Class) for measures.
// It can be thought of as a random variable that is not conditional
// and whose domain is implied (or just indexed by integers from 0 to n-1)

class RVSimple : public StreamObj {
public:
	RVSimple();
	RVSimple(std::istream &is);
	virtual ~RVSimple();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	// virtual copy constructor
	virtual RVSimple *Clone() const = 0;

	virtual std::ostream & Print(std::ostream & out, RV * the_initial_distribution_of_the_markov_process) const;
    virtual std::ostream & PrintJoint(std::ostream & out, RV * the_initial_distribution_of_the_markov_process) const;
    virtual std::ostream & PrintMarginal(std::ostream & out, RV * the_initial_distribution_of_the_markov_process) const;

	// Makes sum to 1 -- returns sum prior to normalization
	virtual double Normalize() = 0;

	// returns the probability of an individual event (index)
	virtual double Prob(int ind, bool log=false) const = 0;

	// Returns the sum (or log thereof)
	virtual double Sum(bool log=false) const = 0;

	// Sets to zero all indexes that are not listed in ind
	// Restrict followed by Normalize is the same as
	// conditioning the RV to the set of ind
	virtual void Restrict(const std::vector<int> &ind) = 0;

	// Marginalizes, reindexes, or expands depending on the argument
	// ind is a list, for each new index, of the list of old indexes
	// which should be summed
	virtual void Reindex(const std::vector<std::vector<int> > &ind) = 0;

	// Multiplies by another measure (point-wise)
	virtual void MultBy(const RVSimple *x) = 0;

	// Adds in another measure
	virtual void Add(const RVSimple *x) = 0;
	virtual void Add(const RVSimple *x, double w) = 0;

	// Multiplies measure by constant
	virtual void Mult(double x) = 0;

	// return suitable sufficient statistics object
	virtual SS *BlankSS() const = 0;

	// add independent draw of x to ss
	virtual void AddSS(int x, SS *ss, double w=1.0) const = 0;
	// add average of independent draws of x to ss
	virtual void AddExpSS(SS *ss, double w=1.0) const = 0;

	virtual void AddSS(const SS *toadd, const RVSimple* rvs,
			const std::vector<std::vector<int> > &mapping, 
			SS *ss, double w=1.0) const = 0;

	// returns a sample
	virtual int Sample(Random &rand = randomizer) const = 0;

	// sets parameters to ML estimate
	virtual void Maximize(const SS *ss) = 0;

	// set measure to be 1 everywhere it has support
	virtual void MakeUniform() = 0;

	virtual void GetDist(vectr &d, double &logfactor) const = 0;
	virtual void SetDist(const vectr &d, double logfactor=0.0) = 0;

	// Draws parameters from Dirichlet distribution
	// with all alphas the same (equal to alpha)
	// If degree equals 1, all parameters are replaced
	// If degree is between 0 and 1, the old parameters are
	// mixed with the new one (linear combination, degree controls
	// the mixing).  
	virtual void Scramble(double alpha=1.0, double degree=1.0, 
			Random &rand=randomizer) = 0; 
	// Returns the log-likelihood of a set of sufficient statistics
	virtual double LLH(const SS *ss) const = 0;

	virtual double GetScore(double numTrans, const SS* ss) const;

	SERIAL_START_V(RVSimple)
	SERIAL_END
};

// This is the same, but with conditioning
// That is, it represents a conditional measure (indexed by 0..n-1)
// conditioned on an integer from 0..m-1
class RVCondSimple : public StreamObj {
public:
	RVCondSimple();
	RVCondSimple(std::istream &is);
	virtual ~RVCondSimple();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	// virtual copy constructor
	virtual RVCondSimple *Clone() const = 0;

    virtual std::ostream & PrintCond(std::ostream & out) {
        throw "Not yet implemented.";
    }

	// returns the probability of an individual conditional event
	virtual double Prob(int ind, int cind, bool log=false) const = 0;

	// Makes sum to 1 for each conditioning value
	virtual void Normalize() = 0;

	// returns an object representing the conditioning
	// on that value
	virtual RVSimple *Condition(int ind) const = 0;

	// Sets to zero all indexes that are not listed in ind and cind
	// That is, for conditions mentioned in cind, the measure is
	// restricted to the indexes in ind.  For conditions not in cind,
	// the measure is made to be zero everywhere.
	virtual void Restrict(const std::vector<int> &ind, 
			const std::vector<int> &cind) = 0;


	// Marginalizes, reindexes, or expands depending on the
	// argument
	// n = # new values  nc = # new conditioning values
	// ind[<i,j>] is a mapping that needs to be applied to old 
	//   conditioning index j and then added to i
	// (see Reindex for RVSimple above to see how a mapping is
	//  applied)
	// (Really, this is probably too confusing to be "simple.")
	typedef std::map<std::pair<int,int>, std::vector<std::vector<int> > > RemapT;
	virtual void Reindex(int n, int nc, const RemapT &ind) = 0;

	// Multiplies by another measure (point-wise)
	virtual void MultBy(const RVCondSimple *x) = 0;

	// return suitable sufficient statistics object
	virtual SS *BlankSS() const = 0;

	// add independent draw of (cind,ind) (conditioning value, value) to ss
	virtual void AddSS(int cind, int ind, SS *ss, double w=1.0) const = 0;
	// add average of independent draws to ss
	virtual void AddExpSS(int cind, SS *ss, double w=1.0) const = 0;
	virtual void AddExpSS(int cind, const RVSimple *rvs, 
			SS *ss, double w) const = 0;

	virtual void AddSS(const SS *toadd, const RVCondSimple* rv, 
			const std::vector<std::vector<int> > &mapping,
			int mycondi, int rvccondi,
			SS *ss, double w=1.0) const = 0;

	// returns a sample
	virtual int Sample(int cx, Random &rand = randomizer) const = 0;

	// sets parameters to ML estimate
	virtual void Maximize(const SS *ss) = 0;

	// assumes that x is over the conditional range
	// multiplies x by *this and then sums out over the conditional
	// range to produce a distribution over the domain
	// x is changed to the result
	// Think of *this as p(x|y) and the input as p(y)
	// This function changes the input to be the resulting
	// margin p(x)
	virtual void Mult(RVSimple *&x) const = 0;
	// this is the reverse (not a generally used operation except
	// in the forward-backward algorithm)
	// it assumes x is over the range of the non-conditioning variable
	// and returns a measure over the conditioning variables
	virtual void RMult(RVSimple *&x) const = 0;
	// Returns a new object that can be multiplied into this
	// object (as per Mult and RMult)
	virtual RVSimple *Convert(const RVSimple *rvs) const = 0;

	// see above class for description
	virtual void Scramble(double alpha=1.0, double degree=1.0, 
			Random &rand=randomizer) = 0; //by Yu
	virtual double LLH(const SS *ss) const = 0;

	virtual double GetScore(double numTrans, const SS *ss) 
								const = 0;

	SERIAL_START_V(RVCondSimple)
	SERIAL_END
};

// This is a template.  It takes in a subclass of RVSimple and
// creates a RVCondSimple from it (by making one copy for each 
// conditioning value)
template<class RVS>
class RVCondSimpleComp : public RVCondSimple {
	SOBJCLASSDECL(RVCondSimpleComp)
public:
	RVCondSimpleComp(int n=0, int nc=0);
	RVCondSimpleComp(std::istream &is);
	virtual ~RVCondSimpleComp();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
	virtual RVCondSimpleComp<RVS> *Clone() const;
	virtual double Prob(int ind, int cind, bool log=false) const;
	virtual void Normalize();
	virtual RVS *Condition(int ind) const;
	virtual void Restrict(const std::vector<int> &ind, const std::vector<int> &cind);
	virtual void Reindex(int n, int nc, const RVCondSimple::RemapT &ind);
	virtual void MultBy(const RVCondSimple *x);
	virtual SS *BlankSS() const;
	virtual void AddSS(int cind, int ind, SS *ss, double w=1.0) const;
	virtual void AddExpSS(int cind, SS *ss, double w=1.0) const;
	virtual void AddExpSS(int cind, const RVSimple *rvs, SS *ss, double w) const;
	virtual void AddSS(const SS *toadd, const RVCondSimple* rv,
				const std::vector<std::vector<int> > &mapping,
				int mycondi, int rvccondi,
				SS* ss, double w=1.0) const;
	virtual int Sample(int cx, Random &rand = randomizer) const;
	virtual void Maximize(const SS *ss);
	virtual void Mult(RVSimple *&x) const;
	virtual void RMult(RVSimple *&x) const;
	virtual RVSimple *Convert(const RVSimple *rvs) const;

	RVS &operator[](int i) { return impl[i]; }
	const RVS &operator[](int i) const { return impl[i]; }

	virtual void Scramble(double alpha=1.0, double degree=1.0, 
			Random &rand=randomizer); //by Yu
	virtual double LLH(const SS *ss) const;
	virtual double GetScore(double numTrans, const SS* ss) const;
private:
	std::vector<RVS> impl;

	SERIAL_START1_V(RVCondSimpleComp,RVS)
		SERIAL_SUPER(RVCondSimple)
		SERIAL_VAR(std::vector<RVS>,impl)
	SERIAL_END
};

// This is the corresponding sufficient statistics for RVCondSimpleComp
// (for any template instantiation)
class RVCSCompSS : public SS {
	SOBJCLASSDECL(RVCSCompSS)
public:
	RVCSCompSS(int nc=0);
	RVCSCompSS(const RVCSCompSS &x);
	RVCSCompSS(std::istream &is);
	virtual ~RVCSCompSS();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
	virtual RVCSCompSS *Clone() const;

	RVCSCompSS &operator=(const RVCSCompSS &x);

	virtual void Scale(double w);
	virtual void AddSS(const SS* nss, double w=1.0);
	// can't make the general set of RVCondSimpleComp<*> friends
	// for some reason in g++
//private:
	std::vector<SS *> ss;

	SERIAL_START_V(RVCSCompSS)
		SERIAL_SUPER(SS)
		SERIAL_VAR(std::vector<SS*>,ss)
	SERIAL_END
};

//#include "rvsimple.tcc"
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

} // end of ctbn namespace

#endif
