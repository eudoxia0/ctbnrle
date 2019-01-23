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
#ifndef CTBNRLE_MARKOVSIMPLE_H
#define CTBNRLE_MARKOVSIMPLE_H

#include "dynsimple.h"
#include "multisimple.h"
#include "matrix.h"
#include "rk.h"



namespace ctbn {

// This main class in this file (see the end) is Markov Simple
// It describes a Markov process over the state space 0..n-1
// with no "conditioning" variables (that is the dynamics do not
// depend on anything else)
class MarkovSimple;

// The sufficient statistics for MarkovSimple:
class MarkovSimpleSS : public SS {
	SOBJCLASSDECL(MarkovSimpleSS)
  public:
	friend class MarkovSimple;
	friend class MarkovSimpleToggle;
	MarkovSimpleSS(int n=0);
	MarkovSimpleSS(std::istream &is);

	virtual ~MarkovSimpleSS();

	virtual MarkovSimpleSS *Clone() const;
	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
	virtual void Scale(double w);
	virtual void AddSS(const SS* nss, double w=1.0);
	virtual double Element(int s1, int s2) const { return c[s1][s2]; }
  private:
	// diagonal elements hold the amount of time in each state
	// off-diagonal element (i,j) is the number of transitions from i to j
	matrix c;

	SERIAL_START_V(MarkovSimpleSS)
		SERIAL_SUPER(SS)
		SERIAL_VAR(matrix,c)
	SERIAL_END
};

// The next 4 classes are special versions of conditional random variables
// (mapping from a domain to the same domain).  They are used by MarkovSimple
// and are not designed to be learned from data or really interact with
// sufficient statistics.  Therefore all of the suffstats methods are
// blank.  For similar reasons they do not have serialization (load/save)
// currently.

// The next 2 classes represent the conditional distribution of the value of a
// Markov process at one time point given its value at a previous time point
// They store the Q matrix and the time interval instead of the stochastic
// matrix explicitly.
// When "Mult" is called, they generate the marginal distribution at the
// later point given the previous point (and the reverse for RMult) using
// numeric integration (see rk.h)

// This one is directly as described above.
template<typename QTYPE>
class CondTransQ : public RVCondSimple {
	//SOBJCLASSDECL(CondTransQ)
  public:
	CondTransQ(const matrix &QQ=matrix(), double deltat=0.0);
	//CondTransQ(std::istream &is);
	virtual ~CondTransQ();
	virtual CondTransQ<QTYPE> *Clone() const;
	virtual double Prob(int ind, int cind, bool log=false) const;
	virtual RVSimple *Condition(int ind) const;
	virtual void Mult(RVSimple *&x) const;
	virtual void RMult(RVSimple *&x) const;
	virtual RVSimple *Convert(const RVSimple *rvs) const;

	virtual void LoadOld(std::istream &is) {}
	virtual void SaveOld(std::ostream &os) const {}
	virtual void Normalize() {}
	virtual void Restrict(const std::vector<int> &ind, 
				const std::vector<int> &cind) {}
	virtual void Reindex(int n, int nc, 
				const RVCondSimple::RemapT &ind) {}
	virtual void MultBy(const RVCondSimple *x) {}
	virtual SS *BlankSS() const { return NULL; }
	virtual void AddSS(int cind, int ind, SS *ss, double w=1.0) const {}
	virtual void AddExpSS(int cind, SS *ss, double w=1.0) const {}
	virtual void AddExpSS(int cind, const RVSimple *rvs, 
				SS *ss, double w=1.0) const {}
	virtual void AddSS(const SS *toadd, const RVCondSimple* rv,
			const std::vector<std::vector<int> > &mapping,
			int mycondi, int rvccondi, 
			SS *ss, double w=1.0) const {}
	virtual int Sample(int cx, Random &rand = randomizer) const;
	virtual void Maximize(const SS *ss) {}
	virtual void Scramble(double alpha=1.0, double degree=1.0, 
						  Random &rand=randomizer) {}
	virtual double LLH(const SS *ss) const {return 0.0;}
	virtual double GetScore(double numTrans, const SS *ss) const {return 0.0;}
  protected:
	QTYPE Q;
	double t;

/*
	SERIAL_START1_V(CondTransQ,QTYPE)
		SERIAL_SUPER(RVCondSimple)
		SERIAL_VAR(double,t)
		SERIAL_VAR(QTYPE,Q)
	SERIAL_END
*/
};

// This one assumes that the system is constrained to a subset of the
// states (ind) and it only represents the submatrix of this state.
// It should be used with SparseMultiZSimple (as the RVSimple for Mult and
// RMult) to also capture this sparseness.  The methods will error if
// used with other versions of RVSimple
template<typename QTYPE>
class SparseCondTransQ : public RVCondSimple {
	//SOBJCLASSDECL(SparseCondTransQ)
  public:
	SparseCondTransQ(const matrix &QQ=matrix(), double deltat=0.0, 
				const std::vector<int> &ind=std::vector<int>(0));
	virtual ~SparseCondTransQ();
	virtual SparseCondTransQ<QTYPE> *Clone() const;
	virtual double Prob(int ind, int cind, bool log=false) const;
	virtual RVSimple *Condition(int ind) const;
	virtual void Mult(RVSimple *&x) const;
	virtual void RMult(RVSimple *&x) const;
	virtual RVSimple *Convert(const RVSimple *rvs) const;

	virtual void LoadOld(std::istream &is) {}
	virtual void SaveOld(std::ostream &os) const {}
	virtual void Normalize() {}
	virtual void Restrict(const std::vector<int> &ind, 
				const std::vector<int> &cind) {}
	virtual void Reindex(int n, int nc, 
				const RVCondSimple::RemapT &ind) {}
	virtual void MultBy(const RVCondSimple *x) {}
	virtual SS *BlankSS() const { return NULL; }
	virtual void AddSS(int cind, int ind, SS *ss, double w=1.0) const {}
	virtual void AddExpSS(int cind, SS *ss, double w=1.0) const {}
	virtual void AddExpSS(int cind, const RVSimple *rvs, 
				SS *ss, double w=1.0) const {}
	virtual void AddSS(const SS *toadd, const RVCondSimple* rv,
			const std::vector<std::vector<int> > &mapping, 
			int mycondi, int rvccondi, 
			SS *ss, double w=1.0) const {}
	virtual int Sample(int cx, Random &rand = randomizer) const {return 0;}
	virtual void Maximize(const SS *ss) {}
	virtual void Scramble(double alpha=1.0, double degree=1.0, 
				Random &rand=randomizer) {}
	virtual double LLH(const SS *ss) const {return 0.0;}
	virtual double GetScore(double numTrans, const SS *ss) const {return 0.0;}
	QTYPE Q;
	double t;
	std::vector<int> ix;
  protected:

/*
	SERIAL_START1_V(SparseCondTransQ,QTYPE)
		SERIAL_SUPER(RVCondSimple)
		SERIAL_VAR(double,t)
		SERIAL_VAR(std::vector<int>,ix)
		SERIAL_VAR(QTYPE,Q)
	SERIAL_END
*/
};


// The next two classes do the same thing (and again, for an "unconstrained"
// system and then for a "constrained" system).  However, they represent
// the conditional density of a transition at a time instance (therefore
// they map from "just before" the transition to "just after" the transition
// instead of across any non-zero interval of time).

// Here is the unconstrained (non-sparse) version.
template<typename QTYPE>
class CondTransQ2 : public RVCondSimple {
	//SOBJCLASSDECL(CondTransQ2)
  public:
	CondTransQ2(const matrix &QQ=matrix());
	virtual ~CondTransQ2();
	virtual CondTransQ2<QTYPE> *Clone() const;
	virtual double Prob(int ind, int cind, bool log=false) const;
	virtual RVSimple *Condition(int ind) const;
	virtual void Mult(RVSimple *&x) const;
	virtual void RMult(RVSimple *&x) const;
	virtual RVSimple *Convert(const RVSimple *rvs) const;

	virtual void LoadOld(std::istream &is) {}
	virtual void SaveOld(std::ostream &os) const {}
	virtual void Normalize() {}
	virtual void Restrict(const std::vector<int> &ind, 
				const std::vector<int> &cind) {}
	virtual void Reindex(int n, int nc, 
				const RVCondSimple::RemapT &ind) {}
	virtual void MultBy(const RVCondSimple *x) {}
	virtual SS *BlankSS() const { return NULL; }
	virtual void AddSS(int cind, int ind, SS *ss, double w=1.0) const {}
	virtual void AddExpSS(int cind, SS *ss, double w=1.0) const {}
	virtual void AddExpSS(int cind, const RVSimple *rvs, 
				SS *ss, double w=1.0) const {}
	virtual void AddSS(const SS *toadd, const RVCondSimple* rv,
			const std::vector<std::vector<int> > &mapping, 
			int mycondi, int rvccondi, 
			SS *ss, double w=1.0) const {}
	virtual int Sample(int cx, Random &rand = randomizer) const {return 0;}
	virtual void Maximize(const SS *ss) {}
	virtual void Scramble(double alpha=1.0, double degree=1.0, 
				Random &rand=randomizer) {}
	virtual double LLH(const SS *ss) const {return 0.0;}
	virtual double GetScore(double numTrans, const SS *ss) const {return 0.0;}
  protected:
	QTYPE Q;

/*
	SERIAL_START1_V(CondTransQ2,QTYPE)
		SERIAL_SUPER(RVCondSimple)
		SERIAL_VAR(QTYPE,Q)
	SERIAL_END
*/
};

// Here is the constrained version.  It is assumed that the state
// before transition is conditioned to be one of the cind set.
// The state after the transition is conditioned to be one of the ind set.
template<typename QTYPE>
class SparseCondTransQ2 : public RVCondSimple {
	//SOBJCLASSDECL(SparseCondTransQ2)
  public:
	SparseCondTransQ2(const matrix &QQ=matrix(),
				const std::vector<int> &cind=std::vector<int>(0), 
				const std::vector<int> &ind=std::vector<int>(0));
	virtual ~SparseCondTransQ2();
	virtual SparseCondTransQ2<QTYPE> *Clone() const;
	virtual double Prob(int ind, int cind, bool log=false) const;
	virtual RVSimple *Condition(int ind) const;
	virtual void Mult(RVSimple *&x) const;
	virtual void RMult(RVSimple *&x) const;
	virtual RVSimple *Convert(const RVSimple *rvs) const;

	virtual void LoadOld(std::istream &is) {}
	virtual void SaveOld(std::ostream &os) const {}
	virtual void Normalize() {}
	virtual void Restrict(const std::vector<int> &ind, 
				const std::vector<int> &cind) {}
	virtual void Reindex(int n, int nc, 
				const RVCondSimple::RemapT &ind) {}
	virtual void MultBy(const RVCondSimple *x) {}
	virtual SS *BlankSS() const { return NULL; }
	virtual void AddSS(int cind, int ind, SS *ss, double w=1.0) const {}
	virtual void AddExpSS(int cind, SS *ss, double w=1.0) const {}
	virtual void AddExpSS(int cind, const RVSimple *rvs, 
				SS *ss, double w=1.0) const {}
	virtual void AddSS(const SS *toadd, const RVCondSimple* rv,
			const std::vector<std::vector<int> > &mapping, 
			int mycondi, int rvccondi, 
			SS *ss, double w=1.0) const {}
	virtual int Sample(int cx, Random &rand = randomizer) const {return 0;}
	virtual void Maximize(const SS *ss) {}
	virtual void Scramble(double alpha=1.0, double degree=1.0, 
				Random &rand=randomizer) {}
	virtual double LLH(const SS *ss) const {return 0.0;}
	virtual double GetScore(double numTrans, const SS *ss) const {return 0.0;}
	QTYPE Q;
	std::vector<int> cix,ix;
  protected:

/*
	SERIAL_START1_V(SparseCondTransQ2,QTYPE)
		SERIAL_SUPER(RVCondSimple)
		SERIAL_VAR(std::vector<int>,cix)
		SERIAL_VAR(std::vector<int>,ix)
		SERIAL_VAR(QTYPE,Q)
	SERIAL_END
*/
};

// see dynsimple.h for a full description of the methods' uses
class MarkovSimple : public DynSimple {
	SOBJCLASSDECL(MarkovSimple)
  public:
	MarkovSimple(int n=1);
	MarkovSimple(std::istream &is);
	virtual ~MarkovSimple();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	virtual MarkovSimple *Clone() const;
	virtual MarkovSimple *MakeNew(int nstates) const;

	virtual void Normalize();
	virtual void Restrict(const std::vector<int> &ind);
	virtual void Expand(const std::vector<std::vector<int> > &ind, 
				int newn);

	virtual CondTransQ<const matrix &> *Cond(double t0, double t1) const;
	
	virtual SparseCondTransQ<matrix> *CondRestrict(double t0, double t1,
			const std::vector<int> &ind) const;
			

	virtual CondTransQ2<matrix> *Cond(double t) const;
	virtual SparseCondTransQ2<matrix> *CondRestrict(double t,
			const std::vector<int> &fromind,
			const std::vector<int> &toind, 
			bool transition=true) const;

	virtual MarkovSimpleSS *BlankSS() const;

	virtual void Mult(const DynSimple *x);

	virtual void AddSS(int x, double t0, double t1,
			SS *ss, double w=1.0) const;
	virtual void AddTransSS(int x1, int x2, double t,
			SS *ss, double w=1.0) const;
	virtual void AddExpSS(const RVSimple *alpha, const RVSimple *beta,
			double t0, double deltat, SS *ss, double w=1.0) const;
	virtual void AddExpTransSS(const RVSimple *x1, const RVSimple *x2,
			double t, SS *ss, double w=1.0) const;
	virtual void AddSS(const SS *toadd, const DynSimple *dyn,
			const std::vector<std::vector<int> > &mapping,
			SS *ss, double w=1.0) const;

	virtual void SampleNextEvent(int ind, 
					double t,
					double &newt, 
					int &newind, 
					Random &rand=randomizer) const;

	virtual void Maximize(const SS *ss);

	inline matrix &Intensity() { return Q; }
	inline const matrix &Intensity() const { return Q; }

	virtual void Scramble(double a=1.0, double b=1.0, 
				double alpha=1.0, double degree=1.0, 
				Random &rand=randomizer);
	virtual double LLH(const SS *ss) const;

	virtual double GetScore(int parentCard, double numTrans, 
			double amtTime, const SS* ss) const;
  protected:
	matrix Q;
	
	SERIAL_START_V(MarkovSimple)
		SERIAL_SUPER(DynSimple)
		SERIAL_VAR(matrix,Q)
	SERIAL_END
};

// MarkovSimpleToggle represents a variable that "toggles" or "flips"
// as an extention to the MarkovSimple, but has parameter tieing;
// e.g Q = [-2 1 1;
//          1 -2 1;
//          1 1 -2];
// a transition of a toggle variable only indicates a state changes;
// like a flip variable, it does not differentiate between changing from
// which to which;
// the probabilities of transit from one state to the other are all the same; 
// for a fixed size of the toggle variable, only 1 free paratmeter:
// the intensity of the state changes
class MarkovSimpleToggle : public MarkovSimple {
  public:
	MarkovSimpleToggle(int n=1);
	MarkovSimpleToggle(std::istream &is);
	virtual ~MarkovSimpleToggle();

	// calls MarkovSimple:MakeNew();
	virtual MarkovSimpleToggle *MakeNew(int nstates) const;
	// calls MarkovSimple:Clone();
	virtual MarkovSimpleToggle *Clone() const;
	// get the ML estimate from sufficient statistics ss; 
	// the state change intensity is estimated by total number of 
	// transitions between all pair of states divided by total number of 
	// time staying at all states
	virtual void Maximize(const SS *ss);

	virtual void Scramble(double a=1.0, double b=1.0, 
				double alpha=1.0, double degree=1.0, 
				Random &rand=randomizer);

	SERIAL_START_V(MarkovSimpleToggle)
		SERIAL_SUPER(MarkovSimple)
	SERIAL_END
};

//-------------------------

template<typename QTYPE>
CondTransQ<QTYPE>::CondTransQ(const matrix &QQ, double deltat) 
 : RVCondSimple(), Q(QQ), t(deltat) {
}

template<typename QTYPE>
CondTransQ<QTYPE>::~CondTransQ() {
}

template<typename QTYPE>
CondTransQ<QTYPE> *CondTransQ<QTYPE>::Clone() const {
	return new CondTransQ<QTYPE>(*this);
}

template<typename QTYPE>
double CondTransQ<QTYPE>::Prob(int ind, int cind, bool log) const {
	vectr a(Q.getm(),0.0);
	a[cind] = 1.0;
	double zlog = vexpmt(a,Q,t);
	if (log) return ::log(a[cind])+zlog;
	else return a[cind]*exp(zlog);
}

template<typename QTYPE>
RVSimple *CondTransQ<QTYPE>::Condition(int ind) const {
	vectr a(Q.getm(),0.0);
	a[ind] = 1.0;
	double logz = vexpmt(a,Q,t);

	return new MultiZSimple(a,logz);
}

template<typename QTYPE>
int CondTransQ<QTYPE>::Sample(int cx, Random &rand) const {
	vectr a(Q.getm(),0.0);
	a[cx] = 1.0;
	vexpmt(a,Q,t);
//	return rand.SampleMultinomial(a,a.getm(),a.sum());
	return rand.SampleMultinomial(a,a.sum());
}

template<typename QTYPE>
void CondTransQ<QTYPE>::Mult(RVSimple *&x) const {
	vectr a(Q.getm(),0.0);
	double logz;
	x->GetDist(a,logz);
	
//	std::cout << "CondTransQ Exp" << std::std::endl;
//	Q.niceprint(std::cout);
	
	logz += vexpmt(a,Q,t);
	x->SetDist(a,logz);
}

template<typename QTYPE>
void CondTransQ<QTYPE>::RMult(RVSimple *&x) const {
	vectr a(Q.getm(),0.0);
	double logz;
	x->GetDist(a,logz);
	logz += expmtv(a,Q,t);
	x->SetDist(a,logz);
}

template<typename QTYPE>
RVSimple *CondTransQ<QTYPE>::Convert(const RVSimple *rvs) const {
	// any will work...
	return rvs->Clone();
}


//-----------

template<typename QTYPE>
SparseCondTransQ<QTYPE>::SparseCondTransQ(const matrix &QQ,
		double deltat, const std::vector<int> &ind) :
			Q(QQ), t(deltat),ix(ind) {
}

template<typename QTYPE>
SparseCondTransQ<QTYPE>::~SparseCondTransQ() {
}

template<typename QTYPE>
SparseCondTransQ<QTYPE> *SparseCondTransQ<QTYPE>::Clone() const {
	return new SparseCondTransQ<QTYPE>(*this);
}

template<typename QTYPE>
double SparseCondTransQ<QTYPE>::Prob(int ind, int cind, bool log) const {
	std::vector<int>::const_iterator cindpos = 
		lower_bound(ix.begin(),ix.end(),cind);
	if (cindpos==ix.end() || *cindpos!=cind) return log?-INFINITY:0.0;
	
	vectr a(Q.getm(),0.0);
	int tmp = cindpos-ix.begin();
	a[tmp] = 1.0;
	double zlog = vexpmt(a,Q,t);
	cindpos = lower_bound(ix.begin(),ix.end(),ind);
	if (cindpos==ix.end() || *cindpos!=ind) return log?-INFINITY:0.0;
	else {
		tmp = cindpos-ix.begin(); 
		return log ? ::log(a[tmp]) + zlog : a[tmp]*exp(zlog);
	}
}

template<typename QTYPE>
RVSimple *SparseCondTransQ<QTYPE>::Condition(int ind) const {
	std::vector<int>::const_iterator cindpos = 
					lower_bound(ix.begin(),ix.end(),ind);
	if (cindpos==ix.end() || *cindpos!=ind) return NULL;
	// probably should do something else
	vectr a(Q.getm(),0.0);
	int tmp = cindpos-ix.begin();
	a[tmp] = 1.0;
	double zlog = vexpmt(a,Q,t);
	return new SparseMultiZSimple(a,ix,zlog);
}

template<typename QTYPE>
void SparseCondTransQ<QTYPE>::Mult(RVSimple *&x) const {
	// at the moment this assumes that x is a SparseMultiZSimple
	// with indexes that match up
	SparseMultiZSimple *xx = dynamic_cast<SparseMultiZSimple *>(x);
	vectr a;
	double logz;
	xx->GetSparseDist(a,logz);
	
	//std::cout << "SparseCondTransQ Exp" << std::endl;
	//Q.niceprint(std::cout);
	
	logz += vexpmt(a,Q,t);
	//a.niceprint(std::cout);
	xx->SetSparseDist(a,logz);
}

template<typename QTYPE>
void SparseCondTransQ<QTYPE>::RMult(RVSimple *&x) const {
	// at the moment this assumes that x is a SparseMultiZSimple
	// with indexes that match up
	SparseMultiZSimple *xx = dynamic_cast<SparseMultiZSimple *>(x);
	vectr a;
	double logz;
	xx->GetSparseDist(a,logz);
	logz += expmtv(a,Q,t);
	xx->SetSparseDist(a,logz);
}

template<typename QTYPE>
RVSimple *SparseCondTransQ<QTYPE>::Convert(const RVSimple *rvs) const {
	if (dynamic_cast<const SparseMultiZSimple *>(rvs)!=NULL)
		return rvs->Clone();
	vectr a;
	double logz;
	rvs->GetDist(a,logz);
	return new SparseMultiZSimple(vectr(a,ix),ix,logz);
}

template<typename QTYPE>
CondTransQ2<QTYPE>::CondTransQ2(const matrix &QQ) : Q(QQ) {
}

template<typename QTYPE>
CondTransQ2<QTYPE>::~CondTransQ2() {
}

template<typename QTYPE>
CondTransQ2<QTYPE> *CondTransQ2<QTYPE>::Clone() const {
	return new CondTransQ2<QTYPE>(*this);
}

template<typename QTYPE>
double CondTransQ2<QTYPE>::Prob(int ind, int cind, bool log) const {
	return log ? ::log(Q[cind][ind]) : Q[cind][ind];
}

template<typename QTYPE>
RVSimple *CondTransQ2<QTYPE>::Condition(int ind) const {
	vectr a(Q.getm(),0.0);
	a[ind] = 1.0;
	a = a*Q;
	return new MultiZSimple(a);
}

template<typename QTYPE>
void CondTransQ2<QTYPE>::Mult(RVSimple *&x) const {
	// at the moment this assumes that x is a MultiZSimple
	MultiZSimple *xx = dynamic_cast<MultiZSimple *>(x);
	vectr a;
	double logz;
	xx->GetDist(a,logz);
	
//	std::cout << "CondTransQ2 Mult" << std::endl;
//	Q.niceprint(std::cout);
	
	a = a*Q;
//	a.niceprint(std::cout);
	xx->SetDist(a,logz);
}

template<typename QTYPE>
void CondTransQ2<QTYPE>::RMult(RVSimple *&x) const {
	// at the moment this assumes that x is a MultiZSimple
	// with indexes that match up
	MultiZSimple *xx = dynamic_cast<MultiZSimple *>(x);
	vectr a;
	double logz;
	xx->GetDist(a,logz);
	a = Q*a;
	xx->SetDist(a,logz);
}

template<typename QTYPE>
RVSimple *CondTransQ2<QTYPE>::Convert(const RVSimple *rvs) const {
	// any will work...
	return rvs->Clone();
}

//-----------

template<typename QTYPE>
SparseCondTransQ2<QTYPE>::SparseCondTransQ2(const matrix &QQ,
		const std::vector<int> &cind, const std::vector<int> &ind) :
			Q(QQ), cix(cind), ix(ind) {
}

template<typename QTYPE>
SparseCondTransQ2<QTYPE>::~SparseCondTransQ2() {
}

template<typename QTYPE>
SparseCondTransQ2<QTYPE> *SparseCondTransQ2<QTYPE>::Clone() const {
	return new SparseCondTransQ2<QTYPE>(*this);
}

template<typename QTYPE>
double SparseCondTransQ2<QTYPE>::Prob(int ind, int cind, bool log) const {
	std::vector<int>::const_iterator cindpos =
		lower_bound(cix.begin(),cix.end(),cind);
	if (cindpos==cix.end() || *cindpos!=cind) return log?-INFINITY:0.0;
	std::vector<int>::const_iterator indpos = lower_bound(ix.begin(),
								ix.end(),ind);
	if (indpos==ix.end() || *indpos!=ind) return log?-INFINITY:0.0;
	return log ? ::log(Q[cindpos-cix.begin()][indpos-ix.begin()]) :
		Q[cindpos-cix.begin()][indpos-ix.begin()];
}

template<typename QTYPE>
RVSimple *SparseCondTransQ2<QTYPE>::Condition(int ind) const {
	std::vector<int>::const_iterator cindpos =
		lower_bound(cix.begin(),cix.end(),ind);
	if (cindpos==cix.end() || *cindpos!=ind) return NULL;

	vectr a(Q.getm(),0.0);
	int tmp = cindpos-ix.begin();
	a[tmp] = 1.0;
	a = a*Q;
	return new SparseMultiZSimple(a,ix);
}

template<typename QTYPE>
void SparseCondTransQ2<QTYPE>::Mult(RVSimple *&x) const {
	// at the moment this assumes that x is a SparseMultiZSimple
	// with indexes that match up
	SparseMultiZSimple *xx = dynamic_cast<SparseMultiZSimple *>(x);
	vectr a;
	double logz;
	xx->GetSparseDist(a,logz);
	//if (Q.sum() > 1e-10) //comment out by Yu
	
	//std::cout << "SparseCondTransQ2 Mult" << std::endl;
	//Q.niceprint(std::cout);
	
	a = a*Q;
	//a.niceprint(std::cout);
	xx->SetSparseDist(a,logz);
	xx->SetSparseIndexes(ix);
}

template<typename QTYPE>
void SparseCondTransQ2<QTYPE>::RMult(RVSimple *&x) const {
	// at the moment this assumes that x is a SparseMultiZSimple
	// with indexes that match up
	SparseMultiZSimple *xx = dynamic_cast<SparseMultiZSimple *>(x);
	vectr a;
	double logz;
	xx->GetSparseDist(a,logz);
	//if (Q.sum() > 1e-10) //by Yu
		a = Q*a;
	xx->SetSparseDist(a,logz);
	xx->SetSparseIndexes(cix);
}

// this assumes that we will do a "Mult" (and not a "RMult") with this
// one
template<typename QTYPE>
RVSimple *SparseCondTransQ2<QTYPE>::Convert(const RVSimple *rvs) const {
	if (dynamic_cast<const SparseMultiZSimple *>(rvs)!=NULL)
		return rvs->Clone();
	vectr a;
	double logz;
	rvs->GetDist(a,logz);
	return new SparseMultiZSimple(vectr(a,cix),cix,logz);
}

} // end of ctbn namespace

#endif
