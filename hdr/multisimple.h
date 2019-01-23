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
#ifndef CTBNRLE_MULTISIMPLE_H
#define CTBNRLE_MULTISIMPLE_H

#include "rvsimple.h"
#include "matrix.h"



namespace ctbn {

// Four different versions of a RVSimple (see rvsimple.h)
// Each one stores the multinomial distribution differently
// However, all assume a general multinomial form (no parameter tieing)

class MultiSimple;
class MultiLogSimple;
class MultiZSimple;
class SparseMultiZSimple;

// This is the Sufficient Statistics for any of the RVSimple
// subclasses in this file
class MultiSimpleSS : public SS {
	SOBJCLASSDECL(MultiSimpleSS)
public:
	friend class MultiSimple;
	friend class MultiLogSimple;
	friend class MultiZSimple;
	friend class SparseMultiZSimple;
	friend class RVSimple;
public:
	MultiSimpleSS(int n=1);
	MultiSimpleSS(std::istream &is);
	MultiSimpleSS(const MultiSimpleSS &msss);
	virtual ~MultiSimpleSS();

	virtual MultiSimpleSS *Clone() const;
	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
        virtual void Scale(double w);
        virtual void AddSS(const SS* nss, double w=1.0);
private:
	vectr c; // vector of counts

	SERIAL_START_V(MultiSimpleSS)
		SERIAL_SUPER(SS)
		SERIAL_VAR(vectr,c)
	SERIAL_END
};

// This class stores the probabilties directly (as theta, below)
class MultiSimple : public RVSimple {
	SOBJCLASSDECL(MultiSimple)
public:
	MultiSimple(int n=1);
	MultiSimple(std::istream &is);
	virtual ~MultiSimple();

	virtual MultiSimple *Clone() const;
	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	virtual double Normalize();
	virtual double Prob(int ind, bool log=false) const;
	virtual double Sum(bool log=false) const;
	virtual void Restrict(const std::vector<int> &ind);
	virtual void Reindex(const std::vector<std::vector<int> > &ind);
	virtual void MultBy(const RVSimple *x);
	virtual void Add(const RVSimple *x);
	virtual void Add(const RVSimple *x, double w);
	virtual void Mult(double x);
	virtual MultiSimpleSS *BlankSS() const;
	virtual void AddSS(int x, SS *ss, double w=1.0) const;
	virtual void AddExpSS(SS *ss, double w=1.0) const;
	virtual void AddSS(const SS *toadd, const RVSimple *rvs,
				const std::vector<std::vector<int> > &mapping,
				SS *ss, double w=1.0) const;
	virtual int Sample(Random &rand = randomizer) const;
	virtual void Maximize(const SS *ss);
	virtual void GetDist(vectr &d, double &logfactor) const;
	virtual void SetDist(const vectr &d, double logfactor=0.0);
	virtual void MakeUniform();

	virtual void Scramble(double alpha=1.0, double degree=1.0, 
			Random &rand = randomizer);
	virtual double LLH(const SS *ss) const;
protected:
	vectr theta; // parameters of the multinomial distribution

	SERIAL_START_V(MultiSimple)
		SERIAL_SUPER(RVSimple)
		SERIAL_VAR(vectr,theta)
	SERIAL_END
};

// same as above, but log of parameters is stored
class MultiLogSimple : public RVSimple {
	SOBJCLASSDECL(MultiLogSimple)
public:
	MultiLogSimple(int n=1);
	MultiLogSimple(std::istream &is);
	virtual ~MultiLogSimple();

	virtual MultiLogSimple *Clone() const;
	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	virtual double Normalize();
	virtual double Prob(int ind, bool log=false) const;
	virtual double Sum(bool log=false) const;
	virtual void Restrict(const std::vector<int> &ind);
	virtual void Reindex(const std::vector<std::vector<int> > &ind);
	virtual void MultBy(const RVSimple *x);
	virtual void Add(const RVSimple *x);
	virtual void Add(const RVSimple *x, double w);
	virtual void Mult(double x);
	virtual MultiSimpleSS *BlankSS() const;
	virtual void AddSS(int x, SS *ss, double w=1.0) const;
	virtual void AddExpSS(SS *ss, double w=1.0) const;
	virtual void AddSS(const SS *toadd, const RVSimple *rvs,
				const std::vector<std::vector<int> > &mapping,
				SS *ss, double w=1.0) const;
	virtual int Sample(Random &rand = randomizer) const;
	virtual void Maximize(const SS *ss);
	virtual void GetDist(vectr &d, double &logfactor) const;
	virtual void SetDist(const vectr &d, double logfactor=0.0);
	virtual void MakeUniform();

	virtual void Scramble(double alpha=1.0, double degree=1.0, 
			Random &rand = randomizer);
	virtual double LLH(const SS *ss) const;
protected:
	vectr logtheta; // parameters of the multinomial distribution

	SERIAL_START_V(MultiLogSimple)
		SERIAL_SUPER(RVSimple)
		SERIAL_VAR(vectr,logtheta)
	SERIAL_END
};

// Same as MultiSimple, but normalization constant is kept separately
// (in log-space)
class MultiZSimple : public RVSimple {
	SOBJCLASSDECL(MultiZSimple)
public:
	MultiZSimple(int n=1);
	MultiZSimple(const vectr &dist, double logfactor=0.0);
	MultiZSimple(std::istream &is);
	virtual ~MultiZSimple();

	virtual MultiZSimple *Clone() const;
	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	virtual double Normalize();
	virtual double Prob(int ind, bool log=false) const;
	virtual double Sum(bool log=false) const;
	virtual void Restrict(const std::vector<int> &ind);
	virtual void Reindex(const std::vector<std::vector<int> > &ind);
	virtual void MultBy(const RVSimple *x);
	virtual void Add(const RVSimple *x);
	virtual void Add(const RVSimple *x, double w);
	virtual void Mult(double x);
	virtual MultiSimpleSS *BlankSS() const;
	virtual void AddSS(int x, SS *ss, double w=1.0) const;
	virtual void AddExpSS(SS *ss, double w=1.0) const;
	virtual void AddSS(const SS *toadd, const RVSimple *rvs,
				const std::vector<std::vector<int> > &mapping,
				SS *ss, double w=1.0) const;
	virtual int Sample(Random &rand = randomizer) const;
	virtual void Maximize(const SS *ss);
	virtual void GetDist(vectr &d, double &logfactor) const;
	virtual void SetDist(const vectr &d, double logfactor=0.0);
	virtual void MakeUniform();

	virtual void Scramble(double alpha=1.0, double degree=1.0, 
			Random &rand = randomizer);
	virtual double LLH(const SS *ss) const;
protected:
	void SetZ();

	vectr theta; // parameters of the multinomial distribution
	double logz; // the normalizing factor: true theta = exp(logz)*theta

	SERIAL_START_V(MultiZSimple)
		SERIAL_SUPER(RVSimple)
		SERIAL_VAR(double,logz)
		SERIAL_VAR(vectr,theta)
	SERIAL_END
};

// Same as MultiZSimple, but only a sparse number of non-zero elements
// are kept
class SparseMultiZSimple : public MultiZSimple {
	SOBJCLASSDECL(SparseMultiZSimple)
public:
	SparseMultiZSimple(int n=1);
	// indexes need to be in sorted order
	SparseMultiZSimple(const vectr &dist, const std::vector<int> &indexes,
			double logfactor=0.0);
	SparseMultiZSimple(std::istream &is);
	virtual ~SparseMultiZSimple();

	virtual SparseMultiZSimple *Clone() const;
	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	// These work as before:
	using MultiZSimple::Normalize;
	using MultiZSimple::Sum;
	using MultiZSimple::Mult;
	using MultiZSimple::MultBy; // this one is a bit dangerous 
			// as it assumes the represented indices are the same
	using MultiZSimple::Add; // as is this one
	using MultiZSimple::BlankSS;
	//using MultiZSimple::AddSS;
	using MultiZSimple::MakeUniform;
	using MultiZSimple::Scramble;  // note that this means that
			// scramble will not "unzero" zero elements (unlike Maximize)

	// These are modified:
	virtual double Prob(int ind, bool log=false) const;
	virtual void Restrict(const std::vector<int> &ind);
	virtual void Reindex(const std::vector<std::vector<int> > &ind);
	virtual void AddSS(int x, SS *ss, double w=1.0) const;
	virtual void AddExpSS(SS *ss, double w=1.0) const;
	virtual void AddSS(const SS *toadd, const RVSimple *rvs,
				const std::vector<std::vector<int> > &mapping,
				SS *ss, double w=1.0) const;
	virtual int Sample(Random &rand = randomizer) const;
	virtual void Maximize(const SS *ss);
	virtual void GetDist(vectr &d, double &logfactor) const;
	virtual void SetDist(const vectr &d, double logfactor=0.0);
	
	//Busra
	virtual void GetSparseDist(vectr &d, vectr &ind, double &logfactor) const;
	virtual void SetSparseDist(const vectr &d, const vectr &ind, double logfactor=0.0);
	//
	
	virtual void SetSparseDist(const vectr &d, double logfactor=0.0);
	virtual void GetSparseDist(vectr &d, double &logfactor) const;
	
	virtual void SetSparseIndexes(const std::vector<int> &ind);
	virtual double LLH(const SS *ss) const;

	inline const std::vector<int> &Indexes() const { return ix; }
	inline double LogFactor() const { return logz; }
	inline const vectr &Dist() const { return theta; }

protected:
	std::vector<int> ix;

	SERIAL_START_V(SparseMultiZSimple)
		SERIAL_SUPER(MultiZSimple)
		SERIAL_VAR(std::vector<int>,ix)
	SERIAL_END
	
};


} // end of ctbn namespace

#endif
