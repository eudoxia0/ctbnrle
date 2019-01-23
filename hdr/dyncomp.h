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
#ifndef CTBNRLE_DYNCOMP_H
#define CTBNRLE_DYNCOMP_H

#include "streamextra.h"
#include "streamserial.h"
#include "dynamics.h"
#include "dynsimple.h"



namespace ctbn {

class DynComp;

// sufficient statistics for a DynComp
class DynCompSS : public SS {
	SOBJCLASSDECL(DynCompSS)
	friend class DynComp;
public:
	DynCompSS(int n=0);
	DynCompSS(const DynCompSS &dcss);
	DynCompSS(std::istream &is);
	DynCompSS &operator=(const DynCompSS &dcss);
	virtual ~DynCompSS();

	virtual DynCompSS *Clone() const;
	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
	virtual void Scale(double w);
	virtual void AddSS(const SS* nss, double w=1.0);
	virtual double Element(std::vector<int> index, int s1, int s2) const;
private:
	std::vector<SS *> ss;
	SERIAL_START_V(DynCompSS)
		SERIAL_SUPER(SS)
		SERIAL_VAR(std::vector<SS*>,ss)
	SERIAL_END
};

// a class that creates a Dynamics object out of a DynSimple class
// see dynamics.h for a description of the methods
// see dynsimple.h for a description of DynSimple
// This class is essentially just a set of DynSimple objects (one
// for each possible value of the conditioning context), with
// wrappers to map from instantiations to the raw indexes used by 
// DynSimple
class DynComp : public Dynamics {
	SOBJCLASSDECL(DynComp)
public:
	// the base argument will only be used for Clone()
	// it is still owned by the caller
	DynComp(const Context &var=Context(), const Context &cvar=Context(), 
			const DynSimple *base=NULL);
	DynComp(std::istream &is);
	DynComp(const DynComp &dc);
	DynComp &operator=(const DynComp &dc);
	virtual ~DynComp();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
	virtual DynComp *Clone() const;
	virtual void Normalize();
	virtual void Restrict(const Instantiation &x);
	virtual void Mult(const Dynamics *x);
	virtual RVCondSimple *Cond(double t0, double t1,
			const Instantiation &x) const;
	virtual RVCondSimple *Cond(double t, const Instantiation &from,
			const Instantiation &to, bool transition=true) const;
	virtual DynCompSS *BlankSS() const;
	virtual void AddSS(const Instantiation &x, double t0,
			double deltat, SS *ss, double w=1.0) const;
	virtual void AddTransSS(const Instantiation &x1, 
				const Instantiation &x2,
				double t, SS *ss, double w=1.0) const;
	virtual void AddExpSS(const RV *alpha, const RV *beta,
			double t0, double deltat, SS *ss, double w=1.0) const;
	virtual void AddExpTransSS(const RV *x1, const RV *x2,
			const Context &changevar,
			double t, SS *ss, double w=1.0) const;
	virtual void AddExpSS(const Instantiation &condi,
			const RVSimple *alpha, const RVSimple *beta,
			double t0, double deltat, SS *ss, double w=1.0) const;
	virtual void AddExpSS(int condi,
			const RVSimple *alpha, const RVSimple *beta,
			double t0, double deltat, SS *ss, double w=1.0) const;
	virtual void AddExpTransSS(const Instantiation &condi,
			const RVSimple *x1, const RVSimple *x2,
			const Context &changevar,
			double t, SS *ss, double w=1.0) const;
	virtual void AddExpTransSS(int condi,
			const RVSimple *x1, const RVSimple *x2,
			const Context &changevar,
			double t, SS *ss, double w=1.0) const;
	virtual void AddSS(const SS *toadd, const Dynamics *dyn,
			SS *ss, double w=1.0) const;

	virtual void SampleNextEvent(const Instantiation &i, double t,
					double &nextt, 
					Instantiation &nexti, 
					Random &rand=randomizer) const;
	virtual void Maximize(const SS *ss);

	inline DynSimple *operator[](const Instantiation &i) {
		return impl[CondDomain().Index(i)];
	}
	inline const DynSimple *operator[](const Instantiation &i) const {
		return impl[CondDomain().Index(i)];
	}

	virtual void Scramble(double a=1.0, double b=1.0, 
						  double alpha=1.0, double degree=1.0, 
						  Random &rand=randomizer);
	virtual double LLH(const SS *ss) const;
	
	virtual double GetScore(double numTrans, double amtTime, 
			const SS* ss) const;

	// It is questionable as to whether this breaks the abstraction...
	// note that the pointer returned is still owned by this object
	inline DynSimple *operator[](int i) {
		return impl[i];
	}
	inline const DynSimple *operator[](int i) const {
		return impl[i];
	}
protected:
	inline int CondIndex(const Instantiation &x) const {
		 if (impl.size()<2) return 0;
		 else return CondDomain().Index(x);
	}


	std::vector<DynSimple *> impl;
	SERIAL_START_V(DynComp)
		SERIAL_SUPER(Dynamics)
		SERIAL_VAR(std::vector<DynSimple*>,impl)
	SERIAL_END
};

} // end of ctbn namespace

#endif
