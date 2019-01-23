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
#ifndef CTBNRLE_RVCOMP_H
#define CTBNRLE_RVCOMP_H

#include "rv.h"
#include "ss.h"
#include <vector>

namespace ctbn {

class RVSimple;
class RVCondSimple;

// a class that creates a RV out of an RVCondSimple
// essentially, it adds mappings from Instantiations to integers
// (Context::Index()) to wrap the RVCondSimple
class RVComp : public RV {
	SOBJCLASSDECL(RVComp)
public:
	// touse becomes the property of RVComp.  Call with touse->Clone()
	// to avoid this
	RVComp(const Context &var=Context(),
		const Context &cvar=Context(), RVCondSimple *touse=NULL);
	RVComp(std::istream &is);
	RVComp(const RVComp &rvc);
	RVComp &operator=(const RVComp &rvc);
	virtual ~RVComp();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
	virtual RVComp *Clone() const;
	virtual void Normalize();
	virtual double Prob(const Instantiation &x, bool log=false) const;
	virtual void Restrict(const Instantiation &x);
	virtual void Project(const Context &c, const Context &cc);
	virtual void MultBy(const RV *x);
	virtual SS *BlankSS() const;
	virtual void AddSS(const Instantiation &x, SS *ss, double w=1.0) const;
	virtual void AddExpSS(const Instantiation &x, SS *ss, double w=1.0) const;
	virtual void AddExpSS(const Instantiation &x, const RVSimple *rvs, 
			SS *ss, double w=1.0) const;
	//This adds the sufficient statistics of a different random variable
	//to this random variable's sufficient statistics
	virtual void AddSS(const SS* toadd, const RV *rv,
			SS *ss, double w=1.0) const;
	virtual void Sample(Instantiation &x, Random &rand = randomizer) const;
	virtual void Maximize(const SS *ss);

	RVCondSimple *Base() { return impl; }
	const RVCondSimple *Base() const { return impl; }

	virtual RVSimple * MakeSimple(Instantiation const & x, bool normalize = true) const;

	virtual void Scramble(double alpha=1.0, double degree=1.0,
		    	Random &rand=randomizer);
	virtual double LLH(const SS *ss) const;

	virtual double GetScore(double numTrans, const SS* ss) const;

private:
	RVCondSimple *impl;

	SERIAL_START_V(RVComp)
		SERIAL_SUPER(RV)
		SERIAL_VAR(RVCondSimple*,impl)
	SERIAL_END

public:
	void serial_preload();
	
};

} // end of ctbn namespace

#endif
