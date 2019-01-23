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
#ifndef CTBNRLE_UNIFORMIZEDFACTOREDINF_H
#define CTBNRLE_UNIFORMIZEDFACTOREDINF_H

#include "markov.h"
#include "ctbndyn.h"
#include "dyncomp.h"
#include "bn.h"
#include "rvcomp.h"
#include "inference.h"
#include "fbinf.h"
#include "factoredvector.h"
#include "factoredmatrix.h"

#include <vector>


namespace ctbn {

// This class performs approximate forward inference (filtering) only
// (currently) by using the uniformized approximation of Celikkaya & Shelton
// (2010)

class UniformizedFactoredInf : public FBInf {
public:
    UniformizedFactoredInf();
	virtual ~UniformizedFactoredInf() throw();
	virtual UniformizedFactoredInf *Clone() const;
	void SetL(int l){ lval = l; }
	void SetTheta(double th) { thval = th; }

private:

	RVSimple * Convert(RVCondSimple *&rvcond, RVSimple *&rv);
	void Restrict(FactoredVector * & rv, Instantiation const & x) const;
	
	
	int lval;
	double thval;
public:

	virtual void Restrict(RVSimple *, Instantiation const & variable_assignment) const;
    
	virtual RVSimple * MakeSimple(RV * prior, Instantiation const &, bool normalize = true) const;
    virtual bool IsAlphaElementValid(RVSimple * const alpha_sub_x);

	virtual FactoredMatrix * IntervalPropagator(Dynamics *, double, double, Instantiation const &) const;
	virtual FactoredMatrix * PointTransitionPropagator(Dynamics *, double, Instantiation const &, Instantiation const &) const;
	virtual FactoredMatrix * PointChangeEvidencePropagator(Dynamics *, double, Instantiation const &, Instantiation const &) const;

	void AddExpectedSufficientStatistics(Dynamics const *, unsigned int, SS *, double) const { throw "Not yet implemented."; }
	void AddExpectedTransitionSufficientStatistics(Dynamics const *, unsigned int, SS *, double) const { throw "Not yet implemented."; }
	void AddExpectedInitialSufficientStatistics(Trajectory const *, RVSimple *, const RV *, SS *, double) const { throw "Not yet implemented."; }
};

} // end of ctbn namespace

#endif

