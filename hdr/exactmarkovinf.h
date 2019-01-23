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
#ifndef CTBNRLE_EXACTMARKOVINF_H
#define CTBNRLE_EXACTMARKOVINF_H

#include "fbinf.h"
#include "markov.h"
#include "ctbndyn.h"
#include "dyncomp.h"
#include "bn.h"
#include "rvcomp.h"
#include "multisimple.h"
#include "inference.h"
#include <vector>



namespace ctbn {

// This class performs exact inference by using the conditional distributions
// supplied by the process's (assumed to be a Markov process) dynamics's
// method "Cond."  For compact structures (like a CTBN) this will probably
// require generation of the full "flattened" dynamics which is exponentially
// (in the number of variables) large.  However, for some processes (like a
// fully factored process), it might be possible to return a conditional
// object that can represent this more compactly.
//
// The dynamics's "cond" methods should return RVCondSimple objects that
// are compatible with the RVSimple object's "MakeSimple" method (in that
// the RVCondSimple can accept the output of the MakeSimple method).

class ExactMarkovInf : public FBInf {
public:
	ExactMarkovInf();
	virtual ~ExactMarkovInf() throw();

	virtual ExactMarkovInf * Clone() const;

	virtual void Restrict(RVSimple *, Instantiation const & variable_assignment) const;

	virtual RVSimple * MakeSimple(RV * prior, Instantiation const &, bool normalize) const;
    virtual bool IsAlphaElementValid(RVSimple * const alpha_sub_x);

	virtual RVCondSimple * IntervalPropagator(Dynamics * const the_markov_process_dynamics, double start_time_point, double past_end_time_point, Instantiation const & current_evidence) const;
	virtual RVCondSimple * PointTransitionPropagator(Dynamics * const the_markov_process_dynamics, double time_point, Instantiation const & prior_evidence, Instantiation const & current_evidence) const;
	virtual RVCondSimple * PointChangeEvidencePropagator(Dynamics * const the_markov_process_dynamics, double time_point, Instantiation const & prior_evidence, Instantiation const & current_evidence) const;

	virtual void AddExpectedSufficientStatistics(Dynamics const * const dynamics, unsigned int const trajectory_index, SS * const sufficient_statistics, double const weight = 1.0) const;
	virtual void AddExpectedTransitionSufficientStatistics(Dynamics const * const dynamics, unsigned int const trajectory_index, SS * const sufficient_statistics, double const weight = 1.0) const;
	virtual void AddExpectedInitialSufficientStatistics(Trajectory const * const trajectory, RVSimple * initial_prob_dist, const RV * ssrv, SS * const sufficient_statistics, double const weight = 1.0) const;

private:
    SparseMultiZSimple * dynamic_cast_to_SparseMultiZSimple(RVSimple * const ptr);

};

} // end of ctbn namespace

#endif
