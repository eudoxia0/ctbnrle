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
#include "exactmarkovinf.h"

namespace ctbn {


	ExactMarkovInf::ExactMarkovInf():
	 FBInf() {
		// Except for invoking the ancestor's constructor, there is nothing
		// else to do.
	}


	ExactMarkovInf::~ExactMarkovInf() throw() {
		// There is no cleanup required above and beyond what the ancestor's
		// class destructor will already clean up for us.
	}


	ExactMarkovInf * ExactMarkovInf::Clone() const {
		ExactMarkovInf * the_clone = new ExactMarkovInf(*this);
		CloneTo(the_clone);

		// This class has no additional data members to be cloned.
		return the_clone;
	}

	void ExactMarkovInf::Restrict(RVSimple * a, Instantiation const & variable_assignment) const {
			vector<int> ind;
			the_markov_process_dynamics()->Domain().ConsistentIndexes(ind, variable_assignment);
			a->Restrict(ind);
	}


    RVSimple * ExactMarkovInf::MakeSimple(RV * prior, Instantiation const &, bool normalize) const {
        Trajectory const * const the_trajectory(getTrajectory());
        Instantiation const trajectory_values(
         the_trajectory->Values(prior->Domain() + prior->CondDomain(),
                                the_trajectory->TimeBegin()));
        RVSimple * ret = prior->MakeSimple(trajectory_values);
        return ret;
    }


	RVCondSimple * ExactMarkovInf::IntervalPropagator(Dynamics * const the_markov_process_dynamics, double start_time_point, double past_end_time_point, Instantiation const & current_evidence) const {
		RVCondSimple * const ret(the_markov_process_dynamics->Cond(start_time_point, past_end_time_point, current_evidence));
		return ret;
	}


	RVCondSimple * ExactMarkovInf::PointTransitionPropagator(Dynamics * const the_markov_process_dynamics, double time_point, Instantiation const & prior_evidence, Instantiation const & current_evidence) const {
		RVCondSimple * const ret(the_markov_process_dynamics->Cond(time_point, prior_evidence, current_evidence));
		return ret;
	}


	RVCondSimple * ExactMarkovInf::PointChangeEvidencePropagator(Dynamics * const the_markov_process_dynamics, double time_point, Instantiation const & prior_evidence, Instantiation const & current_evidence) const {
		RVCondSimple * const ret(the_markov_process_dynamics->Cond(time_point, prior_evidence, current_evidence, false));
		return ret;
	}


	void ExactMarkovInf::AddExpectedSufficientStatistics(Dynamics const * const dynamics, unsigned int const trajectory_index, SS * const sufficient_statistics, double const weight) const {
		unsigned int const & i(trajectory_index);
		dynamics->AddExpSS(condi_at(i), alpha_at(i), beta_at(i+1), time_points_at(i), time_points_at(i+1) - time_points_at(i), sufficient_statistics, weight);
	}


	void ExactMarkovInf::AddExpectedTransitionSufficientStatistics(Dynamics const * const dynamics, unsigned int const trajectory_index, SS * const sufficient_statistics, double const weight) const {
		unsigned int const & i(trajectory_index);
		dynamics->AddExpTransSS(condi_at(i), alpha_at(i), beta_at(i+1), evidence_variables_whose_observed_values_changed_at(i), time_points_at(i), sufficient_statistics, weight);
	}


	void ExactMarkovInf::AddExpectedInitialSufficientStatistics(Trajectory const * const trajectory, RVSimple * initial_prob_dist, const RV *ssrv, SS * const sufficient_statistics, double const weight) const {
		ssrv->AddExpSS(trajectory->Values(ssrv->CondDomain(), trajectory->TimeBegin()), initial_prob_dist, sufficient_statistics, weight);
	}


    bool ExactMarkovInf::IsAlphaElementValid(RVSimple * const alpha_sub_x) {
        SparseMultiZSimple * as_sparse_multi_z_simple = dynamic_cast<SparseMultiZSimple *>(alpha_sub_x);
        bool ret = ((nullptr03 != as_sparse_multi_z_simple) && (as_sparse_multi_z_simple->Dist().isvalid()));
        return ret;
    }

} // end of ctbn namespace
