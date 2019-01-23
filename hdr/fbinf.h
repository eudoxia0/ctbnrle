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
#ifndef CTBNRLE_FBINF_H
#define CTBNRLE_FBINF_H

#include "inference.h"
#include "markov.h"
#include "nullptr03.h"
#include "trajectory.h"

#include <vector>


namespace ctbn {

	using std::vector;

	// This class implements the forward-backward algorithm.

	class FBInf : public Inference {

	public:
		FBInf();
		FBInf(FBInf const * const copy_from);
		virtual ~FBInf() throw();

		// Virtual copy constructor (object clone generator).
		virtual FBInf * Clone() const = 0;

	protected:
		// The next function is NOT a virtual copy "constructor", but rather
		// common code that must be called by any virtual copy "constructors"
		// defined in descendents of this class to initialize fields that are
		// private to this class.  The single argument must be of the
		// appropriate (most specialized) descendant class.
		void CloneTo(FBInf * const copy_to) const;

	public:
		// See inference.h for ownership of the process and trajectory
		// in the following functions
		virtual void SetProcess(const Process * some_process);
		virtual void SetTrajectory(const Trajectory * some_trajectory);

		virtual double Filter(const Instantiation & x, double t,
		 bool log = false);
		virtual double Smooth(const Instantiation & x, double t,
		 bool log = false);

		virtual double Prob(double t, bool log=false);

		virtual void AddExpSuffStats(SS *ss, double w = 1.0);
		virtual void AddExpSuffStats(const Dynamics *dyn, SS *ss,
								double w = 1.0);
		virtual void AddExpSuffStats(const RV *p0, SS *ss,
								double w = 1.0);

		virtual void printdist(std::ostream &os) const {
			the_markov_process->SaveOld(os);
		}

		double CalcQuery(QueryCalculator &calc);

	protected:
		// Function interfaces for Forward-Backward Operations
		
		virtual void Restrict(RVSimple *, Instantiation const & variable_assignment) const = 0;
		
		virtual RVSimple * MakeSimple(RV * prior, Instantiation const &, bool normalize) const = 0;
		virtual bool IsAlphaElementValid(RVSimple * const) = 0;

		// related to P(X_t+delta_t | X_t, E_0:t)
		virtual RVCondSimple * IntervalPropagator(Dynamics * const the_markov_process_dynamics, double start_time_point, double past_end_time_point, Instantiation const & current_evidence) const = 0;
		// related to P(X_t+delta_t | X_t, E_0:t, E_t+delta_t)
		virtual RVCondSimple * PointTransitionPropagator(Dynamics * const the_markov_process_dynamics, double time_point, Instantiation const & prior_evidence, Instantiation const & current_evidence) const = 0;
		// related to lim epsilon->0 P(X_t+epsilon | X_t-epsilon, E_0:t-epsilon, E_t+epsilon)
		virtual RVCondSimple * PointChangeEvidencePropagator(Dynamics * const the_markov_process_dynamics, double time_point, Instantiation const & prior_evidence, Instantiation const & current_evidence) const = 0;

		virtual void AddExpectedSufficientStatistics(Dynamics const * const dynamics, unsigned int const trajectory_index, SS * const sufficient_statistics, double const weight = 1.0) const = 0;
		virtual void AddExpectedTransitionSufficientStatistics(Dynamics const * const dynamics, unsigned int const trajectory_index, SS * const sufficient_statistics, double const weight = 1.0) const = 0;
		virtual void AddExpectedInitialSufficientStatistics(Trajectory const * const trajectory, RVSimple * initial_prob_dist, const RV *ssrv, SS * const sufficient_statistics, double const weight = 1.0) const = 0;

	protected:

		inline RV * the_initial_distribution_of_the_markov_process() const {
			return (nullptr03 == the_markov_process ? nullptr03 : the_markov_process->p0);
		};

		inline Dynamics * the_markov_process_dynamics() const {
			return (nullptr03 == the_markov_process ? nullptr03 : the_markov_process->d);
		};

		inline Trajectory const * getTrajectory() const {
	        return the_trajectory;
	    };

		inline double time_points_at(unsigned int i) const {
			return time_points[i];
		}

		inline RVSimple const * alpha_at(unsigned int i) const {
			return alpha[i];
		};

		inline RVSimple const * beta_at(unsigned int i) const {
			return beta[i];
		};

		inline RVCondSimple const * prop_at(unsigned int i) const {
			return prop[i];
		};

		inline int condi_at(unsigned int i) const {
			return condi[i];
		};

		inline EvidenceChangeEnum transition_kind_at(unsigned int i) const {
			return transition_kind[i];
		};

		inline Context const & evidence_variables_whose_observed_values_changed_at(unsigned int i) const {
			return evidence_variables_whose_observed_values_changed[i];
		}

	private:

		Markov const * the_markov_process;
		Trajectory const * the_trajectory;

		/* Deallocate each element of the vector from the heap. */
		template <typename T>
		inline void delete_each(vector<T*> & the_vector) {
			for (unsigned int i(0); i != the_vector.size(); ++i) {
				delete the_vector[i];
				the_vector[i] = nullptr03;
			}
		}

		/* Ensure that the capacity of the vector becomes zero. */
		template <typename T>
		inline void free_memory(vector<T> & the_vector) {
			vector<T> empty_vector;
			swap(the_vector, empty_vector);
		}

		double FilterFuncImpl(Instantiation const & variable_assignment,
		 double time_point_of_interest, bool use_log_probs = false,
		 bool normalize_probs = true);

		void ClearCache();
		// make up until t
		void EnsurePassesPerformedAndResultsCached();
		void FetchInformationRequiredForPasses();
		void PerformForwardAlphaPass();
		void PerformBackwardBetaPass();

		// the cache
		std::vector<double> time_points;
		std::vector<RVSimple *> alpha;
		std::vector<RVSimple *> beta;
		std::vector<RVCondSimple *> prop;
		std::vector<int> condi;
		std::vector<EvidenceChangeEnum> transition_kind;
		std::vector<Context> evidence_variables_whose_observed_values_changed;

		// Cache the joint SS
		SS* jointP0SS;
		SS* jointDSS;
	};

}

#endif
