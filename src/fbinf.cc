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


#include "fbinf.h"

#include "factoredvector.h"
#include "factoredmatrix.h"
#include "dyncomp.h"
#include "multisimple.h"
#include "nullptr03.h"
#include "rvcomp.h"
#include <algorithm>

using namespace std;


namespace ctbn {

	double ZeroProb(bool use_log_probs);
	SparseMultiZSimple * dynamic_cast_to_SparseMultiZSimple(RVSimple * const ptr);


	FBInf::FBInf():
	 the_markov_process(nullptr03),
	 the_trajectory(nullptr03),
	 jointP0SS(nullptr03),
	 jointDSS(nullptr03) {}


	FBInf::FBInf(FBInf const * const copy_from):
	 the_markov_process(nullptr03 == copy_from->the_markov_process ? nullptr03 : copy_from->the_markov_process->Clone()),
	 the_trajectory(copy_from->the_trajectory),
	 time_points(copy_from->time_points),
	 alpha(copy_from->alpha),
	 beta(copy_from->beta),
	 prop(copy_from->prop),
	 condi(copy_from->condi),
	 transition_kind(copy_from->transition_kind),
	 evidence_variables_whose_observed_values_changed(copy_from->evidence_variables_whose_observed_values_changed),
	 jointP0SS(nullptr03 == copy_from->jointP0SS ? nullptr03 : copy_from->jointP0SS->Clone()),
	 jointDSS(nullptr03 == copy_from->jointDSS ? nullptr03 : copy_from->jointDSS->Clone()) {}


	FBInf::~FBInf() throw() {
		ClearCache();
		// Do not delete the_trajectory: we don't own it (see inference.h).
		// Do not delete the_markov_process: we don't own it (see inference.h).
	}


	void FBInf::CloneTo(FBInf * const copy_to) const {

		copy_to->alpha.clear();
		copy_to->beta.clear();
		copy_to->prop.clear();
		for (unsigned int i(0); i != alpha.size(); ++i)
			copy_to->alpha.push_back(alpha[i]->Clone());
		for (unsigned int i(0); i != beta.size(); ++i)
			copy_to->beta.push_back(beta[i]->Clone());
		for (unsigned int i(0); i != prop.size(); ++i)
			copy_to->prop.push_back(prop[i]->Clone());

		if (jointP0SS)
			copy_to->jointP0SS = jointP0SS->Clone();
		if (jointDSS)
			copy_to->jointDSS = jointDSS->Clone();
	}


	void FBInf::SetProcess(Process const * const some_markov_process) {
		if (this->the_markov_process != some_markov_process) {
			ClearCache();
			this->the_markov_process=nullptr03;
			this->the_markov_process = dynamic_cast<const Markov *>(some_markov_process);
		}
	}


	void FBInf::SetTrajectory(const Trajectory * some_trajectory) {
		if (this->the_trajectory != some_trajectory) {
			ClearCache();
			this->the_trajectory = some_trajectory;
		}
	}


	double FBInf::Prob(double time_point_of_interest, bool use_log_probs) {
		return FilterFuncImpl(Instantiation(), time_point_of_interest, use_log_probs, false);
	}


	double FBInf::Filter(const Instantiation & variable_assignment, double time_point_of_interest, bool use_log_probs) {
		return FilterFuncImpl(variable_assignment, time_point_of_interest, use_log_probs);
	}


	double FBInf::FilterFuncImpl(Instantiation const & variable_assignment, double const time_point_of_interest, bool const use_log_probs, bool const normalize_probs) {
		if (the_markov_process == nullptr03 || the_trajectory == nullptr03) {
			return ZeroProb(use_log_probs);
		}
		EnsurePassesPerformedAndResultsCached();

		unsigned int const first_index_where_time_point_is_not_less_than_time_point_of_interest = lower_bound(time_points.begin(), time_points.end(), time_point_of_interest) - time_points.begin();
		unsigned int const & i = first_index_where_time_point_is_not_less_than_time_point_of_interest;
		RVSimple *a;
		if (time_points[i] == time_point_of_interest) {
			a = alpha[i]->Clone();
		} else {
			if (i == 0) {
				return ZeroProb(use_log_probs);
			}
            
			a = alpha[i-1]->Clone();
			/*RVCondSimple *trans = the_markov_process_dynamics()->Cond(time_points[i-1], time_point_of_interest,
				the_trajectory->Values(the_markov_process_dynamics()->Domain() + the_markov_process_dynamics()->CondDomain(), time_point_of_interest));
				*/
			RVCondSimple *trans = IntervalPropagator(the_markov_process_dynamics(), time_points[i-1], time_point_of_interest,
				the_trajectory->Values(the_markov_process_dynamics()->Domain() + the_markov_process_dynamics()->CondDomain(), time_point_of_interest));
			trans->Mult(a);
			delete trans;
		}

		if (normalize_probs) a->Normalize();
		double ret;
		if (variable_assignment.NumVars() == 0)
			ret = a->Sum(use_log_probs);
		else {
			Restrict(a, variable_assignment);
			/*
			vector<int> ind;
			the_markov_process_dynamics()->Domain().ConsistentIndexes(ind, variable_assignment);
			a->Restrict(ind);
			*/
			ret = a->Sum(use_log_probs);
			
		}
		delete a;
		return ret;
	}


	double FBInf::Smooth(const Instantiation & variable_assignment, double time_point_of_interest, bool use_log_probs) {
		if (the_markov_process == nullptr03 || the_trajectory == nullptr03) {
			return ZeroProb(use_log_probs);
		}
		EnsurePassesPerformedAndResultsCached();

		unsigned int const first_index_where_time_point_is_not_less_than_time_point_of_interest = lower_bound(time_points.begin(), time_points.end(), time_point_of_interest) - time_points.begin();
		unsigned int const & i = first_index_where_time_point_is_not_less_than_time_point_of_interest;
		RVSimple *a;
		if (time_points[i] == time_point_of_interest) {
			a = alpha[i]->Clone();
			a->MultBy(beta[i]);
		} else {
			if (i==0) {
				return ZeroProb(use_log_probs);
			}

			a = alpha[i-1]->Clone();
			Instantiation inst = the_trajectory->Values(the_markov_process_dynamics()->Domain() + the_markov_process_dynamics()->CondDomain(), time_point_of_interest);
			//RVCondSimple *trans = the_markov_process_dynamics()->Cond(time_points[i-1], time_point_of_interest, inst);
			RVCondSimple *trans = IntervalPropagator(the_markov_process_dynamics(),time_points[i-1], time_point_of_interest, inst);
			trans->Mult(a);
			delete trans;

			RVSimple *b = (i==time_points.size() ? beta[i-1]->Clone() : beta[i]->Clone());
			//trans = the_markov_process_dynamics()->Cond(time_point_of_interest, time_points[i], inst);
			trans = IntervalPropagator(the_markov_process_dynamics(), time_point_of_interest, time_points[i], inst);
			trans->RMult(b);
			delete trans;
			a->MultBy(b);
			delete b;
		}

		a->Normalize();
		double ret;
		if (variable_assignment.NumVars()==0)
			ret = a->Sum(use_log_probs);
		else {
			vector<int> ind;
			the_markov_process_dynamics()->Domain().ConsistentIndexes(ind,variable_assignment);
			a->Restrict(ind);
			ret = a->Sum(use_log_probs);
		}
		delete a;
		return ret;
	}


	void FBInf::ClearCache() {
		delete_each(alpha);
		free_memory(alpha);
		delete_each(beta);
		free_memory(beta);
		delete_each(prop);
		free_memory(prop);

		free_memory(time_points);
		free_memory(condi);
		free_memory(transition_kind);
		free_memory(evidence_variables_whose_observed_values_changed);

		if (jointP0SS) {
			delete jointP0SS;
			jointP0SS = nullptr03;
		}
		if (jointDSS) {
			delete jointDSS;
			jointDSS = nullptr03;
		}
	}


	void FBInf::EnsurePassesPerformedAndResultsCached() {
		if (time_points.empty()) {
			FetchInformationRequiredForPasses();
			PerformForwardAlphaPass();
			PerformBackwardBetaPass();
		}
	}


	/* This method retrieves the entire sequence of probability distributions
	 * that are conditioned on the evidence in the trajectory.  This includes
	 * not only transitions in time from one event to the next, but also
	 * transitions where time does not elapse but an event occurs (for example,
	 * the new availability or unavailability of part of the evidence). */
	void FBInf::FetchInformationRequiredForPasses() {
		time_points.push_back(the_trajectory->TimeBegin());
		Context const initial_context(the_markov_process_dynamics()->Domain() + the_markov_process_dynamics()->CondDomain());
		Trajectory::Index i = the_trajectory->Begin(initial_context);
		while (!i.Done()) {
			double const start_time_point(i.Time());
			double const past_end_time_point(start_time_point + i.DeltaT());
			Instantiation const current_evidence(i.Values());
//			RVCondSimple * prob(the_markov_process_dynamics()->Cond(start_time_point, past_end_time_point, current_evidence));
			RVCondSimple * prob(IntervalPropagator(
			 the_markov_process_dynamics(), start_time_point,
			 past_end_time_point, current_evidence));
			prop.push_back(prob);
			transition_kind.push_back(unchanged_with_respect_to_context);
			evidence_variables_whose_observed_values_changed.push_back(Context());
			condi.push_back(the_markov_process_dynamics()->CondDomain().Index(current_evidence));

			EvidenceChangeEnum const varchange(i.TestInc(the_markov_process_dynamics()->Domain()));
			double const double_equality_comparison_epsilon = 1.0e-6;
			assert(i.Time() - past_end_time_point < double_equality_comparison_epsilon);
			time_points.push_back(past_end_time_point);

			if (i.Done()) break;

			Instantiation const new_evidence(i.Values());

			switch (varchange) {
			case unchanged_with_respect_to_context:
				// The conditioning set has changed, but nothing else has changed.
				// Therefore, nothing need be done at this moment (pun intended).
				// We merely want to proceed to the next iteration of the loop.
				continue;
			case value_of_an_observed_variable_changed:
				{
//					prob = the_markov_process_dynamics()->Cond(i.Time(), current_evidence, i.Values());
					prob = PointTransitionPropagator(the_markov_process_dynamics(), past_end_time_point, current_evidence, new_evidence);
					evidence_variables_whose_observed_values_changed.push_back(Context(current_evidence, new_evidence));
				}
				break;
			case only_variable_observability_changed:
				{
//					prob = the_markov_process_dynamics()->Cond(i.Time(), current_evidence, new_evidence, false);
					prob = PointChangeEvidencePropagator(the_markov_process_dynamics(), past_end_time_point, current_evidence, new_evidence);
					evidence_variables_whose_observed_values_changed.push_back(Context());
				}
				break;
			default:
				// Unanticipated circumstance.
				throw;
			}
			prop.push_back(prob);
			condi.push_back(the_markov_process_dynamics()->CondDomain().Index(new_evidence));
			time_points.push_back(past_end_time_point);
			transition_kind.push_back(varchange);
		}
	}


	// This method propagates the classic "alpha pass" starting with the initial
	// distribution indicated by the Markov process's starting distribution and
	// running it forward through the conditional distributions set up by
	// FetchInformationRequiredForPasses.  They are saved for easy access later.
	void FBInf::PerformForwardAlphaPass() {
		unsigned int n = prop.size();
		alpha.resize(n+1);

		RV * prior = the_initial_distribution_of_the_markov_process();
		Instantiation const trajectory_values(
		 the_trajectory->Values(prior->Domain() + prior->CondDomain(),
		                        the_trajectory->TimeBegin()));
		// TODO: replace with descendent's MakeSimple
//		RVSimple *curr = prior->MakeSimple(trajectory_values);
		RVSimple * curr = MakeSimple(prior, trajectory_values, true);
/*
        if (dynamic_cast<FactoredVector *>(curr)) {
		    cout << "yay, curr is a factored vector." << endl;
		} else {
		    cout << "warning, curr is not a factored vector." << endl;
		}

        if (dynamic_cast<FactoredMatrix *>(prop[0])) {
            cout << "yay, prop[0] is a factored matrix." << endl;
        } else {
            cout << "warning, prop[0] is not a factored matrix." << endl;
        }
*/
		if (n > 0) {
		    alpha[0] = prop[0]->Convert(curr);
			delete curr;
			curr = alpha[0]->Clone();
		} else {
			alpha[0] = curr->Clone();
		}

		if (!IsAlphaElementValid(alpha[0])) {
		    cerr << "alpha[0] is invalid" << endl;
		}

		for (unsigned int i = 0; i != n; ++i) {
			prop[i]->Mult(curr);
			alpha[i + 1] = curr->Clone();
			if (!IsAlphaElementValid(alpha[i + 1])) {
				cerr << "alpha[" << i + 1 << "] is invalid" << endl;
			}
		}

		delete curr;
	}


	// This is the same as PerformForwardAlphaPass, but running the
	// distribution (starting with a uniform one) backward through the
	// conditional distributions.
	void FBInf::PerformBackwardBetaPass() {
		int i = prop.size();
		beta.resize(i+1);
		RVSimple *curr = alpha[i]->Clone();
		curr->MakeUniform(); // set to all 1s
		beta[i] = curr->Clone();
		for(i--;i>=0;i--) {
			prop[i]->RMult(curr);
			beta[i] = curr->Clone();
		}
		delete curr;
	}


	void FBInf::AddExpSuffStats(SS *ss, double w) {

		MarkovSS * mss = dynamic_cast<MarkovSS *>(ss);
		EnsurePassesPerformedAndResultsCached();
		int const n(prop.size());
		for(int i = 0; i != n; ++i) {
			// we could construct the correct full RV and call AddExpSS on the
			// markov process dynamics, but that's a lot of extra unneeded
			// work given that we know the underlying structure
			switch (transition_kind[i]) {
			case unchanged_with_respect_to_context:
/*
				the_markov_process_dynamics()->AddExpSS(condi[i], alpha[i],
				 beta[i+1], time_points[i], time_points[i+1] - time_points[i],
				 mss->dss, w);
*/
				AddExpectedSufficientStatistics(
				 the_markov_process_dynamics(), i, mss->dss, w);
				break;
			case value_of_an_observed_variable_changed:
/*
				the_markov_process_dynamics()->AddExpTransSS(condi[i], alpha[i],
				 beta[i+1],	evidence_variables_whose_observed_values_changed[i],
				 time_points[i], mss->dss, w);
*/
				AddExpectedTransitionSufficientStatistics(
				 the_markov_process_dynamics(), i, mss->dss, w);
				break;
			case only_variable_observability_changed:
				break;
			default:
				throw;
			}
		}

		if (!alpha.empty()) {
			// add suffstats for initial distribution
			RVSimple * initial_prob_dist = alpha[0]->Clone();
			initial_prob_dist->MultBy(beta[0]);
			initial_prob_dist->Normalize();
			RV * prior = the_initial_distribution_of_the_markov_process();
/*
			prior->AddExpSS(the_trajectory->Values(prior->CondDomain(),
			 the_trajectory->TimeBegin()), initial_prob_dist, mss->p0ss, w);
*/
			AddExpectedInitialSufficientStatistics(
			 the_trajectory, initial_prob_dist, prior, mss->p0ss, w);
			delete initial_prob_dist;
		}
	}


	void FBInf::AddExpSuffStats(const Dynamics *dyn, SS *ss, double w) {
		CTBNDyn const * const cdyn = dynamic_cast<CTBNDyn const *>(the_markov_process_dynamics());
		DynComp const * const jcdyn = dynamic_cast<DynComp const *>(cdyn->GetJoint());
		if(!jointDSS) {
			DynCompSS * const dcss = dynamic_cast<DynCompSS *>(jcdyn->BlankSS());
			EnsurePassesPerformedAndResultsCached();
			int const n(prop.size());
			for(int i(0); i != n; ++i) {
				if (transition_kind[i]==0)
/*
					jcdyn->AddExpSS(condi[i],alpha[i],beta[i+1],
					 time_points[i],time_points[i+1]-time_points[i],dcss,w);
*/
					AddExpectedSufficientStatistics(jcdyn, i, dcss, w);
				else if (transition_kind[i]==2)
/*
					jcdyn->AddExpTransSS(condi[i],alpha[i],beta[i+1],
					 evidence_variables_whose_observed_values_changed[i],
					 time_points[i],dcss,w);
*/
					AddExpectedTransitionSufficientStatistics(jcdyn, i, dcss, w);
			}
			jointDSS = dcss;
		}
		dyn->AddSS(jointDSS,jcdyn,ss,w);
	}


	void FBInf::AddExpSuffStats(const RV * p0, SS * ss, double w) {
		const BN* bn = dynamic_cast<const BN *>(
		 the_initial_distribution_of_the_markov_process());
		RVComp *jbn = dynamic_cast<RVComp*>(bn->Joint(
				Instantiation(bn->Domain()+bn->CondDomain())));
		if(!jointP0SS) {
			RVCSCompSS* rvcscss = dynamic_cast<RVCSCompSS *>(jbn->BlankSS());
			EnsurePassesPerformedAndResultsCached();
			if (!alpha.empty()) {
				RVSimple * initial_prob_dist = alpha[0]->Clone();
				initial_prob_dist->MultBy(beta[0]);
				initial_prob_dist->Normalize();
/*
				jbn->AddExpSS(the_trajectory->Values(jbn->CondDomain(),
				 the_trajectory->TimeBegin()), init, rvcscss, w);
*/
				AddExpectedInitialSufficientStatistics(
						the_trajectory, initial_prob_dist, jbn, rvcscss, w);
				delete initial_prob_dist;
			}
			jointP0SS = rvcscss;
		}
		p0->AddSS(jointP0SS, jbn, ss, w);
	}


	double FBInf::CalcQuery(QueryCalculator &calc) {
		SS *ss = the_markov_process->BlankSS();
		//ClearCache();
		AddExpSuffStats(ss);
		double ret = calc.Calculate(ss);
		delete ss;
		return ret;
	}



	inline double ZeroProb(bool use_log_probs) {
		return use_log_probs ? -INFINITY : 0.0;
	}

} // end of ctbn namespace
