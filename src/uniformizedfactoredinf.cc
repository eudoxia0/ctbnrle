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
#include "uniformizedfactoredinf.h"

#include <algorithm>



namespace ctbn {

using namespace std;

UniformizedFactoredInf::UniformizedFactoredInf() : FBInf() {
	lval = 10;
	thval = 10;
}

UniformizedFactoredInf::~UniformizedFactoredInf() throw() {}


UniformizedFactoredInf * UniformizedFactoredInf::Clone() const {
    UniformizedFactoredInf * cloneInf = new UniformizedFactoredInf(*this);
	CloneTo(cloneInf);
	cloneInf->lval = lval;
	cloneInf->thval = thval;
	return cloneInf;
}


FactoredMatrix * UniformizedFactoredInf::IntervalPropagator(Dynamics * const the_markov_process_dynamics, double start_time_point, double past_end_time_point, Instantiation const & current_evidence) const {
    FactoredMatrix * m = new FactoredMatrix(the_markov_process_dynamics);
	m->SetL(lval);
	m->SetTheta(thval);
	m->Cond(start_time_point, past_end_time_point, current_evidence);
	return m;
}


FactoredMatrix * UniformizedFactoredInf::PointTransitionPropagator(Dynamics * const the_markov_process_dynamics, double time_point, const Instantiation & prior_evidence, const Instantiation & current_evidence) const {
    FactoredMatrix * m = new FactoredMatrix(the_markov_process_dynamics);
	m->SetL(lval);
	m->SetTheta(thval);
	m->Cond(time_point, prior_evidence, current_evidence, true);
	return m;
}


FactoredMatrix * UniformizedFactoredInf::PointChangeEvidencePropagator(Dynamics * const the_markov_process_dynamics, double time_point, const Instantiation & prior_evidence, const Instantiation & current_evidence) const {
    FactoredMatrix * m = new FactoredMatrix(the_markov_process_dynamics);
    m->SetL(lval);
    m->SetTheta(thval);
    m->Cond(time_point, prior_evidence, current_evidence, false);
    return m;
}

void UniformizedFactoredInf::Restrict(RVSimple * a, Instantiation const & variable_assignment) const {
	vector<int> cind;
	CTBNDyn *cdyn = dynamic_cast<CTBNDyn*>(the_markov_process_dynamics());
	for (int i=0; i<cdyn->NumofNodes(); i++){
		cind.clear();		
		(cdyn->Node(i)->Domain()).ConsistentIndexes(cind, variable_assignment);
		//dynamic_cast<FactoredVector*>(a)->Restrict(i,cind);
		if(cind.size() < (unsigned int)cdyn->Node(i)->Domain().Size()) {
			for(unsigned int j=0; j<cind.size(); j++){
				dynamic_cast<FactoredVector*>(a)->RestrictById(i, cind[j]);
			}
		}
	}

}

//FactoredVector * UniformizedFactoredInf::MakeSimple(Instantiation const & x, bool normalize) {
RVSimple * UniformizedFactoredInf::MakeSimple(RV * prior, Instantiation const & x, bool normalize) const {
//  FactoredVector * ret =  new FactoredVector(the_initial_distribution_of_the_markov_process());
    FactoredVector * ret =  new FactoredVector(prior);
    FactoredVector * & ref_to_ptr_to_approximation_vector(ret);
	Restrict(ref_to_ptr_to_approximation_vector, x);
	if (normalize){
		ret->Normalize();
	}
	return ret;
}


RVSimple * UniformizedFactoredInf::Convert(RVCondSimple *&rvcond, RVSimple *&rv){
	return rv->Clone();
}


void UniformizedFactoredInf::Restrict(FactoredVector * & rv, Instantiation const & x) const {
	CTBNDyn * ctbndyn = dynamic_cast<CTBNDyn*>(the_markov_process_dynamics());
	int n = ctbndyn->NumofNodes();
	for (int i=0; i<n; i++){
		vector<int> ind;
		ctbndyn->NodeByVar(i)->Domain().ConsistentIndexes(ind, x);
		if (ind.size() == 1){
			rv->RestrictById(i, ind.at(0));
		}
	}
	//dynamic_cast<FactoredVector *>(rv)->Print(cout, nullptr03);
}


bool UniformizedFactoredInf::IsAlphaElementValid(RVSimple * const alpha_sub_x) {
    FactoredVector * as_factored_vector(dynamic_cast<FactoredVector *>(alpha_sub_x));
    // TODO: consider if there is a tighter validity condition to be checked here.
    bool ret = (nullptr03 != as_factored_vector);
    return ret;
}


} // end of ctbn namespace
