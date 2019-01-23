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
#include "expmethod.h"
#include "markovdyn.h"
#include "rk.h"



namespace ctbn {

using namespace std;

ExpMethod::ExpMethod() : VarSample() {
}

ExpMethod::~ExpMethod() {
}

double ExpMethod::SampleTime(const matrix &Q, int id, int val,
		const Instantiation &currinst, double currt,
		const Trajectory *evid, double te, int e,
		Random &rand) const {
	double deltat;
	if (e==-1 || val==e) {
		deltat = rand.SampleExp(-Q[val][val]);
	} else {
		deltat = rand.SampleTruncateExp(-Q[val][val], te-currt);
	}
	return deltat;
}

int ExpMethod::SampleTransition(const matrix &Q, int id, int val,
		const Instantiation &inst, double currt,
		const Trajectory *evid, double te, int e,
		double &weight, Random &rand) const {
	int newval = -1;
	//vectr prob = Q.row(val);
	//double sum = -prob[val];
	// unnecessary -- SampleMultinormial ignores negative values
	//prob[val] = 0;
	//newval = rand.SampleMultinomial(prob, sum);
	double sum = -Q[val][val];
	newval = rand.SampleMultinomial(&Q[val][0], Q.getn(), sum);
	return newval;
}

double ExpMethod::TimeWeight(const matrix &Q, const Instantiation &currinst,
		int id, int transitionid, int val,
		double currt, double deltat,
		const Trajectory *evid, double te, int e) const {
	if (e==-1||val==e) return 0;
	double negq = Q[val][val];
	double weight = log(1-exp(negq * (te-currt)));
	if (id!=transitionid) weight -= log(1 - exp(negq*(te-currt-deltat)));
	if (!finite(weight)) exit(0);
	return weight;
}
/*
double ExpMethod::TransitionWeight(const matrix &Q, int transitionid, int val,
		const Instantiation &currinst,
		const Instantiation &nextinst,
		double currt, double deltat,
		const Trajectory *evid, double te, int e,
		bool logscale) const {
	double weight = 0;
	//   if (te-currt<1e-8&&Q.getm()>2)
	//     if (nextinst.Value(transitionid)==e)
	//         weight = Q[val][e];
	return weight;
}

*/

//********************************************************
LookAheadMethod::LookAheadMethod() : ExpMethod() {

}

int LookAheadMethod::SampleTransition(const matrix &Q, int id, int val,
		const Instantiation &inst, double currt,
		const Trajectory *evid, double te, int e,
		double &weight, Random &rand) const {
	if (e==-1) return ExpMethod::SampleTransition(Q,id,val,inst,
						currt,evid,te,e,weight,rand);
	transitionprob = 1;
	vectr prob = Q.row(val);
	double sum = -prob[val];
	prob[val] = 0;
	int newval = -1;
	int size = prob.length();
	if (size<3) {
		newval = rand.SampleMultinomial(prob, sum);
		transitionprob = 1.0;
	} else {
		vectr b(prob.length(), 0.0);
		b[e] = 1;
		expmtv(b, Q, te-currt);
		b = b.dotstar(Q.row(val));
		b[val] = 0.0;
		double bsum = b.sum();
		newval = randomizer.SampleMultinomial(b, bsum);
		transitionprob = prob[newval] * bsum / (b[newval] * sum);
	}
	weight += log(transitionprob);
	return newval;
}

/*
double LookAheadMethod::TransitionWeight(const matrix &Q, int transitionid, int val,
		const Instantiation &currinst,
		const Instantiation &nextinst,
		double currt, double deltat,
		const Trajectory *evid, double te, int e,
		bool logscale) const {
	//   vectr prob = Q.row(val);
	//   if (prob.length()<3) return 0;
	//   else
	return log(transitionprob);
}
*/

} // end of ctbn namespace
