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
#ifndef CTBNRLE_VARSAMPLE_H
#define CTBNRLE_VARSAMPLE_H

#include "dynamics.h"
#include "context.h"
#include "trajectory.h"



namespace ctbn {

//VarSample is abstract class to define sample method
//for Markov process, including how to sample the transition
//time,  how to sample the next state and the corresponding
//sampling weights.
//*weights are all in log scale!*

class VarSample{
public:
	VarSample();
	virtual ~VarSample();
	//sample the next transition time for 
	virtual double SampleTime(const matrix &Q, int id, int val,
			const Instantiation &currinst, double currt, 
			const Trajectory *evid, double te, int e,
			Random &rand=randomizer) const = 0;
	//sample the next state
	//and update the corresponding weight contribution
	virtual int SampleTransition(const matrix &Q, int id, int val,
			const Instantiation &inst, double currt, 
			const Trajectory *evid, double te, int e,
			double &weight, Random &rand=randomizer) const  = 0;
	//calculate the weight contribution of sampling time
	virtual double TimeWeight(const matrix &Q, const Instantiation &currinst, 
			int id, int transitionid, int val, 
			double currt, double deltat, 
			const Trajectory *evid, double te, int e) const = 0;
	//calculate the weight contribution of sampling next state
/*	virtual double TransitionWeight(const matrix &Q, int transitionid, int val, 
			const Instantiation &currinst, 
			const Instantiation &nextinst, 
			double currt, double deltat, 
			const Trajectory *evid, double te, int e,
			bool logscale = 1) const = 0;
*/
};

} // end of ctbn namespace

#endif
