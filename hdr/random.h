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
#ifndef CTBNRLE_RANDOM_H
#define CTBNRLE_RANDOM_H

#include <stdlib.h>
#include <vector>
#include "matrix.h"


namespace ctbn {

// Class for generating (pseduo-)random numbers.  Keeps the state of the
// generator so that different instantiations can be used for different
// parts of the code and be predictable (their pseudo-random draws don't
// affect each other)
// Class also includes method for drawing from non uniform distributions
class Random {
	public:
		typedef unsigned long uint32; // just needs to be at least 32 bits

		Random (uint32=-1);

		void Reset(uint32 in_seed=-1);

		// returns [0.0, 1.0)
		double  RandReal();
		// same, but only uses 32 bits of randomness (instead of 64)
		double  RandRealFast();

		// returns [0, range-1]
		int RandInt (int range);

		// returns ~N(0,1)
		double  SampleNormal (void);

		// returns ~Beta(a,b) where a and b are positive reals
		double SampleBeta (double a, double b);

		// returns ~Gamma(a,b) where a and b are positive reals
		double SampleGamma (double a, double b);

		/// samples from the Dirichlet joint distribution over thetas
		// theta  ~ theta[0]^(alpha[0]-1) * ... * theta[n]^(alpha[n]-1)
		void SampleDirichlet(double const *alpha, int n, double *theta);

		// samples from Exponential Dist for intensity a = 1/mean
		double SampleExp (double a);

		//sample from Exponential dist for intensity q, given value is
		// less than t
		double SampleTruncateExp(double q, double t);

		// samples from a set of additive weights (for multi)
		int SampleAddWeights (const std::vector<double>& additive_weights);
	
        // samples from a weighted distribution
		int SampleMultinomial (const std::vector<double>& weights,
				double sum=1.0);
		int SampleMultinomial (const double *weights, int n, double sum=1.0);
		int SampleMultinomial (const vectr &weights, double sum=1.0);
		//sample n samples a time
		void SampleMultinomial(const std::vector<double> &weights,
				std::vector<int> &samples, int n);

	private:
		double DblGammaLessThanOne(double dblaphla);
		double DblGammaGreaterThanOne(double dblaphla);

		uint32 Rand32();
		void generateMore();

	private:
		// starting seed number
		uint32 seed;
		// state of the generator:
		uint32 MT[624];
		int loc;

		static const double Conv32Double, Conv64Double;

		/*
	SERIAL_START(Random)
		SERIAL_VAR(uint32,seed)
		TODO: add fixed-length array to serialization code
		SERIAL_VAR_ARRAY(uint32,MT,624)
		SERIAL_VAR(int,loc)
	SERIAL_END
		*/
};


extern Random randomizer;

} // end of ctbn namespace

#endif
