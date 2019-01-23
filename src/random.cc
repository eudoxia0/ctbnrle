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
#include "random.h"
#include "defines.h"

#include <algorithm>
#include <cassert>

namespace ctbn {

// Define a constant to be 2 * sqrt(2/e); used in SampleNormal()
const double TWICE_A = 1.7155277699214135;

using namespace std;

const double Random::Conv32Double = 1.0/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2/2;
const double Random::Conv64Double = Random::Conv32Double*Random::Conv32Double;

Random randomizer(GET_TIME % (1 << 24));

Random::Random (Random::uint32 in_seed) {
	static int init_seed = 5489;
	if (in_seed > 0)
		seed = in_seed;
	else
		seed = init_seed;
	init_seed++;
	Reset(seed);
}


void Random::Reset (Random::uint32 in_seed) {
	if (in_seed > 0) seed = in_seed;

	MT[0] = seed;
	for(int i=1;i<624;i++)
		MT[i] = (1812433253 * (MT[i-1] ^ (MT[i-1]>>30))+i) & 0xffffffff;
	loc = 0;
}

Random::uint32 Random::Rand32() {
	if (loc==0) generateMore();
	uint32 ret = MT[loc];
	ret ^= ret>11;
	ret ^= (ret<<7) & 0x9d2c5680;
	ret ^= (ret<<15) & 0xefc60000;
	ret ^= ret>>18;
	loc++;
	if (loc==624) loc=0;
	return ret&0xffffffff;
}

void Random::generateMore() {
	for(int i=0;i<624;i++) {
		int i1 = i==623 ? 0 : i+1;
		int i397 = i>226 ? i-227 : i+397;
		uint32 t = (MT[i]&0x80000000) + (MT[i1]&0x7fffffff);
		MT[i] = MT[i397] ^ (t>>1);
		if (t & 1)
			MT[i] = MT[i] ^ 0x9908b0df;
		MT[i] &= 0xffffffff;
	}
}

double Random::RandRealFast() {
	return Rand32()*Conv32Double;
}

double Random::RandReal() {
	return Rand32()*Conv32Double+Rand32()*Conv64Double;
}

int Random::RandInt(int range) {
	// for small range, this works fine
	return int(RandRealFast()*range);
}

// The algorithm is taken from "Non-Uniform Random Variate Generation" by
// Luc Devroye pages 194-199
double Random::SampleNormal (void) {
	double  u, v, x, sqr_x;

	while (true) {
		u = RandReal();
		if (u == 0.0)
			continue;
		v = (RandReal() - 0.5) * TWICE_A;
		x = v / u;
		sqr_x = x*x;
		if (sqr_x <= 6 - 8*u + 2*u*u) return x;
		if (sqr_x >= 2 / u - 2 * u) continue;
		if (sqr_x <= -4 * ::log(u)) return x;
	}
}

// The algorithm is taken from Cheng, Generating Beta Variates
// with Nonintegral Shape Parameters. Communications of the ACM
// 1978 Vol 21 No 4 p317 - p322
double Random::SampleBeta (double a, double b) {
	double  s, lambda, u, u1, u2, v, y;

	s = a+b;

	double minab = a<b ? a : b;
	if (minab <= 1.0) {
		lambda = minab;
	} else {
		lambda = sqrt((2*a*b-s)/(s-2));
	}

	u = a + lambda;

	while (true) {
		u1 = RandReal();
		u2 = RandReal();

		v = 1/lambda * ::log(u1/(1-u1));
		y = a * exp(v);

		if(s*::log(s/(b+y)) + u*v - ::log(4.0) >= ::log(u1*u1*u2)) {
			return (y/(y+b));
		}
	}
}

double Random::SampleGamma (double a, double b) {
	if( a < 1.0 ) return DblGammaLessThanOne(a)/b;
	else if( a > 1.0 ) return DblGammaGreaterThanOne(a)/b;
	return -::log(RandReal());
}


// Code adopted from Nir who adapted it from David Heckerman
//-----------------------------------------------------------
//	DblGammaGreaterThanOne(dblAlpha)
//
//	routine to generate a gamma random variable with unit scale and
//      alpha > 1
//	reference: Ripley, Stochastic Simulation, p.90
//	Chang and Feast, Appl.Stat. (28) p.290
//-----------------------------------------------------------
double Random::DblGammaGreaterThanOne(double dblAlpha) {
	double rgdbl[6];

	rgdbl[1] = dblAlpha - 1.0;
	rgdbl[2] = (dblAlpha - (1.0 / (6.0 * dblAlpha))) / rgdbl[1];
	rgdbl[3] = 2.0 / rgdbl[1];
	rgdbl[4] = rgdbl[3] + 2.0;
	rgdbl[5] = 1.0 / sqrt(dblAlpha);

	for (;;) {
		double  dblRand1;
		double  dblRand2;

		do {
			dblRand1 = RandReal();
			dblRand2 = RandReal();

			if (dblAlpha > 2.5)
				dblRand1 = dblRand2 + rgdbl[5] * (1.0 - 1.86 * dblRand1);

		} while (!(0.0 < dblRand1 && dblRand1 < 1.0));

		double dblTemp = rgdbl[2] * dblRand2 / dblRand1;

		if (rgdbl[3] * dblRand1 + dblTemp + 1.0 / dblTemp <= rgdbl[4] ||
		        rgdbl[3] * ::log(dblRand1) + dblTemp - ::log(dblTemp) < 1.0) {
			return dblTemp * rgdbl[1];
		}
	}

	assert(false);
	return 0.0;
}

/* routine to generate a gamma random variable with unit scale and
 * alpha < 1 reference: Ripley, Stochastic Simulation, p.88
 */
double Random::DblGammaLessThanOne(double dblAlpha) {
	double dblTemp;
	const double	dblexp = exp(1.0);

	for (;;) {
		double dblRand0 = RandReal();
		double dblRand1 = RandReal();

		if (dblRand0 <= (dblexp / (dblAlpha + dblexp))) {
			dblTemp = pow(((dblAlpha + dblexp) * dblRand0) /
			              dblexp, 1.0 / dblAlpha);
			if (dblRand1 <= exp(-1.0 * dblTemp)) return dblTemp;
		} else {
			dblTemp = -1.0 * ::log((dblAlpha + dblexp) * (1.0 - dblRand0) /
			                     (dblAlpha * dblexp));
			if (dblRand1 <= pow(dblTemp,dblAlpha - 1.0)) return dblTemp;
		}
	}

	assert(false);
	return 0.0;
}  /* DblGammaLessThanOne */

void Random::SampleDirichlet(double const *alpha, int n, double *theta) {
	double sum = 0;

	for(int i =0; i < n; i++)
		sum += (theta[i] = SampleGamma(alpha[i], 1));
	// normalize
	for(int i =0; i < n; i++)
		theta[i] /= sum;
}

double Random::SampleExp (double a) {
	// where a = intensity = 1/mean
	///// ///// Need to fix so it never samples a 0?
  if(isfinite(1/a))
	return (-::log(1.0 - RandReal()) / a);
     else
       return INFINITY;
}

double Random::SampleTruncateExp(double q, double t) {
	double u1 = 1 - exp(-q * t);
	double u = RandReal() * u1;
	double ret = (-::log(1.0 - u)/q);
	return ret;
}

int Random::SampleAddWeights (const vector<double>& additive_weights) {
	double val = RandReal() * additive_weights[additive_weights.size()-1];
	return lower_bound(additive_weights.begin(),additive_weights.end(),val)
				- additive_weights.begin();
}

int Random::SampleMultinomial(const double *weights, int n, double sum) {
	double j = RandReal()*sum;

	for(int i=0;i<n;i++) {
		if (weights[i]>0) { // added so that an intensity row can be used
						// directly
			j -= weights[i];
			if (j<=0) return i;
		}
	}
	return n-1;
}

int Random::SampleMultinomial(const vectr &weights, double sum) {
	double j= RandReal()*sum;
	int n = weights.length();
	for(int i=0;i<n;i++) {
		if (weights[i]>0) { // added so that an intensity row can be used
						// directly
			j -= weights[i];
			if (j<=0) return i;
		}
	}
	return weights.length()-1;
}

int Random::SampleMultinomial(const vector<double> &weights, double sum) {
	double j= RandReal()*sum;
	int n = weights.size();
	for(int i=0;i<n;i++) {
		if (weights[i]>0) { // added so that an intensity row can be used
						// directly
			j -= weights[i];
			if (j<=0) return i;
		}
	}
	return weights.size()-1;
}

//sample n samples from multinomial
void Random::SampleMultinomial(const vector<double> &weights, vector<int> &samples, int n) {
	int length = weights.size();
	double * cdf = new double [length];
	cdf[0] = weights[0];
	for(int i=1; i<length; i++) cdf[i] = cdf[i-1] + weights[i];
	double pos = RandReal() * cdf[length-1] / n;
	double step = cdf[length-1] / n;
	int idx = 0 ;
	for(int i=0; i<n; pos += step, i++)
	{
		while(pos>cdf[idx]) idx++;
		samples.push_back(idx);
	}
	delete [] cdf;
}

} // end of ctbn namespace
