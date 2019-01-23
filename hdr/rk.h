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
#ifndef CTBNRLE_RK_H
#define CTBNRLE_RK_H
#include "matrix.h"
#include "defines.h"
#include "meanfieldinf.h"

namespace ctbn {

// eps is the desired tolerance of the answer
// (EPS below is the default value) which can be overridden
// by the parameter (see params.h) RKEPS
// (the default of -1.0 below forces the code to check for
//  the default parameter)
#define RKEPS 1e-8

// h is a "recommended starting step size"
// If h<0, then it uses 0.1/max(|diag(Q)|) as h

// computes a * exp(Qt), via integration with desired error of eps
// (h, if non-negative, gives a starting estimate of step size)
// a is replaced with the result
// the return value is the log of the factor by which a should
// be multiplied for the true result (the returned result is normalized
//  to sum to 1)
double vexpmt(vectr &a, const matrix &Q, double t,
		double eps=-1.0, double h=-1.0, vector<double> *timesteps=NULL);

// same as above, but calculates exp(Qt) * b
double expmtv(vectr &b, const matrix &Q, double t,
		double eps=-1.0, double h=-1.0);

// calculates  \int_u=0^t a*exp(Q*u)*e_i e_j*exp(Q*(t-u))*b du
//       [this is a matrix indexed by (i,j)]
// the transpose of this calculation is returned in c
//  (it is natural to calculate the transpose... the callee can
//   transpose it back)
// Again, the returned value is the log of the factor by which the
// result should be multiplied to find the "true result"
double suffstat(const vectr &a, const vectr &b, matrix &c,
              const matrix &Q, double t, double eps=-1.0, double h=-1.0);

// And here are two versions for 2-by-2 matrices.  They are
// often faster.  Not currently used, but here for reference
double vexpmt2(vectr &a, const matrix &Q, double t);
double expmtv2(vectr &b, const matrix &Q, double t);
double suffstat2(const vectr &a, const vectr &b, matrix &c,
		 const matrix &Q, double t); 

void mfBackward(MeanFieldInf *inf, int varid, 
	double t1, double t0, double eps=-1.0, double h=-1.0);
void mfForward(MeanFieldInf *inf, int varid, 
	double t0, double t1,
	double eps=-1.0, double h=-1.0);
double mfEnergy(MeanFieldInf *inf, int varid, 
	double t0, double t1, double eps=-1.0, double h=-1.0);
double mfSuffStat(MeanFieldInf *inf, const Instantiation &x,
	double t0, double t1, double eps=-1.0, double h=-1.0);
double mfTransSuffStat(MeanFieldInf *inf, int varid,
	const Instantiation &x1, const Instantiation &x2,
	double t0, double t1, double eps=-1.0, double h=-1.0);

} // end of ctbn namespace

#endif
