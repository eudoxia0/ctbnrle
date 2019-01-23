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
#ifndef CTBNRLE_GIBBSSAMPLER_H
#define CTBNRLE_GIBBSSAMPLER_H

#include <map>
#include <vector>
#include "sampler.h"
#include "gibbsbase.h"



namespace ctbn {

class matrix;
class MultiSimple;
class RV;
class MarkovDyn;

// GibbsSampler is an implementation of the paper:
// Tal El-Hay, Nir Friedman, and Raz Kupferman (2008). "Gibbs Sampling in
// Factorized Continuous-Time Markov Processes."
// This class will take in an evidence trajectory, a markov process, 
// the number of samples needed as parameters.
// It will first burn in for several iterations (number of burnin_iter is
// specified in construction function), then generate a set of sampled 
// trajectories. It is an approximate inference method used for learning 
// and inference on CTBNs.

// Sample usage:
//
// Markov ctbn;
// Trajectory evidence;
// GibbsSampler sampler(&ctbn, &evidence);
// Trajectory sample_traj;
// sampler.SampleTrajectory(sample_traj);

class GibbsSampler : public GibbsBase {
public:
	// the caller of the function owns the process and evidence pointers
	 // but they must remain valid for the life of the GibbsSampler object
	GibbsSampler(const Process *process, const Trajectory *evidence, int burnin);
	virtual ~GibbsSampler();

	// inherited interface functions from base class
	virtual Sampler *Clone() const { return new GibbsSampler(p, evid, numBurninIter); }

protected:
	virtual void SampleVariableInterval(int v, int x0, int xT, double t0, double tT,
					Random &rand = randomizer) const;
private:
	// equivalent to Qss;
	matrix GetRMatrix(int v, const Instantiation &inst) const;
	// search the exact time points between t0 and t1 using binary search;
	double SearchPoint(double ksi, int x0, const matrix &R, 
				double log_p_future0, double log_p_past, 
				MultiSimple p_future, double log_z, 
				double t0, double t1, 
				double epsilon = 1e-6) const;
	// sample the initial bn conditioned on evindence "partial"
	Instantiation SampleConditionalRv(const Instantiation &partial, double t, int v,
					Random &rand = randomizer) const; 
};

} // end of ctbn namespace

#endif
