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
#ifndef CTBNRLE_MEANFIELDINF_H_
#define CTBNRLE_MEANFIELDINF_H_

#include "inference.h"
#include "contfunction.h"
#include "markov.h"

#include <vector>


namespace ctbn {

// MeanFieldInf is an implementation of the paper:
// Ido Cohn, Tal El-Hay, Nir Friedman, and Raz Kupferman (2009). "Mean Field
// Variational Approximation for Continuous-Time Bayesian Networks."
// This class takes in an evidence trajectory and a markov process. 
// It assumes that the process is defined with a BN as its initial 
// distribution and a CTBN for its dynamics. This also performs the mean 
// field method for the underlying BN after each backward integration before 
// integrating forward.

// Sample usage:
//
// Markov ctbn;
// Trajectory evidence;
// MeanFieldInf inf;
// inf.SetProcess(&ctbn);
// inf.SetTrajectory(&evidence);
// inf.SetEPS(desiredEPS); // (optional, default is 1e-5)

// After the process and trajectory are set, the various inference functions
// (Smooth, AddExpSuffStats, CalcQuery) may be called. For the first time any
// one of them are called, the optimization procedure is also called and the 
// results are cached. Subsequent calls will then used these cached results
// unless the process or evidence has changed.

class MeanFieldInf : public Inference {
public:
	MeanFieldInf();
	virtual ~MeanFieldInf();

	virtual MeanFieldInf *Clone() const;
	
	// See inference.h for ownership of p and tr in the following 
	// functions

	virtual void SetProcess(const Process *p);
	virtual void SetTrajectory(const Trajectory *tr);
	
	virtual double Filter(const Instantiation &x, double t,
				bool log=false);
	virtual double Smooth(const Instantiation &x, double t,
				bool log=false);

	virtual double Prob(double t, bool log=false);
	
	virtual void AddExpSuffStats(SS *ss, double w=1.0);
	virtual void AddExpSuffStats(const Dynamics *dyn, SS *ss,
					double w=1.0);
	virtual void AddExpSuffStats(const RV *p0, SS *ss, 
					double w=1.0);

	double CalcQuery(QueryCalculator &calc);

	// Sets the EPS for Runge-Kutta calls. Calls Reset().
	virtual void SetEPS(double neweps);

	// Computes the gradient of the rho at t and the current rho value
	// Result is returned in rhoGrad
	virtual void CalcRhoGrad(int varid, double t, const vectr &rhoVal, 
		vectr &rhoGrad);

	// Computes the gradient of the mu at t and the current mu value
	// Result is returned in muGrad
	virtual void CalcMuGrad(int varid, double t, const matrix &gammaVal,
		vectr &muGrad);

	// Computes gamma at t based on the current mu and rho values
	// Result is returned in gamma
	virtual void CalcGamma(int varid, double t, 
				const vectr &muVal, const vectr &rhoVal,
				matrix &gamma);

	// Computes the integrand of the suffstats integral at t for 
	// transitions based on the gammas and mus
	virtual double PointTransSuffStat(int varid, 
			const Instantiation &x1, const Instantiation &x2,
			double t);

	// Calculates the integrand of the energy for variable varid at 
	// time t
	virtual double CalcPointCompEnergy(int varid, double t);


	// These allow the external RK function to access and 
	// update the temp functions as it integrates
	virtual void GetRhoTempVal(double t, vectr &rhoVal);
	virtual void GetMuTempVal(double t, vectr &muVal);
	virtual void GetGammaTempVal(double t, matrix &gammaVal);

	virtual void AddRhoTempVal(double t, const vectr &rhoVal);
	virtual void AddMuTempVal(double t, const vectr &muVal);
	virtual void AddGammaTempVal(double t, const matrix &gammaVal);

	
	// Print out the marginals and joints
	void PrintMu() const;
	void PrintGamma() const;

	// Print out the times that where the state is hidden
	void PrintHidden() const;

	int GetNumStates(int varid) const {
		const CTBNDyn *cdyn =
			dynamic_cast<const CTBNDyn* >(p->GetDynamics());
		int ret = cdyn->NodeByVar(varid)->Domain().Size();
		return ret;
	}

private:
	// Create a factored CTBN by randomly choosing one of
	// the intensity matrices for each variable
	virtual void TransToFactored();

	// Performs forward-backward passes on each of the variables using
	// the factored process created by TransToFactored() to initialize 
	// the mus and gammas
	virtual void Initialize();

	// Performs the backward-forward ODE passes until convergence 
	// (or specified KMAX steps[default is 30])
	virtual void Optimize();

	// Calculates the expected Q matrix based on the current states
	// of the system.
	// Specifying jvarid and xj makes this Calculate the conditional
	// expected rates
	virtual void CalcExpectQ(int varid, double t, matrix &expectQ, 
					int jvarid = -1, int xj = -1,
					bool log = true);

	virtual void CalcPsi(int varid, double t, vectr &psi);


	// Calculates the energy for the specified varid
	double CalcCompEnergy(int varid);

	// Calculates the energy for the entire system
	double CalcEnergy();

	// Resets the inference object to a state before any 
	// initialization/optimization
	void Reset();

	// Returns the state pair containing the states at the 
	// instants before and after time t for a variable
	std::pair<int,int> GetTrans(int varid, double t) const;

	// Given a variable, returns all the observed transition times of 
	// its children mapped to the varid of its children, respectively
	std::map<double,int> GetChildTransTimes(int varid,
					double start, double end) const;

	// Returns all of the observed transition time of a variable
	std::vector<double> GetObservedTransTimes(int varid, 
					double start, double end) const;

	// Perform mean-field for BNs at the start of the trajectory
	void CalcMuStart(int varid, const vectr &rho0, vectr &mu);
	
	// Calculates the probability of the system being in instantiation x
	// based on the marginals
	double GetProb(const Instantiation &x, double t,
				bool log=false);

	// Add a variable to maps storing the marginals, joints, and hidden
	// times
	void AddVar(int varid);

	// Adds values to the maps storing the marginals, joints, and hidden
	// times
	void AddMuVal(int varid, double t, const vectr &vals);
	void AddGammaVal(int varid, double t, const matrix &vals);
	void AddHiddenInterval(int varid, double start, double end);


	// Not owned by this class
	const Markov *p;
	const Trajectory *tr;

	// Owned by this class
	Markov *factored;

	// The solved ODEs are stored here
	std::map<int,ContFunction<vectr> > mus;
	std::map<int,ContFunction<matrix> > gammas;

	// Determines where backward-forward integration is needed
	typedef std::pair<double,double> double_pair;
	std::map<int,std::vector<double_pair> > hiddenIntervals;

	// Cache the transition times, used to keep the energy calculation
	// more stable by not integrating over a transition in the system,
	// where the function is not differentiable
	std::map<int, std::vector<double> > transTimes;

	// Used during the integration process
	ContFunction<vectr> rhoTemp;
	ContFunction<vectr> muTemp;
	ContFunction<matrix> gammaTemp;

	// Flags to mark whether initialize/optimize have been executed
	bool initialized;
	bool optimized;

	double eps;
};

} // end of ctbn namespace

#endif

