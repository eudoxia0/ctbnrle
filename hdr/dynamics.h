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
#ifndef CTBNRLE_DYNAMICS_H
#define CTBNRLE_DYNAMICS_H

#include "context.h"
#include "trajectory.h"
#include "rvsimple.h"
#include "rv.h"

namespace ctbn {

// Dynamics is a continuous-time semi-Markov process except for the
// starting distribution
// The domain is explicit
class Dynamics : public StreamObj {
public:
	// var is the set of variables (state space of the process)
	// cvar is the set of variables on which the process is
	// conditioned.  Let cvar be a null context (just create a
	// context without any arguments to the constructor) if this
	// process is not dependent on any others.
	Dynamics(const Context &var=Context(), const Context &cvar=Context());
	// load from file..
	Dynamics(std::istream &is);
	// destructor
	virtual ~Dynamics();

	// load and save...
	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	// virtual copy constructor
	virtual Dynamics *Clone() const = 0;

	// Returns the state space
	const Context &Domain() const { return v; }
	// Returns the conditioning variables
	const Context &CondDomain() const {return cv; }

	// Cause the dynamics to be "closed" (no probability of leaving
	// the system)
	virtual void Normalize() = 0;

	// Add evidence that the system stays within the states
	// defined in x
	virtual void Restrict(const Instantiation &x) = 0;

	// the version here just sets up the contexts (v and cv)
	virtual void Mult(const Dynamics *x);

	// The next two return a conditional random variable that
	// describes the distribution at one time point conditioned on 
	// a previous time point.  They are returned as RVCondSimple
	// because the domains are implied

	// the implied domain for RVCondSimple is the Context
	// of this dynamic process
	// x gives the values for the conditional context
	// as well as any subsystem restrictions for the process's context

	virtual RVCondSimple *Cond(double t0, double t1,
			const Instantiation &x) const = 0;

	// same, but instantaneous intensity of transitioning
	// if transition==false, this doesn't correspond to a
	// transition, but rather just a change of support
	virtual RVCondSimple *Cond(double t, const Instantiation &from,
		const Instantiation &to, bool transition=true) const = 0;

	// Return a blank sufficient statistics object suitable for
	// this dynamics
	virtual SS *BlankSS() const = 0;

	// Adds sufficient statistics for being in state x from t0
	// to t0+deltat by adjusting ss.  Addition can be weighted by w
	virtual void AddSS(const Instantiation &x, double t0,
			double deltat, SS *ss, double w=1.0) const = 0;
	// Same, but adds ss for a transition from x1 to x2 at time t
	virtual void AddTransSS(const Instantiation &x1, 
				const Instantiation &x2,
				double t, SS *ss, double w=1.0) const = 0;

	// Adds expected sufficient statistics for the interval t0 to
	// to+deltat if alpha is the distribution at t0 given evidence
	// before t0 and beta is the distribution at t0+deltat given
	// evidence after t0+deltat.
	virtual void AddExpSS(const RV *alpha, const RV *beta,
				double t0, double deltat, 
				SS *ss, double w=1.0) const = 0;
	// Same as above, but for a transition at t.  x1 is the
	// distribution just before the transition (conditioned on evidence
	// before t) and x2 is the same just after t (conditioned on
	// evidence after t)
	// changevar is the variable that changed value
	virtual void AddExpTransSS(const RV *x1, const RV *x2,
			const Context &changevar,
			double t, SS *ss, double w=1.0) const = 0;

	// These are the same as above, but where the context is implied
	// and the conditioning variables (cvar) takes on a known value (not
	// a distribution over them) -- either its index or a full 
	// instantiation is supplied...
	virtual void AddExpSS(const Instantiation &condi,
			const RVSimple *alpha, const RVSimple *beta,
			double t0, double deltat, SS *ss, 
			double w=1.0) const = 0;
	virtual void AddExpSS(int condi, const RVSimple * const alpha,
	 const RVSimple * const beta, double t0, double deltat, SS * const ss,
	 double w = 1.0) const = 0;
	virtual void AddExpTransSS(const Instantiation &condi,
			const RVSimple *x1, const RVSimple *x2,
			const Context &changevar,
			double t, SS *ss, double w=1.0) const = 0;
	virtual void AddExpTransSS(int condi,
			const RVSimple *x1, const RVSimple *x2,
			const Context &changevar,
			double t, SS *ss, double w=1.0) const = 0;

	// This adds the sufficient statistics of a different process
	// to this process's sufficient statistics
	virtual void AddSS(const SS *toadd, const Dynamics *dyn,
			SS *ss, double w=1.0) const = 0;

	// i is the assignment at t.  nextt becomes the next time when
	// a transition occurs.  nexti is the next assignment
	// assumes that the context of nexti equals Domain
	virtual void SampleNextEvent(const  Instantiation &i, double t,
			double &nextt, Instantiation &nexti,
			Random &rand=randomizer) const = 0;
	// completes the trajectory, starting at t and going until 
	// tr.EndTime()
	// the method implemented here repeatedly calls SampleNextEvent
	// override it if you can do something more efficient
	virtual void SampleTrajectory(Trajectory &tr, double t,
			Random &rand=randomizer) const;

	// perform ML by adjusting the parameters, using the sufficient
	// statistics in ss (which you may assume are of the correct type)
	virtual void Maximize(const SS *ss) = 0;

	// Randomize the parameters
	virtual void Scramble(double a=1.0, double b=1.0, 
			double alpha=1.0, double degree=1.0, 
			Random &rand=randomizer) = 0;
	// return the log-likelihood from the sufficient statistics ss
	virtual double LLH(const SS *ss) const = 0;

	// returnn the sufficient statistics from the set of trajectories
	// tr, each weighted by the weight w (w is *not* in log space)
	SS* SuffStats(const std::vector<Trajectory> &tr, 
			const std::vector<double> &w) const;

	// returns a structure search score based on ss
	// numTrans and amtTime are the Bayesian prior parameters
	// for a BDe prior
	virtual double GetScore(double numTrans, double amtTime, 
				const SS *ss) const = 0;
protected:
	Context v,cv;

	SERIAL_START_V(Dynamics)
		SERIAL_VAR(Context,v)
		SERIAL_VAR(Context,cv)
	SERIAL_END

};

} // end of ctbn namespace

#endif
