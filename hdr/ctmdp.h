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
#ifndef CTBN_CTMDP_H_
#define CTBN_CTMDP_H_

#include "clique.h"
#include "linearprogram.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>

namespace ctbn {

// class to store (and solve) a continuous-time Markov Decision process
class CTMDP {
public:
	// a mapping from an assignment to a subset of the variables
	// (those variables in var)
	// val is the value of the mapping, organized in the same sequence
	// that Instantiation will iterate over the assignments to the Context var
	class Factor {
	public:
		Context var;
		std::vector<double> val;
		double operator()(const Instantiation &i) const
			{ return val[var.Index(i)]; }
	};

	// Note that you really should set beta.  It defaults to
	// 1.0, but that's pretty arbitrary (beta=0.0 does not have
	// a solution for many problems)
	CTMDP();
	~CTMDP();

	// beta is the discount rate of the process in that
	// a reward achived at time t is worth exp(-beta*t) that of a reward
	// at time 0
	// (if unset, it defaults to 1.0)
	void SetBeta(const double &beta);
	// removes all actions, basis fns, and rewards that have been defined
	void Clear();

	// forms full exponental rate matrices and solves exactly
	// using lp solver passed in
	// returns value function for each state
	std::vector<double> SolveExact(LinearProgram *lpsolver) const;
	// forms full exponential rate matrices and solves using 
	// basis fns (below) to approximate value function
	// returns weights for each basis function
	std::vector<double> SolveApprox_exp(LinearProgram *lpsolver) const;
	// same as SolveApprox_exp, except without forming the full rate matrix
	// in particular, forms clique tree and adds extra variables for
	// each sepset value-action pair.  For small networks, this takes
	// longer than _exp.  For large networks, it takes less time.
	std::vector<double> SolveApprox(LinearProgram *lpsolver) const;

	
	// returns the optimal policy for each state (this is exponential
	// in number!)  isapprox is true iff sol came from SolveApprox
	// or SolveApprox_exp
	std::vector<int> Policy(const std::vector<double> &sol,
				bool isapprox=false) const;
	// this is the same, but for a given state x
	// (essentially policy above returns the vector of the output
	//  of this method for each value of x)
	int Policy(const Instantiation &x, const std::vector<double> &sol,
			bool isapprox=false) const;

	// adds an action that has the dynamics of dyn
	// and reward rates given by the sum of the factors of r
	// ---> this object owns the dyn pointer after this call <----
	void AddAction(const Dynamics *dyn, const std::vector<Factor> &r);

	// adds a basis function for the value function approximation
	// (not used by SolveExact)
	void AddBasis(const Factor &phi);

protected:

	class ActionInfo {
	public:
		const Dynamics *dyn; // dynamics (rate matrix specification)
		std::vector<Factor> rs; // reward factors (reward is their sum)
	};


	double beta;
	std::vector<Factor> basis;
	std::vector<ActionInfo> actions;

	// expands factor by dynamics
	Factor ExpandFactor(const Factor &h, const Dynamics *dyn) const;

	// write contraints for exact factored LP 
	void writeConstraints(const std::vector<CliqueTree::Node> &cliqueAct,
			const int &a, const int &ci, const int &parent,
			const std::vector<Factor> &gs,
			const std::vector<std::vector<int> > &gloc,
			const std::vector<std::vector<int> > &rloc,
			LinearProgram *lp) const;
	// subroutine of above for writing them for a clique-tree node
	void writeNodeConstraints(const std::vector<CliqueTree::Node> &clique,
			const int &a, const int &ci, const int &parent,
			const std::vector<Factor> &gs,
			const std::vector<std::vector<int> > &gloc,
			const std::vector<std::vector<int> > &rloc,
			LinearProgram *lp) const;
	// set weights for optimization fn for approx LP
	void SetApproxOptWeights(int totalsize, LinearProgram *lpsolver) const;

	// from here down should probably get moved out of here
	// (as they aren't needed between method calls)

	class IDDP{ // to store an action-node-connection-value pair
			// (which has to be mapped onto an LP variable ID
	public:
		int action, node1, node2, value;
		IDDP(int a, int n1, int n2, int v) {
			action=a; value=v;
			if (n1<n2) { node1=n1; node2=n2; }
			else { node1=n2; node2=n1; }
		}
		bool operator<(const IDDP &n) const {
			if (action<n.action) return true;
			if (action>n.action) return false;
			if (node1<n.node1) return true;
			if (node1>n.node1) return false;
			if (node2<n.node2) return true;
			if (node2>n.node2) return false;
			if (value<n.value) return true;
			return false;
		}
	};
	int getLPVar(int action, int node1, int node2, int value) const;
	void resetlpvarmap() const;

	mutable std::map<IDDP,int> lpvarmap;
	mutable int maxvar;
};

}

#endif /* LINEARPROGRAM_H_ */
