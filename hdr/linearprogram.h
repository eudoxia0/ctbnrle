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
#ifndef CTBN_LINEARPROGRAM_H_
#define CTBN_LINEARPROGRAM_H_

#include <vector>
#include <glpk.h>

namespace ctbn {

// an abstract class (interface) for linear program solvers
// solves finding the *minimum* of the objective fn s.t. the constraints
class LinearProgram {
public:
	inline virtual ~LinearProgram() throw() {}

	// You do *not* need to call this before the first problem... only
	// between problems.
	virtual void ClearProblem() = 0;

	virtual void SetObjective(const std::vector<double> &coeff) = 0;
	// same as above, but only non-zero coefficients are mentioned
	// (the corresponding variable numbers are given in the parallel vector ix
	//  with 0 being the "first" variable)
	virtual void SetObjective(const std::vector<int> &ix,
			const std::vector<double> &coeff) = 0;

	// add the constraint that x'*coeff = val
	virtual void AddEqConstraint(const std::vector<double> &coeff,
			double val) = 0;
	// same as above, but coeff only stores the non-zero coefficients
	// their indicies are stored in ix
	// indicies in ix start at 0
	virtual void AddEqConstraint(const std::vector<int> &ix,
			const std::vector<double> &coeff, double val) = 0;

	// add the constraint that x'*coeff >= val
	virtual void AddGrEqConstraint(const std::vector<double> &coeff,
			double val) = 0;
	// same as above, but coeff only stores the non-zero coefficients
	// their indicies are stored in ix
	// indicies in ix start at 0
	virtual void AddGrEqConstraint(const std::vector<int> &ix,
			const std::vector<double> &coeff, double val) = 0;

	virtual std::vector <double> Solve(double &value) = 0;
};

// a subclass of LinearProgram that employs the GLPK libraries
class GLPKSolver: public LinearProgram {
public:
	GLPKSolver();
	GLPKSolver(const GLPKSolver &lp);
	GLPKSolver &operator=(const GLPKSolver &lp);
	virtual ~GLPKSolver() throw();

	virtual void ClearProblem();
	
	virtual void SetObjective(const std::vector<double> &coeff);
	virtual void SetObjective(const std::vector<int> &ix,
			const std::vector<double> &coeff);

	virtual void AddEqConstraint(const std::vector<double> &coeff,
			double val);
	virtual void AddEqConstraint(const std::vector<int> &ix,
			const std::vector<double> &coeff, double val);
	virtual void AddGrEqConstraint(const std::vector<double> &coeff,
			double val);
	virtual void AddGrEqConstraint(const std::vector<int> &ix,
			const std::vector<double> &coeff, double val);

	virtual std::vector <double> Solve(double &value);

protected:
	void AddGenRow(const std::vector<double> &coeff,
			int *ix,double val,int type);
	void ReassertVarBounds(int varid);
	void ResizeVars(int maxvarnum);
	void ResizeVars(const std::vector<int> &ix);

private:
	glp_prob *p;
	int nconstr;
	std::vector<double> mins,maxs;
};

}

#endif /* LINEARPROGRAM_H_ */
