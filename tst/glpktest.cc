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

#include "linearprogram.h"
#include "streamextra.h"

#include <iostream>

using namespace ctbn;
using namespace std;

namespace{
	int testing_result = 0;
}

#define EXPECT(cond) do { \
        if (!(cond)) {                                                  \
               std::cerr << "Test failure at " << __FILE__ << ":" << __LINE__ \
                      << " in " << __func__ \
                      << ", failed condition is " #cond << std::endl; \
            testing_result = 1; \
        } \
    } while (0)


int main(int argc, char **argv) {
	LinearProgram *lp = new GLPKSolver();

	vector<double> w(3);
	w[0] = -10; w[1] = -6; w[2] = -4;
	lp->SetObjective(w);
	
	w[0] = -1; w[1] = -1; w[2] = -1;
	lp->AddGrEqConstraint(w,-100);

	w[0] = -10; w[1] = -4; w[2] = -5;
	lp->AddGrEqConstraint(w,-600);

	w[0] = -2; w[1] = -2; w[2] = -6;
	lp->AddGrEqConstraint(w,-300);

	w.resize(1);
	w[0] = 2;

	vector<int> i(1);
	i[0] = 0;

	lp->AddGrEqConstraint(i,w,-30);
	lp->AddGrEqConstraint(i,w,-20);
	i[0] = 1;
	lp->AddGrEqConstraint(i,w,-20);
	lp->AddGrEqConstraint(i,w,-30);
	i[0] = 2;
	w[0] = 1;
	lp->AddEqConstraint(i,w,-10);


	double val;
	w = lp->Solve(val);

	EXPECT(fabs(val+760)<0.0001);
	EXPECT(w.size()==3);
	double expw[3] = {35,75,-10};
	EXPECT(fabs(expw[0]-w[0])<0.0001);
	EXPECT(fabs(expw[1]-w[1])<0.0001);
	EXPECT(fabs(expw[2]-w[2])<0.0001);

	/* output should be:
3 35 75 -10 
-760
	 */
	//cout << w << endl;
	//cout << val << endl;
	delete lp;
	if (testing_result == 0)
		cout << "PASS exactquery unit test." << endl;
	return testing_result!=0;
}
