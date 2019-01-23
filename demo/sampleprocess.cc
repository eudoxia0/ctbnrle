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

// this file reads in a model file and outputs a trajectory sampled from it
// from time 0 to time as given by the first argument
#include "process.h"
#include <iostream>
#include "params.h"
#include "ensurectbn.h"

using namespace std;
using namespace ctbn;

int main(int argc, char **argv)
{
	InitParams(argc, argv);

	if (argc < 2) {
		cerr << "Not enough arguments! ./sampleprocess ending_time" << endl;
		exit(1);
	}
	addnan(cout);
	addnan(cin);
 
	double endt = atof(argv[1]);
	if (endt<0.0) {
		cerr << "time must be non-negative" << endl;
		exit(1);
	}

	Process *p = Process::LoadPtr(cin);
	if (p==NULL) {
		cerr << "stdin is not a saved process pointer" << endl;
		exit(1);
	}
	Trajectory tr;
	tr.SetBeginTime(0.0);
	tr.SetEndTime(endt);
	p->Sample(tr);

	tr.Save(cout);
	cout << endl;
	return 0;
}
