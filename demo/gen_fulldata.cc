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

//This code generates n fully observed trajectories from
//a given ctbn. 
//
//Usage: gen_fulldata <.ctbnfile> <output file> 

/* Command line parameters:
 *
 * Set the begin time of the trajectory
 * -DBeginTime <double> (default: 0.0)
 *
 * Set the end time of the trajectory
 * -DEndTime <double> (default: 5.0)
 *
 * Set the number of trajectories
 * -DNumTraj <integer> (default: 1)
 *
 */
#include "markov.h"
#include "trajectory.h"
#include "params.h"
#include <fstream>
#include "ensurectbn.h"

using namespace std;
using namespace ctbn;

int main(int argc, char**argv) {
    if (argc<3) {
        cout << "Usage: ./gen_fulldata <.ctbn file> <output file>" << endl;
        exit(0);
    }
	InitParams(argc, argv);
	double begintime = ParamDouble("BeginTime" , 0.0);
	double endtime = ParamDouble("EndTime", 5.0);
    if (endtime <= begintime) {
        cout << "invalid length of trajectory" << endl;
        exit(0);
    }
    int n = ParamInt("NumTraj", 1);
	ifstream fin(argv[1]);
	if (!fin.good()) {
		cout << "bad CTBN file" << endl;
		exit(1);
	}
	//load CTBN model from the file
	Markov ctbn;
	ctbn.Load(fin);
    ofstream fdata(argv[2]);
    if (!fdata.good()) {
        cout << "bad output file" << endl;
        exit(1);
    }
    fdata << n << endl;
	//generate data
	for(int i=0; i<n; i++) {
		Trajectory tr;
		tr.SetBeginTime(begintime);
		tr.SetEndTime(endtime);
		ctbn.Sample(tr);
        tr.Save(fdata);
        fdata << endl;
	}

	return 0;
}
