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

//This code generates n partially observed trajectories from
//a given ctbn by sampling trajectories and randomly remove
//certain fraction of information from the trajectories. 
//
//Usage: gen_partialdata <.ctbnfile> <output file> 

/* Command line parameters:
 *
 * Set input model file
 * -DModelFile <string> (default: - => stdin)
 *
 * Set output trajectories file
 * -DOutputFile <string> (default: - => stdout)
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
 * Set the fraction of information removed in each iteration
 * -DRemoveFrac <double> (default: 0.1)
 *
 * Set the number of iterations 
 * -DRemoveIt <integer> (default 1)
 */

#include "markov.h"
#include "trajectory.h"
#include "params.h"
#include "utils.h"
#include <fstream>
#include "ensurectbn.h"

using namespace std;
using namespace ctbn;

int main(int argc, char**argv) {
	InitParams(argc, argv);
	double begintime = ParamDouble("BeginTime" , 0.0);
	double endtime = ParamDouble("EndTime", 5.0);
    int n = ParamInt("NumTraj", 1);
    double frac = ParamDouble("RemoveFrac", 0.1);
    int nit = ParamInt("RemoveIt", 1);
    if (endtime <= begintime) {
        cout << "invalid length of trajectory" << endl;
        exit(0);
    }
	//load CTBN model from the file
	string instr = ParamStr("ModelFile","-");
	Markov *model;
	if (instr != "-") {
		ifstream fin(instr.c_str());
		if (!fin.good()) {
			cout << "bad CTBN file" << endl;
			exit(1);
		}
		model = Markov::LoadPtr(fin);
	} else {
		model = Markov::LoadPtr(cin);
	}
	string outstr = ParamStr("OutputFile","-");
	ofstream outf;
	if (outstr != "-") {
	    outf.open(outstr.c_str());
	    if (!outf.good()) {
		   cout << "bad output file" << endl;
		   exit(1);
	    }
	}
	ostream &fdata = (outstr=="-" ? cout : outf);
	fdata.precision(10);
    fdata << n << endl;
	//generate partially observed data
    Context context = model->GetDynamics()->Domain();
    int numvar = context.NumVars();
	for(int i=0; i<n; i++) {
		Trajectory tr;
		tr.SetBeginTime(begintime);
		tr.SetEndTime(endtime);
		model->Sample(tr);
        //randomly remove some information of the trajectory
        RemoveNodesInformation(tr, context, numvar, nit*numvar, frac);
        tr.Save(fdata);
        fdata << endl;
	}

	return 0;
}
