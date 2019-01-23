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

//This code shows how to learn the paramters of a CTBN 
//given fully observed data. Data are n trajectories
//saved in datafile. They can be generated using 
//gen_fulldata

#include "markov.h"
#include "trajectory.h"
#include "params.h"
#include <fstream>
#include "ensurectbn.h"

using namespace std;
using namespace ctbn;

int main(int argc, char**argv) {
	InitParams(argc, argv);

    if (argc<3) {
        cout << "Usage: ./learnparams <.ctbn file> <datafile>" << endl;
        exit(0);
    }
	ifstream fin(argv[1]);
	if (!fin.good()) {
		cout << "bad CTBN file" << endl;
		exit(1);
	}
	//load CTBN model from the file
	Markov ctbn;
	ctbn.Load(fin);
    ifstream fdata(argv[2]);
    if (!fdata.good()) {
        cout << "bad datafile" << endl;
        exit(1);
    }
    int n;
    fdata >> n; 
	//load data
    vector<Trajectory> data;
	for(int i=0; i<n; i++) {
		Trajectory tr;
        tr.Load(fdata);
        data.push_back(tr);
	}
 
	vector<double> w(n, 1.0);
	//Calculate the sufficient statistics of the data
	SS *ss = ctbn.SuffStats(data, w);
	//randomly choose parameters
	ctbn.Scramble();
	//learn the parameters from complete data
	ctbn.Maximize(ss);
	ctbn.Save(cout);
	if(ss!=NULL) delete ss;

	return 0;
}
