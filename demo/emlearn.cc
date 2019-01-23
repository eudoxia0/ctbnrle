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
//given partially observed data using EM algorithm. Data 
//are n trajectories saved in datafile. They can be generated 
//using gen_partialdata. Exact inference method is used in EM.

#include "markov.h"
#include "trajectory.h"
#include "params.h"
#include "em.h"
#include "exactmarkovinf.h"
#include <fstream>
#include "ensurectbn.h"

using namespace std;
using namespace ctbn;

int main(int argc, char**argv) {
	InitParams(argc, argv);
	//load CTBN model from the file
	Markov *model;
	string modelstr = ParamStr("ModelFile","-");
	if (modelstr != "-") {
		ifstream fin(modelstr.c_str());
		if (!fin.good()) {
			cerr << "bad CTBN file" << endl;
			exit(1);
		}
		model = Markov::LoadPtr(fin);
	} else {
		model = Markov::LoadPtr(cin);
	}
	string datastr = ParamStr("InputData","-");
	ifstream dataf;
	if (datastr != "-") {
	    dataf.open(datastr.c_str());
	    if (!dataf.good()) {
		   cerr << "bad datafile" << endl;
		   exit(1);
	    }
	}
	istream &fdata = datastr =="-" ? cin : dataf;
	
    int n;
    fdata >> n; 
	//load data
    vector<Trajectory> data;
	for(int i=0; i<n; i++) {
		Trajectory tr;
        tr.Load(fdata);
        data.push_back(tr);
	}
    Inference *inf = new ExactMarkovInf;
    //start from randomly selected parameters
	if (ParamInt("Scramble",0)) model->Scramble();
    EM(data, model, inf);
	string outstr = ParamStr("OutputFile","-");
	if (outstr != "-") {
		ofstream outf(outstr.c_str());
		if (!outf.good()) {
			cerr << "bad datafile" << endl;
			exit(1);
		}
		model->Save(outf);
	} else {
	    model->Save(cout);
	}
    delete inf;
	return 0;
}
