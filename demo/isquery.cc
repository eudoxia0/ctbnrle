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
#include "markovdyn.h"
#include "streamextra.h"
#include "context.h"
#include "multirv.h"
#include "samplinginf.h"
#include "importancesampler.h"
#include "expmethod.h"
#include "params.h"
#include "bn.h"
#include "ctbndyn.h"
#include "em.h"
#include "utils.h"
#include <iostream>
#include "ensurectbn.h"

using namespace std;
using namespace ctbn;

// this file reads in the input model & evidence & true query results from "queryinput.data"
// it calculate the queries using impotrance sampling

// Usage:
// no inputs to run: ./importancequery >output_file
int main(int argc, char **argv)
{
	// load the input model, evidence and query answers in "queryinput.data"
	const char *input_file = "queryinput.data";
	addnan(cout);
	InitParams(argc, argv);
 
	ifstream fin(input_file);
	// load model
	Markov m;
	m.Load(fin);
	const CTBNDyn *ctbndyn = dynamic_cast<const CTBNDyn *>(m.GetDynamics());
	const BN *bn = dynamic_cast<const BN *>(m.GetStartDist());
	Context context = ctbndyn->Domain() + ctbndyn->CondDomain();

    //load evidence
	Trajectory evid;
	evid.Load(fin);
	// load true query answers
	double expect_q1, expect_q2;
	fin >> expect_q1 >> expect_q2;
	fin.close();

	// construct the queries
	Instantiation x(context, -1);
	int querynode = ParamInt("QueryNode", 1);
	int querytimeval = ParamInt("QueryTimeVal", 1);
	Context c = ctbndyn->Node(querynode)->Domain();
	x.SetVal(querynode, querytimeval);
	// query the amount of time node 1 spends on state 1
	QueryTime query1(x);
	// query the number of times node 1 transition from state 0 to state 1 
	QueryTransition query2(c, 0, 1); 

	// set up the approximate inference engine using importance sampler
	SamplingInf iinf;
	// set the number of samples
	int numofsamples = 1000;
	Sampler *sampler;
	ExpMethod method;
	sampler = new ImportanceSampler(&m, &evid, &method);
	iinf.SetSampler(sampler);
	iinf.SetProcess(&m);
	iinf.SetNumSample(numofsamples);

	// output the results
	cout << "the expected query answer is:" << endl;
	cout << expect_q1 <<  " " << expect_q2 << endl;
	cout << "the approximate inference answer using importance sampler to the query is:" << endl;
	cout << iinf.CalcQuery(query1) << " " << iinf.CalcQuery(query2) << endl;
	return 0;
}
