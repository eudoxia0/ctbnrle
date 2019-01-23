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
#include <iostream>
#include "graphic.h"
#include "markov.h"
#include "params.h"
#include "ensurectbn.h"

#include "bn.h"
#include "ctbn.h"
#include "ensurectbn.h"

using namespace std;
using namespace ctbn;

int main(int argc, char **argv)
{

	//to initialize or use the params.h class
	InitParams(argc, argv);
	//to see if the file opens or not
	string File;
	ifstream ls("drug.ctbn");
	if(!ls.is_open())
	{
		cout << "drug.ctbn is not presented "
			 << "in current directory. " << endl;
		return 1;
	}

	ifstream output1;
	ifstream output2;
	//building the Markov to use in graphic program
	Markov ctbn;
	ctbn.Load(ls);

	//several different file to be written
	//dot extention for graphviz and eps extention for post-script
	ofstream file("output_test_ctbn.dot");

	//creating the Graphic class name test
	Graphic test;

	//setting the font for the trajectory graph
	test.SetFontTrajectory(16, 10);

	//drawing the CTBN graph with rank system intergrated
	test.DrawCTBN(ctbn, file);
	file.close();
	char compare1, compare2;
	output1.open("unit_test_ctbn.dot");
	output2.open("output_test_ctbn.dot");
	if(output1.is_open() and output2.is_open())
	{

		while(output1 >> compare1 and output2 >> compare2)
		{
			if(compare1 != compare2)
			{
				cout << "CTBN graphic it "
				     << "generated is wrong. " << endl;
				return 2;
			}
		}

		output1.close();
		output2.close();
	}

	ofstream testoutput1("output.dot");
	ofstream testoutput2("output2.eps");
	string filename = "output.dot";
	test.DrawCTBNRank(ctbn,testoutput1, filename, testoutput2);

	const CTBNDyn * B = dynamic_cast<const CTBNDyn *>(ctbn.GetDynamics());
	Context context = Context();
	for (int i = 0; i < B->NumofNodes(); ++i)
		context = context + B->Node(i)->Domain();

	//to create the Trajectory::Index
	Trajectory tr;
	ls.close();
	ls.open("unit_test_traj.graph");
	if(!ls.is_open())
	{
		cout << "unit_test_traj.graph is"
			 << "not presented in current directory. " << endl;
		return 3;
	}
	tr.Load(ls);

	ofstream file2("output_unit_test.eps");
	//to draw the Index into trajectory graph
	test.DrawTrajectory(tr, context, file2);
	file2.close();
	//end of the file
	output1.open("unit_test_traj.eps");
	output2.open("output_unit_test.eps");

	if(output1.is_open() and output2.is_open())
	{
		string buffer;
		while(output1 >> compare1 and output2 >> compare2)
		{
			if(compare1 != compare2)
			{
				cout << compare1 << " " << compare2 << endl;
				cout << "Trajectory graphic "
				     << "it generated is wrong. " << endl;
				return 4;
			}
			if(compare1 == '%')
			{
				getline(output1, buffer);
				getline(output2, buffer);
			}
		}
		output1.close();
		output2.close();
	}
	else
	{
		cout << "Error on trajectory graphic test."
		     << "\ncouldn't find all the necessary file." << endl;
		return 5;
	}
	cout << "PASS" << endl;
	return 0;
}
