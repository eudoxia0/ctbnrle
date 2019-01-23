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
#include "utils.h"
#include "ensurectbn.h"

using namespace std;
using namespace ctbn;

/*
 * This is code for generating the illustration of CTBN and trajectory
 * It takes in the MARKOV file (.ctbn) to extract the information.
 *
 */

/* Command line parameters:
 *
 * Pass in CTBN file
 * -DFile <char *> (default: drug.ctbn)
 *
 * Set the window x coordinate size in trajectory 
 * -Dtraj_win_size <integer> (default: 800)
 *
 *
 * Set the font on trajectory graph 
 * -Dtraj_font <integer> (default: 15)
 *
 *
 * Set the font on ctbn graph 
 * -Dctbn_font <integer> (default: 16)
 *
 *
 * Set font on trajectory graph ruler 
 * -Dtraj_ruler <integer> (default: 5)
 */

int main(int argc, char **argv)
{

	//to initialize or use the params.h class
	InitParams(argc, argv);
	//to see if the file opens or not
	string File;
	//Read in the -DFile 
	File = ParamStr("File", "drug.ctbn");
	ifstream ls(File.c_str());
	if(!ls.is_open()) return 0;

	//building the Markov to use in graphic program
	Markov ctbn;
	ctbn.Load(ls);

	//several different file to be written
	//dot extention for graphviz and eps extention for post-script
	//do not use temp.dot as the file name
	ofstream file("TrajRank.eps");
	ofstream file2("CTBNRank.dot");

	ofstream file3("TrajWName.eps");

	ofstream file4("TrajWoName.eps");

	ofstream file5("CTBNMatrix.dot");

	//creating the Graphic class name test
	Graphic test;

	//setting the font for the trajectory graph
	test.SetFontTrajectory(16, 10);

	//to change the size of page in X use the below function
	//test.SetSizeTrajectory(const int &);

	//setting the font size for CTBN graphic
	test.SetFontCTBN(20);

	//string_ is the vector of name which will replace the number system
	//automatically assigned when graphic.cc has envoked.
	vector <string> string_;
	string_.push_back("Eating");
	string_.push_back("Full Stomach");
	string_.push_back("Hungry");
	string_.push_back("Uptake");
	string_.push_back("Concentration");
	string_.push_back("Barometer");
	string_.push_back("Joint pain");
	string_.push_back("Drowsy");
	//changing the name of the CTBN and trajectory graphics
	test.Changes(string_);

	// use DrawBoth for specific trajectory need to be draw 
	// DrawBoth take in
	// const Markov &, Trajectory & 
	// std::ofstream &,std::string &, std::ofstream &
	
	//drawing the CTBN graph with rank system intergrated
	//it need to know what the file2 name is
	string filename = "CTBNRank.dot";
	test.DrawCTBNRank(ctbn,file2, filename, file);
	
	//use the DrawCTBN if CTBN without matrix is needed
	//drawing the CTBN with matrixs for each node
	test.DrawCTBNMatrix(ctbn, file5);


	/*
	/ Generate the trajectory using the CTBNDyn
	*/
	//to create the Trajectory::Index
	const CTBNDyn * B = dynamic_cast<const CTBNDyn *>(ctbn.GetDynamics());
	Context context = Context();
	for (int i = 0; i < B->NumofNodes(); ++i)
		context = context + B->Node(i)->Domain();
	Trajectory tr;
	//generating the trajectory starting from 0 to 3 second 
	tr.SetBeginTime(0);
	tr.SetEndTime(3);
	ctbn.Sample(tr);
	//just save the trajectory which is created for later use
	ofstream os("using_graphic.traj");
	//remove the part of the Information frmo the trajectory
	for(int i = 0; i < B->NumofNodes(); i ++)
	{
		RemoveInformation(tr,context, i, 1.0, 1.5);
	}
	tr.Save(os);


	//to draw the Index into trajectory graph
	test.DrawTrajectory(tr, context ,file3);
	//clear the changed name
	test.ClearName();
	//draw the trajectory with name unchanged
	test.DrawTrajectory(tr,context ,file4);
	
	//example on DrawBoth
	//test.DrawBoth(ctbn,tr, file2, filename, file);

	//end of the file
	cout << "done" << endl;
	file.close();
	file2.close();
	file3.close();
	file4.close();
	file5.close();
	return 0;
}
