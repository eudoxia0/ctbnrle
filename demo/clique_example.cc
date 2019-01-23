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
#include <fstream>
#include "clique.h"
#include "markov.h"

using namespace std;
using namespace ctbn;

/*
 * Pass in ctbn file
 * <char *> 
 */

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		cout << " you need to include a ctbn file" << endl;
		return 0;
	}

	//to initialize the CTBN
	ifstream ls(argv[1]);
	Markov ctbn(ls);

	//initializing the CliqueTree class
	CliqueTree check;

	//build the clique tree
	check.BuildClique(ctbn);

	cout << "Outputting the clique tree" << endl << endl;

	//print out the clique tree in shell
    	check.CliquePrint(cout);

   	//recieve the clique that been created by vector of <clique::cliques>
    	vector <CliqueTree::Node> test2 = check.ReturnCliques();
    	cout << "received item" << endl;

    	//this will output clique::cliques
    	//conn_context is whatthis contextes is connected to
    	//contextes is contexts
    	for(unsigned int i = 0; i < test2.size(); i++)
    	{
	    cout << "clique " << i << ": " << endl;
	    cout << "   variable ids(size): ";
	    vector<int> varlist = test2[i].vars.VarList();
	    for (unsigned int j=0;j<varlist.size();j++)
		    cout << varlist[j] << '(' << test2[i].vars.Cardinality(varlist[j]) << ") ";
	    cout << endl << "   adjacent cliques: ";
	    for(unsigned int j=0;j<test2[i].adj.size();j++)
		    cout << test2[i].adj[j] << ' ';
	    cout << endl;
    	}


	ls.close();
	//delete Dyn;


	//to save the clique tree that has been built
	ofstream os ("test.clique");
	check.Save(os);
	cout << "End" << endl;
	return 0;
}


