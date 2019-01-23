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

#include "clique.h"
#include "ensurectbn.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace ctbn;

const int unit_test_check [4][3] = {{1}, {0, 2, 3} , {1} , {1}};
const int unit_test_check1 [4][3] = {{0,1,2}, {1, 3, 4}, {4,5,6}, {4,7}};
const int unit_test_save [3][80] = {{0, 2, -1, 0, 1, -1, 1, 2,
				   -1, 3, -1, 1, 3, 4, -1, 5, -1, 4, 5, 6,
				   -1, 4, 7, -1, -100, 0, 0, 1, 2, -1, 1,
				    1, 3, 4, -1, 2, 4, 5, 6, -1, 2, 4, 7, -1},
				   {0, 2, -1, 0, 1, -1, 1, 2, -1, 3, -1, 1,
				    3, 4, -1, 5, -1, 4, 5, 6, -1, 4, 7, -1,
				    0, -1, 0, 1, -1, 2, -1, 2, 3, -1, 4, -1,
				    5, -1, 4, 6, -1, 7, -1, -100, 0, 1, 2, 3,
				   -1, 1, 1, 3, 4, -1, 2, 4, 5, 6, -1, 2, 4,
				    7, -1, 1, 0, 1, 2, -1 }};

bool checkingSAVE (ifstream &, const int &);

int main(int argc, char **argv)
{

	//to initialize the CTBNDyn
	ifstream ls("drug.ctbn");
	if(!ls.is_open())
	{
		cout << "error: can not find the drug.ctbn file";
		return 1;
	}

	Markov ctbn; 
	try {
		ctbn.Load(ls);
	} catch (const serial::streamexception &e) {
		cerr << e.what() << endl;
		exit(-1);
	}
	CliqueTree check;

	//set the random generated value to always be 0
	check.IgnoreRand(true);

	const Dynamics * d = ctbn.GetDynamics();
	const CTBNDyn * Dyn = dynamic_cast<const CTBNDyn *>(d);

	//build the clique tree
	check.BuildClique(Dyn);
	//recieve the clique that been created by
	//vector of <clique::cliques>
	vector <CliqueTree::Node> test2 = check.ReturnCliques();

	//use the receive clique to test the correctness

	if(test2.size() != 4)
	{
		cout << "error: wrong number of clique item" << endl;
		cout << test2.size() << endl;
		return 2;
	}
	set <int> ::iterator it;
	int count;
	for(unsigned int i = 0; i < test2.size(); i++)
	{
		vector<int> varlist = test2[i].vars.VarList();
		for (unsigned int j=0;j<varlist.size();j++)
		{
			if(varlist[j] != unit_test_check1[i][j])
			{
				cout << "error: wrong tree been generated" << endl;
				return 3;
			}
		}

		for(unsigned int j=0;j<test2[i].adj.size();j++)
		{
			if(test2[i].adj[j] != unit_test_check[i][j])
			{
				cout << "error: wrong tree been generated" << endl;
				return 4;
			}
		}
	}
	ls.close();

	ofstream os ("testing_drug.clique");
	check.Save(os);
	os.close();
	ls.open("testing_drug.clique");
	if(!checkingSAVE(ls, 0))
	{
		cout << "error: failed to created saved file from CTBNDyn" << endl;
		return 5;
	}
	ls.close();
	CliqueTree check2;
	ls.open("testing_drug.clique");
	check2.Load(ls);
	test2 = check.ReturnCliques();
	//this will output clique::cliques
	//conn_context is whatthis contextes is connected to
	//contextes is contexts
	if(test2.size() != 4)
	{
		cout << "error: fail to correctly generate "
		     << "tree from saved data" << endl;
		return 6;
	}

	for(unsigned int i = 0; i < test2.size(); i++)
	{
		vector<int> varlist = test2[i].vars.VarList();
		for (unsigned int j=0;j<varlist.size();j++)
		{
			if(varlist[j] != unit_test_check1[i][j])
			{
				cout << "error: wrong tree been generated" << endl;
				return 7;
			}
		}

		for(unsigned int j=0;j<test2[i].adj.size();j++)
		{
			if(test2[i].adj[j] != unit_test_check[i][j])
			{
				cout << "error: wrong tree been generated" << endl;
				return 8;
			}
		}
	}
		//build the clique tree using ctbn markov class
	check.BuildClique(ctbn);
	os.open("testing_drug.clique");
	check.Save(os);
	os.close();
	ls.open("testing_drug.clique");
	//test the saved data
	if(!checkingSAVE(ls, 1))
	{
		cout << "error: failed to created saved "
                     << "file from markov" << endl;
		return 9;
	}


	cout << "PASS the unit test" << endl;
	return 0;
}


bool checkingSAVE (ifstream & ls, const int & j)
{
	if(!ls.is_open())
		return false;
	int save_testing;
	int count_ = 0;
	while(ls >> save_testing)
	{
		if(save_testing != unit_test_save[j][count_])
		{
			cout << "error: wrong saved data been created" << endl;
			return false;
		}
		count_ ++;
	}
	return true;
}


