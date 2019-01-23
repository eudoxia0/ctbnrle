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
#include "params.h"
#include "bn.h"
#include "ctbndyn.h"
#include "em.h"
#include "utils.h"
#include "ensurectbn.h"

#include <ios>
#include <iostream>

using namespace std;
using namespace ctbn;

class myBN : public BN {
public:
	myBN(const BN &b) : BN(b) { }
	RV *J() { return Joint(Context()); }
};

int main(int argc, char **argv) {
	const char *input_file = "queryinput.data";
 
	ifstream fin(input_file);
	Markov m;
	try {
		m.Load(fin);
	} catch (const serial::streamexception &e) {
		cerr << e.what() << endl;
		exit(-1);
	}

	const CTBNDyn *ctbndyn = dynamic_cast<const CTBNDyn *>(m.GetDynamics());
	const BN *bn = dynamic_cast<const BN *>(m.GetStartDist());
	Trajectory evid;
	evid.Load(fin);

	myBN mybn(*bn);
	vectr p0;
	double lf;
	const MultiZSimple &dist = (dynamic_cast<MultiRV *>(mybn.J()))->operator[](0);
	dist.GetDist(p0,lf);
	p0 *= exp(lf);
	cout << p0 << endl;

	cout << ctbndyn->JointMatrix() << endl;

	int c = 0;
	for(Trajectory::Index i=evid.Begin(ctbndyn->Domain());!i.Done();++i) c++;
	cout << c << endl;
	Trajectory::Index i=evid.Begin(ctbndyn->Domain());
	int type =0;
	for(;!i.Done();type=i.TestInc(ctbndyn->Domain())) {
		cout << i.Time() << ' ' << (type>1) << ' ';
		vector<int> l;
		ctbndyn->Domain().ConsistentIndexes(l,i.Values());
		cout << l << endl;
	}
	cout << i.Time() << endl; // one last time at the end

	while(1) {
		int qtype;
		fin >> qtype;
		if (qtype==-1 || fin.eof()) break;
		int var, val;
		double t1;
		fin >> var >> val;
		double compans;
		Instantiation x(ctbndyn->Domain());
		x.SetVal(var,val);
		vector<int> l;
		ctbndyn->Domain().ConsistentIndexes(l,x);
		switch(qtype) {
			case 0: // filter distribution at point
				fin >> t1;
				cout << "0 " << l << t1 << endl;
			break;
			case 1: // smooth distribution at point
				fin >> t1;
				cout << "1 " << l << t1 << endl;
			break;
			case 2: // expected time query
				cout << "2 " << l << endl;
			break;
			case 3: // expected # trans query
				{ int val2; fin >> val2;
				vector<int> l2;
				Instantiation x2(x);
				x2.SetVal(var,val2);
				ctbndyn->Domain().ConsistentIndexes(l2,x2);
				cout << "3 " << l << l2 << endl;
				}
			break;
		}
		double ans;
		fin >> ans;
	}
	fin.close();
	return 0;
}
