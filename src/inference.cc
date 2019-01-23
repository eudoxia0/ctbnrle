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
#include "inference.h"
#include "markov.h"
#include <vector>



namespace ctbn {

using namespace std;

QueryTime::QueryTime(const Instantiation &c) {
	condi = Instantiation(c.KnownVars());
	vector<int> varlist = condi.VarList();
	for (unsigned int i=0; i<varlist.size(); i++)
		condi.SetVal(varlist[i], c.Value(varlist[i]));
	queryindex = condi.Index();
	//   cout << "Querytime" << endl;
	//   cout << condi.VarList() << endl;
	//   cout << queryindex << endl;
}

double QueryTime::Calculate(const Trajectory &tr) {
	Trajectory::Index i = tr.Begin(condi);
	//tr.Draw(cout);
	double retval = 0.0;
	for (; !i.Done(); ++i) {
		if (queryindex ==i.Values().Index()) {
			//cout << "Add deltat : " << i.DeltaT() << endl;
			retval += i.DeltaT();
		}
	}
	return retval;
}

double QueryTime::Calculate(const SS *ss) {
	vector<int> varlist = condi.VarList();
	int id = varlist[0];
	return dynamic_cast<const MarkovSS *>(ss)->NodeSS(id, 
							 condi.Value(id), 
							 condi.Value(id), 
							 condi);
}

QueryTransition::QueryTransition(const Context &c, int from, int to) {
	queryfrom = from;
	queryto = to;
	querycontext = c;
}

double QueryTransition::Calculate(const SS *ss) {
	vector<int> varlist = querycontext.VarList();
	int id = varlist[0];
	// condi is not used here. Just a place-holder to call ss.NodeSS.
	Instantiation condi = Instantiation(querycontext);
	return dynamic_cast<const MarkovSS *>(ss)->NodeSS(id, 
							 queryfrom, 
							 queryto, 
							 condi);
}

double QueryTransition::Calculate(const Trajectory &tr) {
	Trajectory::Index i = tr.Begin(querycontext);
	double retval = 0.0;
	int previndex = i.Values().Index();
	++i;
	//tr.Draw(cout);
	for(; !i.Done(); ++i) {
		int currindex = i.Values().Index();
		//cout << currindex << ' ' << previndex << " retval=" 
		//     << retval << endl;
		if( (currindex==queryto && previndex==queryfrom) ||
				(currindex==queryto && queryfrom==-1) ||
				(queryto==-1 && queryfrom==previndex))
			retval += 1.0;
		previndex = currindex;
	}
	return retval;
}

QueryProb::QueryProb(const Instantiation &val, double time) : x(val.KnownVars()), t(time) {
	xind = x.Index(val);
}

double QueryProb::Calculate(const Trajectory &tr) {
	return tr.Values(x,t).Index()==xind ? 1.0 : 0.0;
}

double QueryProb::Calculate(const SS *ss) {
	assert(0); // not implementable -- see .h file
}

//SOBJCLASSDEF(Inference)
Inference::~Inference() {
}

} // end of ctbn namespace
