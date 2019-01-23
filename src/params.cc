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
#include "params.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <map>
#include <algorithm>

namespace ctbn {

using namespace std;

map<string,string> Params;

void ReadParams(istream &pf) {
	while(!pf.eof()) {
		string line;
		getline(pf,line);
		if (line.empty() || line[0] == ';' || line[0] == '#') continue;
		string::size_type nstart = line.find_first_not_of(" \t\n\r",0);

		if (nstart==string::npos) continue;
		string::size_type nend = line.find_first_of(" \t\n\r",nstart);
		if (nend==string::npos) continue;
		string::size_type valstart = line.find_first_not_of(" \t\n\r",nend);
		if (valstart==string::npos) continue;
		string::size_type valend = line.size();
		if (valend==string::npos) continue;
		Params[line.substr(nstart,nend-nstart)] = 
			line.substr(valstart,valend-valstart);
	}
}

void SetParam(const string &name, const string &val) {
	Params[name] = val;
}

void InitParams(int argc, char **argv) {

	ifstream pf("params");
	if (pf.good()) {
		ReadParams(pf);
		pf.close();
	}
	pf.open((string(argv[0]) + string(".params")).c_str());
	if (pf.good()) {
		ReadParams(pf);
		pf.close();
	}
	bool listem = false;
	for(int i=1;i<argc;i++) {
		if (argv[i][0] == '-' && argv[i][1] == 'D') {
			if (i+1>=argc) {
				cerr << "invalid parameters:  missing value for " << argv[i]+2 << endl;
				exit(-1);
			}
			Params[argv[i]+2] = argv[i+1];
			i++;
		} else if (strcmp(argv[i],"-params")==0) {
			pf.open(argv[i+1]);
			if (pf.good()) {
				ReadParams(pf);
				i++;
				pf.close();
			}
		} else if (strcmp(argv[i],"-listparams")==0) listem = true;
	}
	if (listem) {
		for(map<string,string>::iterator i = Params.begin();
				i != Params.end(); ++i) {
			cout << i->first << " => " << i->second << endl;
		}
		exit(0);
	}
}

string ParamStr(const string &name, const string &dflt) {
	map<string,string>::iterator i = Params.find(name);
	if (i==Params.end()) return dflt;
	return i->second;
}

int ParamInt(const string &name, int dflt) {
	map<string,string>::iterator i = Params.find(name);
	if (i==Params.end()) return dflt;
	return atoi(i->second.c_str());
}

double ParamDouble(const string &name, double dflt) {
	map<string,string>::iterator i = Params.find(name);
	if (i==Params.end()) return dflt;
	return atof(i->second.c_str());
}

} // end of ctbn namespace
