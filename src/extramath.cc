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
#include "extramath.h"

#include <algorithm>
#include <cmath>
#include <vector>


namespace ctbn {

using namespace std;

double lngamma(double x) {
	double y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	                      24.01409824083091,-1.231739572450155,
	                      0.1208650973866179e-2,-0.5395239384953e-5};
	y=x;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
						
	for (int j=0;j<6;j++) ser += cof[j]/++y;
							
	return -tmp+log(2.5066282746310005*ser/x);
}

vector <vector <int> > combo(int nItems, int nChoose) {
	vector <vector <int> > ret(0);
	ret.push_back(vector <int> (0));

	vector <int> vec;
	for(int i=0;i<nItems;i++)
		vec.push_back(i);
	if(nChoose >= 1) {
		vector <vector <int> > temp = ret;

		for(int i=0; i<nChoose; i++) {
			vector<vector <int> > current(0);
			for(unsigned int j=0; j<temp.size();j++) {
				vector<int> element = temp[j];
				int end = element.size()-1;
				for(int k=end<0 ? 0 : element[end]+1;k<nItems;k++) {
					vector<int> temp2 = element;
					temp2.push_back(k);
					current.push_back(temp2);
				}
			}
			temp = current;

			ret.insert(ret.end(),current.begin(),current.end());
		}
	}
	return ret;
}

int listDiff(const vector<int> &listA, const vector<int> &listB) {
	int ret = 0;
	for( vector<int>::const_iterator it = listB.begin();
			it != listB.end(); it++) {
		if(find(listA.begin(),listA.end(),*it) != listA.end())
			ret++;
	}
	return listA.size() + listB.size() - (2*ret);
}

} // end of ctbn namespace
