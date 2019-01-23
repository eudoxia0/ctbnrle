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
#include "brutestructuresearch.h"
#include "extramath.h"
#include <vector>


namespace ctbn {

using namespace std;

BruteStructureSearch::BruteStructureSearch(const Context &c, FamScore* fs,
								int mp) 
 : StructureSearch(c, fs), maxParents(mp) {
}

Structure BruteStructureSearch::LearnStructure() {
        int numNodes = vars.size();
        Context nullc;

        vector <vector<int> > exhaustIdx = combo(numNodes-1,maxParents);
        
	unsigned int numTests = exhaustIdx.size();

        for(int i=0; i<numNodes; i++) {
                vector<Context> parentSet;
                parentSet.insert(parentSet.end(),
				 vars.begin(),vars.begin() + i);
                parentSet.insert(parentSet.end(),
				 vars.begin() + (i+1),vars.end());

                double score;
                double bestScore = 1;

		Node* bestnode = NULL;

                for(unsigned int j=0; j<numTests; j++) {
                        Context cvTest;
                        for(unsigned int k=0;k<exhaustIdx[j].size();k++)
                                cvTest = cvTest + parentSet[exhaustIdx[j][k]];
			
			score = fScore->GetScore(vars[i],cvTest);
                       	if(bestScore == 1 || score > bestScore) {
				if(bestScore != 1 && bestnode!=NULL) delete bestnode;
                               	bestScore = score;
                                bestnode = new Node(vars[i],cvTest);
			}
		}
		nodes.push_back(bestnode->Clone()); 
		vector<int> tmp = bestnode->CondDomain().VarList();
		for (unsigned int k=0; k < tmp.size(); k++)
			cout << tmp[k] << " ";
		cout << endl;
		delete bestnode;
        }
        Structure s;
	GetStructure(s);
        return s;
}

} // end of ctbn namespace
