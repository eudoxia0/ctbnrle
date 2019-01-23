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
#include "nullptr03.h"
#include "params.h"

#include <algorithm>


namespace ctbn {

using namespace std;

CliqueTree::CliqueTree() :
	Head_compare(nullptr03),
	numrandrestarts(ParamInt("CliqueTreeRestarts",100)),
	Head(nullptr03),
	unit_test(false) {
}


CliqueTree::CliqueNode::CliqueNode() {
	prev = nullptr03;
}

CliqueTree::CliqueNode::~CliqueNode() {
}

void CliqueTree::CliqueNode::AddVar(const int &varid) {
	connected.insert(varid);
}


void CliqueTree::CliqueNode:: print(ostream &os) {
	set <int> ::iterator it;
	for(it = connected.begin(); it != connected.end(); it ++)
	{
		os << *it << " ";
	}
	os << endl;
}

void CliqueTree::CliqueNode::serial_preload() {
	prev = nullptr03;
}

void CliqueTree::CliqueNode::serial_postload() {
	for(unsigned int i=0;i<next.size();i++)
		next[i]->prev = this;
}

//destructor
CliqueTree::~CliqueTree() {
	Destroy(Head);
	Destroy(Head_compare);
}

void CliqueTree::BuildClique(const Markov &process, Random &rand) {
    vars = Context();
	data.clear();
	Destroy(Head_compare);
	Head_compare = nullptr03;

	AddFromBN(dynamic_cast<const BN *>(process.GetStartDist()));
	AddFromCTBNDyn(dynamic_cast<const CTBNDyn *>(process.GetDynamics()));
	
	BuildIt(rand);
}

void CliqueTree::BuildClique(const BN *bn, Random &rand) {
    vars = Context();
	data.clear();
	Destroy(Head_compare);
	Head_compare = nullptr03;

	AddFromBN(bn);
	BuildIt(rand);
}

void CliqueTree::BuildClique(const CTBNDyn *dyn, Random &rand) {
    vars = Context();
	data.clear();
	Destroy(Head_compare);
	Head_compare = nullptr03;

	AddFromCTBNDyn(dyn);
	BuildIt(rand);
}

void CliqueTree::BuildClique(const vector <Context> & Basis, Random & rand) {
    vars = Context();
    data.clear();
    Destroy(Head_compare);
    Head_compare = nullptr03;

    AddFromContext(Basis);
    BuildIt(rand);
}

void CliqueTree::AddFromContext(const vector <Context> & Basis) {
    for(unsigned int j = 0; j < Basis.size(); j ++) {
        CliqueNode newitem;
        vector<int> varids = Basis[j].VarList();
        if (varids.size()==0) continue;

        for(vector<int>::iterator i=varids.begin(); i!=varids.end();++i) {
            newitem.AddVar(*i);
            vars.AddVarCheck(*i,Basis[j].Cardinality(*i));
        }
        data.push_back(newitem);
    }
}

void CliqueTree::AddFromBN(const BN *bn) {
	int n = bn->NumofNodes();

	for(int i = 0; i < n; i++) {
		CliqueNode newitem;

		Context allvar = bn->Node(i)->Domain() + bn->Node(i)->CondDomain();
		vector<int> varids = allvar.VarList();
        if (varids.empty()) continue;

		for(vector<int>::iterator i=varids.begin(); i!=varids.end();++i) {
			newitem.AddVar(*i);
//			AddVarSize(*i,allvar.Cardinality(*i));
			vars.AddVarCheck(*i,allvar.Cardinality(*i));
		}
		data.push_back(newitem);
	}
}

void CliqueTree::AddFromCTBNDyn(const CTBNDyn *dyn) {
	int n = dyn->NumofNodes();
	for(int i = 0; i < n; i++) {
		CliqueNode newitem;

		Context allvar = dyn->Node(i)->Domain() + dyn->Node(i)->CondDomain();
		vector<int> varids = allvar.VarList();

		for(vector<int>::iterator i=varids.begin(); i!=varids.end();++i) {
			newitem.AddVar(*i);
//			AddVarSize(*i,allvar.Cardinality(*i));
			vars.AddVarCheck(*i,allvar.Cardinality(*i));
		}
		data.push_back(newitem);
	}
}

/*
int CliqueTree::VarSize(int varid) const {
	map<int,int>::const_iterator loc = context_size.find(varid);
	if (loc==context_size.end()) return 0;
	return loc->second;
}

void CliqueTree::AddVarSize(int varid, int size) {
	map<int,int>::iterator loc = context_size.lower_bound(varid);
	if (loc != context_size.end() && loc->first==varid) {
		if (loc->second==size) return;
		cerr << "inconsistent model input to CliqueTree: variable sizes do not agree" << endl;
		exit(1);
	}
	context_size.insert(loc,make_pair(varid,size));
}
*/

void CliqueTree::BuildIt(Random &rand) {
	save_data = data;
	double max_compare = 0;
	double maximum;
	//run the clique tree build on random arrangement as it set
	for(int i = 0; i < numrandrestarts; i++) {
		CliqueTrees(rand);
		maximum = MaxSize();
		//get the better clique tree
		if(maximum < max_compare || i == 0) {
			Destroy(Head_compare);
			Head_compare = Head;
			Head = nullptr03;
			max_compare = maximum;
		}
		else {
			Destroy(Head);
			Head = nullptr03;
		}
	}
}

//remove all the database;
void CliqueTree::Destroy(CliqueNode * place) {
	if(place) {
		for(unsigned int i = 0; i < place->next.size(); i++) {
			Destroy(place->next[i]);
		}
		delete place;
	}
}


//remove the item from the vector while generating
void CliqueTree::RemoveElem(vector <int> & random_clique, int remove) {
	if((random_clique.size()-1) != (unsigned)remove ) {
		random_clique[remove] = 
			random_clique[random_clique.size()-1];
	}
	random_clique.pop_back();
}

//building a clique_tree base on the node with same 
//context connected and combine into the tree
//eliminate the used node from the list
//run it til all the node in the data vector is gone or 
//all the number has been called
void CliqueTree::CliqueTrees(Random &rand) {
	data = save_data;
	clique_built.clear();
	int i = 0;

    vector <int> random_clique = vars.VarList();
/*
	vector <int> random_clique;
	for(map<int,int>::iterator i=context_size.begin();
			i!=context_size.end();++i)  {
		random_clique.push_back(i->first);
	}
*/

	while (data.size() > 0) {
		i =  unit_test ? 0 : rand.RandInt(random_clique.size());
		CliqueTrees(random_clique[i]);
		RemoveElem(random_clique, i);
	}
	Eliminate();
}

void CliqueTree::CliqueTrees(int search) {
	vector <int> node_exist;
	set <int> extend;
	set <int> extend_out;
	extend.insert(search);
	set <int> ::iterator it;

	for(unsigned int j = 0; j < data.size(); j++ ) {
		if(data[j].connected.find(search) != data[j].connected.end()) {
			node_exist.push_back(j);
		}
	}
	if(node_exist.size() > 0) {
		for(unsigned int j = 0; j < node_exist.size(); j ++) {
			set <int> ::iterator it;
			for(it = data[node_exist[j]].connected.begin();
					it != data[node_exist[j]].connected.end(); it ++) {
				if(extend.find(*it) == extend.end()) {
					extend.insert(*it);
					extend_out.insert(*it);
				}
			}
		}
		if(extend_out.size() != 0) {
			CliqueNode * newitem = new CliqueNode();
			newitem->connected = extend;
			Dequeue(node_exist, newitem);
			if(data.size() > 0) {
				CliqueNode newitem2;
				newitem2.connected = extend_out;
				newitem2.prev = newitem;
				data.push_back(newitem2);
			}
			clique_built.push_back(newitem);
		}
		else {
			CliqueNode * newitem = new CliqueNode();
			newitem->connected = extend;
			Dequeue(node_exist, newitem);
		}
	}
}

//set up the head and remove all the duplicate node within the tree
void CliqueTree::Eliminate() {
	if(clique_built.size() > 0) {
		CliqueNode * find_head = clique_built[0];
		while(find_head->prev != nullptr03) {
			find_head = find_head->prev;
		}
		Head = find_head;
		EliminateDups(Head);
		clique_built.clear();
		data = save_data;
	}
}

//eliminate the duplicate with in clique tree
// returns whether the previous node was eliminated
//  (because if it was, then the calling function needs to know!)
bool CliqueTree::EliminateDups(CliqueNode * checking) {
	bool ret = false;
	if(checking != Head) {
		if(checking->connected.size() >=
				checking->prev->connected.size()) {
			if(includes(checking->connected.begin(),
						checking->connected.end(),
						checking->prev->connected.begin(),
						checking->prev->connected.end())) {
				EliminatePrev(checking);
				ret = true;
			}
		}
	}
	for(unsigned int i = 0; i < checking->next.size(); i++) {
		// if checking was removed, then I cannot continue
		// (after removal, other "nexts" will have been checked
		//  by recursive call, so it is okay)
		if (EliminateDups(checking->next[i])) return ret;
	}
	return ret;
}

void CliqueTree::EliminatePrev(CliqueNode * current) {
	CliqueNode * delet = current->prev;
	if(delet == Head) {
		Head = current;
		for(unsigned int i = 0;i < delet->next.size(); i++) {
			if(delet->next[i] != current) {
				delet->next[i]->prev = current;
				current->next.push_back(delet->next[i]);
			}
		}
		current->prev = nullptr03;
		delete delet;
	}
	else {
		for(unsigned int i = 0;i < delet->next.size(); i++) {
			if(delet->next[i] != current) {
				current->next.push_back(delet->next[i]);
				delet->next[i]->prev = current;
			}
		}
		current->prev = delet->prev;
		//current->prev->next.push_back(current);
		for(unsigned int i = 0;i < current->prev->next.size(); i++) {
			if(current->prev->next[i] == delet) {
				current->prev->next[i] = current;
			}
		}
		delete delet;
	}
}

//remove the unnecessary node from the clique tree
void CliqueTree::Dequeue(vector <int> remove, CliqueNode * newitem) {
	for(int i = remove.size()-1; i > -1; i--) {
		if(data[remove[i]].prev != nullptr03) {
			newitem->next.push_back(data[remove[i]].prev);
			data[remove[i]].prev->prev = newitem;
		}
		if(&data[remove[i]].connected != 
		   &data[data.size()-1].connected) {
			data[remove[i]] = data[data.size()-1];
		}
		data.pop_back();
	}
}

//return the largest value to get the best tree
double CliqueTree::MaxSize() const {
	return MaxSize(Head);
}

//return the maximum size in the tree by looking into the each node
double CliqueTree::MaxSize(CliqueNode * place) const {
	double ret = 0;
	if(place) {
		int compare_size = 1;
		for(set <int>::iterator i = place->connected.begin();
				i != place->connected.end(); ++i) {
//			compare_size = compare_size * VarSize(*i);
			compare_size = compare_size * vars.Cardinality(*i);
		}
		ret = compare_size;
		for(unsigned int i = 0; i < place->next.size(); ++i) {
			double temp = MaxSize(place->next[i]);
			if (temp>ret) ret = temp;
		}
	}
	return ret;
}

//print the clique tree to stream os
void CliqueTree::CliquePrint(ostream &os) {
	CliqueNodePrint(Head_compare, 0, os);
}

void CliqueTree::CliqueNodePrint(CliqueNode * print, int size, ostream &os) {
	int a = size;
	if(print) {
		size ++;
		for(unsigned int i = 0; i < print->next.size(); ++i) {
			CliqueNodePrint(print->next[i], size, os);
		}
		while(a > 0) {
			os << "     ";
			a--;
		}
		if(print->prev != nullptr03) {
			prints(print->connected,os);
			os << "(";
			prints(print->prev->connected,os);
			os << ")" << endl << endl;
		}
		else {
			prints(print->connected,os);
			os << endl << endl;
		}
	}
}


void CliqueTree::prints(set <int> print_out, ostream &os) const {
	for(set <int>::iterator i = print_out.begin();
		i != print_out.end(); ++i) {
		os << *i << " ";
	}
}

//return the CliqueTree tree by CliqueTree::cliques;
vector<CliqueTree::Node> CliqueTree::ReturnCliques() const {
	vector<Node> ret;
	ReturnClique(Head_compare, 0, ret);
	return ret;
}


void CliqueTree::ReturnClique(CliqueNode * place,
		int parent, vector<CliqueTree::Node> &ret) const {
	if(place) {
		InsertCliqueNode(place, parent, ret);
		int myspot = ret.size();
		for(unsigned int i = 0; i < place->next.size(); ++i) {
			ReturnClique(place->next[i], myspot, ret);
		}
	}
}

void CliqueTree::InsertCliqueNode(CliqueNode * place, int parent,
			vector<CliqueTree::Node> &ret) const {
	Node new_node;
	new_node.vars = Clique2Context(place);
	if(parent>0) {
		new_node.adj.push_back(parent-1);
		ret[parent-1].adj.push_back(ret.size());
	}
	ret.push_back(new_node);
}

Context CliqueTree::Clique2Context(CliqueNode *node) const {
	Context ret;
	set<int>::const_iterator e = node->connected.end();
	for(set<int>::const_iterator i=node->connected.begin();i!=e;++i)
//		ret.AddVar(*i,VarSize(*i));
		ret.AddVar(*i,vars.Cardinality(*i));
	return ret;
}


void CliqueTree::SaveOld(std::ostream &os) const {
	set<int>::const_iterator it;
	for(unsigned int i = 0; i < data.size(); ++i) {
		for(it = data[i].connected.begin();
			it != data[i].connected.end(); it++) {
			os << *it << os.fill();
		}
		os << -1 <<os.fill();
	}
	os << -100 << os.fill();
	CliqueSave(os, Head_compare, 0);
}

void CliqueTree::CliqueSave(std::ostream &os,
		CliqueNode * print, int size) const {
	if(print) {
		set<int>::iterator it;
		os << size << os.fill();
		for(it = print->connected.begin();
			it != print->connected.end(); it++) {
			os << *it << os.fill();
		}
		os << -1 << os.fill();

		size ++;
		for(unsigned int i = 0; i < print->next.size(); ++i) {
			CliqueSave(os, print->next[i], size);
		}
	}
}

void CliqueTree::LoadOld(std::istream &is) {
	int get_out;
	data.clear();
	is >> get_out;
	while (get_out != -100) {
		CliqueNode newitem;
		while(get_out != -1 && get_out != -100) {
			newitem.AddVar(get_out);
			is >> get_out;
		}

		if(get_out != -100) {
			data.push_back(newitem);
			is >> get_out;
		}
	}
	CliqueNode * curr_pos = NULL;
	int get_before = 0;
	Destroy(Head_compare);
	int move_back;
	while(is >> get_out) {
		if(get_out != -1) {
			CliqueNode * newitem = new CliqueNode();
			if(get_out == 0) {
				get_before = get_out;
				while(is >> get_out) {
					if(get_out == -1) break;
					newitem->AddVar(get_out);
				}
				Head_compare = newitem;
				curr_pos = Head_compare;
			}
			else {
				if (curr_pos==NULL) exit(1); // file format error
				if(get_before >= get_out) {
					move_back = get_before - get_out + 1;
					for(int i = 0; i < move_back; i++) {
						curr_pos = curr_pos ->prev;
					}
				}
				get_before = get_out;
				while(is >> get_out) {
					if(get_out == -1) break;
					newitem->AddVar(get_out);
				}
				curr_pos->next.push_back(newitem);
				newitem ->prev = curr_pos;
				curr_pos = newitem;
			}
		}
	}
}

void CliqueTree::serial_preload() {
	Destroy(Head_compare);
	clique_built.clear();
	save_data.clear();
}

} // end of ctbn namespace
