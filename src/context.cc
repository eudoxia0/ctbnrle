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
#include "context.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>

namespace ctbn {

using namespace std;

ostream& operator<< (ostream& os, const Instantiation& inst) {
	return(inst.SaveOld(os));
}

istream& operator>> (istream& is, Instantiation& inst) {
	return(inst.LoadOld(is));
}

ostream& operator<< (ostream& os, const Context& con) {
	return(con.SaveOld(os));
}

istream& operator>> (istream& is, Context& con) {
	return(con.LoadOld(is));
}

Context::Context() {
	size = 1;
}

Context::Context(Context const & other): card(other.card), size(other.size) {}

Context::Context(istream &is) {
	LoadOld(is);
}

Context::~Context() {
}

Context::Context(const Context &c1, const Context &c2, setop op) {
	switch(op) {
		case UNION: set_union(c1.card.begin(),c1.card.end(),
						  c2.card.begin(),c2.card.end(),
						  inserter(card,card.begin()));
			break;
		case INTERSECTION: set_intersection(c1.card.begin(),c1.card.end(),
						  c2.card.begin(),c2.card.end(),
						  inserter(card,card.begin()));
			break;
		case DIFFERENCE: set_difference(c1.card.begin(),c1.card.end(),
						  c2.card.begin(),c2.card.end(),
						  inserter(card,card.begin()));
			break;
	}
	size = 1;
	for(map<int,int>::iterator i=card.begin();i!=card.end();++i)
		size *= i->second;
}

Context::Context(const Instantiation &i1, const Instantiation &i2) {
	set_intersection(i1.card.begin(),i1.card.end(),
			i2.card.begin(),i2.card.end(),
			inserter(card,card.begin()));

	size = 1;
	for(map<int,int>::iterator i=card.begin();i!=card.end();) {
		int v1 = i1.Value(i->first);
		int v2 = i2.Value(i->first);
		if (v1 == -1 || v2==-1 || v1==v2)
			card.erase(i++);
		else {
			size *= i->second;
			++i;
		}
	}
}

vector<int> Context::VarList() const {
	vector<int> ret;
	for(map<int,int>::const_iterator i=card.begin();i!=card.end();++i)
		ret.push_back(i->first);
	return ret;
}

int Context::Index(const Instantiation &ind) const {
	int ret = 0;
	int inc = 1;
	for(map<int,int>::const_iterator i=card.begin();i!=card.end();++i) {
		int v = ind.Value(i->first);
		if (v==-1) return -1;
		ret += inc*v;
		inc *= i->second;
	}
	return ret;
}

Instantiation Context::Index(int i) const {
	Instantiation ret(*this);
	ret.SetIndex(i);
	return ret;
}

istream &Context::LoadOld(istream &is) {
        return is >> size >> card;
}

void Context::serial_postload() {
	size = 1;
	for(map<int,int>::iterator i=card.begin();i!=card.end();++i)
		size *= i->second;
}

ostream &Context::SaveOld(ostream &os) const {
	return os << size << os.fill() << card;
}

void Context::ConsistentIndexes(vector<int> &l, const Instantiation &v) const {
	l.clear();

	vector<int> maxind,incind;
	int iinc = 1;
	int index = 0;
	for(map<int,int>::const_iterator j=card.begin();j!=card.end();++j) {
		int vv = v.Value(j->first);
		if (vv == -1) {
			maxind.push_back(j->second);
			incind.push_back(iinc);
		} else {
			index += vv*iinc;
		}
		iinc *= j->second;
	}
	CConsistentIndexes(l,maxind,incind,index);
}

void Context::CConsistentIndexes(vector<int> &l, const vector<int> &maxind, const vector<int> &incind, int startindex) const {
	int n = maxind.size();
	vector<int> v(n,0);

	int i=0;
	do {
		l.push_back(startindex);
		for(i=0;i<n;i++) {
			v[i]++;
			if (v[i]==maxind[i]) {
				v[i] = 0;
				startindex -= incind[i]*(maxind[i]-1);
			} else {
				startindex += incind[i];
				break;
			}
		}
	} while(i<n);
}

void Context::InconsistentIndexes(vector<int> &l, const Instantiation &v) const {
	l.clear();
	vector<int> lneg;
	ConsistentIndexes(lneg,v);
	CInconsistentIndexes(l,lneg);
}

void Context::CInconsistentIndexes(vector<int> &l, const vector<int> &lneg) const {
	int s = Size();
	vector<int>::const_iterator negi = lneg.begin(),negend = lneg.end();
	int i = 0;
	while(negi!=negend && i<s) {
		if (*negi == i) { ++negi; ++i; }
		else l.push_back(i++);
	}
	while(i<s) l.push_back(i++);
}

bool Context::IsSubset(const Context &c) const {
	for(map<int,int>::const_iterator i=c.card.begin();i!=c.card.end();++i)
		if (card.find(i->first)==card.end()) return false;
	return true;
}

bool Context::IsOverlap(const Context &c) const {
	for(map<int,int>::const_iterator i=c.card.begin();i!=c.card.end();++i)
		if (card.find(i->first)!=card.end()) return true;
	return false;
}

// ----------------------------------------------------

Instantiation::Instantiation() : Context() {
	index = 0;
	numunknown = 0;
}

Instantiation::Instantiation(const Context &c, int defaultval) :
			Context(c) {
	index = 0;
	numunknown = defaultval == -1 ? c.NumVars() : 0;
	int accum = 1;
	for(map<int,int>::iterator i=card.begin();i!=card.end();++i) {
		vals.insert(make_pair(i->first,defaultval));
		indexinc.insert(make_pair(i->first,accum));
		if (defaultval > 0) index += accum*defaultval;
		accum *= i->second;
	}
}

Instantiation::Instantiation(const Context &c, const Instantiation &v) :
			Context(c) {
	index = 0;
	numunknown = 0;
	int accum = 1;
	for(map<int,int>::iterator i=card.begin();i!=card.end();++i) {
		int vv = v.Value(i->first);
		vals.insert(make_pair(i->first,vv));
		indexinc.insert(make_pair(i->first,accum));
		if (vv!=-1) index += accum*vv;
		else numunknown++;
		accum *= i->second;
	}
}

Instantiation::~Instantiation() {
}

bool Instantiation::Inc(int varid) {
	map<int,int>::iterator i=vals.find(varid);
	map<int,int>::iterator j=card.find(varid);
	map<int,int>::iterator k=indexinc.find(varid);
	if (i!=vals.end() && j!=card.end() && k!=indexinc.end()
			&& i->second != -1) {
		i->second++;
		index += k->second;
		if (i->second==j->second) {
			index -= k->second*j->second;
			i->second = 0;
			return false;
		} else {
			return true;
		}
	}
	return false; // I guess...
}

bool Instantiation::Inc() {
	index++;
	for(map<int,int>::iterator i=vals.begin(),j=card.begin();
			i!=vals.end();++i,++j) {
		i->second++;
		if (i->second == j->second) i->second = 0;
		else return true;
	}
	index = 0;
	return false;
}

bool Instantiation::Inc(const Context &c) {
	for(map<int,int>::iterator
					i=vals.begin(),j=card.begin(),k=indexinc.begin();
			i!=vals.end();++i,++j,++k) {
		if (!c.HasId(i->first)) continue;
		i->second++;
		index += k->second;
		if (i->second == j->second) {
			index -= k->second*j->second;
			i->second = 0;
		}
		else return true;
	}
	return false;
}

void Instantiation::AddVar(int varid, int cardinality, int val) {
	if (HasId(varid)) SetVal(varid,val);
	else {
		map<int,int>::iterator ii = indexinc.lower_bound(varid);
		if (ii==indexinc.end()) { // if no other varids after me...
			indexinc.insert(make_pair(varid,size));
			index += size*(val==-1 ? 0 : val);
		} else { // otherwise we need to adjust the index values
			int iinc = ii->second;
			if (val!=-1) index += iinc*val;
			// for all (lexographically) later indexs,
			//  adjust the index accordinly and adjust the indexinc
			map<int,int>::iterator vi = vals.find(ii->first);
			for(;ii != indexinc.end(); ++ii,++vi) {
				index += (cardinality-1)*ii->second*vi->second;
				ii->second *= cardinality;
			}
			indexinc.insert(make_pair(varid,iinc));
		}
		Context::AddVar(varid,cardinality);
		vals.insert(make_pair(varid,val));
	}
}

bool Instantiation::operator==(const Instantiation &i) const {
	for(iimap::const_iterator vi = vals.begin();vi!=vals.end();++vi)
		if (vi->second != i.Value(vi->first)) return false;
	for(iimap::const_iterator vi = i.vals.begin();vi!=i.vals.end();++vi)
		if (vi->second != Value(vi->first)) return false;
	return true;
}

bool Instantiation::operator!=(const Instantiation &i) const {
	return !(*this==i);
}

ostream& Instantiation::SaveOld(ostream &os) const {
	return Context::SaveOld(os) << os.fill() << vals;
}

istream& Instantiation::LoadOld(istream &is) {
	Context::LoadOld(is) >> vals;
	Recalc();
	return is;
}

void Instantiation::serial_postload() {
	Recalc();
}

void Instantiation::SetIndex(int i) {
	numunknown = 0;
	index = i;
	for(map<int,int>::iterator ci=card.begin();
			ci!=card.end();++ci) {
		vals[ci->first] = i%(ci->second);
		i /= ci->second;
	}
}

void Instantiation::Recalc() {
	indexinc.clear(); index = 0; numunknown = 0;
	int acc = 1;
	for(map<int,int>::iterator i=card.begin(),j=vals.begin();
			i!=card.end();++i,++j) {
		indexinc.insert(make_pair(i->first,acc));
		if (j->second != -1) index += acc*j->second;
		else numunknown++;
		acc *= i->second;
	}
}


Context Instantiation::MissingVars() const {
	Context c;
	for(map<int,int>::const_iterator i=vals.begin(),j=card.begin();
			i!=vals.end();++i,++j)
		if (i->second==-1) c.AddVar(j->first,j->second);
	return c;
}

Context Instantiation::KnownVars() const {
	Context c;
	for(map<int,int>::const_iterator i=vals.begin(),j=card.begin();
			i!=vals.end();++i,++j)
		if (i->second!=-1) c.AddVar(j->first,j->second);
	return c;
}

void Instantiation::ConsistentIndexes(vector<int> &l) const {
	l.clear();

	vector<int> maxind,incind,v;
	int n = 0;
	for(map<int,int>::const_iterator
				i=vals.begin(),j=card.begin(),k=indexinc.begin();
			i!=vals.end();++i,++j,++k)
		if (i->second == -1) {
			v.push_back(0);
			maxind.push_back(j->second);
			incind.push_back(k->second);
			n++;
		}
	CConsistentIndexes(l,maxind,incind,index);
}


void Instantiation::InconsistentIndexes(vector<int> &l) const {
	l.clear();
	vector<int> lneg;
	CInconsistentIndexes(l,lneg);
}

void Instantiation::SetAllVal(int val) {
	for(map<int,int>::iterator i=vals.begin(),j=indexinc.begin();
			i!=vals.end();++i,++j) {
		index += ((val==-1?0:val)-(i->second==-1?0:i->second))*j->second;
		i->second = val;
	}
	numunknown = val==-1 ? vals.size() : 0;
}


void Instantiation::PrintVal(ostream &os) const {
	for(map<int,int>::const_iterator i=vals.begin();i!=vals.end();++i)
		os << i->second << ' ';
	os << endl;
}

// ---------------------------------------------

} // end of ctbn namespace
