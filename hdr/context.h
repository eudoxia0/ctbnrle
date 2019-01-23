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
#ifndef CTBNRLE_CONTEXT_H
#define CTBNRLE_CONTEXT_H

#include "streamserial.h"

#include <vector>
#include <map>
#include <numeric>
#include <iosfwd>

namespace ctbn {

// The two classes here (context and instantiation) are coupled and
// represent the event space of a set of random variables and
// an assignment to the same, respectively

class Instantiation;

// Represents the event space of a set of random variables
// That is, it represents the cross product of a set of discrete sets.
// Each event space (for a random variable) is assumed to be the
// integers from 0 through n-1 (for a variable with n possible values).
// The variables are also indexed by integers (varid), which need not
// be consecutive.
class Context {
public:
	typedef std::map<int,int> iimap; // general int->int map
	Context(); // assumes null context
	Context(Context const &);

	enum setop {UNION, INTERSECTION, DIFFERENCE};
	Context(const Context &c1, const Context &c2, setop op=UNION);
	// creates a context of the variables that are difference between
	//  the two instantiations (variables that are "unmentioned"
	//  in one or the other are not included)
	Context(const Instantiation &i1, const Instantiation &i2);
	Context(std::istream &is);

	virtual ~Context();

	// Maps an assignment to all variables onto an integer from
	// 0 to Size().  If Instantiation is missing some values, returns -1
	int Index(const Instantiation &ind) const;
	// Reverse mapping from above
	Instantiation Index(int i) const;

	// assumes the varid is not already there...
	// if not true, future behavior of the object is undefined
	void AddVar(int varid, int cardinality) {
		size*=cardinality;
		card[varid]=cardinality;
	}

    // this one checks first...
    void AddVarCheck(int varid, int cardinality) {
        if (!HasId(varid)) AddVar(varid,cardinality);
    }

	// Lists the set of variable IDs
	std::vector<int> VarList() const;

	// Total number of possible assignments to these variables
	// (product of the cardinalities of each constituent variable)
	int Size() const {
		return size;
	}

	int NumVars() const { return card.size(); }
	// smallest varid
	int MinVar() const { return card.empty() ? -1 : card.begin()->first; }
	// largest varid
	int MaxVar() const { return card.empty() ? -1 : card.rbegin()->first; }
	// does this varid exist in this context?
	bool HasId(int varid) const { return card.find(varid) != card.end(); }
	// total number of values varid can take on (from 0 - (n-1))
	int Cardinality(int varid) const {
		iimap::const_iterator l = card.find(varid);
		return l != card.end() ? l->second : 0;
	}

	// returns, in l, all of the indexes (like those that would be found
	// from calling Index()) that are consistant with the (partial) 
	// assignment v
	void ConsistentIndexes(std::vector<int> &l, const Instantiation &v) const;
	// Like above, but returns all of the indexes that are *not*
	// consistant.
	void InconsistentIndexes(std::vector<int> &l, const Instantiation &v) const;

	virtual std::ostream &SaveOld(std::ostream &os) const;
	virtual std::istream &LoadOld(std::istream &is);

	// returns whether c is a subset (not proper) of *this
	bool IsSubset(const Context &c) const;
	// Does *this and c share any variables in common?
	bool IsOverlap(const Context &c) const;

protected:
	iimap card; // maps id to max value
	int size;

	// These "complete" the calculation (and are reuseable in other
	//  similar calculations)
	void CConsistentIndexes(std::vector<int> &l, const std::vector<int> &maxind, const std::vector<int> &incind, int startindex) const;
	void CInconsistentIndexes(std::vector<int> &l, const std::vector<int> &lneg) const;
	SERIAL_START(Context)
		SERIAL_VAR(iimap,card)
	SERIAL_END
public:
	void serial_postload();
};

// forms union:
inline Context operator+(const Context &c1, const Context &c2) {
	return Context(c1,c2,Context::UNION);
}
// forms difference
inline Context operator-(const Context &c1, const Context &c2) {
	return Context(c1,c2,Context::DIFFERENCE);
}
// forms intersection
inline Context operator/(const Context &c1, const Context &c2) {
	return Context(c1,c2,Context::INTERSECTION);
}

// This class is a context *plus* a, possibly partial, assignment
// to the variables in the context
// -1 is used to indicate a non-assigned variable
class Instantiation : public Context {
friend class Context;
public:
	Instantiation();  // assumes null context
	Instantiation(const Context &c, int defaultval=-1);
	Instantiation(const Context &c, const Instantiation &v);
	virtual ~Instantiation();

	// returns the value of the variable varid
	// (from 0 to Cardinality(varid)-1), or -1 if the variable is unassigned
	int Value(int varid) const {
		iimap::const_iterator l = vals.find(varid);
		return l != vals.end() ? l->second : -1;
	}

	// the following three Inc methods increment the instantiation
	// they return false iff the increment has "wrapped"

	// increments only varid -- simple
	bool Inc(int varid);
	// increments all variables lexographically
	bool Inc();
	// increments only those variables in *this that are also in c
	// (again, lexographically)
	bool Inc(const Context &c);

	// val needs to be either -1 (to indicate the value is not assigned)
	// or in [0,n-1] where n is the cardinality of varid
	bool SetVal(int varid, int val) { // returns whether it was added
		iimap::iterator l = vals.find(varid);
		if (l != vals.end()) {
			ChangeVar(varid,l->second,val);
			l->second = val;
			return true;
		} else return false;
	}

	// Set the value of variables mentioned in val.  If add is
	// true, then variables mentioned in val but not in *this are
	// added to *this.  If add is false, they are not.
	void SetVal(const Instantiation &val, bool add=false) {
		for(iimap::const_iterator l = val.vals.begin(),
				ll = val.card.begin();
				l!=val.vals.end();++l,++ll)
			if (!SetVal(l->first,l->second) && add)
				AddVar(l->first,ll->second,l->second);
	}

	// Set every variable to value val
	void SetAllVal(int val);

	// Add variable varid with cardinality cardinality to the set
	// initialize the value to val (-1 = unassigned)
	void AddVar(int varid, int cardinality, int val=-1);

	virtual std::ostream &SaveOld(std::ostream &os) const;
	virtual std::istream &LoadOld(std::istream &is);

	// Returns the index (see Context::Index()) for the current
	// instantiation (like calling Index(*this), but faster)
	// Returns -1 if any variables are unassigned.
	int Index() const { return numunknown ? -1 : index; }
	// Makes the index equal to i.  i may *not* be -1
	void SetIndex(int i);

	// used with ZeroOut from IMatrix (and ZeroOut above)
	// See same method names from Context
	void InconsistentIndexes(std::vector<int> &l) const;
	void ConsistentIndexes(std::vector<int> &l) const;

	// Returns a Context of all the variables that are assigned
	Context KnownVars() const;
	// Same for variables that are not assigned
	Context MissingVars() const;
	int NumMissingVars() const { return numunknown; }
        
	void PrintVal(std::ostream &os) const; //debug use only
	bool operator==(const Instantiation &i) const;
	bool operator!=(const Instantiation &i) const;

protected:
	void ChangeVar(int varid, int oldval, int newval) {
		if (oldval==-1) { numunknown--; oldval = 0; }
		if (newval==-1) { numunknown++; newval = 0; }
		index += indexinc[varid]*(newval-oldval);
	}
	void Recalc();

	iimap vals;
	iimap indexinc;
	int index, numunknown;
	SERIAL_START(Instantiation)
		SERIAL_SUPER(Context)
		SERIAL_VAR(iimap,vals)
	SERIAL_END
public:
	void serial_postload();
};

std::ostream& operator<< (std::ostream& os, const Instantiation& inst);
std::istream& operator>> (std::istream& is, Instantiation& inst);

std::ostream& operator<< (std::ostream& os, const Context& con);
std::istream& operator>> (std::istream& is, Context& con);

// range of all assign to context
class InstantRange {
public:
	class iterator {
	public:
		iterator(const Instantiation &si, bool done=false) : i(si), d(done) {}
		iterator &operator++() { d = !i.Inc(); return *this; }
		iterator operator++(int){ iterator ret(*this); ++(*this); return ret;}
		Instantiation &operator*() { return i; }
		Instantiation *operator->() { return &i; }
		bool operator==(const iterator &ii)
			{ return (ii.d&&d) || (ii.d==d && ii.i==i); }
		bool operator!=(const iterator &ii) { return !(*this==ii); }
		bool isdone() { return d; }
	protected:
		Instantiation i;
		bool d;
		
	};

	InstantRange(const Context &con) : c(con) {}
	iterator begin() const { return iterator(Instantiation(c,0)); }
	iterator end() const { return iterator(Instantiation(),true); }

protected:
	const Context &c;
};

class InstantSubRange : public InstantRange {
public:
	class iterator : public InstantRange::iterator {
	public:
		iterator(const Instantiation &si, const Context &incc,
			bool done=false) : InstantRange::iterator(si,done), c(incc) {}
		iterator &operator++() { d= !i.Inc(c); return *this; }
		iterator operator++(int){ iterator ret(*this); ++(*this); return ret;}
	protected:
		const Context c;
	};
	friend class iterator;

	// range of all assign to context, keeping other assignments in start
	// unchanged
	InstantSubRange(const Instantiation &start, const Context &con)
			: InstantRange(con), si(start) {}
	iterator begin() const {
		Instantiation s(si); s.SetVal(Instantiation(c,0),true);
		return iterator(s,c);
	}
	iterator end() const { return iterator(Instantiation(),c,true); }

protected:
	const Instantiation &si;
};

class DiffByOneRange : public InstantRange {
public:
	class iterator : public InstantRange::iterator {
	public:
		iterator(const Instantiation &si, const Context &c,
				bool includesame=false,bool done=false)
				: InstantRange::iterator(si,done) {
			incsame = includesame;
			if (!done) {
				vars = c.VarList();
				chvar = 0;
				nextvar();
			}
		}
		iterator &operator++() {
			if (chvar==vars.size()) d = true;
			else {
				int v = i.Value(vars[chvar]);
				int mval = i.Cardinality(vars[chvar]);
				if (v==mval-1 || (v==mval-2 && origval==mval-1)) {
					i.SetVal(vars[chvar],origval);
					chvar++;
					nextvar();
				} else {
					v++;
					if (v==origval) v++;
					i.SetVal(vars[chvar],v);
				}
			}
			return *this;
		}
		iterator operator++(int){ iterator ret(*this); ++(*this); return ret;}
		int diffvar() const
			{ return (d||chvar==vars.size()) ? -1 : vars[chvar]; }
	protected:
		std::vector<int>::size_type chvar;
		int origval;
		std::vector<int> vars;
		bool incsame;

		void nextvar() {
			while(chvar<vars.size() && i.Cardinality(vars[chvar])<2) chvar++;
			if (chvar==vars.size()) {
				if (!incsame) d = true;
				else incsame = false;
			} else {
				origval = i.Value(vars[chvar]);
				if (origval!=0) i.SetVal(vars[chvar],0);
				else i.SetVal(vars[chvar],1);
			}
		}

	};

	DiffByOneRange(const Instantiation &start, const Context &con,
			bool includesame=false) : InstantRange(con), si(start) {
		incsame = includesame;
	}
	iterator begin() const { return iterator(si,c,incsame); }
	iterator end() const { return iterator(si,c,incsame,true); } 

protected:
	const Instantiation &si;
	bool incsame;
};

// some macros to make the above a little easier to read
// 
// forallassign(Instantiation &i,c)
//      <loop-body>
// will iterate over all instantiations to c (in i)
// Can declare i to not be ref (but then copy each loop)
// Can declare i outside of loop and use as forallassign(i,c)
//   (but then cannot be reference)
#define forallassign(insvar,context) \
	for(InstantRange::iterator _contextit = InstantRange(context).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true)

// same as above, but instantiation (insvar) cycles through all
// assignments to context, but it also holds values for variables in 
// startins (those not mentiond in context)  startins need not have assignment
// to variables in context
#define forallassign_subset(insvar,context,startins) \
	for(InstantSubRange::iterator _contextit = \
			InstantSubRange(startins,context).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true)

// same as above, but the variable to change (context) are take to be those
// unassigned in startins
#define forallassign_missing(insvar,startins) \
	for(InstantSubRange::iterator _contextit = \
			InstantSubRange(startins,startins.MissingVars()).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true)

// similar to above, but total variable set is specified by context
// and startins specifies the values for variables not to change
#define forallassign_except(insvar,context,startins) \
	for(InstantSubRange::iterator _contextit = \
			InstantSubRange(startins,context-startins).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true)


// similar to those above, this iterates over all instantiations that
// differ only by one variable from startins
#define foralldiffone(insvar,startins) \
	for(DiffByOneRange::iterator _contextit = \
			DiffByOneRange(startins,startins).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true)

// same as above, but also interates (as last step) over startins itself
#define foralldiffone_incself(insvar,startins) \
	for(DiffByOneRange::iterator _contextit = \
			DiffByOneRange(startins,startins,true).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true)

// next two same as above two, except that the variables that can change
// are limited to those in context
#define foralldiffone_subset(insvar,startins,context) \
	for(DiffByOneRange::iterator _contextit = \
			DiffByOneRange(startins,context).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true)

#define foralldiffone_incself_subset(insvar,startins,context) \
	for(DiffByOneRange::iterator _contextit = \
			DiffByOneRange(startins,context,true).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true)

// next four are the same as the above four, except that the 
// second argument is a variable (or declaration of one) of type int
// that will be set to the variable number of the variable that
// differs (or -1 if no variable differs)
#define foralldiffone_diff(insvar,dvar,startins) \
	for(DiffByOneRange::iterator _contextit = \
			DiffByOneRange(startins,startins).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true) \
		if (bool _contextbool2 = false) {} \
		else for(dvar = _contextit.diffvar(); !_contextbool2; \
			_contextbool2=true)

#define foralldiffone_incself_diff(insvar,dvar,startins) \
	for(DiffByOneRange::iterator _contextit = \
			DiffByOneRange(startins,startins,true).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true) \
		if (bool _contextbool2 = false) {} \
		else for(dvar = _contextit.diffvar(); !_contextbool2; \
			_contextbool2=true)

#define foralldiffone_subset_diff(insvar,dvar,startins,context) \
	for(DiffByOneRange::iterator _contextit = \
			DiffByOneRange(startins,context).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true) \
		if (bool _contextbool2 = false) {} \
		else for(dvar = _contextit.diffvar(); !_contextbool2; \
			_contextbool2=true)

#define foralldiffone_incself_subset_diff(insvar,dvar,startins,context) \
	for(DiffByOneRange::iterator _contextit = \
			DiffByOneRange(startins,context,true).begin(); \
				!_contextit.isdone();++_contextit) \
		if (bool _contextbool = false) {} \
		else for(insvar = *_contextit; !_contextbool; _contextbool=true) \
		if (bool _contextbool2 = false) {} \
		else for(dvar = _contextit.diffvar(); !_contextbool2; \
			_contextbool2=true)

		
} // end of ctbn namespace

#endif
