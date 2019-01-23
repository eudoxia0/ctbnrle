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
#ifndef CTBNRLE_BN_H
#define CTBNRLE_BN_H

#include "rv.h"
#include "streamserial.h"
#include "structure.h"


namespace ctbn {

class BN;
class RVSimple;

// sufficient statistics object for a Bayesian network
class BNSS : public SS {
	friend class BN;
	SOBJCLASSDECL(BNSS)
public:
	BNSS(std::istream &is);
	BNSS(const BNSS &bnss);

	BNSS &operator=(const BNSS &bnss);
	virtual BNSS *Clone() const;

	virtual ~BNSS();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
        virtual void Scale(double w);
        virtual void AddSS(const SS* nss, double w=1.0);

	BNSS();
private:
	std::vector<SS *> rvss;

	SERIAL_START_V(BNSS)
		SERIAL_SUPER(SS)
		SERIAL_VAR(std::vector<SS*>,rvss)
	SERIAL_END
public:
	void serial_preload();
};


// BN is a class to represent a Bayesian Network
// See RV for description of most of the methods
class BN : public RV {
	SOBJCLASSDECL(BN)

friend class FBInf;
public:
	BN();
	BN(const BN &bn);
	// load from stream
	BN(std::istream &is);
	// build with structure and variables given
	BN(const Structure &s, const Context &vars);

	BN &operator=(const BN &bn);
	virtual BN *Clone() const;

	virtual ~BN();

	// next set: standard RV methods:

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	virtual void Normalize();

	virtual double Prob(const Instantiation &x, bool log=false) const;

	virtual void Restrict(const Instantiation &x);

	virtual void Project(const Context &c, const Context &cc);

	virtual void MultBy(const RV *x);

	virtual BNSS *BlankSS() const;

	virtual void AddSS(const Instantiation &x, 
			SS *ss, double w=1.0) const;
	virtual void AddExpSS(const Instantiation &x, 
			SS *ss, double w=1.0) const;
	virtual void AddExpSS(const Instantiation &x, const RVSimple *rvs,
			SS *ss, double w=1.0) const;

	virtual void AddSS(const SS* toadd, const RV *rv,
			SS *ss, double w=1.0) const {}

	virtual void Sample(Instantiation &x, Random &rand=randomizer) const;

	virtual void Maximize(const SS *ss);

	virtual RVSimple * MakeSimple(Instantiation const & x, bool normalize = true) const;

	virtual void Scramble(double alpha=1.0, double degree=1.0, 
			Random &rand=randomizer);

	virtual double LLH(const SS *ss) const;

	virtual double GetScore(double numTrans, const SS* ss) const;

	// these are specific to a BN:

	// the BN now owns the pointer rv (call with rv->Clone()) if the
	// caller doesn't want to give up the variable
	void AddNode(RV *rv);

	// implements importance sampling and returns the weight of the
	// sample
	double ImportanceSample(Instantiation &x, 
			bool log=false, Random &rand=randomizer) const;

	// get parameters from a BN that might not have the same
	// structure (say to initialize the E step inside SEM)
	void FillParams(const BN &src);

	// i is the index of the node (might not be variable "i")
	const RV *Node(int i) const {return nodes[i]; }
	// i is the "name" of a variable
	const RV *NodeByVar(int i) const;
	int NumofNodes() const { return nodes.size(); }

	// Get a list of children by varid; The children list are also varids
	std::vector<int> GetChildrenByVar(int varid) const;

	void GetStructure(Structure &s) const;


protected:
	void SetSampleOrder() const;
	std::vector<RV *> nodes;
	mutable std::vector<int> sampleorder; // as a cache
	// generates the items below
	// call before use -- will do nothing if they already exist
	void PrepNode2IdMaps() const;
	// clears them (call if nodes changes)
	void ClearNode2IdMaps() const;

	// returns joint distribution restricted to x
	RV *Joint(const Instantiation &x) const;

	RVSimple *MakeSimpleSparse(const Instantiation &x, bool normalize=true) const;

	// cache of mapping from index into "nodes" to the variable id represented
	mutable std::vector<int> node2var;
	// cache of mapping from variable id to the index into "nodes"
	mutable std::map<int,int> var2node;
	// cache of the indexes into "nodes" of the children for a particular
	// variable id (Note this one is slightly odd as it maps a variable id into
	// a list of node-indexes (not back into varids)
	mutable std::map<int,std::vector<int> > children;

	SERIAL_START_V(BN)
		SERIAL_SUPER(RV)
		SERIAL_VAR(std::vector<RV *>,nodes)
	SERIAL_END
public:
	void serial_preload();
};

} // end of ctbn namespace
#endif
