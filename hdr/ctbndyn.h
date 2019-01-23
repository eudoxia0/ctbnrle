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
#ifndef CTBNRLE_CTBNDYN_H
#define CTBNRLE_CTBNDYN_H

#include "dynamics.h"
#include "structure.h"
#include "suffstatsquery.h"



namespace ctbn {

class CTBNDyn;
class Structure;

// These are the sufficient statistics for a CTBN Dynamics
// (that is everything in a CTBN, except the initial distribution)
// It consists of a vector (for each node in the graph) of its
// sufficient statistics.
class CTBNDynSS : public SS {
	friend class CTBNDyn;
	friend class MeanFieldInf;
	SOBJCLASSDECL(CTBNDynSS)
public:
	CTBNDynSS(std::istream &is);
	CTBNDynSS(const CTBNDynSS &ctbnss);

	CTBNDynSS &operator=(const CTBNDynSS &ctbnss);
	virtual CTBNDynSS *Clone() const;

	virtual ~CTBNDynSS();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
	virtual void Scale(double w);
	virtual void AddSS(const SS* nss, double w=1.0);
	virtual double NodeSS(int id, int val1, int val2, 
				const Instantiation &cond) const;
	const SS* GetNodeSS(int id) const {return dynss[id];}   


	CTBNDynSS(const CTBNDyn *d=NULL);
private:
	void Compact() const;
	mutable std::vector<SS *> dynss;
	mutable SS *jss;
	const CTBNDyn *dyn;

	SERIAL_START_V(CTBNDynSS)
		SERIAL_SUPER(SS)
		SERIAL_VAR(std::vector<SS*>,dynss)
	SERIAL_END
public:
	void serial_presave() const;
	void serial_preload();
};

// This is a CTBN Dynamics (that is, everything except the
// initial distribution).  All methods work just as described in
// dynamics.h
class CTBNDyn : public Dynamics {
	friend class CTBNDynSS;
	SOBJCLASSDECL(CTBNDyn)
public:
	CTBNDyn();
	CTBNDyn(std::istream &is);
	CTBNDyn(const CTBNDyn &ctbn);
	CTBNDyn(const Structure &s, const Context& vars);

	CTBNDyn &operator=(const CTBNDyn &ctbn);
	virtual CTBNDyn *Clone() const;

	virtual ~CTBNDyn();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	virtual void Normalize();
	virtual void Restrict(const Instantiation &x);
	virtual void Mult(const Dynamics *x);
	virtual RVCondSimple *Cond(double t0, double t1,
			const Instantiation &x) const;
	virtual RVCondSimple *Cond(double t, const Instantiation &from,
			const Instantiation &to, bool transition=true) const;

  	virtual CTBNDynSS *BlankSS() const;
	virtual void AddSS(const Instantiation &x, double t0,
			double deltat, SS *ss, double w=1.0) const;
	virtual void AddTransSS(const Instantiation &x1, 
				const Instantiation &x2,
				double t, SS *ss, double w=1.0) const;
	virtual void AddExpSS(const RV *alpha, const RV *beta,
			double t0, double deltat, SS *ss, double w=1.0) const;
	virtual void AddExpTransSS(const RV *x1, const RV *x2,
			const Context &changevar,
			double t, SS *ss, double w=1.0) const;
	virtual void AddExpSS(const Instantiation &condi,
			const RVSimple *alpha, const RVSimple *beta,
			double t0, double deltat, SS *ss, double w=1.0) const;
	virtual void AddExpSS(int condi,
			const RVSimple *alpha, const RVSimple *beta,
			double t0, double deltat, SS *ss, double w=1.0) const;
	virtual void AddExpTransSS(const Instantiation &condi,
			const RVSimple *x1, const RVSimple *x2,
			const Context &changevar,
			double t, SS *ss, double w=1.0) const;
	virtual void AddExpTransSS(int condi,
			const RVSimple *x1, const RVSimple *x2,
			const Context &changevar,
			double t, SS *ss, double w=1.0) const;
	virtual void AddSS(const SS *toadd, const Dynamics *dyn,
			SS *ss, double w=1.0) const;
	virtual void SampleNextEvent(const Instantiation &i, double t,
			double &nextt, Instantiation &nexti,
			Random &rand=randomizer) const;
	virtual void SampleTrajectory(Trajectory &tr, double t,
			Random &rand=randomizer) const;
	virtual void Maximize(const SS *ss);
	virtual double LLH(const SS *ss) const;

	// These methods are new to a CTBNDyn:

	// CTBNDyn now owns the pointer (call ->Clone() to give it a copy)
	void AddNode(Dynamics *d);

	//Places the structure of the graph in s
	void GetStructure(Structure &s) const;
	// returns the number of nodes in the graph
	int NumofNodes() const {return nodes.size();}
	// returns the ith node
	const Dynamics *Node(int i) const {return nodes[i];}
	// returns the node for the variable with ID i
	const Dynamics *NodeByVar(int i) const;
	// randomizes the parameters of the network
	virtual void Scramble(double a=1.0, double b=1.0, 
			double alpha=1.0, double degree=1.0, 
			Random &rand=randomizer);
	// returns a joint intensity matrix for the entire
	// CTBN (not recommended if your CTBN is large!)
	// the rows (and columns) are ordered as by the Inc() method
	// for Instantiation (that is, construct an instantiation
	// from the context of this CTBNDyn and increment it to iterate
	// through the assignments for the rows).
	matrix JointMatrix() const;
	// returns a dynamics that represents the entire process
	// (flattened) -- again not recommended for large CTBNs
	const Dynamics *GetJoint() const {return Joint();}

	// returns the Bayesian structure score for the CTBN (given
	// the sufficient statistics) with the BDe prior with numTrans
	// and amtTime as the two meta-parameters
	virtual double GetScore(double numTrans, double amtTime, 
			const SS* ss) const;


	// fills in the parameters from a different CTBNDyn.  Used
	// in structural EM to initial the parameters after the
	// structure has been changed
	void FillParams(const CTBNDyn &src, bool useRandom=false);
	
	// returns the sufficient statistics just for variable ID i
	const SS *NodeSSByVar(int i, SS* ss) const;

	// Get a list of children by varid. The children list are also varids
	std::vector<int> GetChildrenByVar(int varid) const;
protected:
	// cannot delete returned value
	const Dynamics *Joint() const;

	std::vector<Dynamics *> nodes;

	mutable Dynamics *joint; // cache of joint RV


	// generates the items below
	// call before use -- will do nothing if they already exist
	void PrepNode2IdMaps() const;
	// clears them (call if nodes changes)
	void ClearNode2IdMaps() const;

	// cache of mapping from index into "nodes" to the variable id represented
	mutable std::vector<int> node2var;
	// cache of mapping from variable id to the index into "nodes"
	mutable std::map<int,int> var2node;
	// cache of the indexes into "nodes" of the children for a particular
	// variable id (Note this one is slightly odd as it maps a variable id into
	// a list of node-indexes (not back into varids)
	mutable std::map<int,std::vector<int> > children;

	mutable std::vector<CTBNDynSS *> assocss;

	SERIAL_START_V(CTBNDyn)
		SERIAL_SUPER(Dynamics)
		SERIAL_VAR(std::vector<Dynamics*>,nodes)
	SERIAL_END
public:
	void serial_preload();
};

} // end of ctbn namespace

#endif
