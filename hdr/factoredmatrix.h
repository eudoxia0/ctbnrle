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
#ifndef CTBNRLE_FACTORED_MATRIX_H
#define CTBNRLE_FACTORED_MATRIX_H

#include "ctbndyn.h"
#include "factoredvector.h"
// cshelton: added
#include "expmcache.h"
#include "rk.h"

#include <vector>

namespace ctbn {
	
#define RKEPS 1e-8

class FactoredMatrix :public RVCondSimple {

	friend class FactoredVector;
	friend class TOPETree;

public:
	FactoredMatrix();
	FactoredMatrix(Dynamics *dyn);
	FactoredMatrix(const FactoredMatrix &m);
	~FactoredMatrix();
	
	RVCondSimple *Clone() const;
	
	double FindMinVal(Dynamics *node) const;
	double FindMinRate() const;
	void Cond(double t0, double t1, const Instantiation &x);
	void Cond(double t0, const Instantiation &from, const Instantiation &to, bool transition=true);
	void CondMultByVec(FactoredVector &v, bool transpose, double &maxval) const;
	void MultByVec(FactoredVector &v, bool transpose, double &maxval) const;
	void TransMultByVec(FactoredVector &v, bool transpose) const;
	void TransMult(RVSimple *&x, bool transpose) const;
	void Mult(RVSimple *&x) const;
	void Mult(RVSimple *&x, bool transpose) const;
	double FindAvgRate(FactoredVector &v, const Dynamics *node, int from, int to,
		const std::vector<int> &cind, bool transpose, double &maxval) const;

	void RMult(RVSimple *&x) const;
	double Unifvexpmt(FactoredVector &v, bool conditioned, bool transition) const;
	
	//tope
	// cshelton: added
	void initRKcache(double maxtime, double eps=1e-12) const;
	double RKvexpmtA(FactoredVector & v, double time, bool conditioned, bool transpose,
			vector< vector<double> > *timesteps=NULL) const;
	void BMultByVec(FactoredVector &v, int varid, int bno, bool transpose) const;
	matrix *MakeAB(Dynamics *bnode);
	//

	void SetL(int terms_in_taylor_expansion);
	void SetTheta(double theta);

	void Load(std::istream &is);
	void Save(std::ostream &os) const;


	// returns the probability of an individual conditional event
	 double Prob(int ind, int cind, bool log=false) const;

	// Makes sum to 1 for each conditioning value
	 void Normalize();

	// returns an object representing the conditioning
	// on that value
	 RVSimple *Condition(int ind) const;

	// Sets to zero all indexes that are not listed in ind and cind
	// That is, for conditions mentioned in cind, the measure is
	// restricted to the indexes in ind.  For conditions not in cind,
	// the measure is made to be zero everywhere.
	 void Restrict(const std::vector<int> &ind, 
			const std::vector<int> &cind);


	// Marginalizes, reindexes, or expands depending on the
	// argument
	// n = # new values  nc = # new conditioning values
	// ind[<i,j>] is a mapping that needs to be applied to old 
	//   conditioning index j and then added to i
	// (see Reindex for RVSimple above to see how a mapping is
	//  applied)
	// (Really, this is probably too confusing to be "simple.")
	typedef std::map<std::pair<int,int>, std::vector<std::vector<int> > > RemapT;
	void Reindex(int n, int nc, const RemapT &ind);

	// Multiplies by another measure (point-wise)
	void MultBy(const RVCondSimple *x);

	// return suitable sufficient statistics object
	SS *BlankSS() const;

	// add independent draw of (cind,ind) (conditioning value, value) to ss
	void AddSS(int cind, int ind, SS *ss, double w=1.0) const;
	// add average of independent draws to ss
	void AddExpSS(int cind, SS *ss, double w=1.0) const;
	void AddExpSS(int cind, const RVSimple *rvs,
			SS *ss, double w) const;

	void AddSS(const SS *toadd, const RVCondSimple* rv,
			const std::vector<std::vector<int> > &mapping,
			int mycondi, int rvccondi,
			SS *ss, double w=1.0) const;

	// returns a sample
	int Sample(int cx, Random &rand = randomizer) const;

	// sets parameters to ML estimate
	void Maximize(const SS *ss);

	// Returns a new object that can be multiplied into this
	// object (as per Mult and RMult)
	 RVSimple * Convert(const RVSimple *rvs) const;

	// see above class for description
	void Scramble(double alpha=1.0, double degree=1.0, Random &rand=randomizer); //by Yu
	double LLH(const SS *ss) const;
	double GetScore(double numTrans, const SS *ss)const;
///////

private:

	void Uniformize(Dynamics *node);
	
	double minrate;
	double t;
	bool trans;

	int taylor_expansion_term_count;
	double theta_value;

    //holds nodes by varids
    std::vector<Dynamics *> nodes;
    
    //for tope
    std::vector<Dynamics *> bnodes;
	std::vector<matrix *> amatrices;

	// cshelton: added
	mutable double cachemaxt;
	mutable std::vector<ExpMCache> rkcache;
	//
	
	//consistent parent indices of the variables that are not in the evidence
	std::map<int, std::vector <int> > condpars;
	std::map<int, std::vector <int> > children_in_evidence;
	///consistent parent indices of the variables that are in the evidence
	std::map<int, std::vector <int> > condpars_evid;
	std::map<int, int> evid;
	std::map<int, std::vector<int> > trans_evid; //for from to conditioning
	std::map<int, int> notrans_evid; //for evidence that went from uninst to inst (when there is no transition)	

};


}
#endif
