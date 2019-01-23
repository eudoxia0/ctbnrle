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

#include "ctmdp.h"
#include "params.h"

namespace ctbn {

using namespace std;

CTMDP::CTMDP() {
	beta = 1.0;
	resetlpvarmap();
}

CTMDP::~CTMDP() {
	for(unsigned int i=0;i<actions.size();i++)
		delete actions[i].dyn;
}

void CTMDP::Clear() {
	beta = 1.0;
	basis.clear();
	for(unsigned int i=0;i<actions.size();i++)
		delete actions[i].dyn;
	actions.clear();
	resetlpvarmap();
}

//adding the CTBNDyn with the associating Factor
void CTMDP::AddAction(const Dynamics *dyn, const vector<Factor> &r) {
	ActionInfo a;
	a.dyn = dyn;
	a.rs = r;
	actions.push_back(a);
}

void CTMDP::AddBasis(const Factor &phi) {
	basis.push_back(phi);
}

// expands factor h by the dynamics dyn (see paper)
CTMDP::Factor CTMDP::ExpandFactor(const CTMDP::Factor &h,
		const Dynamics *dyn) const {
	Factor newfactor;
	const vector<int> &basisVar = h.var.VarList();
	// get the union of all the contexts in ctbndyn
	// which associate with factor h
	for(unsigned int i = 0; i < basisVar.size(); i++) {
		const MarkovDyn *node = dynamic_cast<const MarkovDyn *>
				(dynamic_cast<const CTBNDyn *>(dyn)
					->NodeByVar(basisVar[i]));
		newfactor.var = newfactor.var + node->Domain()+node->CondDomain();
	}

	//initiate the instantiation for union of the context
	vector<double> basisAdj = h.val;
	forallassign(Instantiation &inst,newfactor.var) {
		double valueG = 0.0;
		int placeH, placeH2 = 0;
		//h_i (x) place
		placeH = h.var.Index(inst);
		//sum all the Q which x_j to x_j^' and h_i(x <- x_j^')
		for(unsigned int k = 0; k < basisVar.size(); k ++) {
			const MarkovDyn *node = dynamic_cast<const MarkovDyn *>
					(dynamic_cast<const CTBNDyn *>(dyn)
						->NodeByVar(basisVar[k]));
			const matrix &Q = (*node)(inst)->Intensity();
			for (int m = 0; m < Q.getn(); m ++) {
				if(inst.Value(basisVar[k]) != m)
				{
					Instantiation newinst(inst);
					newinst.SetVal(basisVar[k],m);
					placeH2 = h.var.Index(newinst);
					int i1 = node->Domain().Index(inst);
					int i2 = node->Domain().Index(newinst);
					valueG += Q[i1][i2]*(h.val[placeH2]-h.val[placeH]);
				}
			}
		}
		valueG -= beta*h.val[placeH];
		newfactor.val.push_back(valueG);
		valueG = 0.0;
	}
	return newfactor;
}

int findClique(const vector<CliqueTree::Node> &c, const Context &v) {
	for(unsigned int i=0;i<c.size();i++)
		if (c[i].vars.IsSubset(v)) return i;
	return -1;
}

vector<double> CTMDP::SolveApprox(LinearProgram *lpsolver) const {
	lpsolver->ClearProblem();
	bool norandom = ParamInt("CTMDP_NoRandom",0);

	// first first set are reserved for the basis function weights:
	maxvar = basis.size();
	set <int> DoneBF;
	for(unsigned int a = 0; a < actions.size(); a++ ) {
		// build clique tree for expanded basis factors and reward factors:
		vector<Context> fcon;
		vector<Factor> g;
		for(unsigned int i=0;i<basis.size();i++) {
			Factor gg = ExpandFactor(basis[i],actions[a].dyn);
			g.push_back(gg);
			fcon.push_back(gg.var);
		}
		for(unsigned int i=0;i<actions[a].rs.size();i++)
			fcon.push_back(actions[a].rs[i].var);
		CliqueTree clique_build;
		clique_build.IgnoreRand(norandom);
		clique_build.BuildClique(fcon);
		vector<CliqueTree::Node> cliqueAct = clique_build.ReturnCliques();
		/*
		for(int i=0;i<cliqueAct.size();i++) {
			cout << "clique " << i << ":" << endl;
			cout << "\tvars: " << cliqueAct[i].vars.VarList() << endl;
			cout << "\tadj: " << cliqueAct[i].adj << endl;
		}
		*/
		// find location of g and r factors in clique tree
		vector<vector<int> > gloc(cliqueAct.size());
		vector<vector<int> > rloc(cliqueAct.size());
		for(unsigned int i=0;i<g.size();i++) {
			int ci = findClique(cliqueAct,fcon[i]);
			assert(ci!=-1);
			gloc[ci].push_back(i);
		}
		for(unsigned int i=0;i<actions[a].rs.size();i++) {
			int ci = findClique(cliqueAct,actions[a].rs[i].var);
			assert(ci!=-1);
			rloc[ci].push_back(i);
		}
		// write out constraints for each clique:
		writeConstraints(cliqueAct,a,0,-1,g,gloc,rloc,lpsolver);
	}
	// set optimization fn weights:
	SetApproxOptWeights(actions[0].dyn->Domain().Size(),lpsolver);

	double value;
	vector<double> sol = lpsolver->Solve(value);
	vector<double> ret;
	for(unsigned int i = 0; i < basis.size(); i++)
		ret.push_back(sol[i]);
	
	resetlpvarmap();
	return ret;
}

void CTMDP::SetApproxOptWeights(int totalsize,LinearProgram *lpsolver) const {
	vector<int> obji;
	vector<double> objw;
	for(unsigned int k=0; k < basis.size(); k++) {
		obji.push_back(k);
		double alpha = 0.0;
		for(unsigned int i=0;i<basis[k].val.size();i++)
			alpha += basis[k].val[i];
		alpha *= totalsize / basis[k].val.size();
		objw.push_back(alpha);
	}
	lpsolver->SetObjective(obji,objw);
}

void CTMDP::writeConstraints(const vector<CliqueTree::Node> & cliqueAct,
		const int &a, const int &ci, const int &parent,
		const std::vector<Factor> &gs,
		const std::vector<std::vector<int> > &gloc,
		const std::vector<std::vector<int> > &rloc,
		LinearProgram *lp) const {

	writeNodeConstraints(cliqueAct,a,ci,parent,gs,gloc,rloc,lp);
	for(unsigned int i=0;i<cliqueAct[ci].adj.size();i++) {
		int conn = cliqueAct[ci].adj[i];
		if (conn!=parent)
			writeConstraints(cliqueAct,a,conn,ci,gs,gloc,rloc,lp);
	}
}

void CTMDP::writeNodeConstraints(const vector<CliqueTree::Node> & cliqueAct,
		const int &a, const int &ci, const int &parent,
		const std::vector<Factor> &gs,
		const std::vector<std::vector<int> > &gloc,
		const std::vector<std::vector<int> > &rloc,
		LinearProgram *lp) const {

	forallassign(Instantiation &inst,cliqueAct[ci].vars) {
		vector<int> ix;
		vector<double> coeff;
		double cons = 0.0;
		for(unsigned int j=0;j<rloc[ci].size();j++)
			cons += actions[a].rs[rloc[ci][j]](inst);
		for(unsigned int j=0;j<gloc[ci].size();j++) {
			ix.push_back(gloc[ci][j]);
			coeff.push_back(-gs[gloc[ci][j]](inst));
		}
		for(unsigned int j=0;j<cliqueAct[ci].adj.size();j++) {
			int adj = cliqueAct[ci].adj[j];
			Context sepset = cliqueAct[ci].vars / cliqueAct[adj].vars;
			ix.push_back(getLPVar(a,ci,adj,sepset.Index(inst)));
			coeff.push_back(adj==parent ? 1 : -1);
		}
		lp->AddGrEqConstraint(ix,coeff,cons);
	}
}

void CTMDP::SetBeta(const double & beta) {
	this -> beta = beta;
}

vector<double> CTMDP::SolveExact(LinearProgram *lpsolver) const {
	const Context &allvars = actions[0].dyn->Domain();
	const vector<int> varlist = allvars.VarList();
	lpsolver->SetObjective(vector<double>(allvars.Size(),1.0));

	for(unsigned int a=0;a<actions.size();a++) {
		forallassign(Instantiation &inst,allvars) {
			int i = inst.Index();
			double cons = 0;
			for(unsigned int j=0;j<actions[a].rs.size();j++)
				cons += actions[a].rs[j](inst);
			vector<double> ws;
			vector<int> xs;
			double qa = 0.0;
			// we assume that only one variable can change at a
			// time (if not true, then we need to change this to a
			// complete loop over all newinst).
			// We also assume it is a CTBN with MarkovDyn nodes, so
			// that would also need to be changed
			foralldiffone_diff(Instantiation &newinst,int chvar,inst) {
				const MarkovDyn *node = dynamic_cast<const MarkovDyn *>
						(dynamic_cast<const CTBNDyn *>(actions[a].dyn)
							->NodeByVar(chvar));
				const matrix &Q = (*node)(inst)->Intensity();
				xs.push_back(newinst.Index());
				double negq = -Q[inst.Value(chvar)][newinst.Value(chvar)];
				ws.push_back(negq);
				qa -= negq;
			}
			xs.push_back(i);
			ws.push_back(beta+qa);
			lpsolver->AddGrEqConstraint(xs,ws,cons);
		}
	}
	double value;
	return lpsolver->Solve(value);
}

vector<int> CTMDP::Policy(const vector<double> &sol, bool isapprox) const {
	vector<int> ret;
	forallassign(Instantiation &x,actions[0].dyn->Domain())
		ret.push_back(Policy(x,sol,isapprox));
	return ret;
}

int CTMDP::Policy(const Instantiation &x, const vector<double> &sol,
		bool isapprox) const {
	int ret=-1;
	double val=-INFINITY;
	const Context &allvars = actions[0].dyn->Domain();
	const vector<int> varlist = allvars.VarList();
	for(unsigned int a=0;a<actions.size();a++) {
		double cval = 0;
		for(unsigned int j=0;j<actions[a].rs.size();j++)
			cval += actions[a].rs[j](x);
		double qa = 0;
		// we assume that only one variable can change at a
		// time (if not true, then we need to change this to a
		// complete loop over all xp).
		// We also assume it is a CTBN with MarkovDyn nodes, so
		// that would also need to be changed
		foralldiffone_diff(Instantiation &xp,int chvar,x) {
		const MarkovDyn *node = dynamic_cast<const MarkovDyn *>
					(dynamic_cast<const CTBNDyn *>(actions[a].dyn)
						->NodeByVar(chvar));
			const matrix &Q = (*node)(x)->Intensity();
			int oldval = x.Value(chvar), newval=xp.Value(chvar);
			double q = Q[oldval][newval];
			qa += q;
			double vfn = 0;
			if (!isapprox) vfn = sol[xp.Index()];
			else for(unsigned int j=0;j<basis.size();j++)
					vfn += basis[j](xp)*sol[j];
			cval += Q[oldval][newval]*vfn;
		}
		cval /= (beta + qa);
		if (cval>val) {
			val = cval;
			ret = a;
		}
	}
	return ret;
}

vector<double> CTMDP::SolveApprox_exp(LinearProgram * lpsolver) const {
	lpsolver->ClearProblem();

	const Context &allvars = actions[0].dyn->Domain();
	const vector<int> varlist = allvars.VarList();

	for(unsigned int a=0;a<actions.size();a++) {
		forallassign(Instantiation &inst,allvars) {
			int i = inst.Index();
			double cons = 0;
			for(unsigned int j=0;j<actions[a].rs.size();j++)
				cons += actions[a].rs[j](inst);
			vector<double> ws(basis.size(),0.0);
			double qa = 0;
			// we assume that only one variable can change at a
			// time (if not true, then we need to change this to a
			// complete loop over all newinst).
			// We also assume it is a CTBN with MarkovDyn nodes, so
			// that would also need to be changed
			foralldiffone_diff(Instantiation &newinst,int chvar,inst) {
				const MarkovDyn *node = dynamic_cast<const MarkovDyn *>
						(dynamic_cast<const CTBNDyn *>(actions[a].dyn)
							->NodeByVar(chvar));
				const matrix &Q = (*node)(inst)->Intensity();
				double negq = -Q[inst.Value(chvar)][newinst.Value(chvar)];
				qa -= negq;
				for(unsigned int j=0;j<basis.size();j++)
					ws[j] += basis[j](newinst)*negq;
			}
			for(unsigned int j=0;j<basis.size();j++)
				ws[j] += basis[j](inst)*(beta+qa);
			lpsolver->AddGrEqConstraint(ws,cons);
		}
	}
	SetApproxOptWeights(allvars.Size(),lpsolver);
	double value;
	return lpsolver->Solve(value);
}

int CTMDP::getLPVar(int action, int node, int connect, int value) const {
	IDDP x(action,node,connect,value);
	map<IDDP,int>::iterator i = lpvarmap.find(x);
	if (i==lpvarmap.end())
		i = lpvarmap.insert(make_pair(x,maxvar++)).first;
	return i->second;
}

void CTMDP::resetlpvarmap() const {
	lpvarmap.clear();
	maxvar = 0;
}
}
