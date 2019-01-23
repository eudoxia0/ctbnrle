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

#include "linearprogram.h"
#include "defines.h"

#include <cmath>
#include <iostream>


#define USE_SIMPLEX


// could be duplicated from defines.h, but this makes this package
// more independent

#ifndef INFINITY
#define INFINITY (std::numeric_limits<double>::infinity())
#endif

using namespace std;

namespace ctbn {

GLPKSolver::GLPKSolver() : LinearProgram() {
	p = glp_create_prob();
	nconstr = 0;
}

GLPKSolver::GLPKSolver(const GLPKSolver &lp) : mins(lp.mins), maxs(lp.maxs) {
	p = glp_create_prob();
	glp_copy_prob(p,lp.p,GLP_OFF);
	nconstr = lp.nconstr;
}

GLPKSolver::~GLPKSolver() throw() {
	glp_delete_prob(p);
}

void GLPKSolver::ClearProblem() {
	glp_delete_prob(p);
	p = glp_create_prob();
	nconstr = 0;
	mins.clear(); maxs.clear();
}

GLPKSolver &GLPKSolver::operator=(const GLPKSolver &lp) {
	if (this == &lp) return *this;
	glp_delete_prob(p);
	glp_copy_prob(p,lp.p,GLP_OFF);
	nconstr = lp.nconstr;
	mins = lp.mins;
	maxs = lp.maxs;
	return *this;
}

void GLPKSolver::SetObjective(const vector<double> &coeff) {
	ResizeVars(coeff.size()-1);
	for(unsigned int i=0;i<coeff.size();i++)
		glp_set_obj_coef(p,i+1,coeff[i]);
}

void GLPKSolver::SetObjective(const vector<int> &ix, const vector<double> &coeff) {
	ResizeVars(ix);
	for(unsigned int i=0;i<coeff.size();i++)
		glp_set_obj_coef(p,ix[i]+1,coeff[i]);
}

void GLPKSolver::AddEqConstraint(const vector<double> &coeff, double val) {
	int ncoeff = coeff.size();
	ResizeVars(ncoeff-1);
	// indices in GLPK are 1-indexed (hence the +1 everywhere)
	int *ixs = new int[ncoeff+1];
	for(int i=1;i<=ncoeff;i++)
		ixs[i] = i;
	AddGenRow(coeff,ixs,val,GLP_FX);
	delete []ixs;
}

void GLPKSolver::AddEqConstraint(const vector<int> &ix,
		const vector<double> &coeff, double val) {
	int ncoeff = coeff.size();
	ResizeVars(ix);
	// indices in GLPK are 1-indexed (hence the +1 everywhere)
	int *ixs = new int[ncoeff+1];
	for(int i=0;i<ncoeff;i++)
		ixs[i+1] = ix[i]+1;
	AddGenRow(coeff,ixs,val,GLP_FX);
	delete []ixs;
}

void GLPKSolver::AddGrEqConstraint(const vector<double> &coeff, double val) {
	int ncoeff = coeff.size();
	ResizeVars(ncoeff-1);
	// indices in GLPK are 1-indexed (hence the +1 everywhere)
	int *ixs = new int[ncoeff+1];
	for(int i=1;i<=ncoeff;i++)
		ixs[i] = i;
	AddGenRow(coeff,ixs,val,GLP_LO);
	delete []ixs;
}

void GLPKSolver::AddGrEqConstraint(const vector<int> &ix,
		const vector<double> &coeff, double val) {
	int ncoeff = coeff.size();
	ResizeVars(ix);
	// indices in GLPK are 1-indexed (hence the +1 everywhere)
	int *ixs = new int[ncoeff+1];
	for(int i=0;i<ncoeff;i++)
		ixs[i+1] = ix[i]+1;
	AddGenRow(coeff,ixs,val,GLP_LO);
	delete []ixs;
}

void GLPKSolver::ReassertVarBounds(int varid) {
	if (isinfinite(mins[varid])) {
		if (isinfinite(maxs[varid]))
			glp_set_col_bnds(p,varid+1,GLP_FR,0,0);
		else glp_set_col_bnds(p,varid+1,GLP_UP,0,maxs[varid]);
	} else {
		if (isinfinite(maxs[varid]))
			glp_set_col_bnds(p,varid+1,GLP_LO,mins[varid],0);
		else if (maxs[varid]==mins[varid]) 
			glp_set_col_bnds(p,varid+1,GLP_FX,mins[varid],maxs[varid]);
		else glp_set_col_bnds(p,varid+1,GLP_DB,mins[varid],maxs[varid]);
	}
}

void GLPKSolver::AddGenRow(const std::vector<double> &coeff,
		int *ix, double val, int type) {
	if (coeff.size()==1) { // just place bounds directly on variable
		int varid = ix[1]-1;
		double realval = val/coeff[0];
		int newtype = type;
		if (coeff[0]<0.0) {
                        if (type==GLP_LO) {
                                newtype = GLP_UP;
			} else {
                                if (type==GLP_UP) newtype = GLP_LO;
			}
		}
		switch (newtype) {
			case GLP_LO:
				if (realval > mins[varid]) {
					if (realval < maxs[varid]) {
						mins[varid] = realval;
						ReassertVarBounds(varid);
						return;
					} else break; // we have inconsistent... so add as normal
						// so that the solver will report it (same below)
				} return; // bound is subsumed by current bound (same below)
			case GLP_UP:
				if (realval < maxs[varid]) {
					if (realval > mins[varid]) {
						maxs[varid] = realval;
						ReassertVarBounds(varid);
						return;
					} else break;
				} return;
			case GLP_FX:
				if (realval < maxs[varid] && realval > mins[varid]) {
					mins[varid] = maxs[varid] = realval;
					ReassertVarBounds(varid);
					return;
				} else break;
		}
	}
	glp_add_rows(p,1);
	nconstr++;
	int ncoeff = coeff.size();
	// indices in GLPK are 1-indexed (hence the +1 everywhere)
	double *c = new double[ncoeff+1];
	for(int i=0;i<ncoeff;i++)
		c[i+1] = coeff[i];
	glp_set_mat_row(p,nconstr,coeff.size(),ix,c);
	glp_set_row_bnds(p,nconstr,type,val,val);
	delete []c;
}

void GLPKSolver::ResizeVars(int maxvarnum) {
	int toadd = maxvarnum+1-mins.size();
	if (toadd>0) {
		glp_add_cols(p,toadd);
		for(int i=0;i<toadd;i++) {
			glp_set_col_bnds(p,mins.size()+1,GLP_FR,0,0);
			mins.push_back(-INFINITY);
			maxs.push_back(+INFINITY);
		}
	}
}

void GLPKSolver::ResizeVars(const vector<int> &ix) {
	int maxv = -1;
	for(vector<int>::const_iterator i=ix.begin();i!=ix.end();++i)
		if (*i>maxv) maxv = *i;
	ResizeVars(maxv);
}


vector <double> GLPKSolver::Solve(double &val) {
	vector<double> ret;
	glp_set_obj_dir(p,GLP_MIN);
	//glp_write_lp(p,NULL,"debug.lp");
#ifdef USE_SIMPLEX
	glp_smcp param;
	glp_init_smcp(&param);
	param.meth = GLP_PRIMAL; // GLP_DUAL;  // GLP_DUALP;
	param.presolve = GLP_OFF; // GLP_ON;
	param.msg_lev = GLP_MSG_OFF; // GLP_MSG_ON;
	int retcode = glp_simplex(p,&param);
	int status = glp_get_status(p);
	if (retcode || status==GLP_INFEAS) {
		val = nan("");
		return ret;
	}
	if (status==GLP_UNBND || status==GLP_NOFEAS) {
		val = -INFINITY;
		return ret;
	}
	ret.resize(mins.size());
	for(unsigned int i=0;i<mins.size();i++)
		ret[i] = glp_get_col_prim(p,i+1);
	val = glp_get_obj_val(p);
#else // use interior-point / barrier
	glp_iptcp param;
	glp_init_iptcp(&param);
	param.msg_lev = GLP_MSG_OFF; // GLP_MSG_ON;
	int retcode = glp_interior(p,&param);
	int status = glp_ipt_status(p);
	if (retcode!=GLP_EFAIL && (retcode!=0
			|| status==GLP_INFEAS || status==GLP_UNDEF)) {
		val = nan("");
		return ret;
	}
	if (retcode==GLP_EFAIL || status==GLP_UNBND || status==GLP_NOFEAS) {
		val = -INFINITY;
		return ret;
	}
	ret.resize(mins.size());
	for(size_t i=0;i<mins.size();i++)
		ret[i] = glp_ipt_col_prim(p,i+1);
	val = glp_ipt_obj_val(p);
#endif

	return ret;
}

} // end ctbn namespace
