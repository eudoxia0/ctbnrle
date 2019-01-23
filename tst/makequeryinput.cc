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
#include "markovdyn.h"
#include "streamextra.h"
#include "context.h"
#include "multirv.h"
#include "params.h"
#include "bn.h"
#include "ctbndyn.h"
#include "em.h"
#include "utils.h"
#include <iostream>
#include "ensurectbn.h"
#include <cstdarg>

using namespace std;
using namespace ctbn;


// this file is used to create an example CTBN model and a partial trajectory sampled from the model
// the model, the evidence, and the true results of queries will be saved to a file "queryinput.data"
// "queryinput.data" will be feeded as an input to:
// exactquery.cc, importancequery.cc and gibbsquery.cc
// the above three .cc file will do the queries, and compare the results with the true ones in "queryinput.data"
// the queries will be specified in these three files

// In case you want to change the example model, you can modify this file and rewrite the output to "queryinput.data"
// and then do the exactquery.cc/importancequery.cc/gibbsquery.cc 
/*
 * Network of variable A-G (with non-consecutive numbers to test code)
 *
 * D<-\
 * |  |
 * v  |
 * A->C<-B
 *    |
 *    v
 * F->E->G
 *
 * Initial BN looks like
 * D
 * |
 * v
 * A->C<-B
 *    |
 *    v
 * F->E->G
 *
 */

const int avar=11,bvar=14,cvar=32,dvar=33,evar=34,fvar=45,gvar=50;
//const int avar=0,bvar=1,cvar=2,dvar=3,evar=4,fvar=5,gvar=6;

Instantiation assign(const Context &c, ...) {
	Instantiation ret(c);
	va_list args;
	va_start(args,c);
	int vnum;
	while((vnum=va_arg(args,int))>=0) {
		int val = va_arg(args,int);
		ret.SetVal(vnum,val);
	}
	va_end(args);
	return ret;
}

vectr v(int nel, ...) {
	vectr ret(nel);
	va_list args;
	va_start(args,nel);
	for(int i=0;i<nel;i++)
		ret[i] = va_arg(args,double);
	va_end(args);
	return ret;
}

int main (int argc, char **argv) {
	Context A; A.AddVar(avar,2);
	Context B; B.AddVar(bvar,3);
	Context C; C.AddVar(cvar,3);
	Context D; D.AddVar(dvar,2);
	Context E; E.AddVar(evar,4);
	Context F; F.AddVar(fvar,2);
	Context G; G.AddVar(gvar,2);
	Context NullC;

	MarkovDyn ADyn(A,D), BDyn(B,NullC), CDyn(C,A+B), DDyn(D,C),
		EDyn(E,C+F), FDyn(F,NullC), GDyn(G,E);

	MultiRV Ap0(A,D), Bp0(B,NullC), Cp0(C,A+B), Dp0(D,NullC),
		Ep0(E,C+F), Fp0(F,NullC), Gp0(G,E);

	// below works as long as compiled with Eigen matrices
	ADyn(assign(D,dvar,0,-1))->Intensity() << -2,  2,
	                                           1, -1;
	ADyn(assign(D,dvar,1,-1))->Intensity() << -0.5, 0.5,
	                                           3,  -3;

	BDyn(0)->Intensity() << -5,   2,    3,
	                         1,  -2,    1,
	                         0.5, 0.5, -1;

	CDyn(assign(A+B,avar,0,bvar,0,-1))->Intensity() << -2,    0.5,  1.5,
	                                                    0.2, -0.3,  0.1,
	                                                    0.1,  0.3, -0.4;
	CDyn(assign(A+B,avar,1,bvar,0,-1))->Intensity() << -3,    0.5,  2.5,
	                                                    0.2, -1.3,  1.1,
	                                                    1.1,  0.3, -1.4;
	CDyn(assign(A+B,avar,0,bvar,1,-1))->Intensity() << -0.5,  0.2,  0.3,
	                                                    0.1, -0.2,  0.1,
	                                                    .1,  0.1, -0.2;
	CDyn(assign(A+B,avar,1,bvar,1,-1))->Intensity() << -1.5,  0.2,  1.3,
	                                                    0.1, -1.2,  1.1,
	                                                    1.1,  0.1, -1.2;
	CDyn(assign(A+B,avar,0,bvar,2,-1))->Intensity() << -1,    0.5,  0.5,
	                                                    0.1, -0.2,  0.1,
	                                                    1.0,  1.0, -2.0;
	CDyn(assign(A+B,avar,1,bvar,2,-1))->Intensity() << -2,    0.5,  1.5,
	                                                    0.1, -1.2,  1.1,
	                                                    2.0,  1.0, -3.0;

	DDyn(assign(C,cvar,0,-1))->Intensity() << -4,    4,
	                                           0.1, -0.1;
	DDyn(assign(C,cvar,1,-1))->Intensity() << -2,    2,
	                                           1.5, -1.5;
	DDyn(assign(C,cvar,2,-1))->Intensity() << -0.1,  0.1,
	                                           4.5, -4.5;

	EDyn(assign(C+F,cvar,0,fvar,0,-1))->Intensity() << -0.3,  0.1,  0.1,  0.1,
	                                                    0.1, -0.3,  0.1,  0.1,
	                                                    0.1,  0.1, -0.3,  0.1,
	                                                    0.1,  0.1,  0.1, -0.3;
	EDyn(assign(C+F,cvar,1,fvar,0,-1))->Intensity() << -0.7,  0.3,  0.3,  0.1,
	                                                    0.1, -0.6,  0.4,  0.1,
	                                                    0.1,  0.4, -0.6,  0.1,
	                                                    0.1,  0.3,  0.3, -0.7;
	EDyn(assign(C+F,cvar,2,fvar,0,-1))->Intensity() << -0.6,  0.1,  0.1,  0.4,
	                                                    0.3, -0.7,  0.1,  0.3,
	                                                    0.3,  0.1, -0.7,  0.3,
	                                                    0.4,  0.1,  0.1, -0.6;
	EDyn(assign(C+F,cvar,0,fvar,1,-1))->Intensity() << -0.6,  0.1,  0.1,  0.4,
	                                                    0.1, -0.6,  0.1,  0.4,
	                                                    0.1,  0.1, -0.6,  0.4,
	                                                    0.1,  0.1,  0.1, -0.3;
	EDyn(assign(C+F,cvar,1,fvar,1,-1))->Intensity() << -0.7,  0.3,  0.3,  0.1,
	                                                    0.4, -0.9,  0.4,  0.1,
	                                                    0.4,  0.4, -0.9,  0.1,
	                                                    0.4,  0.3,  0.3, -1.0;
	EDyn(assign(C+F,cvar,2,fvar,1,-1))->Intensity() << -1.0,  0.3,  0.3,  0.4,
	                                                    0.3, -0.7,  0.1,  0.3,
	                                                    0.3,  0.1, -0.7,  0.3,
	                                                    0.4,  0.3,  0.3, -1.0;

	FDyn(0)->Intensity() << -0.05, 0.05,
	                         0.07, -0.07;

	GDyn(assign(E,evar,0,-1))->Intensity() << -3,    3,
	                                           0.5, -0.5;
	GDyn(assign(E,evar,1,-1))->Intensity() << -2,   2,
	                                           1,  -1;
	GDyn(assign(E,evar,2,-1))->Intensity() << -1,   1,
	                                           2,  -2;
	GDyn(assign(E,evar,3,-1))->Intensity() << -0.5, 0.5,
	                                           3,  -3;

	Ap0[assign(D,dvar,0,-1)].SetDist(v(2,0.9,0.1));
	Ap0[assign(D,dvar,1,-1)].SetDist(v(2,0.1,0.9));

	Bp0[0].SetDist(v(3,0.2,0.6,0.2));
	
	Cp0[assign(A+B,avar,0,bvar,0,-1)].SetDist(v(3,0.8,0.1,0.1));
	Cp0[assign(A+B,avar,1,bvar,0,-1)].SetDist(v(3,0.2,0.5,0.3));
	Cp0[assign(A+B,avar,0,bvar,1,-1)].SetDist(v(3,0.3,0.6,0.1));
	Cp0[assign(A+B,avar,1,bvar,1,-1)].SetDist(v(3,0.1,0.5,0.4));
	Cp0[assign(A+B,avar,0,bvar,2,-1)].SetDist(v(3,0.2,0.3,0.5));
	Cp0[assign(A+B,avar,1,bvar,2,-1)].SetDist(v(3,0.1,0.1,0.8));

	Dp0[0].SetDist(v(2,0.7,0.3));

	Ep0[assign(C+F,cvar,0,fvar,0,-1)].SetDist(v(4,0.7,0.1,0.1,0.1));
	Ep0[assign(C+F,cvar,1,fvar,0,-1)].SetDist(v(4,0.3,0.4,0.2,0.1));
	Ep0[assign(C+F,cvar,2,fvar,0,-1)].SetDist(v(4,0.1,0.2,0.3,0.4));
	Ep0[assign(C+F,cvar,0,fvar,1,-1)].SetDist(v(4,0.1,0.2,0.5,0.2));
	Ep0[assign(C+F,cvar,1,fvar,1,-1)].SetDist(v(4,0.1,0.2,0.4,0.3));
	Ep0[assign(C+F,cvar,2,fvar,1,-1)].SetDist(v(4,0.1,0.1,0.1,0.7));

	Fp0[0].SetDist(v(2,0.2,0.8));

	Gp0[assign(E,evar,0,-1)].SetDist(v(2,0.9,0.1));
	Gp0[assign(E,evar,1,-1)].SetDist(v(2,0.7,0.3));
	Gp0[assign(E,evar,2,-1)].SetDist(v(2,0.3,0.7));
	Gp0[assign(E,evar,3,-1)].SetDist(v(2,0.1,0.9));
	
	// set up CTBN:
	CTBNDyn ctbndyn;
	ctbndyn.AddNode(ADyn.Clone());
	ctbndyn.AddNode(BDyn.Clone());
	ctbndyn.AddNode(CDyn.Clone());
	ctbndyn.AddNode(DDyn.Clone());
	ctbndyn.AddNode(EDyn.Clone());
	ctbndyn.AddNode(FDyn.Clone());
	//ctbndyn.AddNode(GDyn.Clone());

	BN bn;
	bn.AddNode(Ap0.Clone());
	bn.AddNode(Bp0.Clone());
	bn.AddNode(Cp0.Clone());
	bn.AddNode(Dp0.Clone());
	bn.AddNode(Ep0.Clone());
	bn.AddNode(Fp0.Clone());
	//bn.AddNode(Gp0.Clone());

	Markov ctbn(bn.Clone(),ctbndyn.Clone());

	// set trajectory to something we know the answer (via Matlab script)
	Trajectory evid;
	evid.SetBeginTime(0.0);
	evid.SetEndTime(2.0);
	evid.AddTransition(avar,0.0,-1);
	evid.AddTransition(bvar,0.0,2);
	evid.AddTransition(cvar,0.0,1);
	evid.AddTransition(dvar,0.0,-1);
	evid.AddTransition(evar,0.0,0);
	evid.AddTransition(fvar,0.0,-1);

	evid.AddTransition(evar,0.05,1);
	evid.AddTransition(cvar,0.3,0);
	evid.AddTransition(bvar,0.4,1);
	evid.AddTransition(cvar,0.5,1);
	evid.AddTransition(avar,0.7,1);
	evid.AddTransition(bvar,1.0,-1);
	evid.AddTransition(fvar,1.0,0);
	evid.AddTransition(avar,1.1,0);
	evid.AddTransition(dvar,1.5,1);
	evid.AddTransition(dvar,1.6,0);
	
	ctbn.Save(cout);
	evid.Save(cout);

	// "true" values (those at the end of the lines)
	// were calculated using octave from calcquery.m & queries.txt
	// (queries.txt is generated from the output of this program
	//  using genquerytxt) with a delta=0.00001
	// query filtered prob of D=0 at time 0.9
	cout << "0 " << dvar << " 0 0.9 0.546618" << endl;
	// query smoothed prob of D=0 at time 0.9
	cout << "1 " << dvar << " 0 0.9 0.432647" << endl;
	// query expected total amount of time D spent in state 1
	cout << "2 " << dvar << " 1 0.940300" << endl;
	// query expected total # of trans of D from 2 to 0
	cout << "3 " << dvar << " 1 0 1.76552" << endl;
	return 0;
}
