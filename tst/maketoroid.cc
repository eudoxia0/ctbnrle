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
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>

#include "bn.h"
#include "em.h"
#include "context.h"
#include "ctbn.h"
#include "ctbndyn.h"
#include "markov.h"
#include "markovdyn.h"
#include "multirv.h"
#include "trajectory.h"

#include "params.h"

using namespace std;
using namespace ctbn;

double tao = 1;
double beta = 0.5;

void FillIntensityIsing(MarkovDyn &md, Instantiation &i, int pa1, int pa2){
	int pv1, pv2;
	for (int v1=0; v1<2; v1++){
		i.SetVal(pa1, v1);
		for (int v2=0; v2<2; v2++){
			i.SetVal(pa2, v2);
			pv1 = (v1==0)?-1:1;
			pv2 = (v2==0)?-1:1;

			md(i)->Intensity()[1][0] = tao*pow((1+exp(2*beta*(pv1+pv2))), -1.0); 
			md(i)->Intensity()[0][1] = tao*pow((1+exp(-2*beta*(pv1+pv2))), -1.0);	
			md(i)->Intensity()[0][0] = -1*(md(i)->Intensity()[0][1]);
			md(i)->Intensity()[1][1] = -1*(md(i)->Intensity()[1][0]);
		}
	}

}

void FillIntensityIsingOp(MarkovDyn &md, Instantiation &i, int pa1, int pa2){
	int pv1, pv2;
	for (int v1=0; v1<2; v1++){
		i.SetVal(pa1, v1);
		for (int v2=0; v2<2; v2++){
			i.SetVal(pa2, v2);
			pv1 = (v1==0)?-1:1;
			pv2 = (v2==0)?-1:1;

			md(i)->Intensity()[1][0] = tao*pow((1+exp(-2*beta*(pv1+pv2))), -1.0); 
			md(i)->Intensity()[0][1] = tao*pow((1+exp(2*beta*(pv1+pv2))), -1.0);	
			md(i)->Intensity()[0][0] = -1*(md(i)->Intensity()[0][1]);
			md(i)->Intensity()[1][1] = -1*(md(i)->Intensity()[1][0]);
		}
	}
}

void FillIntensityIsing(MarkovDyn &md, Instantiation &i, int pa1){

	int pv1;
	for (int v1=0; v1<2; v1++){
		i.SetVal(pa1, v1);
			pv1 = (v1==0)?-1:1;
			
			md(i)->Intensity()[1][0] = tao*pow((1+exp(2*beta*(pv1))), -1.0); 
			md(i)->Intensity()[0][1] = tao*pow((1+exp(-2*beta*(pv1))), -1.0);	
			md(i)->Intensity()[0][0] = -1*(md(i)->Intensity()[0][1]);
			md(i)->Intensity()[1][1] = -1*(md(i)->Intensity()[1][0]);

	}

}


void FillIntensityIsing(MarkovDyn &md){
	double pv1 = 0;
		//xo and pv = {-1, +1}
		//tao*(1+exp(-2*xo*beta*sum(pv)))^-1
		md(0)->Intensity()[1][0] = tao*pow((1+exp(2*beta*(pv1))), -1.0); 
		md(0)->Intensity()[0][1] = tao*pow((1+exp(-2*beta*(pv1))), -1.0);	
		md(0)->Intensity()[0][0] = -1*(md(0)->Intensity()[0][1]);
		md(0)->Intensity()[1][1] = -1*(md(0)->Intensity()[1][0]);
		
	

}

int main (int argc, char **argv) {
	if (argc > 2){
		tao = atof(argv[1]);
		beta = atof(argv[2]);
	}

	Context x0, x1, x2, x3, x4,x5,x6,x7,x8,x9,x10, x11;
	Context Null;
	x0.AddVar(0,2);
	x1.AddVar(1,2);
	x2.AddVar(2,2);
	x3.AddVar(3,2);
	x4.AddVar(4,2);
	x5.AddVar(5,2);
	x6.AddVar(6,2);
	x7.AddVar(7,2);
	x8.AddVar(8,2);
	x9.AddVar(9,2);
	x10.AddVar(10,2);
	x11.AddVar(11,2);

	MarkovDyn PX0(x0,Context(x3,x8)); // process of X
	MarkovDyn PX1(x1,Context(x0,x9)); 
	MarkovDyn PX2(x2,Context(x1,x10));
	MarkovDyn PX3(x3,Context(x2,x11));
	MarkovDyn PX4(x4,Context(x0,x7));
	MarkovDyn PX5(x5,Context(x1,x4)); 
	MarkovDyn PX6(x6,Context(x2, x5));
	MarkovDyn PX7(x7,Context(x3, x6)); 
	MarkovDyn PX8(x8,Context(x4, x11));
	MarkovDyn PX9(x9,Context(x5, x8));
	MarkovDyn PX10(x10,Context(x6, x9));
	MarkovDyn PX11(x11,Context(x7, x10)); 
	
	Instantiation i0(Context(x3,x8));
	FillIntensityIsing(PX0, i0, 3, 8);
	
	Instantiation i1(Context(x0,x9));
	FillIntensityIsing(PX1, i1, 0, 9);
	
	Instantiation i2(Context(x1,x10));
	FillIntensityIsing(PX2, i2, 1, 10);
	
	Instantiation i3(Context(x2, x11));
	FillIntensityIsing(PX3, i3, 2, 11);
	
	Instantiation i4(Context(x0, x7));
	FillIntensityIsing(PX4, i4, 0, 7);
	
	Instantiation i5(Context(x1, x4));
	FillIntensityIsing(PX5, i5, 1, 4);
	
	Instantiation i6(Context(x2, x5));
	FillIntensityIsing(PX6, i6, 2, 5);
	
	Instantiation i7(Context(x3, x6));
	FillIntensityIsing(PX7, i7, 3, 6);
	
	Instantiation i8(Context(x4, x11));
	FillIntensityIsing(PX8, i8, 4, 11);
	
	Instantiation i9(Context(x5, x8));
	FillIntensityIsing(PX9, i9, 5, 8);
	
	Instantiation i10(Context(x6, x9));
	FillIntensityIsing(PX10, i10, 6, 9);
	
	Instantiation i11(Context(x7, x10));
	FillIntensityIsing(PX11, i11, 7, 10);
	
	

	//Independent initial dist.
	MultiRV P0x0(x0,Null);
	MultiRV P0x1(x1,Null);
	MultiRV P0x2(x2,Null);
	MultiRV P0x3(x3,Null);
	MultiRV P0x4(x4,Null);
	MultiRV P0x5(x5,Null);
	MultiRV P0x6(x6,Null);
	MultiRV P0x7(x7,Null);
	MultiRV P0x8(x8,Null);
	MultiRV P0x9(x9,Null);
	MultiRV P0x10(x10,Null);
	MultiRV P0x11(x11,Null);
	
	vectr p0(2,0.0);
	p0[0] = 1;  p0[1] = 0;
	P0x0[0].SetDist(p0);

	p0[0] = 1; p0[1] = 0;
	P0x1[0].SetDist(p0);
	
	p0[0] = 1; p0[1] = 0;
	P0x2[0].SetDist(p0);

	p0[0] = 1; p0[1] = 0;
	P0x3[0].SetDist(p0);
	
	p0[0] = 0; p0[1] = 1;
	P0x4[0].SetDist(p0);

	p0[0] = 0; p0[1] = 1;
	P0x5[0].SetDist(p0);

	p0[0] = 1; p0[1] = 0;
	P0x6[0].SetDist(p0);

	p0[0] = 1; p0[1] = 0;
	P0x7[0].SetDist(p0);

	p0[0] = 1; p0[1] = 0;
	P0x8[0].SetDist(p0);

	p0[0] = 1; p0[1] = 0;
	P0x9[0].SetDist(p0);	
	
	p0[0] = 1; p0[1] = 0;
	P0x10[0].SetDist(p0);

	p0[0] = 1; p0[1] = 0;
	P0x11[0].SetDist(p0);	

	// set up CTBN:
	CTBNDyn ctbndyn;
	ctbndyn.AddNode(PX0.Clone());
	ctbndyn.AddNode(PX1.Clone());
	ctbndyn.AddNode(PX2.Clone());
	ctbndyn.AddNode(PX3.Clone());
	ctbndyn.AddNode(PX4.Clone());
	ctbndyn.AddNode(PX5.Clone());
	ctbndyn.AddNode(PX6.Clone());
	ctbndyn.AddNode(PX7.Clone());
	ctbndyn.AddNode(PX8.Clone());
	ctbndyn.AddNode(PX9.Clone());
	ctbndyn.AddNode(PX10.Clone());
	ctbndyn.AddNode(PX11.Clone());

	BN bn;
	bn.AddNode(P0x0.Clone());
	bn.AddNode(P0x1.Clone());
	bn.AddNode(P0x2.Clone());
	bn.AddNode(P0x3.Clone());
	bn.AddNode(P0x4.Clone());
	bn.AddNode(P0x5.Clone());
	bn.AddNode(P0x6.Clone());
	bn.AddNode(P0x7.Clone());
	bn.AddNode(P0x8.Clone());
	bn.AddNode(P0x9.Clone());
	bn.AddNode(P0x10.Clone());
	bn.AddNode(P0x11.Clone());

	Markov ctbn(bn.Clone(),ctbndyn.Clone());
	Context context = ctbndyn.Domain() + ctbndyn.CondDomain();
	ctbn.Save(cout);

	// trajectory
	Trajectory tr;
	tr.SetBeginTime(0.0);
	tr.SetEndTime(1.0);
	tr.AddTransition(0,0.5,1);
	tr.AddTransition(1,0.5,0);
	tr.Save(cout);

	// query
	cout << "0 2 0 1.0 0.8983018202" << endl;

	return 0;
}
