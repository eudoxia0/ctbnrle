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

#include "expmcache.h"
#include "params.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

namespace ctbn {
void vexpmt(vectr &a, const matrix &Q, double t, double eps, double h, 
		ContFunction<vectr> &output);
void vexpmt(matrix &a, const matrix &Q, double t, double eps, double h, 
		ContFunction<matrix> &output);

ExpMCache::ExpMCache(const matrix &Q, double maxt, double eps)
			//: expm(Q.getm())
			{
	int m = Q.getm();
/* version 1:
	for(int i=0;i<m;i++) {
		vectr v(m,0.0);
		v[i] = 1.0;
		ctbn::vexpmt(v,Q,maxt,eps,-1.0,expm[i]);
	}
*/
/* version 2:
 */
	matrix eye(m,m,0.0);
	for(int i=0;i<m;i++) eye[i][i] = 1.0;
	ctbn::vexpmt(eye,Q,maxt,eps,-1.0,expm);
}

void ExpMCache::vexpmt(vectr &a, double t) const {
/* version 1:
	int m = a.getm();
	if (m==0) return;
	vectr ret = expm[0].GetVal(t)*a[0];
	for(int i=1;i<m;i++)
		ret += expm[i].GetVal(t)*a[i];
	a = ret;
*/
/* version 2:
 */
	a = a*expm.GetVal(t);
}

// Perhaps below should be in rk.h/.cc instead?

// used by all RKF methods, placed here for easy initialization
static double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
	b21=0.2,
	b31=3.0/40.0, b32=9.0/40.0,
	b41=0.3, b42=-0.9, b43=1.2,
	b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0,
	b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0,
	b64=44275.0/110592.0, b65=253.0/4096.0,
	c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
	dc5=-277.00/14336.0;
static double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
	dc4=c4-13525.0/55296.0, dc6=c6-0.25;

#define RKEPS 1e-8

void vexpmt(vectr &a, const matrix &Q, double t, double eps, double h, 
		ContFunction<vectr> &output) {
	if (t<=0.0) return;

	static double defaulteps = ParamDouble("RKEPS",RKEPS);
	if (eps<0.0) eps = defaulteps;

	int p = Q.getm();

	if (h<=0) {
		double dd = Q.diag().absmax();
		if (dd>0) h = 0.1/Q.diag().absmax();
		else h = 0.1;
	}

	vectr ka1(p),ka2(p),ka3(p),ka4(p),ka5(p),ka6(p);

	double currt = 0;
	while(currt<t) {
		output.AddVal(currt,a);
		if (currt+h>=t) h = t-currt; // don't overshoot!

		double err;
		for(;;) { // find stepsize
			// take 5th-order step:

			vectr temp;

			Q.leftmult(a,ka1);
			ka1 *= h;
			//ka1 = (a*Q)*h;

			temp = a;
			temp.add(ka1,b21);
			//temp = ka1*b21+a;
			Q.leftmult(temp,ka2);
			ka2 *= h;
			//ka2 = (temp*Q)*h;

			temp = a;
			temp.add(ka1,b31);
			temp.add(ka2,b32);
			//temp = ka1*b31+ka2*b32+a;
			Q.leftmult(temp,ka3);
			ka3 *= h;
			//ka3 = (temp*Q)*h;

			temp = a;
			temp.add(ka1,b41);
			temp.add(ka2,b42);
			temp.add(ka3,b43);
			//temp = ka1*b41+ka2*b42+ka3*b43+a;
			Q.leftmult(temp,ka4);
			ka4 *= h;
			//ka4 = (temp*Q)*h;

			temp = a;
			temp.add(ka1,b51);
			temp.add(ka2,b52);
			temp.add(ka3,b53);
			temp.add(ka4,b54);
			//temp = ka1*b51+ka2*b52+ka3*b53+ka4*b54+a;
			Q.leftmult(temp,ka5);
			ka5 *= h;
			//ka5 = (temp*Q)*h;

			temp = a;
			temp.add(ka1,b61);
			temp.add(ka2,b62);
			temp.add(ka3,b63);
			temp.add(ka4,b64);
			temp.add(ka5,b65);
			//temp = ka1*b61+ka2*b62+ka3*b63+ka4*b64+ka5*b65+a;
			Q.leftmult(temp,ka6);
			ka6 *= h;
			//ka6 = (temp*Q)*h;

			vectr erra = ka1;
			erra *= dc1;
			erra.add(ka3,dc3);
			erra.add(ka4,dc4);
			erra.add(ka5,dc5);
			erra.add(ka6,dc6);
			//vectr erra = dc1*ka1+dc3*ka3+dc4*ka4+dc5*ka5+dc6*ka6;

			// check error:
			/*err = erra.absmax();
			  err /= eps;
			  if (err <= 1.0) break;*/
			double amin = a.absmin();
			if(amin<1.0) amin = 1.0;
			err = erra.absmax()/amin;
			err /= eps;
			if(err<=1.0) break;

			// adjust step:
			double htemp = 0.9*h*pow(err,-0.25);
			if (htemp>h/10) h = htemp;
			else h /= 10;
			if (currt+h == currt) return; // should throw exception?
		}
		currt += h;
		if (err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		a.add(ka1,c1);
		a.add(ka3,c3);
		a.add(ka4,c4);
		a.add(ka6,c6);
		//a += c1*ka1+c3*ka3+c4*ka4+c6*ka6;
	}
	output.AddVal(currt,a);
}

// could probably be made a template function (matrix or vectr)
// if some way of initializing ka1, ka2, ... could be arranged
void vexpmt(matrix &a, const matrix &Q, double t, double eps, double h, 
		ContFunction<matrix> &output) {
	if (t<=0.0) return;

	static double defaulteps = ParamDouble("RKEPS",RKEPS);
	if (eps<0.0) eps = defaulteps;

	int p = Q.getm();

	if (h<=0) {
		double dd = Q.diag().absmax();
		if (dd>0) h = 0.1/Q.diag().absmax();
		else h = 0.1;
	}

	matrix ka1(p,p),ka2(p,p),ka3(p,p),ka4(p,p),ka5(p,p),ka6(p,p);

	double currt = 0;
	while(currt<t) {
		output.AddVal(currt,a);
		if (currt+h>=t) h = t-currt; // don't overshoot!

		double err;
		for(;;) { // find stepsize
			// take 5th-order step:

			matrix temp;

			//Q.leftmult(a,ka1);
			ka1 = a*Q;
			ka1 *= h;
			//ka1 = (a*Q)*h;

			temp = a;
			temp.add(ka1,b21);
			//temp = ka1*b21+a;
			//Q.leftmult(temp,ka2);
			ka2 = temp*Q;
			ka2 *= h;
			//ka2 = (temp*Q)*h;

			temp = a;
			temp.add(ka1,b31);
			temp.add(ka2,b32);
			//temp = ka1*b31+ka2*b32+a;
			//Q.leftmult(temp,ka3);
			ka3 = temp*Q;
			ka3 *= h;
			//ka3 = (temp*Q)*h;

			temp = a;
			temp.add(ka1,b41);
			temp.add(ka2,b42);
			temp.add(ka3,b43);
			//temp = ka1*b41+ka2*b42+ka3*b43+a;
			//Q.leftmult(temp,ka4);
			ka4 = temp*Q;
			ka4 *= h;
			//ka4 = (temp*Q)*h;

			temp = a;
			temp.add(ka1,b51);
			temp.add(ka2,b52);
			temp.add(ka3,b53);
			temp.add(ka4,b54);
			//temp = ka1*b51+ka2*b52+ka3*b53+ka4*b54+a;
			//Q.leftmult(temp,ka5);
			ka5 = temp*Q;
			ka5 *= h;
			//ka5 = (temp*Q)*h;

			temp = a;
			temp.add(ka1,b61);
			temp.add(ka2,b62);
			temp.add(ka3,b63);
			temp.add(ka4,b64);
			temp.add(ka5,b65);
			//temp = ka1*b61+ka2*b62+ka3*b63+ka4*b64+ka5*b65+a;
			//Q.leftmult(temp,ka6);
			ka6 = temp*Q;
			ka6 *= h;
			//ka6 = (temp*Q)*h;

			matrix erra = ka1;
			erra *= dc1;
			erra.add(ka3,dc3);
			erra.add(ka4,dc4);
			erra.add(ka5,dc5);
			erra.add(ka6,dc6);
			//matrix erra = dc1*ka1+dc3*ka3+dc4*ka4+dc5*ka5+dc6*ka6;

			// check error:
			/*err = erra.absmax();
			  err /= eps;
			  if (err <= 1.0) break;*/
			double amin = a.absmin();
			if(amin<1.0) amin = 1.0;
			err = erra.absmax()/amin;
			err /= eps;
			if(err<=1.0) break;

			// adjust step:
			double htemp = 0.9*h*pow(err,-0.25);
			if (htemp>h/10) h = htemp;
			else h /= 10;
			if (currt+h == currt) return; // should throw exception?
		}
		currt += h;
		if (err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		a.add(ka1,c1);
		a.add(ka3,c3);
		a.add(ka4,c4);
		a.add(ka6,c6);
		//a += c1*ka1+c3*ka3+c4*ka4+c6*ka6;
	}
	output.AddVal(currt,a);
}

}
