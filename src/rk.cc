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

#include "rk.h"
#include "matrix.h"
#include "params.h"
#include "utils.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>


using namespace std;

namespace ctbn {
const double MU_END_PRECISION = 1e-2;

#define USE_SPARSE_SS

// This version has a calculation cache
// comment out the USECACHE definition to remove the caching
//#define USECACHE
#define MAXCACHE 1000


#ifdef USECACHE
static int reuse = 0;
static int tries = 0;
#endif

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
//static double defaulteps = ParamDouble("RKEPS",RKEPS);

// computes a * exp(Qt), via integration with desired error of eps
// (h, if non-negative, gives a starting estimate of step size)
// return value replaces a
// This is a version based on adaptive step sizes (and therefore
// a fifth-order Runge-Kutta method with an embedded 4th order method
// -- see Numerical Recipes in C)

double vexpmt(vectr &a, const matrix &Q, double t, double eps, double h, 
		vector<double> *timesteps) {
//double vexpmt(vectr &a, const matrix &Q, double t, double eps, double h) {
	double ret = 0.0;
	if (t<=0.0) return ret;

	static double defaulteps = ParamDouble("RKEPS",RKEPS);
	if (eps<0.0) eps = defaulteps;

#ifdef USECACHE
	static vector<vectr> cachea(MAXCACHE);
	static vector<matrix> cacheQ(MAXCACHE);
	static vector<double> cachet(MAXCACHE);
	static vector<vectr> cacheret(MAXCACHE);
	static vector<double> cachellh(MAXCACHE);
	static int cachepos = 0;
	tries++;
	//if (tries%1000==0) cout << "cache: " << reuse << '/' << tries << endl;
	for(int i=0;i<cachea.size();i++)
		if (cachea[i]==a && cacheQ[i]==Q && cachet[i]==t) {
			reuse++;
			a = cacheret[i];
			return cachellh[i];
		}
	cachea[cachepos] = a;
	cacheQ[cachepos] = Q;
	cachet[cachepos] = t;
#endif

	int p = Q.getm();

	if (h<=0) {
		double dd = Q.diag().absmax();
		if (dd>0) h = 0.1/Q.diag().absmax();
		else h = 0.1;
	}

	vectr ka1(p),ka2(p),ka3(p),ka4(p),ka5(p),ka6(p);

	//static double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
	/*
	static double
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
		*/
		
	if(timesteps)	timesteps->clear();
	
	double currt = 0;
	while(currt<t) {
                if (currt+h>=t) h = t-currt; // don't overshoot!
        
        if(timesteps) timesteps->push_back(currt);

		double err;
		for(;;) { // find stepsize
			// take 5th-order step:
			//iter++;
                    
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
			if (currt+h == currt) {
                            /*cerr << "step size underflow in rk2" << endl;
				cerr << "in alpha pass:" << endl;
				cerr << Q << endl;
				cerr << a << endl;
				cerr << t << ' ' << currt << endl;
				exit(1);*/
                            return ret;
			}
		}
		currt += h;
		if (err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		a.add(ka1,c1);
		a.add(ka3,c3);
		a.add(ka4,c4);
		a.add(ka6,c6);
		//a += c1*ka1+c3*ka3+c4*ka4+c6*ka6;
		double ss = a.sum();
		if (ss>0.0) {
			ret += log(ss);
			a /= ss;
		}	
		//ret += log(a.normalize());
	}
#ifdef USECACHE
	cacheret[cachepos] = a;
	cachellh[cachepos] = ret;
	cachepos++;
	//if (cachepos>MAXCACHE) { cerr << "cache overflow" << endl; cachepos=0;}
	if (cachepos>=MAXCACHE) cachepos=0;
#endif
	return ret;
}


// computes exp(Qt) * b, via integration with desired error of eps
// (h, if non-negative, gives a starting estimate of step size)
// return value replaces b
// This is a version based on adaptive step sizes (and therefore
// a fifth-order Runge-Kutta method with an embedded 4th order-method
// -- see Numerical Recipes in C)
double expmtv(vectr &b, const matrix &Q, double t, double eps, double h) {
    double ret = 0.0;
	if (t<=0.0) return ret;

	static double defaulteps = ParamDouble("RKEPS",RKEPS);
	if (eps<0.0) eps = defaulteps;

#ifdef USECACHE
	static vector<vectr> cacheb(MAXCACHE);
	static vector<matrix> cacheQ(MAXCACHE);
	static vector<double> cachet(MAXCACHE);
	static vector<vectr> cacheret(MAXCACHE);
	static vector<double> cachellh(MAXCACHE);
	static int cachepos = 0;
	tries++;
	//if (tries%1000==0) cout << "cache: " << reuse << '/' << tries << endl;
	for(int i=0;i<cacheb.size();i++)
		if (cacheb[i]==b && cacheQ[i]==Q && cachet[i]==t) {
			reuse++;
			b = cacheret[i];
			return cachellh[i];
		}
	cacheb[cachepos] = b;
	cacheQ[cachepos] = Q;
	cachet[cachepos] = t;
#endif

	int p = Q.getm();

	if (h<=0) {
		double dd = Q.diag().absmax();
		if (dd>0) h = 0.1/Q.diag().absmax();
		else h = 0.1;
	}

	vectr kb1(p),kb2(p),kb3(p),kb4(p),kb5(p),kb6(p);

	//static double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
	/*
	static double 
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
		*/

	double currt = 0;
	while(currt<t) {
		if (currt+h>=t) h = t-currt; // don't overshoot!

		double err;
		for(;;) { // find stepsize
			// take 5th-order step:
			vectr temp;

			Q.rightmult(b,kb1);
			kb1 *= h;
			//kb1 = (Q*b)*h;

			temp = b;
			temp.add(kb1,b21);
			//temp = kb1*b21+b;
			Q.rightmult(temp,kb2);
			kb2 *= h;
			//kb2 = (Q*temp)*h;

			temp = b;
			temp.add(kb1,b31);
			temp.add(kb2,b32);
			//temp = kb1*b31+kb2*b32+b;
			Q.rightmult(temp,kb3);
			kb3 *= h;
			//kb3 = (Q*temp)*h;

			temp = b;
			temp.add(kb1,b41);
			temp.add(kb2,b42);
			temp.add(kb3,b43);
			//temp = kb1*b41+kb2*b42+kb3*b43+b;
			Q.rightmult(temp,kb4);
			kb4 *= h;
			//kb4 = (Q*temp)*h;

			temp = b;
			temp.add(kb1,b51);
			temp.add(kb2,b52);
			temp.add(kb3,b53);
			temp.add(kb4,b54);
			//temp = kb1*b51+kb2*b52+kb3*b53+kb4*b54+b;
			Q.rightmult(temp,kb5);
			kb5 *= h;
			//kb5 = (Q*temp)*h;


			temp = b;
			temp.add(kb1,b61);
			temp.add(kb2,b62);
			temp.add(kb3,b63);
			temp.add(kb4,b64);
			temp.add(kb5,b65);
			//temp = kb1*b61+kb2*b62+kb3*b63+kb4*b64+kb5*b65+b;
			Q.rightmult(temp,kb6);
			kb6 *= h;
			//kb6 = (Q*temp)*h;

			vectr errb = kb1;
			errb *= dc1;
			errb.add(kb3,dc3);
			errb.add(kb4,dc4);
			errb.add(kb5,dc5);
			errb.add(kb6,dc6);
			//vectr errb = dc1*kb1+dc3*kb3+dc4*kb4+dc5*kb5+dc6*kb6;


			// check error:
			/*err = errb.absmax();
			err /= eps;
			if (err <= 1.0) break;*/
                        double bmin = b.absmin();
                        if(bmin<1.0) bmin = 1.0;
                        err = errb.absmax()/bmin;
                        err /= eps;
                        if(err<=1.0) break;

			// adjust step:
			double htemp = 0.9*h*pow(err,-0.25);
			if (htemp>h/10) h = htemp;
			else h /= 10;
			if (currt+h == currt) {
                            return ret;
                            /*cerr << "step size underflow in rk2" << endl;
				cerr << "in beta pass:" << endl;
				cerr << Q << endl;
				cerr << b << endl;
				cerr << t << ' ' << currt << endl;
				exit(1);*/
			}
		}
		currt += h;
		if (err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		b.add(kb1,c1);
		b.add(kb3,c3);
		b.add(kb4,c4);
		b.add(kb6,c6);
		//b += c1*kb1+c3*kb3+c4*kb4+c6*kb6;
		double ss = b.sum();
		if (ss>0.0) {
			ret += log(ss);
			b /= ss;
		}
		//ret += log(b.normalize());
	}
#ifdef USECACHE
	cacheret[cachepos] = b;
	cachellh[cachepos] = ret;
	cachepos++;
	//if (cachepos>MAXCACHE) { cerr << "cache overflow" << endl; cachepos=0;}
	if (cachepos>=MAXCACHE) cachepos=0;
#endif
	return ret;
}


// calculates  \int_u=0^t a*exp(Q*u)*e_i e_j*exp(Q*(t-u))*b du
//       [this is a matrix indexed by (i,j)]
// the output is returned on c
// (actually, it returns the transpose...)
// This is a version based on adaptive step sizes (and therefore
// a fifth-order Runge-Kutta method with an embedded 4th order method
// -- see Numerical Recipes in C)
double suffstat(const vectr &a, const vectr &b, matrix &c,
              const matrix &Q, double t, double eps, double h) {
              
    double bgtime = getcputime();
	double ret = 0.0;
	if (t<=0.0) return ret;

	static double defaulteps = ParamDouble("RKEPS",RKEPS);
	if (eps<0.0) eps = defaulteps;

#ifdef USECACHE
	static vector<vectr> cachea(MAXCACHE);
	static vector<vectr> cacheb(MAXCACHE);
	static vector<matrix> cacheQ(MAXCACHE);
	static vector<double> cachet(MAXCACHE);
	static vector<matrix> cacheret(MAXCACHE);
	static vector<double> cachellh(MAXCACHE);
	static int cachepos = 0;
	tries++;
	//if (tries%1000==0) cout << "cache: " << reuse << '/' << tries << endl;
	for(int i=0;i<cacheb.size();i++)
		if (cachea[i]==a && cacheb[i]==b && cacheQ[i]==Q && cachet[i]==t) {
			reuse++;
			c = cacheret[i];
			return cachellh[i];
		}
	cachea[cachepos] = a;
	cacheb[cachepos] = b;
	cacheQ[cachepos] = Q;
	cachet[cachepos] = t;
#endif

	int p = Q.getm();
	
	if (h<=0) {
		double dd = Q.diag().absmax();
		if (dd>0) h = 0.1/Q.diag().absmax();
		else h = 0.1;
	}

	double maxh = h, minh = h;
	//int evals = 0;
	
#ifdef USE_SPARSE_SS
	spmatrix spQ(Q);
#else
	const matrix &spQ = Q;
#endif

	vectr mya(a);
	c = matrix(p,p,0.0);
	vectr ka1(p),ka2(p),ka3(p),ka4(p),ka5(p),ka6(p);
	matrix kc1(p,p),kc2(p,p),kc3(p,p),kc4(p,p),kc5(p,p),kc6(p,p);

	
	
	//static double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
	/*
	static double
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
		*/

	double currt = 0;

	
	while(currt<t) {
		if (h>maxh) maxh = h;
		if (h<minh) minh = h;
		if (currt+h>=t) h = t-currt; // don't overshoot!

		double err;

		for(;;) { 
		
			// find stepsize
			//evals++;
			//if(evals>MAX_INT) return;
			// take 5th-order step:
			vectr temp;

			//ka1 = (mya*Q)*h;
			spQ.leftmult(mya, ka1);
			ka1 *= h;

				
								
			//kc1 = (matrix(b,mya)+ Q*c)*h;
			kc1 = spQ.rightmult(c);
			kc1.add(b,mya);
			kc1 *= h;
				
	
			//temp = ka1*b21+mya;
			temp = mya;
			temp.add(ka1, b21);
			
			//ka2 = (temp*Q)*h;
			spQ.leftmult(temp,ka2);
			ka2 *= h;			

		
						
			//kc2 = (matrix(b,temp)	+ Q*(kc1*b21+c))*h;
			kc2 = c;
			kc2.add(kc1, b21);
			kc2 = spQ.rightmult(kc2);
			kc2.add(b, temp);
			kc2 *= h;

		
			//temp = ka1*b31+ka2*b32+mya;
			temp = mya;
			temp.add(ka1, b31);
			temp.add(ka2, b32);
			
			//ka3 = (temp*Q)*h;
			spQ.leftmult(temp, ka3);
			ka3 *= h;

			
						
			//kc3 = (matrix(b,temp) + Q*(kc1*b31+kc2*b32+c))*h;
			kc3 = c;
			kc3.add(kc1, b31);
			kc3.add(kc2, b32);
			kc3 = spQ.rightmult(kc3);
			kc3.add(b,temp);
			kc3 *= h;
		
			
			//temp = ka1*b41+ka2*b42+ka3*b43+mya;
			temp = mya;
			temp.add(ka1, b41);
			temp.add(ka2, b42);
			temp.add(ka3, b43);
			
			//ka4 = (temp*Q)*h;
			spQ.leftmult(temp, ka4);
			ka4 *= h;

	
			//kc4 = (matrix(b,temp) + Q*(kc1*b41+kc2*b42+kc3*b43+c))*h;
			kc4 = c;
			kc4.add(kc1, b41);
			kc4.add(kc2, b42);
			kc4.add(kc3, b43);
			kc4 = spQ.rightmult(kc4);
			kc4.add(b, temp);
			kc4 *= h;

		
			//temp = ka1*b51+ka2*b52+ka3*b53+ka4*b54+mya;
			temp = mya;
			temp.add(ka1, b51);
			temp.add(ka2, b52);
			temp.add(ka3, b53);
			temp.add(ka4, b54);
			
			//ka5 = (temp*Q)*h;
			spQ.leftmult(temp,ka5);
			ka5 *= h;


			//kc5 = (matrix(b,temp) + Q*(kc1*b51+kc2*b52+kc3*b53+kc4*b54+c))*h;
			kc5 = c;
			kc5.add(kc1, b51);
			kc5.add(kc2, b52);
			kc5.add(kc3, b53);
			kc5.add(kc4, b54);
			kc5 = spQ.rightmult(kc5);
			kc5.add(b,temp);
			kc5 *= h;

		//temp = ka1*b61+ka2*b62+ka3*b63+ka4*b64+ka5*b65+mya;
			temp = mya;
			temp.add(ka1, b61);
			temp.add(ka2, b62);
			temp.add(ka3, b63);
			temp.add(ka4, b64);
			temp.add(ka5, b65);

			
			//ka6 = (temp*Q)*h;
			spQ.leftmult(temp, ka6);
			ka6 *= h;
	
					
			//kc6 = (matrix(b,temp)	+ Q*(kc1*b61+kc2*b62+kc3*b63+kc4*b64+kc5*b65+c))*h;
			kc6 = c;
			kc6.add(kc1, b61);
			kc6.add(kc2, b62);
			kc6.add(kc3, b63);
			kc6.add(kc4, b64);
			kc6.add(kc5, b65);
			kc6 = spQ.rightmult(kc6);
			kc6.add(b, temp);
			kc6 *= h;

			
			//vectr erra = dc1*ka1+dc3*ka3+dc4*ka4+dc5*ka5+dc6*ka6;
			vectr erra = ka1;
			erra *= dc1;
			erra.add(ka3, dc3);
			erra.add(ka4, dc4);
			erra.add(ka5, dc5);
			erra.add(ka6, dc6);

						
			//matrix errc = dc1*kc1+dc3*kc3+dc4*kc4+dc5*kc5+dc6*kc6;
			matrix errc = kc1;
			errc *= dc1;
			errc.add(kc3, dc3);
			errc.add(kc4, dc4);
			errc.add(kc5, dc5);
			errc.add(kc6, dc6);
			
			// check error:
			/*err = erra.absmax();
			double terr = errc.absmax();
			if (err<terr) err = terr;
			err /= eps;
			if (err <= 1.0) break;*/

						
                        double amin = a.absmin();
                        double bmin = b.absmin();
                        if(amin<1.0) amin = 1.0;
                        if(bmin<1.0) bmin = 1.0;
                        if(bmin>amin) amin = bmin;
                        err = erra.absmax()/amin;
                        double terr = errc.absmax()/amin;
                        if(err<terr) err = terr;
                        err /= eps;
                        if(err <= 1.0) break;
                        
                        /*
                        double amax = erra.absmax();
                        if (amax > 1.0) amax = 1.0;
                        
                        double cmax = errc.absmax();
                        if (cmax > 1.0) cmax = 1.0;
                        
                        err = amax > cmax ? amax : cmax;
                        
                        err /= eps;
                        if(err <= 1.0) break;
						*/
			// adjust step:
			double htemp = 0.9*h*pow(err,-0.25);
			if (htemp>h/10) h = htemp;
			else h /= 10;
			if (currt+h == currt) {
                            return ret;
                            /*cerr << "step size underflow in rk2" << endl;
				cerr << "in ss pass:" << endl;
				cerr << Q << endl;
				cerr << a << endl;
				cerr << b << endl;
				cerr << t << ' ' << currt << endl;
				exit(1);*/
				
		
			
			}
		}
		
		currt += h;
		if (err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		//mya += c1*ka1+c3*ka3+c4*ka4+c6*ka6;
		mya.add(ka1, c1);
		mya.add(ka3, c3);
		mya.add(ka4, c4);
		mya.add(ka6, c6);
		
		//c += c1*kc1+c3*kc3+c4*kc4+c6*kc6;
		c.add(kc1, c1);
		c.add(kc3, c3);
		c.add(kc4, c4);
		c.add(kc6, c6);
		
		double ss = mya.sum();
		if (ss>0.0) {
			ret += log(ss);
			c /= ss;
			mya /= ss;
		}
		//double z=mya.normalize();
		//ret += log(z);
		//c/=z;
                
 	           
	}

#ifdef USECACHE
	cacheret[cachepos] = c;
	cachellh[cachepos] = ret;
	cachepos++;
	//if (cachepos>MAXCACHE) { cerr << "cache overflow" << endl; cachepos=0;}
	if (cachepos>=MAXCACHE) cachepos=0;
#endif
	return ret;
}


// these are versions for 2-by-2 matrices that are often faster
double vexpmt2(vectr &a, const matrix &Q, double t) {
     double del = Q[0][0] - Q[1][1];
     double det = del*del + 4*Q[0][1]*Q[1][0];

     double a0 = a[0];
     double a1 = a[1];

	double ch,sh;
	double logfactor = t*(Q[0][0]+Q[1][1])/2;
	if (det>=0.0) {
		det = sqrt(det);
		double ex1 = exp(-det*t);
		logfactor += det*t/2.0;
		ch = (1.0 + ex1)/2.0;
		sh = (1.0 - ex1)/2.0/det;
		if (isnan(sh)) sh = 0.5*exp(-det*t/2.0);
	} else { // should not need for rate matrices
		det = sqrt(-det);
		ch = cos(det*t/2.0);
		sh = sin(det*t/2.0)/det;
	}

     double delsh = del*sh;

     a[0] = a0*(ch + delsh) + 2.0*a1*Q[1][0]*sh;
     a[1] = a1*(ch - delsh) + 2.0*a0*Q[0][1]*sh;
	
     //if (a[0]<0) a[0] = 0.0;
     //if (a[1]<0) a[1] = 0.0;
     double s = a[0]+a[1];
	if (s<0.0) s = 1.0;
     a[0] /= s; a[1] /= s;
     return log(s)+logfactor; //t*(Q[0][0]+Q[1][1])/2;
}

double expmtv2(vectr &b, const matrix &Q, double t) {
     double del = Q[0][0] - Q[1][1];
     double det = del*del + 4*Q[0][1]*Q[1][0];

     double b0 = b[0];
     double b1 = b[1];

	double ch,sh;
	double logfactor = t*(Q[0][0]+Q[1][1])/2;
	if (det>=0.0) {
		det = sqrt(det);
		double ex1 = exp(-det*t);
		logfactor += det*t/2.0;
		ch = (1.0 + ex1)/2.0;
		sh = (1.0 - ex1)/2.0/det;
		if (isnan(sh)) sh = 0.5*exp(-det*t/2.0);
	} else { // should not need for rate matrices
		det = sqrt(-det);
		ch = cos(det*t/2.0);
		sh = sin(det*t/2.0)/det;
	}
	double delsh = del*sh;

     b[0] = b0*(ch + delsh) + 2.0*b1*Q[0][1]*sh;
     b[1] = b1*(ch - delsh) + 2.0*b0*Q[1][0]*sh;

     //if (b[0]<0) b[0] = 0.0;
     //if (b[1]<0) b[1] = 0.0;
     double s = b[0]+b[1];
	if (s<0.0) s = 1.0;
     b[0] /= s; b[1] /= s;
     return log(s)+logfactor; //t*(Q[0][0]+Q[1][1])/2;
}


double suffstat2(const vectr &a, const vectr &b, matrix &c,
		 const matrix &Q, double t) {
  double ret = 0.0;
  
  double del = Q[0][0] - Q[1][1];
  double det = sqrt(del*del + 4*Q[0][1]*Q[1][0]);

  double sa = a.sum();
  double sb = b.sum();
  double a0 = a[0]/sa;
  double a1 = a[1]/sa;
  double b0 = b[0]/sb;
  double b1 = b[1]/sb;

  double ex = exp(det*t/2.0);

  double ch = (ex + 1.0/ex)/2.0;
  double sh = (ex - 1.0/ex)/2.0/det;

  double integral_shsh = (-sh +  t/2.0 * ch)/(det*det);
  if (isnan(sh)) {
    sh = 0.5;
    integral_shsh = 0.0;
  }

  double integral_chch = sh + t/2.0 * ch;
  double integral_chsh = t/2.0 * sh;
  double integral_shch = t/2.0 * sh;

  double s00, s01, s10, s11;
  double tmp;
  tmp = integral_chch;
  tmp += del * 2 * integral_chsh;
  tmp += del * del * integral_shsh;
  s00 = tmp * a0 * b0;
  //----------------
  tmp = integral_chsh + del * integral_shsh;
  s00 += tmp * 2 * Q[0][1] * a0 * b1;
  //----------------
  s00 += tmp * 2 * Q[1][0] * a1 * b0;
  //----------------
  s00 += 4 * Q[0][1] * Q[1][0] * integral_shsh * a1 * b1; 

  //=========================================================

  tmp = integral_chsh + del * integral_shsh;
  s01 = tmp * 2 * Q[1][0] * a0 * b0;
  //----------------
  tmp = integral_chch - del * del * integral_shsh;
  s01 += tmp * a0 * b1;
  //----------------
  s01 += 4 * Q[1][0] * Q[1][0] * integral_shsh * a1 * b0;
  //----------------
  tmp = integral_shch - del * integral_shsh;
  s01 += tmp * 2 * Q[1][0] * a1 * b1;
  //=========================================================

  tmp = integral_shch + del * integral_shsh;
  s10 = tmp * 2 * Q[0][1] * a0 * b0;
  //----------------
  s10 += 4 * Q[0][1] * Q[0][1] * integral_shsh * a0 * b1;
  //----------------
  tmp = integral_chch - del * del * integral_shsh;
  s10 += tmp * a1 * b0;
  //----------------
  tmp = integral_chsh - del * integral_shsh;
  s10 += tmp * 2 * Q[0][1] * a1 * b1;

  //=========================================================

  s11 = 4 * Q[0][1] * Q[1][0] * integral_shsh * a0 * b0;
  //----------------
  tmp = integral_shch - del * integral_shsh;
  s11 += tmp * 2 * Q[0][1] * a0 * b1;
  //----------------
  s11 += tmp * 2 * Q[1][0] * a1 * b0;
  //----------------
  tmp = integral_chch - 2 * del * integral_shch + del * del * integral_shsh;
  s11 += tmp * a1 * b1;

  c[0][0] = s00;
  c[0][1] = s10;
  c[1][0] = s01;
  c[1][1] = s11;

  c *= exp(t*(Q[0][0]+Q[1][1])/2) * sb;
  ret = log(sa);

  return ret;
}

void mfBackward(MeanFieldInf *inf, int varid, double t1, double t0,
		double eps, double h) {

	static double defaulteps = ParamDouble("RKEPS",RKEPS);
	if(eps<0.0) eps = defaulteps;

	if(h<=0.0) h = 0.05;

	int card = inf->GetNumStates(varid);

	vectr kd2(card),kd3(card),kd4(card),kd5(card),kd6(card),
		stepres(card), temp(card), reserr(card);

	/*
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
		*/

	double currt = t1;
	double endt = t0;
	vectr curRho(card,0.0);

	while(currt>endt) {
		if(currt-h<=endt) h = currt-endt;
		double err = 0;
		double hstep = -h;
		for( ; ; ) {
			hstep = -h;
			inf->GetRhoTempVal(currt,curRho);
			vectr rhoGrad(card,0.0);
			inf->CalcRhoGrad(varid,currt,curRho,rhoGrad);
			rhoGrad *= hstep;

			temp = curRho;
			temp.add(rhoGrad,b21);
//			temp = curRho+b21*rhoGrad;
			inf->CalcRhoGrad(varid,currt+a2*hstep,temp,kd2);
			kd2 *= hstep;

			temp = curRho;
			temp.add(rhoGrad,b31);
			temp.add(kd2,b32);
//			temp = curRho+(b31*rhoGrad+b32*kd2);
			inf->CalcRhoGrad(varid,currt+a3*hstep,temp,kd3);
			kd3 *= hstep;

			temp = curRho;
			temp.add(rhoGrad,b41);
			temp.add(kd2,b42);
			temp.add(kd3,b43);
//			temp = curRho+(b41*rhoGrad+b42*kd2+b43*kd3);
			inf->CalcRhoGrad(varid,currt+a4*hstep,temp,kd4);
			kd4 *= hstep;

			temp = curRho;
			temp.add(rhoGrad,b51);
			temp.add(kd2,b52);
			temp.add(kd3,b53);
			temp.add(kd4,b54);
//			temp = curRho+(b51*rhoGrad+b52*kd2+b53*kd3+b54*kd4);
			inf->CalcRhoGrad(varid,currt+a5*hstep,temp,kd5);
			kd5 *= hstep;

			temp = curRho;
			temp.add(rhoGrad,b61);
			temp.add(kd2,b62);
			temp.add(kd3,b63);
			temp.add(kd4,b64);
			temp.add(kd5,b65);
//			temp = curRho+(b61*rhoGrad+b62*kd2+b63*kd3+b64*kd4
//						+b65*kd5);
			inf->CalcRhoGrad(varid,currt+a6*hstep,temp,kd6);
			kd6 *= hstep;

			stepres = curRho;
			stepres.add(rhoGrad,c1);
			stepres.add(kd3,c3);
			stepres.add(kd4,c4);
			stepres.add(kd6,c6);
//			stepres = curRho+(c1*rhoGrad+c3*kd3+c4*kd4+c6*kd6);

			for(int i=0; i<card; i++) {
				if(stepres[i] < 0.0) stepres[i] = 0.0;
			}

			reserr = vectr(card,0.0);
			reserr.add(rhoGrad,dc1);
			reserr.add(kd3,dc3);
			reserr.add(kd4,dc4);
			reserr.add(kd5,dc5);
			reserr.add(kd6,dc6);
//			reserr = dc1*rhoGrad+dc3*kd3+dc4*kd4+dc5*kd+dc6*kd6);


			double mmin = curRho.absmin();
			if(mmin<1.0) mmin = 1.0;
			err = reserr.absmax()/mmin;
			err /= eps;

/*
			cout << "currt: " << currt << endl;
			cout << "endt: " << endt << endl;
			cout << "hstep: " << hstep << endl;
			cout << "rhoGrad: " << rhoGrad << endl;
			cout << "kd2: " << kd2 << endl;
			cout << "kd3: " << kd3 << endl;
			cout << "kd4: " << kd4 << endl;
			cout << "kd5: " << kd5 << endl;
			cout << "kd6: " << kd6 << endl;
			cout << "stepres: " << stepres << endl;
			cout << "reserr: " << reserr << endl;
			cout << "err: " << err << endl;
			*/
			if(err<=1.0) break;
			if(fabs(currt-endt) < eps) break;

			double htemp = 0.9*h*pow(err,-0.25);
			htemp = abs(htemp);
			if(htemp>h/10) h = htemp;
			else h /= 10;
			if(currt+h == currt) {
				cout << "Bad exit" << endl;
				return;
			}
		}
		hstep = -h;
		currt += hstep;
		if(err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		curRho = stepres;
		inf->AddRhoTempVal(currt, curRho);
	}
}

void mfForward(MeanFieldInf *inf, int varid, double t0, double t1,
		double eps, double h) {

	static double defaulteps = ParamDouble("RKEPS",RKEPS);

	if(eps<0.0) eps = defaulteps;

	if(h<=0.0) h = 0.05;

	int card = inf->GetNumStates(varid);

	vectr kd2(card),kd3(card),kd4(card),kd5(card),kd6(card),
		stepres(card), temp(card), rhoTemp(card), reserr(card);
		
	matrix gammaTemp(card,card,0.0);

	/*
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
		*/

	double currt = t0;
	double endt = t1-1e-11;
	vectr curMu(card,0.0);
	matrix curGamma(card,card,0.0);

	bool reachedEnd=false;

	while(currt<endt) {
		if(currt+h>=endt) {
			h = endt-currt;
		}
		if(abs(endt-currt) < MU_END_PRECISION) {
			h = endt-currt;
			vectr endTemp(card,0.0);
			matrix endGamma(card,card,0.0);
			inf->GetRhoTempVal(endt,endTemp);
			int endState = -1;
			if(endTemp.sum() > 1) {
				inf->GetMuTempVal(currt,endTemp);
			}
			else {
				for(int i=0;i<card;i++) {
					if(abs(endTemp[i]-1) < 1e-3) {
						endState = i;
						break;
					}
				}
			}
			if(endState == -1) {
				reachedEnd = true;
			}
			else {
				for(int i=0;i<card;i++) {
					if(i != endState) {
						double gamma = curMu[i]/h;
						endGamma[i][endState] += 
							gamma;
						endGamma[i][i] -= gamma;
					}
				}
				// cout << endGamma << endl;
				inf->AddMuTempVal(endt,endTemp);
				inf->AddGammaTempVal(endt,endGamma);
				return;
			}
		}
		double err = 0;
		double hstep;
		for( ; ; ) {
			hstep = h;

			inf->GetMuTempVal(currt,curMu);
			inf->GetGammaTempVal(currt,curGamma);
			vectr muGrad(card,0.0);
			inf->CalcMuGrad(varid,currt,
					curGamma,muGrad);
			muGrad *= hstep;

			temp = curMu;
			temp.add(muGrad,b21);
//			temp = curMu+b21*muGrad;
			inf->GetRhoTempVal(currt+a2*hstep,rhoTemp);
			inf->CalcGamma(varid,currt+a2*hstep,
					temp,rhoTemp,gammaTemp);

			inf->CalcMuGrad(varid,currt+a2*hstep,
					gammaTemp,kd2);
			kd2 *= hstep;

			temp = curMu;
			temp.add(muGrad,b31);
			temp.add(kd2,b32);
//			temp = curMu+(b31*muGrad+b32*kd2);
			inf->GetRhoTempVal(currt+a3*hstep,rhoTemp);
			inf->CalcGamma(varid,currt+a3*hstep,
					temp,rhoTemp,gammaTemp);
			inf->CalcMuGrad(varid,currt+a3*hstep,
					gammaTemp,kd3);
			kd3 *= hstep;

			temp = curMu;
			temp.add(muGrad,b41);
			temp.add(kd2,b42);
			temp.add(kd3,b43);
//			temp = curMu+(b41*muGrad+b42*kd2+b43*kd3);
			inf->GetRhoTempVal(currt+a4*hstep,rhoTemp);
			inf->CalcGamma(varid,currt+a4*hstep,
					temp,rhoTemp,gammaTemp);
			inf->CalcMuGrad(varid,currt+a4*hstep,
					gammaTemp,kd4);
			kd4 *= hstep;

			temp = curMu;
			temp.add(muGrad,b51);
			temp.add(kd2,b52);
			temp.add(kd3,b53);
			temp.add(kd4,b54);
//			temp = curMu+(b51*muGrad+b52*kd2+b53*kd3+b54*kd4);
			inf->GetRhoTempVal(currt+a5*hstep,rhoTemp);
			inf->CalcGamma(varid,currt+a5*hstep,
					temp,rhoTemp,gammaTemp);
			inf->CalcMuGrad(varid,currt+a5*hstep,
					gammaTemp,kd5);
			kd5 *= hstep;

			temp = curMu;
			temp.add(muGrad,b61);
			temp.add(kd2,b62);
			temp.add(kd3,b63);
			temp.add(kd4,b64);
			temp.add(kd5,b65);
//			temp = curMu+(b61*muGrad+b62*kd2+b63*kd3+b64*kd4
//						+b65*kd5);
			inf->GetRhoTempVal(currt+a6*hstep,rhoTemp);
			inf->CalcGamma(varid,currt+a6*hstep,
					temp,rhoTemp,gammaTemp);
			inf->CalcMuGrad(varid,currt+a6*hstep,
					gammaTemp,kd6);
			kd6 *= hstep;

			stepres = curMu;
			stepres.add(muGrad,c1);
			stepres.add(kd3,c3);
			stepres.add(kd4,c4);
			stepres.add(kd6,c6);
//			stepres = curMu+(c1*muGrad+c3*kd3+c4*kd4+c6*kd6);

			for(int i=0; i<card; i++) {
				if(stepres[i] < 0.0) stepres[i] = 0.0;
				if(stepres[i] > 1.0) stepres[i] = 1.0;
			}

			reserr = vectr(card,0.0);
			reserr.add(muGrad,dc1);
			reserr.add(kd3,dc3);
			reserr.add(kd4,dc4);
			reserr.add(kd5,dc5);
			reserr.add(kd6,dc6);
//			reserr = dc1*muGrad+dc3*kd3+dc4*kd4+dc5*kd+dc6*kd6);


			double mmin = curMu.absmin();
			if(mmin<1.0) mmin = 1.0;
			err = reserr.absmax()/mmin;
			err /= eps;
//			cout << "err=" << err << endl;
			if(err<=1.0) break;

			double htemp = 0.9*h*pow(err,-0.25);
			htemp = abs(htemp);
			if(htemp>h/10) h = htemp;
			else h /= 10;
			if(currt+h == currt) {
				return;
			}
		}
		hstep = h;
		currt += hstep;
		if(err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		curMu = stepres;
		curMu.normalize();
		inf->GetRhoTempVal(currt,rhoTemp);
		inf->CalcGamma(varid,currt,curMu,rhoTemp,curGamma);
		inf->AddMuTempVal(currt, curMu);
		inf->AddGammaTempVal(currt, curGamma);
		if(reachedEnd && currt >= endt) break;
	}
}

double mfEnergy(MeanFieldInf *inf, int varid, double t0, double t1,
		double eps, double h) {

	static double defaulteps = ParamDouble("RKEPS",RKEPS);
	if(abs(t1-t0)< defaulteps) return 0;
	if(eps<0.0) eps = defaulteps;

	if(h<=0.0) h = 0.05;

	int card = inf->GetNumStates(varid);

	double /*kd2(0.0),*/kd3(0.0),kd4(0.0),kd5(0.0),kd6(0.0),
		stepres(0.0), temp(0.0), reserr(0.0);

	double k1(0.0), k2(0.0), k3(0.0), k4(0.0);
		
	/*
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
		*/

	double currt = t0;
	double endt = t1;

	double curEn = 0.0;

	while(currt<endt) {
		if(currt+h>=endt) h = endt-currt;
		double err = 0;
		double hstep;
		for( ; ; ) {
			hstep = h;
			double enIntegrand = 
				inf->CalcPointCompEnergy(varid,currt);

//			kd2 = inf->CalcPointCompEnergy(varid,currt+a2*hstep);
			kd3 = inf->CalcPointCompEnergy(varid,currt+a3*hstep);
			kd4 = inf->CalcPointCompEnergy(varid,currt+a4*hstep);
			kd5 = inf->CalcPointCompEnergy(varid,currt+a5*hstep);
			kd6 = inf->CalcPointCompEnergy(varid,currt+a6*hstep);

			stepres = curEn + hstep*(c1*enIntegrand
						+c3*kd3
						+c4*kd4
						+c6*kd6);

			reserr = hstep*(dc1*enIntegrand
					+dc3*kd3
					+dc4*kd4
					+dc5*kd5
					+dc6*kd6);

			double mmin = abs(curEn);
			if(mmin<1.0) mmin = 1.0;
			err = reserr/mmin;
			err = reserr;
			err /= eps;
			if(err<=1.0) break;

			double htemp = 0.9*h*pow(err,-0.25);
			htemp = htemp > 0 ? htemp : -htemp;
			if(htemp>h/10) h = htemp;
			else h /= 10;
			if(currt+h == currt) {
				return 0;
			}
		}
		hstep = h;
		currt += hstep;
		if(err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		curEn = stepres;
	}
	return curEn;
}

double mfSuffStat(MeanFieldInf *inf, const Instantiation &x, 
		double t0, double t1,
		double eps, double h) {

	static double defaulteps = ParamDouble("RKEPS",RKEPS);
	if(abs(t1-t0)< defaulteps) {
		return 0;
	}
	if(eps<0.0) eps = defaulteps;

	if(h<=0.0) h = 0.05;

//	int card = inf->GetNumStates(varid);

	double /*kd2(0.0),*/kd3(0.0),kd4(0.0),kd5(0.0),kd6(0.0),
		stepres(0.0), temp(0.0), reserr(0.0);

	double k1(0.0), k2(0.0), k3(0.0), k4(0.0);
		
	double currt = t0;
	double endt = t1;

	double curTimeSS = 0.0;

	while(currt<endt) {
		if(currt+h>=endt) h = endt-currt;
		double err = 0;
		double hstep;
		for( ; ; ) {
			hstep = h;
			double integrand = inf->Smooth(x,currt);

//			kd2 = inf->Smooth(x,currt+a2*hstep);
			kd3 = inf->Smooth(x,currt+a3*hstep);
			kd4 = inf->Smooth(x,currt+a4*hstep);
			kd5 = inf->Smooth(x,currt+a5*hstep);
			kd6 = inf->Smooth(x,currt+a6*hstep);

			stepres = curTimeSS + hstep*(c1*integrand
						+c3*kd3
						+c4*kd4
						+c6*kd6);

			reserr = hstep*(dc1*integrand
					+dc3*kd3
					+dc4*kd4
					+dc5*kd5
					+dc6*kd6);

			double mmin = abs(curTimeSS);
			if(mmin<1.0) mmin = 1.0;
			err = reserr/mmin;
			err = reserr;
			err /= eps;
			if(err<=1.0) break;

			double htemp = 0.9*h*pow(err,-0.25);
			htemp = htemp > 0 ? htemp : -htemp;
			if(htemp>h/10) h = htemp;
			else h /= 10;
			if(currt+h == currt) {
				return 0;
			}
		}
		hstep = h;
		currt += hstep;
		if(err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		curTimeSS = stepres;
	}
	return curTimeSS;
}

double mfTransSuffStat(MeanFieldInf *inf, int varid,
		const Instantiation &x1, const Instantiation &x2,
		double t0, double t1,
		double eps, double h) {

	static double defaulteps = ParamDouble("RKEPS",RKEPS);
	if(abs(t1-t0)< defaulteps) return 0;
	if(eps<0.0) eps = defaulteps;

	if(h<=0.0) h = 0.05;

//	int card = inf->GetNumStates(varid);

	double /*kd2(0.0),*/kd3(0.0),kd4(0.0),kd5(0.0),kd6(0.0),
		stepres(0.0), temp(0.0), reserr(0.0);

	double k1(0.0), k2(0.0), k3(0.0), k4(0.0);
		
	/*
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
		*/

	double currt = t0;
	double endt = t1;

	double curTimeSS = 0.0;

	while(currt<endt) {
		if(currt+h>=endt) h = endt-currt;
		double err = 0;
		double hstep;
		for( ; ; ) {
			hstep = h;
			double integrand = 
				inf->PointTransSuffStat(varid,x1,x2,currt);

//			kd2 = inf->PointTransSuffStat(varid,x1,x2,currt+a2*hstep);
			kd3 = inf->PointTransSuffStat(varid,x1,x2,currt+a3*hstep);
			kd4 = inf->PointTransSuffStat(varid,x1,x2,currt+a4*hstep);
			kd5 = inf->PointTransSuffStat(varid,x1,x2,currt+a5*hstep);
			kd6 = inf->PointTransSuffStat(varid,x1,x2,currt+a6*hstep);

			stepres = curTimeSS + hstep*(c1*integrand
						+c3*kd3
						+c4*kd4
						+c6*kd6);

			reserr = hstep*(dc1*integrand
					+dc3*kd3
					+dc4*kd4
					+dc5*kd5
					+dc6*kd6);

			double mmin = abs(curTimeSS);
			if(mmin<1.0) mmin = 1.0;
			err = reserr/mmin;
			err = reserr;
			err /= eps;
			if(err<=1.0) break;

			double htemp = 0.9*h*pow(err,-0.25);
			htemp = htemp > 0 ? htemp : -htemp;
			if(htemp>h/10) h = htemp;
			else h /= 10;
			if(currt+h == currt) {
				return 0;
			}
		}
		hstep = h;
		currt += hstep;
		if(err > 1.89e-4) h = 0.9*h*pow(err,-0.2);
		else h *= 5;

		curTimeSS = stepres;
	}
	return curTimeSS;
}

} // end of ctbn namespace
