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
#include "matrix.h"
#include <iostream>
#include <stdlib.h>

namespace ctbn {
	
using namespace std;

// adopted from NRiC, pg 45
matrix matrix::inv(bool &worked) const {
	if (m!=n) {
		cerr << "invalid matrix inverse" << endl;
		exit(1);
	}

	matrix temp(n,n);
	int *ix = new int[n];

	if (!LUdecomp(temp,ix)) {
		worked = false;
		delete []ix;
		return matrix(n,n,0,false);
	}

	worked = true;

	matrix ret(n,n);
	int i,j;
	double *col = new double[n];

	for(j=0;j<n;j++) {
		for(i=0;i<n;i++) col[i] = 0;

		col[j] = 1;

		temp.LUbacksub(ix,col);

		for(i=0;i<n;i++) ret.x[i*n+j] = col[i];
	}

	delete []col;
	delete []ix;
	return ret;
}

matrix matrix::inv(double &logdet, int &detsign) const {
	if (m!=n) {
		cerr << "invalid matrix inverse" << endl;
		exit(1);
	}

	matrix temp(n,n);
	int *ix = new int[n];

	bool dodet = (detsign>=0);

	if ((detsign=LUdecomp(temp,ix))==0) {
		delete []ix;
		return matrix(n,n,0,false);
	}

	logdet = 0.0;

	if (dodet) {
		for(int i=0;i<n;i++) {
			double v = temp[i][i];

			if (v<0.0) {
				detsign = -detsign;
				v = -v;
			}

			logdet += log(v);
		}
	}

	matrix ret(n,n);
	int i,j;
	double *col = new double[n];

	for(j=0;j<n;j++) {
		for(i=0;i<n;i++) col[i] = 0;

		col[j] = 1;

		temp.LUbacksub(ix,col);

		for(i=0;i<n;i++) ret.x[i*n+j] = col[i];
	}

	delete []ix;
	delete []col;
	return ret;
}


//Algorithm from Golub,van Loan, p. 558; code based on TBCI project
matrix matrix::expm() const {
	if (m!=n) {
		cerr << "invalid matrix exponential" << endl;
		exit(1);
	}

	int e,p,q,k;
	matrix E(n,n);
	matrix temp(*this);
	matrix eye(n,n,1.0,true); //identity matrix of size n
	frexp(norm_inf(),&e);
	int s = (0>(e+1)?0:(e+1));
	temp /= (pow(2.0,s));

	matrix X(temp);
	matrix cX(n,n);
	double c = 0.5;
	matrix D(eye-temp*c);

	E=eye+temp*c;
	q = 6;
	p = 1;

	for (k=2;k<=q;k++) {
		c    = c*(q-k+1)/(k*(2*q-k+1));
		X    = temp*X;
		cX   = X*c;
		E   += cX;

		if(p)
			D += cX;
		else
			D -= cX;

		p = !p;

	}

	//solve for F in DF = E
	matrix F = E * D.inv();

	for (k=1;k<=s;k++) {
		F = F*F;
	}

	return F;
}

void matrix::sparsepremult(vectr &v) const {
	vectr c(v);
	v = 0.0;
	for(int i=0,s=0;i<m;i++,s+=n)
		if (c.x[i]!=0.0)
			for(int j=0;j<n;j++)
				v.x[j] += c.x[i]*x[s+j];
}

// adopted from NRiC, pg 43
int matrix::LUdecomp(matrix &LU, int *ix) const {
	if (m!=n||LU.m!=LU.n||LU.n!=n) {
		cerr << "invalid matrix in LUdecomp" << endl;
		exit(1);
	}

	int d=1,i,j,k;
	LU = *this;
	double *vv = new double[n];
	double dum;

	k = 0;

	for(i=0;i<n;i++) {
		vv[i] = fabs(x[k]); k++;

		for(j=1;j<n;j++,k++) if(vv[i]<(dum=fabs(x[k]))) vv[i]=dum;

		if (vv[i]==(double)0.0) {
			delete []vv;
			return 0;
		}

		vv[i] = 1/vv[i];
	}

	double sum,big;
	int imax;

	for(j=0;j<n;j++) {
		for(i=0;i<j;i++) {
			sum = LU.x[i*n+j];

			for(k=0;k<i;k++) sum -= LU.x[i*n+k]*LU.x[k*n+j];

			LU.x[i*n+j] = sum;
		}

		big = 0;

		for(i=j;i<n;i++) {
			sum = LU.x[i*n+j];

			for(k=0;k<j;k++) sum -= LU.x[i*n+k]*LU.x[k*n+j];

			LU.x[i*n+j] = sum;

			if ((dum=vv[i]*fabs(sum))>=big) {
				big = dum;
				imax = i;
			}
		}

		if (j!=imax) {
			for(k=0;k<n;k++) {
				dum = LU.x[imax*n+k];
				LU.x[imax*n+k] = LU.x[j*n+k];
				LU.x[j*n+k] = dum;
			}

			d = -d;
			vv[imax] = vv[j];
		}

		ix[j] = imax;

		if (LU.x[j*n+j] == 0) {
			LU.x[j*n+j] = (double)1.0e-20;
		}

		if (j!=n-1) {
			dum = 1/LU.x[j*n+j];

			for(i=j+1;i<n;i++) LU.x[i*n+j] *= dum;
		}
	}

	delete []vv;
	return d;
}

void matrix::LUbacksub(int *ix, double *b) const {
	if (n!=m) {
		cerr << "invalid matrix in LUbacksub" << endl;
		exit(1);
	}

	int ip,ii=-1,i;
	double sum;

	for(i=0;i<n;i++) {
		ip = ix[i];
		sum = b[ip];
		b[ip] = b[i];

		if (ii!=-1)
			for(int j=ii;j<=i-1;j++) sum -= x[i*n+j]*b[j];
		else if (sum!=0) ii=i;

		b[i] = sum;
	}

	for(i=n-1;i>=0;i--) {
		sum = b[i];

		for(int j=i+1;j<=n-1;j++) sum -= x[i*n+j]*b[j];

		b[i] = sum/x[i*n+i];
	}
}

double *matrix::solve(const double *b, bool &worked) const {
	if (m!=n) {
		cerr << "invalid matrix in solve" << endl;
		exit(1);
	}

	double *ret = new double[n];

	for(int i=0;i<n;i++) ret[i] = b[i];

	int *ix = new int[n];

	matrix a(n,n);

	if (!LUdecomp(a,ix)) {
		worked = false;
		delete []ix;
		return ret;
	}

	worked=true;
	a.LUbacksub(ix,ret);
	delete []ix;
	return ret;
}

double matrix::pythag(double a, double b) {
	double absa,absb;
	absa = fabs(a);
	absb = fabs(b);

	if (absa>absb) {
		double sqr = absb/absa;
		return absa*sqrt(1.0+sqr*sqr);
	} else {
		if (absb==0.0) return 0;

		double sqr = absa/absb;

		return absb*sqrt(1.0+sqr*sqr);
	}
}

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))

// adopted from pg 67 of NRiC
void matrix::svd(matrix &u, matrix &v, vectr &w) const {

	u = *this;

	if (v.n!=n || v.m!=n) {
		delete []v.x;
		v.n = n; v.m = n;
		v.s = n*n;
		v.x = new double[v.s];
	}

	if (w.m!=n) {
		delete []w.x;
		w.m = n;
		w.x = new double[n];
	}

	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1 = new double[n];
	g=scale=anorm=(double)0;

	for(i=0;i<n;i++) {
		l = i+1;
		rv1[i] = scale*g;
		g=s=scale=(double)0.0;

		if (i<m) {
			for(k=i;k<m;k++) scale += fabs(u[k][i]);

			if (scale) {
				for(k=i;k<m;k++) {
					u[k][i] /= scale;
					s += u[k][i]*u[k][i];
				}

				f = u[i][i];
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				u[i][i] = f-g;

				for(j=l;j<n;j++) {
					for(s=0.0,k=i;k<m;k++)
						s += u[k][i]*u[k][j];

					f=s/h;

					for(k=i;k<m;k++) u[k][j] += f*u[k][i];
				}

				for(k=i;k<m;k++) u[k][i] *= scale;
			}
		}

		w[i] = scale *g;
		g=s=scale=(double)0;

		if (i<m && i!=n-1) {
			for(k=l;k<n;k++) scale += fabs(u[i][k]);

			if (scale) {
				for(k=l;k<n;k++) {
					u[i][k] /= scale;
					s += u[i][k]*u[i][k];
				}

				f = u[i][l];
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				u[i][l] = f-g;

				for(k=l;k<n;k++) rv1[k] = u[i][k]/h;

				for(j=l;j<m;j++) {
					for(s=(double)0,k=l;k<n;k++)
						s+=u[j][k]*u[i][k];

					for(k=l;k<n;k++) u[j][k] += s*rv1[k];
				}

				for(k=l;k<n;k++) u[i][k] *= scale;
			}
		}

		double temp = fabs(w[i])+fabs(rv1[i]);
		anorm = anorm>temp ? anorm : temp;
	}

	for(i=n-1;i>=0;i--) {
		if (i<n-1) {
			if (g) {
				for(j=l;j<n;j++)
					v[j][i] = (u[i][j]/u[i][l])/g;

				for(j=l;j<n;j++) {
					for(s=(double)0,k=l;k<n;k++)
						s += u[i][k]*v[k][j];

					for(k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}

			for(j=l;j<n;j++) v[i][j]=v[j][i]=(double)0;
		}

		v[i][i] = (double)1;
		g=rv1[i];
		l=i;
	}

	for(i=m>n?n-1:m-1;i>=0;i--) {
		l=i+1;
		g=w[i];

		for(j=l;j<n;j++) u[i][j] = (double)0;

		if (g) {
			g = 1/g;

			for(j=l;j<n;j++) {
				for(s=(double)0,k=l;k<m;k++)
					s += u[k][i]*u[k][j];

				f = (s/u[i][i])*g;

				for(k=i;k<m;k++) u[k][j] += f*u[k][i];
			}

			for(j=i;j<m;j++) u[j][i] *= g;
		} else for (j=i;j<m;j++) u[j][i] = (double)0;

		++u[i][i];
	}

	for(k=n-1;k>=0;k--) {
		for(its=1;its<=30;its++) {
			flag = 1;

			for(l=k;l>=0;l--) {
				nm = l-1;

				if ((double)(fabs(rv1[l])+anorm)==anorm) {
					flag = 0;
					break;
				}

				if ((double)(fabs(w[nm])+anorm)==anorm) break;
			}

			if (flag) {
				c = (double)0;
				s = (double)1;

				for(i=l;i<=k;i++) {
					f = s*rv1[i];
					rv1[i] = c*rv1[i];

					if ((double)(fabs(f)+anorm)==anorm) break;

					g = w[i];

					h = pythag(f,g);

					w[i] = h;

					h = 1/h;

					c = g*h;

					s = -f*h;

					for(j=0;j<m;j++) {
						y = u[j][nm];
						z = u[j][i];
						u[j][nm] = y*c+z*s;
						u[j][i] = z*c-y*s;
					}
				}
			}

			z=w[k];

			if (l==k) {
				if (z<0.0) {
					w[k] = -z;

					for(j=0;j<n;j++) v[j][k] = -v[j][k];
				}

				break;
			}

			if (its==30) {
				// some other method (like an error return
				// value should be put here)
				cerr << "no convergence after 30 svdcmp "
				"iterations" << endl;
				exit(1);
			}

			x = w[l];
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1;

			for(j=l;j<=nm;j++) {
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;

				for(jj=0;jj<n;jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x*c+z*s;
					v[jj][i] = z*c-x*s;
				}

				z = pythag(f,h);
				w[j] = z;

				if (z) {
					z = 1/z;
					c = f*z;
					s = h*z;
				}

				f = c*g+s*y;
				x = c*y-s*g;

				for(jj=0;jj<m;jj++) {
					y=u[jj][j];
					z=u[jj][i];
					u[jj][j]=y*c+z*s;
					u[jj][i]=z*c-y*s;
				}
			}

			rv1[l] = 0;
			rv1[k] = f;
			w[k] = x;
		}
	}

	delete []rv1;
}

// returns det
double matrix::adjoint() {
	int i, j, i0, j0, maxi, lastj = -1;
	double max, pivot;
	int *r = new int [m];
	int *r2 = new int [m];
	int *c = new int [m];
	double D = 1.0;
	matrix retval(m,m);
	//double retval[m][m];

	for(i = 0; i < m; i++)
		for(j = 0; j < m; j++)
			retval[i][j] = x[i*n+j];

	for(i= 0; i < m; i++) { //pivot permutation
		r[i] = -1;
	}

	for(j = 0; j < m; j++) {
		if(D == 0.0) // det=0, degree 2 singular
			return 0;

		max = -1.0;

		maxi = -1;

		for(i = 0; i < m; i++) { // find pivot (just in column)

			if(r[i] < 0 && fabs(retval[i][j]) > max) {
				max = fabs(retval[i][j]);
				maxi = i;
			}
		}

		if(j != lastj && max == 0.0) { // leave this column for last
			lastj = j;

			if(j != m-1)
				continue;
		}

		if(maxi == -1) {
			cout << "oops";  // NaN input
			return 0;
		}

		i = maxi;
		pivot = retval[i][j];

		for(i0 = 0; i0 < m; i0++) {
			if(i0 != i) {
				for(j0 = 0; j0 < m; j0++) {
					if(j0 != j) {
						retval[i0][j0] *= pivot;
						retval[i0][j0] -= retval[i0][j] * retval[i][j0];
						retval[i0][j0] /= D;
					}
				}
			}
		}

		for(i0 = 0; i0 < m; i0++) {
			retval[i0][j] = -retval[i0][j];
		}

		retval[i][j] = D;
		D = pivot;
		r[i] = j;
		c[j] = i;

		if(j == lastj)
			break;

		if(j == m-1 && lastj >= 0) // go back and do the column we left for last
			j = lastj - 1;
	}

	int s=0;
	i = 0;
	memcpy(r2, r, m * sizeof(int));

	while(i < m-1) { // find sign of permutation (i just hacked this up, might be done quicker)
		j = r2[i];

		if(i != j) {
			s++;
			r2[i] = r2[j];
			r2[j] = j;
		} else
			i++;
	}

	//unpermute rows and columns
	for(i = 0; i < m; i++)
		for(j = 0; j < m; j++)
			x[i*n+j] = retval[c[i]][r[j]];

	if(s%2 == 1) {
		negate();
		D = -D;
	}
	
	delete [] r;
	delete [] r2;
	delete [] c;
	return D;
}

double matrix::trace() {
	assert(n == m);
	double sum = 0.0;

	for(int i = 0; i < n; i++) {
		sum += x[i*n+i];
	}

	return sum;
}

// from NRiC, pg 474

void matrix::tridiag(vectr &d, vectr &e) {
	// matrix must also be symmetric, but we aren't going to check that

	if (n!=m) {
		cerr << "must be square matrix in tridiag" << endl;
		exit(1);
	}

	if (d.m != n) {
		if (d.x) delete []d.x;

		d.x = new double[n];

		d.m = n;
	}

	if (e.m != n) {
		if (e.x) delete []e.x;

		e.x = new double[n];

		e.m = n;
	}

	double *xi = x+(n-1)*n;

	for (int i=n-1;i>=1;i--,xi-=n) {
		int l = i-1;
		double h=0.0,scale=0.0;

		if (l>0) {
			for(int k=0;k<=l;k++)
				scale += fabs(xi[k]);

			if (scale==0.0)
				e.x[i] = xi[l];
			else {
				for (int k=0;k<=l;k++) {
					xi[k] /= scale;
					h += xi[k]*xi[k];
				}

				double f = xi[l];
				double g = (f>=0.0 ? -sqrt(h) : sqrt(h));
				e.x[i] = scale*g;
				h -= f*g;
				xi[l] = f-g;
				f = 0.0;
				double *xj = x;

				for(int j=0;j<=l;j++,xj+=n) {
					// next line can be omitted if evectors not needed
					xj[i] = xi[j]/h;
					g = 0.0;

					for(int k=0;k<=j;k++)
						g += xj[k]*xi[k];

					double *xk = x+(j+1)*n;

					for(int k=j+1;k<=l;k++,xk+=n)
						g += xk[j]*xi[k];

					e.x[j] = g/h;

					f += e.x[j]*xi[j];
				}

				double hh = f/(h+h);
				xj = x;

				for(int j=0;j<=l;j++,xj+=n) {
					f = xi[j];
					g = e.x[j]-hh*f;
					e.x[j]=g;

					for(int k=0;k<=j;k++)
						xj[k] -= (f*e.x[k]+g*xi[k]);
				}
			}
		} else
			e.x[i] = xi[l];

		d.x[i] = h;
	}

	// next line can be omitted if evectors not needed
	d.x[0]=0.0;

	e.x[0]=0.0;

	xi = x;

	// the whole loop (excepting noted line) can be omitted
	//    if evectors not needed
	for(int i=0;i<n;i++,xi+=n) {
		int l = i-1;

		if (d.x[i]) {
			for(int j=0;j<=l;j++) {
				double g=0.0;
				double *xk=x;

				for(int k=0;k<=l;k++,xk+=n)
					g += xi[k]*xk[j];

				xk=x;

				for(int k=0;k<=l;k++,xk+=n)
					xk[j] -= g*xk[i];
			}
		}

		d.x[i] = xi[i]; // excepted line
		xi[i] = 1.0;
		double *xj = x;

		for(int j=0;j<=l;j++,xj+=n) xj[i]=xi[j]=0.0;
	}
}

void matrix::tridiagev(vectr &d, vectr &e) {
	if (n!=m) {
		cerr << "must be square matrix in tridiagev" << endl;
		exit(1);
	}

	if (d.getm()!=n || e.getm()!=n) {
		cerr << "invalid input arguments to tridiagev" << endl;
		exit(1);
	}

	matrix save(*this);
	vectr dsave(d),esave(e);

	for(int i=1;i<n;i++) e[i-1] = e[i];

	e[n-1] = 0.0;

	for(int l=0;l<n;l++) {
		int iter = 0;
		int mm;

		do {
			for (mm=l;mm<n-1;mm++) {
				double dd = fabs(d.x[mm])+fabs(d.x[mm+1]);

				if ((double)(fabs(e.x[mm])+dd) == dd) break;
			}

			if (mm != l) {
				if (iter++ == 30) {
					cerr << "too many iterations in tridiagev (" << l << ',' << m << ')' << endl;
					ofstream data("tridiagev.data");
					data << save << endl << dsave << endl << esave << endl;
					data.close();
					exit(1);
				}

				double g = (d.x[l+1]-d.x[l])/(2.0*e.x[l]);
				double r = pythag(g,1.0);
				g = d.x[mm]-d.x[l]+e.x[l]/(g+SIGN(r,g));
				double s=1.0,c=1.0;
				double p=0.0;
				int i;

				for (i=mm-1;i>=l;i--) {
					double f = s*e.x[i];
					double b = c*e.x[i];
					e.x[i+1] = (r=pythag(f,g));

					if (r==0.0) {
						d.x[i+1] -= p;
						e.x[mm] = 0.0;
						break;
					}

					s = f/r;
					c = g/r;
					g = d.x[i+1]-p;
					r =(d.x[i]-g)*s+2.0*c*b;
					d.x[i+1] = g+(p=s*r);
					g = c*r-b;
					// this loop may be omitted if evectors not needed
					double *xk = x;

					for (int k=0;k<n;k++,xk+=n) {
						f = xk[i+1];
						xk[i+1] = s*xk[i]+c*f;
						xk[i] = c*xk[i]-s*f;
					}
				}

				if (r==0.0 && i>=l) continue;

				d.x[l] -= p;

				e.x[l] = g;

				e.x[mm] = 0.0;
			}
		} while(mm!=l);

		//cout << l << ": " << iter << endl;
	}
}

void matrix::permutebydiag(vector<int> &revperm, matrix &out) const {
	if (n!=m) {
		cerr << "permutebydiag currently requires a square matrix" << endl;
		exit(1);
	}

	vector<int> perm(n);
	revperm.resize(n);
	double *diag = new double[n];
	int *dind = new int[n];
	double *xii = x;

	for(int i=0;i<n;i++,xii+=(n+1)) {
		diag[i] = *xii;
		dind[i] = i;
	}

	for(int i=0;i<n;i++) {
		double bv = diag[i];
		int bi = i;

		for(int j=i+1;j<n;j++)
			if (diag[j]<diag[bi]) {
				bi=j;
				bv = diag[j];
			}

		perm[i] = dind[bi];
		revperm[dind[bi]] = i;

		if (bi!=i) {
			double t = diag[i];
			diag[i] = diag[bi];
			diag[bi] = t;
			int td = dind[i];
			dind[i] = dind[bi];
			dind[bi] = td;
		}
	}

	delete []dind;
	delete []diag;
	out = matrix(*this,perm);
}

vectr matrix::eigensymsmart(matrix &evectr) const {
	vector<int> rperm;
	matrix permcpy;
	permutebydiag(rperm,permcpy);
	vectr d,e;
	permcpy.tridiag(d,e);
	permcpy.tridiagev(d,e);
	evectr = matrix(permcpy,rperm);
	return vectr(d,rperm);
}

} // end of ctbn namespace
