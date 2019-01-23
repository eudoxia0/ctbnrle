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
#ifndef CTBNRLE_MATRIX_NOEIGEN_H
#define CTBNRLE_MATRIX_NOEIGEN_H

#include <math.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <iomanip>
#include <vector>
#include "defines.h"

namespace ctbn {

// undefine to cause code to check for OOB conditions
//#define NO_ARRAY_CHECK

// This code defines a vector (class name vectr)
// and a matrix (class name matrix) that have most of the
// common vector-matrix operations.

// Code is kept inline for speed for most of the operations

class matrix;

// Vector class
class vectr {
	friend class matrix;
	friend class spmatrix;
	SERIAL_USESHIFT(vectr)
public:
	// default vector is length 1 (a scalar)
	inline vectr() {
		m = 1;
		x = new double[1];
	}

	// new vector (contents not set)
	inline vectr(int m) {
		this->m = m;
		x = new double[m];
	}

	inline ~vectr() {
		delete []x;
	}

	// copy constructor
	inline vectr(const vectr &v) {
		m = v.m;
		x = new double[m];
		//for(int i=0;i<m;i++) x[i] = v.x[i];
		memcpy(x,v.x,m*sizeof(double));
	}

	// new vector (length m, all elements = a)
	inline vectr(int m, const double &a) {
		this->m = m;
		x = new double[m];
		for(int i=0;i<m;i++) x[i] = a;
	}

	// new vector initialized to array v (of length m)
	// if keep=true, then the memory is not copied and
	// the calling function no longer owns v.
	inline vectr(double *v, int m, bool keep=false) {
		this->m = m;

		if (keep) x = v;
		else {
			x = new double[m];
			//for(int i=0;i<m;i++) x[i] = v[i];
			memcpy(x,v,m*sizeof(double));
		}
	}

	// constructor from a std::vector
	inline vectr(const std::vector<double> &v) {
		m = v.size();
		x = new double[m];
		for(int i=0;i<m;i++) x[i] = v[i];
	}

	// constructor from a "sparce vector"
	// v are the non-zero elements.  ind are their corresponding
	// indexes in the newly constructed vector
	inline vectr(const vectr &v, const std::vector<int> &ind) {
		m = ind.size();
		x = new double[m];
		for(int i=0;i<m;i++)
			x[i] = v.x[ind[i]];
	}

	inline void swap(vectr &v) {
		int t = v.m; v.m = m; m = t;
		double *tp = v.x; v.x = x; x = tp;
	}

	// unary subtraction
	inline vectr operator-() const {
		vectr ret(m);
		for(int i=0;i<m;i++) ret.x[i] = -x[i];
		return ret;
	}

	// assignment to a scalar
	inline vectr& operator=(double a) {
		for(int i=0;i<m;i++) x[i] = a;
		return *this;
	}

	// assignment
	inline vectr& operator=(const vectr &v) {
		if (&v==this) return *this;
		if (v.m != m) {
			delete []x;
			m = v.m;
			x = new double[m];
		}
		//for(int i=0;i<m;i++) x[i] = v.x[i];
		memcpy(x,v.x,m*sizeof(double));
		return *this;
	}

	// returns whether all elements are finite
	inline bool isvalid() const {
		for(int i=0;i<m;i++) if (!finite(x[i])) return false;
		return true;
	}

/*
	// conversion to double*
	// returns the array of values (for read only!)
	// returned value is not valid after non-const method on this
	inline operator const double *() const {
		return x;
	}
*/

	// dot product
	inline double operator*(const vectr &v) const {
#ifndef NO_ARRAY_CHECK
		if (m!=v.m) {
			std::cerr << "invalid vectr dot product" << std::endl;
			assert(0);
		}
#endif
		double ret = 0.0;
		for(int i=0;i<m;i++) ret += x[i]*v.x[i];
		return ret;
	}

	// dot product (assumed v is of correct length)
	inline double operator*(const double *v) const {
		double ret = 0.0;
		for(int i=0;i<m;i++) ret += x[i]*v[i];
		return ret;
	}

/*  Commented out, as it is not used (and not available in the Eigen version
	// outer product, but flattened as a vector
	inline vectr outer(const vectr &v) const {
		vectr ret(m*v.m);
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<v.m;j++,c++)
				ret.x[c] = x[i]*v.x[j];
		return ret;
	}
*/

	// the ".*" operator from Matlab
	inline vectr dotstar(const vectr &v) const {
		return vectr(*this).multby(v);
	}

	// returns m, the number of elements
	inline int length() const{
		return m;
	}

	// the "./" operator from Matlab
	inline vectr dotdiv(const vectr &v) const {
		return vectr(*this).divby(v);
	}

	// an inplace ".*" operation (much faster than dotstar, but overwrites
	// this vector)
	inline vectr &multby(const vectr &v) {
#ifndef NO_ARRAY_CHECK
		if (m!=v.m) {
			std::cerr << "invalid vectr component-wise multiplication" << std::endl;
			assert(0);
		}
#endif
		for(int i=0;i<m;i++)
			x[i] = x[i]*v.x[i];
		return *this;
	}

	// like multby
	inline vectr &divby(const vectr &v) {
#ifndef NO_ARRAY_CHECK
		if (m!=v.m) {
			std::cerr << "invalid vectr component-wise division" << std::endl;
			assert(0);
		}
#endif
		for(int i=0;i<m;i++)
			x[i] = x[i]/v.x[i];
		return *this;
	}

	// access operators
	inline double operator[](int i) const {
#ifndef NO_ARRAY_CHECK
		if (i<0 || i>=m) {
			std::cerr << "invalid array index" << std::endl;
			assert(0);
		}
#endif
		return x[i];
	}
	inline double &operator[](int i) {
#ifndef NO_ARRAY_CHECK
		if (i<0 || i>=m) {
			std::cerr << "invalid array index" << std::endl;
			assert(0);
		}
#endif
		return x[i];
	}

	// These are much faster than using operator+ (but they overwrite this)
	inline vectr &operator+=(const vectr &v) {
#ifndef NO_ARRAY_CHECK
		if (v.m!=m) {
			std::cerr << "invalid vectr addition" << std::endl;
			assert(0);
		}
#endif
		for(int i=0;i<m;i++) x[i] += v.x[i];
		return *this;
	}

	// Same as above, but weight the vector v by weight before
	// multiplying:  this = this + v*weight
	inline vectr &add(const vectr &v, const double &weight) {
#ifndef NO_ARRAY_CHECK
		if (v.m!=m) {
			std::cerr << "invalid vectr addition" << std::endl;
			assert(0);
		}
#endif
		for(int i=0;i<m;i++) x[i] += v.x[i]*weight;
		return *this;
	}

	// faster than operator-, but overwrites this
	inline vectr &operator-=(const vectr &v) {
#ifndef NO_ARRAY_CHECK
		if (v.m!=m) {
			std::cerr << "invalid vectr subtraction" << std::endl;
			assert(0);
		}
#endif
		for(int i=0;i<m;i++) x[i] -= v.x[i];
		return *this;
	}

	// faster than operator+, but overwrites this
	inline vectr &operator*=(const double &a) {
		for(int i=0;i<m;i++) x[i] *= a;
		return *this;
	}

	// faster than operator+, but overwrites this
	inline vectr &operator+=(const double &a) {
		for(int i=0;i<m;i++) x[i] += a;
		return *this;
	}

	// faster than operator-, but overwrites this
	inline vectr &operator-=(const double &a) {
		for(int i=0;i<m;i++) x[i] -= a;
		return *this;
	}

	// faster than operator/, but overwrites this
	inline vectr &operator/=(const double &a) {
		for(int i=0;i<m;i++) x[i] /= a;
		return *this;
	}

	// returns maximum value in vector
	inline double max() const {
		double t,ma = x[0];
		for(int i=1;i<m;i++)
			if((t=x[i])>ma) ma = t;
		return ma;
	}

	// returns minimum value in vector
	inline double min() const {
		double t,mi = x[0];
		for(int i=1;i<m;i++)
			if ((t=x[i])<mi) mi = t;
		return mi;
	}

	// returns maximum absolute value in vector
	inline double absmax() const {
		double t,ma = x[0]>0?x[0]:-x[0];
		for(int i=1;i<m;i++)
			if ((t=(x[i]>0?x[i]:-x[i]))>ma) ma = t;
		return ma;
	}

	// returns minimum absolute value in vector
	inline double absmin() const {
		double t,mi = x[0]>0?x[0]:-x[0];
		for(int i=1;i<m;i++)
			if ((t=(x[i]>0?x[i]:-x[i]))<mi) mi = t;
		return mi;
	}

	// checks to see if two vectors are equal (if lengths are
	// not equal, vectors are not equal)
	inline bool operator==(const vectr &v) const {
		if (m!=v.m) return false;
		return memcmp(x,v.x,m*sizeof(double))==0.0;
		//for(int i=0;i<m;i++) if (v.x[i]!=x[i]) return false;
		//return true;
	}

	// checks to see if all elements equal a
	inline bool operator==(const double &a) const {
		for(int i=0;i<m;i++) if (a!=x[i]) return false;
		return true;
	}

	// checked to see if vectors are not equal (negation of == above)
	inline bool operator!=(const vectr &v) const {
		if (m!=v.m) return true;
		return memcmp(x,v.x,m*sizeof(double))!=0.0;
		//for(int i=0;i<m;i++) if (v.x[i]!=x[i]) return true;
		//return false;
	}

	// checked to see if any element does not equal a (negation of == above)
	inline bool operator!=(const double &a) const {
		for(int i=0;i<m;i++) if (a!=x[i]) return true;
		return false;
	}

	// the square of the L-2 norm of the vector
	inline double norm2() const {
		double ret=x[0]*x[0];
		for(int i=1;i<m;i++) ret += x[i]*x[i];
		return ret;
	}

	// the L-2 norm of the vector (sum of square elements)
	inline double norm() const {
		return sqrt(norm2());
	}

	// divide the vector by its sum, to make it sum to 1
	// (useful usually only if vector has all non-negative elements)
	inline double normalize() {
		double s = sum();
		(*this) /= s;
		return s;
	}

	// save and load
	friend std::ostream& operator<<(std::ostream &s, const vectr &v);
	friend std::istream& operator>>(std::istream &s, vectr &v);

	// same as length, but added 
	inline int getm() const {
		return m;
	}

	// a nicer way of displaying the vector (than <<) for humans
	// cannot be loaded (size not explicitly given)
	inline std::ostream &niceprint(std::ostream &s) const {
		for(int i=0;i<m;i++)
			s << x[i] << ' ';
		s << std::endl;
		return s;
	}

	// sets to zero all elements within fuzz of zero
	inline void unfuzz(double fuzz) {
		for(int i=0; i < m; i++)
			if(x[i] < fuzz) x[i] = 0.0;
	}

	// returns sum
	inline double sum() const {
		double total = 0.0;
		for(int i = 0; i < m; i++)
			total += x[i];
		return total;
	}

	// sets to 0 all values for which corresponding index in s is zero
	inline void support(int *s) {
		for(int i = 0; i < m; i++)
			if(!s[i]) x[i] = 0.0;
	}

	// inline version of unary -
	inline void negate() {
		for(int i = 0; i < m; i++)
			x[i] = -x[i];
	}

	// change the length of the vector (keep as many elements as
	// possible; any added elements will be uninitialized)
	inline void resize(int newm) {
		double *newx = new double[newm];
		memcpy(newx,x,(newm>m?m:newm)*sizeof(double));
		delete []x;
		x = newx;
		m = newm;
	}

private:
	int m;
	double *x;
};

// implementations of binary operators in terms of inplace operators
// (see Effective C++ for why this is a good idea)
inline vectr operator+(const vectr &a, const vectr &b) {
	return vectr(a)+=b;
}

inline vectr operator-(const vectr &a, const vectr &b) {
	return vectr(a)-=b;
}

inline vectr operator+(const vectr &a, const double &b) {
	return vectr(a)+=b;
}

inline vectr operator-(const vectr &a, const double &b) {
	return vectr(a)-=b;
}

inline vectr operator+(const double &a, const vectr &b) {
	return vectr(b)+=a;
}

inline vectr operator-(const double &a, const vectr &b) {
	return vectr(b.getm(),a)-=b;
}

inline vectr operator*(const vectr &a, const double &b) {
	return vectr(a)*=b;
}

inline vectr operator*(const double &a, const vectr &b) {
	return vectr(b)*=a;
}

inline vectr operator/(const vectr &a, const double &b) {
	return vectr(a)/=b;
}

// saving a loading uses the format: <m> <element1> <element2> ... <elementm>
inline std::ostream &operator<<(std::ostream &s, const vectr& v) {
	s << v.m;
	for(int i=0;i<v.m;i++)
		s << ' ' << v.x[i];
	return s;
}

inline std::istream &operator>>(std::istream &s, vectr& v) {
	int tm;
	s >> tm;
	if (tm!=v.m) {
		delete []v.x;
		v.m = tm;
		v.x = new double[tm];
	}
	for(int i=0;i<tm;i++) s >> v.x[i];
	return s;
}

class matrix {
	friend class vectr;
	friend class spmatrix;
	SERIAL_USESHIFT(matrix)
public:

	// default matrix is 1-by-1
	inline matrix(int m=1, int n=1) {
		this->m = m; this->n = n;
		s = m*n;
		x = new double[s];
	}

	inline ~matrix() {
		delete []x;
	}

	// copy constructor and ability to construct transpose
	// matrix a;
	// ...
	// matrix b(a,true); // b is now transpose of a
	inline matrix(const matrix &ma, bool transpose=false) {
		s = ma.m*ma.n;
		x = new double[s];
		if (transpose) {
			int i,j,c;
			n = ma.m; m = ma.n;
			c = 0;
			for(i=0;i<m;i++) for(j=0;j<n;j++,c++)
				x[c] = ma.x[i+j*m];
		} else {
			n = ma.n; m = ma.m;
			int i;
			for(i=0;i<s;i++) x[i] = ma.x[i];
		}
	}

	// create m-by-n matrix.
	// If diaonly=true, all diagional elements=a and all others=0
	// If diaonly=false (default), all elements=a
	inline matrix(int m, int n,const double &a,bool diaonly=false) {
		this->m = m;
		this->n = n;
		s = m*n;
		x = new double[s];
		if (diaonly) {
			//for(i=0;i<s;i++) x[i] = 0;
			memset(x,0,s*sizeof(double));
			if (n>=m) for(int i=0;i<m;i++) x[i*n+i] = a;
			else      for(int i=0;i<n;i++) x[i*n+i] = a;
		} else {
			if (a==0.0) memset(x,0,s*sizeof(double));
			else for(int i=0;i<s;i++) x[i] = a;
		}
	}

	// put v on the diagonal
	inline matrix(int m, int n,const vectr &v) {
		this->m = m;
		this->n = n;
		s = m*n;
		x = new double[s];
		//for(int i=0;i<s;i++) x[i] = 0;
		memset(x,0,s*sizeof(double));
		int l = m;

		if (n<l) l = n;

		if (v.m<l) l = v.m;

		for(int i=0,c=0;i<l;i++,c+=n+1) x[c] = v.x[i];
	}

	// make an m-by-1 matrix of v
	inline matrix(const vectr &v) {
		m = v.m;
		n = 1;
		s = m;
		x = new double[s];
		//for(int i=0;i<s;i++) x[i] = v.x[i];
		memcpy(x,v.x,s*sizeof(double));
	}

	// use the elements of v to make an m-by-n
	// matrix (assumes v is row-major order)
	inline matrix(double *v,int m, int n) {
		this->m = m;
		this->n = n;
		s = m*n;
		//int i;
		x = new double[s];
		//for(i=0;i<s;i++) x[i] = v[i];
		memcpy(x,v,s*sizeof(double));
	}

	// forms a matrix of the outer product (ie v1*v2') -- v2 is
	// "transposed" temporarily for this operation
	inline matrix(const matrix &v1, const matrix &v2) {
		if (v1.n!=v2.n) { // this is technically an error
			s = 1;
			m=1; n=1; x = new double[1];
			x[0] = 0;
		} else {
			n = v2.m; m = v1.m;
			s = n*m;
			x = new double[s];
			int i,j,k,c=0;
			for(i=0;i<m;i++) for(j=0;j<n;j++,c++) {
				x[c] = 0;
				for(k=0;k<v1.n;k++)
					x[c] += v1.x[i*v1.n+k]*
						v2.x[j*v1.n+k];
			}
		}
	}

	// form outer product of v1*v2'
	inline matrix(const vectr &v1, const vectr &v2) {
		n = v2.m; m = v1.m;
		s = n*m;
		x = new double[s];
		int i,j,c=0;
		for(i=0;i<m;i++) for(j=0;j<n;j++,c++)
			x[c] = v1.x[i]*v2.x[j];
	}

	// form out product of v1*v2'*w
	inline matrix(const vectr &v1, const vectr &v2, const double &w) {
		n = v2.m; m = v1.m;
		s = n*m;
		x = new double[s];
		int i,j,c=0;
		for(i=0;i<m;i++) for(j=0;j<n;j++,c++)
			x[c] = v1.x[i]*v2.x[j]*w;
	}

	// forms the submatrix of v of just the indices in ind
	inline matrix(const matrix &v, const std::vector<int> &ind) {
		n = m = ind.size();
		s = n*m;
		x = new double[s];
		int c= 0;
		for(int i=0;i<m;i++) for(int j=0;j<n;j++,c++)
			x[c] = v.x[ind[i]*v.n+ind[j]];
	}

	// forms the submatrix of v of just the indicies
	// in ind1 (rows) and ind2 (cols)
	inline matrix(const matrix &v, const std::vector<int> &ind1,
			const std::vector<int> &ind2) {
		m = ind1.size();
		n = ind2.size();
		s = m*n;
		x = new double[s];
		int c = 0;
		for(int i=0;i<m;i++) for(int j=0;j<n;j++,c++)
			x[c] = v.x[ind1[i]*v.n+ind2[j]];
	}

	inline void swap(matrix &ma) {
		int t = ma.m; ma.m = m; m = t;
		t = ma.n; ma.n = n; n = t;
		t = ma.s; ma.s = s; s = t;
		double *tp = ma.x; ma.x = x; x = tp;
	}

	// unary -
	inline matrix operator-() const {
		matrix ret(m,n);
		for(int i=0;i<s;i++) ret.x[i] = -x[i];
		return ret;
	}

	// set all elements to a
	inline matrix& operator=(double a) {
		if (a==0) memset(x,0,s*sizeof(double));
		else for(int i=0;i<s;i++) x[i] = a;
		return *this;
	}

	// assignment operator
	inline matrix& operator=(const matrix &ma) {
		if (&ma==this) return *this;
		if (ma.n != n || ma.m != m) {
			s = ma.s; m = ma.m; n = ma.n;
			delete []x;
			x = new double[s];
		}
		//for(int i=0;i<s;i++) x[i] = ma.x[i];
		memcpy(x,ma.x,s*sizeof(double));
		return *this;
	}
	
	// return whether all elements are finite
	inline bool isvalid() const {
		for(int i=0;i<s;i++) if(!finite(x[i])) return false;
		return true;
	}

	// matrix multiplication
	inline matrix operator*(const matrix &ma) const {
#ifndef NO_ARRAY_CHECK
		if (n!=ma.m) {
			std::cerr << "invalid matrix multiply" << std::endl;
			assert(0);
		}
#endif
		matrix ret(m,ma.n);
		int c=0,c2=0;
		for(int i=0;i<m;i++,c2+=n) for(int j=0;j<ma.n;j++,c++) {
			ret.x[c] = 0;
			for(int k=0,c3=0;k<n;k++,c3+=ma.n)
				ret.x[c] += x[c2+k] * ma.x[c3+j];
		}
		return ret;
	}
	
	// matrix by matrix multiplication out = this * in
	inline void rightmult(const matrix &in, matrix &out) const {
#ifndef NO_ARRAY_CHECK
		if (n!=in.m || m!=out.m || in.n!=out.n) {
			std::cerr << "invalid matrix multiply" << std::endl;
			assert(0);
		}
#endif
		int c=0,c2=0;
		for(int i=0;i<m;i++,c2+=n) for(int j=0;j<in.n;j++,c++) {
			out.x[c] = 0;
			for(int k=0,c3=0;k<n;k++,c3+=in.n)
				out.x[c] += x[c2+k] * in.x[c3+j];
		}
	}

	inline matrix rightmult(const matrix &in) const {
		matrix ret(m,in.n, 0.0);
		rightmult(in,ret);
		return ret;
	}

	// out = this*in (standard multiplication by a vector on the right)
	// also available as an operator, but this version doesn't
	// need to allocate a new vectr
	inline void rightmult(const vectr &in, vectr &out) const {
#ifndef NO_ARRAY_CHECK
		if (n!=in.m || m!=out.m) {
			std::cerr << "invalid vectr-matrix multiply" << std::endl;
			assert(0);
		}
#endif
		for(int i=0,c=0;i<m;i++) {
			out.x[i] = 0;
			for(int j=0;j<n;j++,c++)
				out.x[i] += x[c] * in.x[j];
		}
	}

	// standard multiplication of a vector on the right
	inline vectr operator*(const vectr &v) const {
		vectr ret(m);
		rightmult(v,ret);
		return ret;
	}

	// same as above, but out = in*this (left multiplication)
	inline void leftmult(const vectr &in, vectr &out) const {
#ifndef NO_ARRAY_CHECK
		if (m!=in.m || n!=out.m) {
			std::cerr << "invalid matrix multiplication" << std::endl;
			assert(0);
		}
#endif
		//std::cerr << "in old leftmult" << std::endl;
		out = 0.0;
		for(int i=0;i<n;i++)
			for(int j=0,c=i;j<m;j++,c+=n)
				out.x[i] += x[c]*in.x[j];
	}

	// same as operator*, but for left multiplication
	// does not call leftmult above, because it is faster not to
	// (in terms of setting elements in ret to 0)
	inline vectr transmult(const vectr &v) const {
#ifndef NO_ARRAY_CHECK
		if (m!=v.m) {
			std::cerr << "invalid matrix multiplication" << std::endl;
			assert(0);
		}
#endif
		vectr ret(n,0.0);
		for(int i=0;i<n;i++)
			for(int j=0,c=i;j<m;j++,c+=n)
				ret.x[i] += x[c]*v.x[j];
		return ret;
	}

	// returns the diagonal as a vector
	inline vectr diag() const {
		int sz = (m<n ? m:n);
		double *d = new double[sz];
		vectr ret(sz);
		for(int i=0,c=0;i<sz;i++,c+=n+1)
			d[i] = x[c];
		return vectr(d,sz,true);
	}

	// returns the ith row vector (indexed from 0)
	// (or rather, its transpose)
	inline vectr row(int i) const {
		return vectr(x+(i*n),n,false);
	}

	// returns the ith column vector (indexed from 0)
	inline vectr column(int j) const {
		double *d = new double[m];
		for(int i=0,c=j;i<m;i++,c+=n) d[i] = x[c];
		return vectr(d,m,true);
	}

	// report the sum of multiplying each element in row r
	// by the corresponding element in v
	inline double rowmult(int r, const vectr &v) const {
#ifndef NO_ARRAY_CHECK
		if (n!=v.m) {
			std::cerr << "invalid matrix-vectr multiply" << std::endl;
			assert(0);
		}
#endif
		double ret = 0.0;
		for(int j=0,c=n*r;j<n;j++,c++)
			ret += x[c]*v[j];
		return ret;
	}

	// report the sum of multiplying each element in this
	// by the corresponding element from ma
	inline double dot(const matrix &ma) const {
#ifndef NO_ARRAY_CHECK
		if (n!=ma.n || m!=ma.m) {
			std::cerr << "invalid matrix dot-product" << std::endl;
			assert(0);
		}
#endif
		double ret = 0.0;
		for(int i=0,c=0;i<m;i++,c+=n)
			for(int j=0;j<n;j++)
				ret += x[c+j]*ma.x[c+j];
		return ret;
	}

	// construct the outer product this * ma'
	// place it in ret
	// assumes ret is already the correct size
	inline void outer(const matrix &ma, matrix &ret) const {
#ifndef NO_ARRAY_CHECK
		if (n!=ma.n || ret.m!=m || ret.n!=ma.m) {
			std::cerr << "invalid matrix outer multiply" << std::endl;
			assert(0);
		}
#endif
		for(int i=0,c=0;i<m;i++) for(int j=0;j<ma.m;j++,c++) {
			ret.x[c] = 0;
			for(int k=0;k<n;k++)
				ret.x[c] += x[i*n+k] * ma.x[j*ma.m+k];
		}
	}

	// construct the inner product this*ma
	// place it in ret
	// assumes ret is already the correct size
	inline void inner(const matrix &ma, matrix &ret) const {
#ifndef NO_ARRAY_CHECK
		if (m!=ma.m || ret.m!=n || ret.n!=ma.n) {
			std::cerr << "invalid matrix inner multiply" << std::endl;
			assert(0);
		}
#endif
		int c=0;
		for(int i=0;i<n;i++) for(int j=0;j<ma.n;j++,c++) {
			ret.x[c] = 0;
			for(int k=0;k<m;k++)
				ret.x[c] += x[k*n+i] * ma.x[k*m+j];
		}
	}

	// inplace ".*" operator from Matlab
	inline matrix &multby(const matrix &ma) {
#ifndef NO_ARRAY_CHECK
		if (ma.m!=m || ma.n!=n) {
			std::cerr << "invalid matrix point-wise multiplication" << std::endl;
			assert(0);
		}
#endif
		for(int i=0;i<s;i++)
			x[i] *= ma.x[i];
		return *this;
	}

	// returns this.*m (the non-inplace version of above)
	inline matrix dotstar(const matrix &m) const {
		return matrix(*this).multby(m);
	}

	// inplace, multiply each row's elements by the corresponding
	// elements of v
	// (v should be of size n if this is of size m-by-n)
	inline matrix &multbyrow(const double *v) {
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<n;j++,c++)
				x[c] *= v[j];
		return *this;
	}

	// same as above, but for columns (so v should be of size m)
	inline matrix &multbycol(const double *v) {
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<n;j++,c++)
				x[c] *= v[i];
		return *this;
	}

	// same as above, but with vectr
	inline matrix &multbyrow(const vectr &v) {
#ifndef NO_ARRAY_CHECK
		if (n!=v.m) {
			std::cerr << "invalid multbyrow" << std::endl;
			assert(0);
		}
#endif
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<n;j++,c++)
				x[c] *= v[j];
		return *this;
	}

	// same as above, but with vectr
	inline matrix &multbycol(const vectr &v) {
#ifndef NO_ARRAY_CHECK
		if (m!=v.m) {
			std::cerr << "invalid multbycol" << std::endl;
			assert(0);
		}
#endif
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<n;j++,c++)
				x[c] *= v[i];
		return *this;
	}

	// next four are the same as multby... but for division
	inline matrix &dividebyrow(const double *v) {
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<n;j++,c++)
				x[c] /= v[j];
		return *this;
	}

	inline matrix &dividebycol(const double *v) {
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<n;j++,c++)
				x[c] /= v[i];
		return *this;
	}

	inline matrix &dividebyrow(const vectr &v) {
#ifndef NO_ARRAY_CHECK
		if (n!=v.m) {
			std::cerr << "invalid dividebycol" << std::endl;
			assert(0);
		}
#endif
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<n;j++,c++)
				x[c] /= v[j];
		return *this;
	}

	inline matrix &dividebycol(const vectr &v) {
#ifndef NO_ARRAY_CHECK
		if (m!=v.m) {
			std::cerr << "invalid dividebycol" << std::endl;
			assert(0);
		}
#endif
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<n;j++,c++)
				x[c] /= v[i];
		return *this;
	}

	// access operators:
	// matrix m;
	// m(1,2) // only for r-value
	// m[1][2] // r-value or l-value
	// [] doesn't use proxy class for speed
	inline double operator()(int i, int j) const {
#ifndef NO_ARRAY_CHECK
		if (i<0 || i>=m || j<0 || j>=n) {
			std::cerr << "invalid array index" << std::endl;
			assert(0);
		}
#endif
		return x[i*n+j];
	}

	inline double &operator()(int i, int j) {
#ifndef NO_ARRAY_CHECK
		if (i<0 || i>=m || j<0 || j>=n) {
			std::cerr << "invalid array index" << std::endl;
			assert(0);
		}
#endif
		return x[i*n+j];
	}

	inline const double *operator[](int i) const {
#ifndef NO_ARRAY_CHECK
		if (i<0 || i>=m) {
			std::cerr << "invalid array index" << std::endl;
			assert(0);
		}
#endif
		return x+(i*n);
	}

	inline double *operator[](int i) {
#ifndef NO_ARRAY_CHECK
		if (i<0 || i>=m) {
			std::cerr << "invalid array index" << std::endl;
			assert(0);
		}
#endif
		return x+(i*n);
	}

	// return transpose of the matrix
	inline matrix t() const {
		return matrix(*this,true);
	}
	
	// add ma*weight to this (inplace) 
	inline matrix &add(const matrix &ma, const double &weight) {
#ifndef NO_ARRAY_CHECK
		if (m!=ma.m||n!=ma.n) {
			std::cerr << "invalid matrix addition" << std::endl;
			assert(0);
		}
#endif
		for(int i=0;i<s;i++) x[i] += ma.x[i]*weight;
		return *this;
	}

	// adds outer product of v with itself multiplied by weight...
	// ...pretty specialized but much faster than other methods
	inline matrix &add(const vectr &v, double weight) {
#ifndef NO_ARRAY_CHECK
		if (v.m!=m||n!=m) {
			std::cerr << "invalid matrix addition" << std::endl;
			assert(0);
		}
#endif
		for(int i=0,c=0;i<m;i++)
			for(int j=0;j<n;j++) x[c++] += v.x[i]*v.x[j]*weight;
		return *this;
	}

	// inplace additon of outer product of v1*v2' 
	inline matrix &add(const vectr &v1, const vectr &v2) {
#ifndef NO_ARRAY_CHECK
		if (n!=v2.m || m!=v1.m) {
			std::cerr << "invalid matrix addition" << std::endl;
			assert(0);
		}
#endif
		//n = v2.m; m = v1.m;
		//s = n*m;
		//x = new double[s];
		int i,j,c=0;
		for(i=0;i<m;i++) for(j=0;j<n;j++,c++)
			x[c] += v1.x[i]*v2.x[j];
		return *this;
	}
	

	// inplace addition
	inline matrix &operator+=(const matrix &ma) {
#ifndef NO_ARRAY_CHECK
		if (m!=ma.m||n!=ma.n) {
			std::cerr << "invalid matrix addition" << std::endl;
			assert(0);
		}
#endif
		for(int i=0;i<s;i++) x[i] += ma.x[i];
		return *this;
	}
	
	// inplace subtraction
	inline matrix &operator-=(const matrix &ma) {
#ifndef NO_ARRAY_CHECK
		if (m!=ma.m||n!=ma.n) {
			std::cerr << "invalid matrix addition" << std::endl;
			assert(0);
		}
#endif
		for(int i=0;i<s;i++) x[i] -= ma.x[i];
		return *this;
	}

	// inplace matrix multiplication
	inline matrix &operator*=(const matrix &ma) {
#ifndef NO_ARRAY_CHECK
		if (n!=ma.m || n != ma.n) {
			std::cerr << "invalid matrix multiply" << std::endl;
			assert(0);
		}
#endif
		double *newrow = new double[n];
		for(int i=0,c=0;i<m;i++) {
			for(int j=0;j<n;j++) {
				newrow[j] = 0;
				for(int k=0;k<n;k++)
					newrow[j] += x[c+k] * ma.x[k*n+j];
			}
			for(int j=0; j < n; j++,c++)
				x[c] = newrow[j];
		}
		delete []newrow;
		return *this;
	}

	// inplace addition of scalar
	inline matrix &operator+=(const double &a) {
		for(int i=0;i<s;i++) x[i] += a;
		return *this;
	}

	// inplace subtraction of scalar
	inline matrix &operator-=(const double &a) {
		for(int i=0;i<s;i++) x[i] -= a;
		return *this;
	}

	// inplace multiplication by scalar
	inline matrix &operator*=(const double &a) {
		for(int i=0;i<s;i++) x[i] *= a;
		return *this;
	}

	// inplace division by scalar
	inline matrix &operator/=(const double &a) {
		for(int i=0;i<s;i++) x[i] /= a;
		return *this;
	}

	// returns the maximal element
	inline double max() const {
		double t,ma = x[0];
		for(int i=1;i<s;i++)
			if ((t=x[i])>ma) ma=t;
		return ma;
	}

	// returns the minimal element
	inline double min() const {
		double t,mi = x[0];
		for(int i=1;i<s;i++)
			if ((t=x[i])<mi) mi=t;
		return mi;
	}

	// returns the element with the maximal absolute value
	inline double absmin() const {
		double t,mi = x[0]>0?x[0]:-x[0];
		for(int i=1;i<s;i++)
			if ((t=(x[i]>0?x[i]:-x[i]))<mi) mi=t;
		return mi;
	}

	// returns the element with the minimal absolute value
	inline double absmax() const {
		double t,ma = x[0]>0?x[0]:-x[0];
		for(int i=1;i<s;i++)
			if ((t=(x[i]>0?x[i]:-x[i]))>ma) ma=t;
		return ma;
	}

	// returns the sum of all elements
	inline double sum() const {
		double ret = 0.0;
		for(int i=0;i<s;i++)
			ret += x[i];
		return ret;
	}

	// returns whether the two matrices are equal
	// (if sizes different, they are not equal)
	inline bool operator==(const matrix &ma) const {
		if (ma.n!=n||ma.m!=m) return false;
		return memcmp(ma.x,x,s*sizeof(double))==0.0;
		//for(int i=0;i<s;i++) if (ma.x[i]!=x[i]) return false;
		//return true;
	}

	// returns whether all elements equal a
	inline bool operator==(const double &a) const {
		for(int i=0;i<s;i++) if (x[i]!=a) return false;
		return true;
	}

	// inverse of operator==
	inline bool operator!=(const matrix &ma) const {
		if (ma.n!=n||ma.m!=m) return true;
		return memcmp(ma.x,x,s*sizeof(double))!=0.0;
		//for(int i=0;i<s;i++) if (ma.x[i]!=x[i]) return true;
		//return false;
	}

	// inverse of operator==
	inline bool operator!=(const double &a) const {
		for(int i=0;i<s;i++) if (x[i]!=a) return true;
		return false;
	}

	// returns the square of the frobenius norm
	inline double norm2() const {
		double ret=x[0]*x[0];
		for(int i=1;i<s;i++) ret += x[i]*x[i];
		return ret;
	}

	// returns the frobenius norm
	inline double norm() const {
		return sqrt(norm2());
	}

	// returns the L-inf norm of the matrix
	// (largest L1 row norm)
	inline double norm_inf() const {
		double ret = 0.0;
		for(int i=0,c=0; i<m;i++) {
			double rowsum = 0.0;
			for(int j=0; j<n;j++,c++) {
				rowsum += (x[c]<0?-x[c]:x[c]);
			}
			if (ret<rowsum) ret = rowsum;
		}
		return ret;
	}

	friend std::ostream& operator<<(std::ostream& s, const matrix& ma);
	friend std::istream& operator>>(std::istream& s, matrix& ma);

	// LU decomposition -- ix is the row permutations
	int LUdecomp(matrix &LU, int *ix) const;
	// LU back substitution --
	//    ix from above fn call (this should be an LU combination)
	void LUbacksub(int *ix, double *col) const;

	// solves equation Ax=b (A is this, x is the returned value)
	double *solve(const double *b, bool &worked) const;
	inline double *solve(const double *b) const {
		bool w; return solve(b,w);
	}

	// inplace version of unary -
	inline void negate() {
		for(int i = 0; i < s; i++) x[i] = -x[i];
	}

	// returns inverse of the matrix
	matrix inv(bool &worked) const;
	// dets is negative to skip determinant calculation
	// upon return dets \in {-1,0,1} (sign of determinant)
	//    where 0 => could not perform inverse (upon return)
	// if det calc done, logdet is the log of the absolute value
	//   of the determinant
	matrix inv(double &logdet, int &dets) const;
	matrix expm() const;

	// same as inv, but don't report whether it worked
	inline matrix inv() const {
		bool w; return inv(w);
	}

	// inplace, calculates the adjoint (and returns the determinant)
	double adjoint();
	// return the trace of the matrix (sum of diagonal elements)
	inline double trace();

	// dest = this*source (but we do not need to allocate memory for dest
	inline void multiply(const vectr &source, vectr &dest) {
#ifndef NO_ARRAY_CHECK
		assert(n == source.m && m == dest.m);
#endif
		for(int i = 0,c=0; i < m; i++) {
			dest[i] = 0;
			for(int j = 0; j < n; j++,c++)
				dest[i] += x[c] * source[j];
		}
	}

	// inplace pre-mult for sparse vector v
	// v = v*this (is faster if many elements of v are zero)
	void sparsepremult(vectr &v) const;

	// singular value decomposition of this
	// if this is m-by-n, then upon return
	//  u is m-by-n
	//  v is n-by-n
	//  w is n-by-1
	// and (as much as possible), this = u*W*v' (where W is
	// a diagonal matrix with elements from w)
	void svd(matrix &u, matrix &v, vectr &w) const;

	// makes the matrix tri-diagonal (if symmetric, otherwise ?)
	// The matrix is replaced by Q (the pre-multipler), and the
	// vectors d and e hold the diagonal and off-diagonal
	// elements respectively
	void tridiag(vectr &d, vectr &e);
	// converts the Q matrix from tridiag into the eigenvectors
	// d and e are the diagonal and off-diagonal elements from tridiag
	// and the eigenvalues are left in d
	void tridiagev(vectr &d, vectr &e);

	// uses the above two methods to compute eigen decomp of symmetric
	// matrix (destroys matrix and replaces with eigenvectors)
	vectr eigensym() {
		vectr d,e; tridiag(d,e); tridiagev(d,e); return d;
	}

	// same as above, but doesn't destroy current matrix
	// also does reordering (and the unreordering) of rows/cols
	// for better numeric stability
	vectr eigensymsmart(matrix &evectr) const;

	// prints out the matrix is a some-what nicer form
	// keeps m and n at the top so that it can be read by operator>>
	inline std::ostream &niceprint(std::ostream& s) const {
		s << m << ' ' << n << std::endl;
		for(int i=0;i<m;i++) {
			for(int j=0;j<n;j++)
				s << x[i*n+j] << ' ';
			s << std::endl;
		}
		return s;
	}

	inline int getm() const {
		return m;
	}

	inline int getn() const {
		return n;
	}

protected:
	void permutebydiag(std::vector<int> &revperm, matrix &out) const;
	static double pythag(double a, double b);

private:
	int m,n,s;
	double *x;
};

class spmatrix {
	friend class matrix;
public:
	inline spmatrix(const matrix &ma){
		m = ma.m;
		n = ma.n;
		
		val.reserve(m+n);
		colind.reserve(m+n);
		rowptr.reserve(m);
		s = 0;
		int c=0;
		for(int i=0; i<m; i++){
			rowptr.push_back(s);
			for(int j=0; j<n; j++){
				if (ma.x[i*n+j] > 1e-16 || ma.x[i*n+j] < -1e-16){
					val.push_back(ma.x[i*n+j]);
					colind.push_back(j);
					s++;
				}
			}
		}
	}
	
	
	// matrix multiplication
	inline matrix rightmult(const matrix &ma) const {
#ifndef NO_ARRAY_CHECK
		if (n!=ma.m) {
			std::cerr << "invalid matrix multiply" << std::endl;
			assert(0);
		}
#endif
		matrix ret(m,ma.n, 0.0);
		int c=0;
		for(int i=0; i<m; i++){
			int endind = (i < m-1)? rowptr[i+1] : s;
			
			for(int c1=0; c1<ma.n; c1++, c++){
				
				for(int j=rowptr[i]; j<endind; j++){
					ret.x[c] += val[j] * ma.x[colind[j]*ma.n + c1];
				}
			}
		}
		return ret;
	
	}

	// same as above, but out = in*this (left multiplication)
	inline void leftmult(const vectr &in, vectr &out) const {
#ifndef NO_ARRAY_CHECK
		if (m!=in.m || n!=out.m) {
			std::cerr << "invalid matrix multiplication" << std::endl;
			assert(0);
		}
#endif
		out = 0.0;
		for(int i=0;i<m;i++) {
			int endind = (i<m-1)?rowptr[i+1] : s;
			for(int j=rowptr[i]; j<endind;j++)
				out.x[colind[j]] += val[j] * in.x[i];
		}
	}
	
	
	inline std::ostream &niceprint(std::ostream &s) const {
		for(int i=0;i<this->s;i++)
			s << val[i] << ' ';
		
		s << std::endl;
		
		for(int i=0;i<this->s;i++)
			s << colind[i] << ' ';
		s << std::endl;
		
		for(int i=0;i<this->m;i++)
			s << rowptr[i] << ' ';
		s << std::endl;
		return s;
	}
	
private:
	int m, n, s;
	std::vector<double> val;
	std::vector<int> colind;
	std::vector<int> rowptr;	
};

// binary operators implemented in terms of inplace operators
inline matrix operator+(const matrix &a, const matrix &b) {
	return matrix(a)+=b;
}

inline matrix operator-(const matrix &a, const matrix &b) {
	return matrix(a)-=b;
}

inline matrix operator+(const matrix &a, const double &b) {
	return matrix(a)+=b;
}

inline matrix operator-(const matrix &a, const double &b) {
	return matrix(a)-=b;
}

inline matrix operator+(const double &a, const matrix &b) {
	return matrix(b)+=a;
}

inline matrix operator-(const double &a, const matrix &b) {
	return matrix(b.getn(),b.getm(),a)-=b;
}

inline matrix operator*(const matrix &a, const double &b) {
	return matrix(a)*=b;
}

inline matrix operator*(const double &b, const matrix &a) {
	return matrix(a)*=b;
}

inline vectr operator*(const vectr &v, const matrix &m) {
	return m.transmult(v);
}

inline matrix operator/(const matrix &a, const double &b) {
	return matrix(a)/=b;
}

inline matrix outer(const vectr &v1, const vectr &v2) {
	return matrix(v1,v2);
}


// loading and saving.  Format is <m> <n> <[0][0]> <[0][1]> ... <[m][n]>
inline std::ostream& operator<<(std::ostream& s, const matrix& ma) {
	s << ma.m << ' ' << ma.n << ' ';
	for(int i=0;i<ma.s;i++) {
		if (i%ma.n==0) s << std::endl; s << std::setw(4) << ma.x[i]; if (i!=ma.s) s << ' ';
	}
	return s;
}

inline std::istream& operator>>(std::istream& s, matrix& ma) {
	int tn,tm;
	s >> tm >> tn;
	if (tm!=ma.m || tn!=ma.n) {
		delete []ma.x;
		ma.s = tm*tn;
		ma.x = new double[ma.s];
		ma.m = tm; ma.n = tn;
	}
	for(int i=0;i<ma.s;i++) {
		s >> ma.x[i];
	}
	return s;
}

} // end of ctbn namespace

#endif
