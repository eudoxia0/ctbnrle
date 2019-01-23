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
#ifndef CTBNRLE_MATRIX_EIGEN_H
#define CTBNRLE_MATRIX_EIGEN_H

#define EIGEN_MATRIXBASE_PLUGIN "matrix-eigenext.h"

#include "streamserial.h"
#include "Eigen/Core"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET


#include "Eigen/Sparse"

#include <iostream>
#include <vector>

namespace ctbn {
	namespace ctbninternal {
		template<typename EClass, typename CClass> class mbase;

		typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
			RowMatrix;
	}
	class vectr;
	class matrix;
}


namespace Eigen {
	namespace internal {

		template<typename EClass, typename CClass>
		struct traits<ctbn::ctbninternal::mbase<EClass,CClass> > {
			typedef typename traits<EClass>::Scalar Scalar;
			typedef typename traits<EClass>::StorageKind StorageKind;
			typedef typename traits<EClass>::Index Index;
			typedef typename traits<EClass>::XprKind XprKind;
			enum {
				RowsAtCompileTime
				   = traits<EClass>::RowsAtCompileTime,
				ColsAtCompileTime
				   = traits<EClass>::ColsAtCompileTime,
				MaxColsAtCompileTime
				   = traits<EClass>::MaxColsAtCompileTime,
				MaxRowsAtCompileTime
				   = traits<EClass>::MaxRowsAtCompileTime,
				Flags
				   = traits<EClass>::Flags,
				CoeffReadCost
				   = traits<EClass>::CoeffReadCost,
				Options
				   = traits<EClass>::Options,
				InnerStrideAtCompileTime
				= traits<EClass>::InnerStrideAtCompileTime,
				OuterStrideAtCompileTime
				= traits<EClass>::OuterStrideAtCompileTime
			};
		};
		template<>
		struct traits<ctbn::vectr> {
			typedef traits<Eigen::VectorXd>::Scalar Scalar;
			typedef traits<Eigen::VectorXd>::StorageKind StorageKind;
			typedef traits<Eigen::VectorXd>::Index Index;
			typedef traits<Eigen::VectorXd>::XprKind XprKind;
			enum {
				RowsAtCompileTime
				   = traits<Eigen::VectorXd>::RowsAtCompileTime,
				ColsAtCompileTime
				   = traits<Eigen::VectorXd>::ColsAtCompileTime,
				MaxColsAtCompileTime
				   = traits<Eigen::VectorXd>::MaxColsAtCompileTime,
				MaxRowsAtCompileTime
				   = traits<Eigen::VectorXd>::MaxRowsAtCompileTime,
				Flags
				   = traits<Eigen::VectorXd>::Flags,
				CoeffReadCost
				   = traits<Eigen::VectorXd>::CoeffReadCost,
				Options
				   = traits<Eigen::VectorXd>::Options,
				InnerStrideAtCompileTime
				   = traits<Eigen::VectorXd>::InnerStrideAtCompileTime,
				OuterStrideAtCompileTime
				   = traits<Eigen::VectorXd>::OuterStrideAtCompileTime
			};
		};
		template<>
		struct traits<ctbn::matrix> {
			typedef traits<ctbn::ctbninternal::RowMatrix>::Scalar Scalar;
			typedef traits<ctbn::ctbninternal::RowMatrix>::StorageKind StorageKind;
			typedef traits<ctbn::ctbninternal::RowMatrix>::Index Index;
			typedef traits<ctbn::ctbninternal::RowMatrix>::XprKind XprKind;
			enum {
				RowsAtCompileTime
				   = traits<ctbn::ctbninternal::RowMatrix>::RowsAtCompileTime,
				ColsAtCompileTime
				   = traits<ctbn::ctbninternal::RowMatrix>::ColsAtCompileTime,
				MaxColsAtCompileTime
				   = traits<ctbn::ctbninternal::RowMatrix>::MaxColsAtCompileTime,
				MaxRowsAtCompileTime
				   = traits<ctbn::ctbninternal::RowMatrix>::MaxRowsAtCompileTime,
				Flags
				   = traits<ctbn::ctbninternal::RowMatrix>::Flags,
				CoeffReadCost
				   = traits<ctbn::ctbninternal::RowMatrix>::CoeffReadCost,
				Options
				   = traits<ctbn::ctbninternal::RowMatrix>::Options,
				InnerStrideAtCompileTime
				   = traits<ctbn::ctbninternal::RowMatrix>::InnerStrideAtCompileTime,
				OuterStrideAtCompileTime
				   = traits<ctbn::ctbninternal::RowMatrix>::OuterStrideAtCompileTime
			};
		};
	}
}


namespace ctbn{

namespace ctbninternal {

	template<typename EClass, typename CClass>
	class mbase : public EClass {
	public:
		typedef mbase<EClass,CClass> MyType;

		mbase<EClass,CClass>() : EClass() { }
		mbase<EClass,CClass>(int m) : EClass(m) { }
		mbase<EClass,CClass>(int m, int n) : EClass(m,n) { }
		mbase<EClass,CClass>(const EClass &m) : EClass(m) { }
		template<typename OD>
		mbase<EClass,CClass>(const Eigen::MatrixBase<OD> &other)
			: EClass(other) { }


		typedef EClass Base;

		template<typename OD>
		mbase<EClass,CClass> &operator=(const Eigen::MatrixBase<OD> &other) {
			this->Base::operator=(other);
			return *this;
		}


		typedef typename EClass::Scalar Scalar;
		typedef typename EClass::Index Index;

		using EClass::rows;
		using EClass::cols;
		using EClass::fill;
		using EClass::sum;
		using EClass::norm;

		/*  now included in MatrixBase<.>
		struct valid_visitor {
			bool valid;
			inline void init(const Scalar &value, Index i, Index j) {
				valid = (isfinite(value) ? true : false);
			}
			inline void operator() (const Scalar &value, Index i, Index j) {
				if (!isfinite(value)) valid = false;
			}
		};

		inline bool isvalid() const {
			valid_visitor vv;
			this->visit(vv);
			return vv.valid;
		}

		inline Scalar min() const {
			return this->minCoeff();
		}
		inline Scalar max() const {
			return this->maxCoeff();
		}
		inline Scalar absmin() const {
			return this->array().abs().minCoeff();
		}
		inline Scalar absmax() const {
			return this->array().abs().maxCoeff();
		}
		*/

		inline bool operator==(const MyType &o) const {
			return o.rows()==rows() && o.cols()==cols() && EClass::operator==(o);
		}
		inline bool operator!=(const MyType &o) const {
			return o.rows()!=rows() || o.cols()!=cols() || EClass::operator!=(o);
		}
		inline Scalar norm2() const {
			return EClass::squaredNorm();
		}

		template<typename OtherDerived>
		inline const
			typename Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<typename Eigen::internal::traits<EClass>::Scalar,typename Eigen::internal::traits<OtherDerived>::Scalar>,
				const EClass, const OtherDerived >
		dotstar(const Eigen::MatrixBase<OtherDerived> &other) const {
			return this->cwiseProduct(other);
		}

		template<typename OtherDerived>
		inline const
			typename Eigen::CwiseBinaryOp<Eigen::internal::scalar_quotient_op<typename Eigen::internal::traits<EClass>::Scalar>,
				const EClass, const OtherDerived >
		dotdiv(const Eigen::MatrixBase<OtherDerived> &other) const {
			return this->cwiseQuotient(other);
		}

		template<typename OtherDerived>
		inline typename Eigen::MatrixWrapper<
				Eigen::ArrayWrapper<EClass> >
		multby(const Eigen::MatrixBase<OtherDerived> &other) {
				return (this->array() *= other.array()).matrix();
		}

		inline CClass& operator=(const Scalar &a) {
			if (a==0.0) memset(this->data(),0,sizeof(Scalar)*cols()*rows());
			else setConstant(a);
			return *(static_cast<CClass *>(this));
		}

		inline void unfuzz(const Scalar &fuzz) {
			for(int i=0;i<rows();i++)
				for(int j=0;j<cols();j++)
					if ((*this)(i,j)<fuzz) (*this)(i,j) = 0.0;
		}

		inline int getm() const {
			return rows();
		}
		inline int getn() const {
			return cols();
		}
		// version in Eigen does not return sum
		inline Scalar normalize() {
			Scalar s= sum();
			(*this) /= s;
			return s;
		}

		inline CClass &add(const MyType &m, const Scalar &s) {
			*this += m*s;
			return *(static_cast<CClass *>(this));
		}

		using Base::operator+;
		using Base::operator+=;
		using Base::operator-;
		using Base::operator-=;

		inline CClass &operator+=(const Scalar &s) {
			(this->array()) += s;
			return *(static_cast<CClass *>(this));
		}
		inline CClass &operator-=(const Scalar &s) {
			(this->array()) -= s;
			return *(static_cast<CClass *>(this));
		}


	};
}

template<typename D>
inline const
typename Eigen::CwiseUnaryOp<
		typename Eigen::internal::scalar_add_op<
			typename Eigen::internal::traits<D>::Scalar>,
		const D >
operator+(const typename Eigen::MatrixBase<D> &m,
		const typename Eigen::internal::traits<D>::Scalar &s) {
	return Eigen::CwiseUnaryOp<
		typename Eigen::internal::scalar_add_op<
			typename Eigen::internal::traits<D>::Scalar>,
		const D >(m.derived(),Eigen::internal::scalar_add_op<
				typename Eigen::internal::traits<D>::Scalar>(s));
}

template<typename D>
inline const
typename Eigen::CwiseUnaryOp<
		typename Eigen::internal::scalar_add_op<
			typename Eigen::internal::traits<D>::Scalar>,
		const D >
operator-(const typename Eigen::MatrixBase<D> &m,
		const typename Eigen::internal::traits<D>::Scalar &s) {
	return Eigen::CwiseUnaryOp<
		typename Eigen::internal::scalar_add_op<
			typename Eigen::internal::traits<D>::Scalar>,
		const D >(m.derived(),Eigen::internal::scalar_add_op<
				typename Eigen::internal::traits<D>::Scalar>(-s));
}

class vectr
   : public ctbninternal::mbase<Eigen::VectorXd,vectr> {
	   SERIAL_USESHIFT(vectr)
public:
	typedef ctbninternal::mbase<Eigen::VectorXd,vectr> RealBase;
	typedef RealBase Base;
	//typedef Eigen::Matrix<double,1,Eigen::Dynamic> Base;
	typedef Base::Scalar Scalar;
	typedef Base::Index Index;

	vectr() : RealBase() { }
	vectr(int n) : RealBase(n) { }
	vectr(int n, const double &a) : RealBase(n) { fill(a); }
	vectr(const vectr &v, const std::vector<int> &ind) : RealBase(ind.size()) {
		int s = ind.size();
		for(int i=0;i<s;i++) coeffRef(i) = v.coeff(ind[i]);
	}
	vectr(const Eigen::VectorXd &v) : RealBase(v) { }
	template<typename OD>
	vectr(const Eigen::MatrixBase<OD> &other) : RealBase(other) { }

	inline vectr &operator=(const vectr &m) {
		resizeLike(m);
		memcpy(this->data(),m.data(),sizeof(Scalar)*rows()*cols());
		return *this;
	}

	template<typename OD>
	inline vectr &operator=(const Eigen::MatrixBase<OD> &other) {
		this->Base::operator=(other);
		return *this;
	}

	inline vectr& operator=(Scalar a) {
		if (a==0.0) memset(this->data(),0,sizeof(Scalar)*cols()*rows());
		else setConstant(a);
		return *this;
	}

	inline void swap(vectr &m) {
		Base::swap(m);
	}

	inline Scalar operator*(const vectr &v) const { return this->dot(v); }

	inline const Eigen::CwiseUnaryOp<Eigen::internal::scalar_multiple_op<Scalar>, const Eigen::VectorXd> operator*(const Scalar &s) const {
		return Eigen::VectorXd::operator*(s);
	}

	inline int length() const { return getm(); }

	//busra: changed it to row_count from col_count. why was it cols?
	inline std::ostream & niceprint(std::ostream &s) const {
	  unsigned int row_count = rows();
		for (unsigned int i = 0;i != row_count; ++i) 
			s << Eigen::VectorXd::coeff(i) << ' ' ;
		s << std::endl;
		return s;
	}

	template<typename OD>
	inline const typename Eigen::ProductReturnType<Eigen::Transpose<const OD>, Eigen::VectorXd>::Type operator*(const Eigen::MatrixBase<OD> &m) const {
		return (m.transpose()*(*this));
	}



	inline const Eigen::SparseDenseProductReturnType<Eigen::Transpose<const Eigen::SparseMatrix<double,Eigen::RowMajor> >, Eigen::VectorXd>::Type operator*(const Eigen::SparseMatrix<double,Eigen::RowMajor> &m) const {
		return (m.transpose()*(*this));
	}

};

inline std::ostream &operator<<(std::ostream &s, const vectr &v) {
	s << v.getm() << std::endl;
	for(int i=0;i<v.getm();i++) (i==0 ? s : s << ' ') << v.Eigen::VectorXd::coeff(i);
	return s;
}

inline std::istream &operator>>(std::istream &s, vectr &v) {
	int tm;
	s >> tm;
	v.Eigen::VectorXd::resize(tm);
	double t;
	for(int i=0;i<tm;i++)
		s >> v.Eigen::VectorXd::coeffRef(i);
	return s;
}

template<typename D,typename OD>
inline const typename Eigen::ProductReturnType<D,Eigen::Transpose<const OD> >::Type
outer(const Eigen::MatrixBase<D> &v1,const Eigen::MatrixBase<OD> &v2) {
	return v1*(v2.transpose());
}


class matrix : public ctbninternal::mbase<ctbninternal::RowMatrix,matrix> {
	   //SERIAL_USESHIFT(matrix)
public:
	typedef ctbninternal::RowMatrix EBase;
	typedef ctbninternal::mbase<EBase,matrix> RealBase;
	typedef RealBase Base;
	typedef Base::Scalar Scalar;
	typedef Base::Index Index;

	matrix() : RealBase() { }
	matrix(int m,int n) : RealBase(m,n) { }
	matrix(int m, int n, const Scalar &a) : RealBase(m,n) { 
		if (a==0.0) memset(this->data(),0,sizeof(Scalar)*m*n);
		else setConstant(a);
	}
	matrix(const ctbninternal::RowMatrix &m) : RealBase(m) { }
	template<typename OD>
	matrix(const Eigen::MatrixBase<OD> &other) : RealBase(other) { }

	matrix(const matrix &m, const std::vector<int> &ind)
				: RealBase(ind.size(),ind.size()) {
		int s = ind.size();
		for(int i=0;i<s;i++)
			for(int j=0;j<s;j++)
				coeffRef(i,j) = m.coeff(ind[i],ind[j]);
	}
	matrix(const matrix &m, const std::vector<int> &ind1,
			const std::vector<int> &ind2)
				: RealBase(ind1.size(),ind2.size()) {
		int s1=ind1.size(), s2=ind2.size();
		for(int i=0;i<s1;i++)
			for(int j=0;j<s2;j++)
				coeffRef(i,j) = m.coeff(ind1[i],ind2[j]);
	}

	matrix(int m, int n, const vectr &diag) : RealBase(m,n) {
		memset(this->data(),0,sizeof(Scalar)*m*n);
		for(int i=0;i<m && i<n;i++) coeffRef(i,i) = diag[i];
	
	}

	matrix(const vectr &v1, const vectr &v2) : RealBase(v1.getm(),v2.getm()) {
		*this = v1.Eigen::VectorXd::operator*(v2.transpose());
	}


	inline matrix &operator=(const matrix &m) {
		resizeLike(m);
		memcpy(this->data(),m.data(),sizeof(Scalar)*rows()*cols());
		return *this;
	}
	template<typename OD>
	inline matrix &operator=(const Eigen::MatrixBase<OD> &other) {
		this->Base::operator=(other);
		return *this;
	}
	inline matrix& operator=(Scalar a) {
		if (a==0.0) memset(this->data(),0,sizeof(Scalar)*rows()*cols());
		else setConstant(a);
		return *this;
	}

	inline void swap(matrix &m) {
		Base::swap(m);
	}

	struct IndexProxy {
		Index i;
		matrix &m;
		IndexProxy(Index ii, matrix &mm) : i(ii), m(mm) {}
		inline Scalar &operator[](Index index)
			{ return m(i,index); }
		inline const Scalar &operator[](Index index) const
			{ return m(i,index); }
	};

	struct CIndexProxy {
		Index i;
		const matrix &m;
		CIndexProxy(Index ii, const matrix &mm) : i(ii), m(mm) {}
		inline const Scalar &operator[](Index index) const
			{ return m(i,index); }
	};

	inline IndexProxy operator[](Index i)
		{ return IndexProxy(i,*this); }
	inline const CIndexProxy operator[](Index i) const
		{ return CIndexProxy(i,*this); }

	inline std::ostream &niceprint(std::ostream& s) const {
		s << getm() << ' ' << getn() << std::endl;
		for(int i=0;i<getm();i++) {
			for(int j=0;j<getn();j++)
				s << (*this)(i,j) << ' ';
			s << std::endl;
		}
		return s;
	}

	inline matrix &multbyrow(const vectr &v) {
		for(Index j=0;j<getm();j++)
			Base::RowXpr(derived(),j).array() *= v.array();
		return *this;
	}
	inline matrix &multbycol(const vectr &v) {
		for(Index j=0;j<getn();j++)
			Base::ColXpr(derived(),j).array() *= v.array();
		return *this;
	}

	inline void leftmult(const vectr &in, vectr &out) const {
		out = in*(*this);
	}

	inline void rightmult(const vectr &in, vectr &out) const {
		out = (*this)*in;
	}

	inline const Eigen::ProductReturnType<EBase,EBase>::Type rightmult(const matrix &ma) const {
		return (*this)*ma;
	}
	
	inline Eigen::Diagonal<ctbninternal::RowMatrix> diag()
		{ return diagonal(); }
	inline const Eigen::Diagonal<const ctbninternal::RowMatrix> diag() const
		{ return diagonal(); }

	using Base::add;

	inline matrix &add(const vectr &v, const Scalar &w) {
		(*this) += v.Eigen::MatrixBase<Eigen::VectorXd>::operator*(v.transpose())*w;
		return *this;
	}
	inline matrix &add(const vectr &v1, const vectr &v2) {
		(*this) += v1.Eigen::MatrixBase<Eigen::VectorXd>::operator*(v2.transpose());
		return *this;
	}

};

inline std::ostream &operator<<(std::ostream &s, const matrix &m) {
	s << m.getm() << ' ' << m.getn() << std::endl;
	for(int i=0;i<m.getm();i++) {
		for(int j=0;j<m.getn();j++)
			(j==0 ? s : s << ' ') << m.ctbninternal::RowMatrix::coeff(i,j);
		s << std::endl;
	}
	return s;
}

inline std::istream &operator>>(std::istream &s, matrix &m) {
	int tm,tn;
	s >> tm >> tn;
	m.ctbninternal::RowMatrix::resize(tm,tn);
	for(int i=0;i<tm;i++)
		for(int j=0;j<tn;j++)
			s >> m.ctbninternal::RowMatrix::coeffRef(i,j);
	return s;
}

} // end of ctbn namespace

namespace ctbn {

class spmatrix : public Eigen::SparseMatrix<double,Eigen::RowMajor> {
public:
	typedef Eigen::SparseMatrix<double,Eigen::RowMajor> Base;

	spmatrix(const matrix &ma)
		: Eigen::SparseMatrix<double,Eigen::RowMajor>(ma.getm(),ma.getn()) {
		int m = ma.getm(), n = ma.getn();
		for(int i=0;i<m;i++) {
			startVec(i);
			for(int j=0;j<n;j++)
				if (ma.coeff(i,j) > 1e-16 || ma.coeff(i,j) < -1e-16)
					insertBack(i,j) = ma.coeff(i,j);
		}
		finalize();
	}

	inline const Eigen::SparseDenseProductReturnType<Base,ctbninternal::RowMatrix>::Type rightmult(const matrix &ma) const {
		return (*this)*ma;
	}

	inline void leftmult(const vectr &in, vectr &out) const {
		out = in*(*this);
	}
};

} // end of ctbn namespace

#endif
