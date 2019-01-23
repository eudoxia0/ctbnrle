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
#ifndef CTBNRLE_MATRIX_STREAM_H
#define CTBNRLE_MATRIX_STREAM_H

#include "matrix-eigen.h"

namespace SERIALNAMESPACE {

template<>
struct TypeInfo<ctbn::matrix,void> {
	inline static const char *namestr() { return "matrix"; }
	inline static void writeotherattr(std::ostream &os, const ctbn::matrix &m) {
		os << " m=\"" << m.getm() << "\" n=\"" << m.getn() << "\"";
	}
	inline static bool isshort(const ctbn::matrix &) { return false; }
	inline static bool isinline(const ctbn::matrix &) { return false; }
	inline static void save(const ctbn::matrix &m, std::ostream &os, int indent) {
		os << std::endl;
		for(int i=0;i<m.getm();i++) {
			Indent(os,indent+1);
			for(int j=0;j<m.getn();j++) 
				(j==0 ? os : os << ' ')
					<< m.ctbn::ctbninternal::RowMatrix::coeff(i,j);
			os << std::endl;
		}
		Indent(os,indent);
			
	}
	inline static void load(ctbn::matrix &m, const XMLTagInfo &info,
			std::istream &is) {
		int tm,tn;
		std::map<std::string,std::string>::const_iterator vi
			 = info.attr.find("m");
		if (vi==info.attr.end())
			throw streamexception("missing m-parameter in matrix tag");
		std::istringstream convm(vi->second);
		convm >> tm;
		vi = info.attr.find("n");
		if (vi==info.attr.end())
			throw streamexception("missing n-parameter in matrix tag");
		std::istringstream convn(vi->second);
		convn >> tn;
		m.ctbn::ctbninternal::RowMatrix::resize(tm,tn);
		for(int i=0;i<tm;i++) for(int j=0;j<tn;j++)
			is >> m.ctbn::ctbninternal::RowMatrix::coeffRef(i,j);
		XMLTagInfo einfo;
		ReadTag(is,einfo);
		if (einfo.name!=info.name || !einfo.isend || einfo.isstart)
			throw streamexception("missing end tag for matrix");
	}
};


template<>
struct TypeInfo<ctbn::vectr,void> {
	inline static const char *namestr() { return "vectr"; }
	inline static void writeotherattr(std::ostream &os, const ctbn::vectr &v) {
		os << " m=\"" << v.getm() << "\"";
	}
	inline static bool isshort(const ctbn::vectr &) { return false; }
	inline static bool isinline(const ctbn::vectr &) { return false; }
	inline static void save(const ctbn::vectr &v, std::ostream &os, int indent) {
		os << std::endl;
		Indent(os,indent+1);
		for(int i=0;i<v.getm();i++)
			(i==0 ? os : os << ' ')
				<< v.Eigen::VectorXd::coeff(i);
		os << std::endl;
		Indent(os,indent);
			
	}
	inline static void load(ctbn::vectr &v, const XMLTagInfo &info,
			std::istream &is) {
		int tm;
		std::map<std::string,std::string>::const_iterator vi
			 = info.attr.find("m");
		if (vi==info.attr.end())
			throw streamexception("missing m-parameter in vectr tag");
		std::istringstream convm(vi->second);
		convm >> tm;
		v.Eigen::VectorXd::resize(tm);
		for(int i=0;i<tm;i++) 
			is >> v.Eigen::VectorXd::coeffRef(i);
		XMLTagInfo einfo;
		ReadTag(is,einfo);
		if (einfo.name!=info.name || !einfo.isend || einfo.isstart)
			throw streamexception("missing end tag for vectr");
	}
};

}

#endif
