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

// This file is added into the "middle" of the class definition of
// MatrixBase from the Eigen package (done by the #def and #inc in
// matrix-eigen.h -- see Eigen matrix package

struct valid_visitor {
	bool valid;
	inline void init(const Scalar &, Index i, Index j) { valid = true; }
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

inline Scalar norm2() const {
	return squaredNorm();
}

