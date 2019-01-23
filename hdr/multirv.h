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
#ifndef CTBNRLE_MULTIRV_H
#define CTBNRLE_MULTIRV_H

#include "rvcomp.h"
#include "rvsimple.h"
#include "multisimple.h"

namespace ctbn {

// just a wrapper to make things simpler
// A general multinomial random variable can be described as
// a composition (a set) of RVSimple objects, one for each value
// of the conditioning set.  Each RVSimple is just a standard
// multinomial random variable
class MultiRV : public RVComp {
	SOBJCLASSDECL(MultiRV)
public:
	MultiRV(const Context &var=Context(), const Context &cvar=Context());
	MultiRV(std::istream &is);
	virtual ~MultiRV();

	virtual MultiRV *Clone() const;

	MultiZSimple &operator[](int i) {
		RVCondSimpleComp<MultiZSimple> *ptr =
			dynamic_cast<RVCondSimpleComp<MultiZSimple> *>(Base());
		return ptr->operator[](i);
	}
	const MultiZSimple &operator[](int i) const {
		const RVCondSimpleComp<MultiZSimple> *ptr =
			dynamic_cast<const RVCondSimpleComp<MultiZSimple> *>(Base());
		return ptr->operator[](i);
	}

	MultiZSimple &operator[](const Instantiation &i) {
		RVCondSimpleComp<MultiZSimple> *ptr =
			dynamic_cast<RVCondSimpleComp<MultiZSimple> *>(Base());
		return ptr->operator[](CondDomain().Index(i));
	}
	const MultiZSimple &operator[](const Instantiation &i) const {
		const RVCondSimpleComp<MultiZSimple> *ptr =
			dynamic_cast<const RVCondSimpleComp<MultiZSimple> *>(Base());
		return ptr->operator[](CondDomain().Index(i));
	}

	SERIAL_START_V(MultiRV)
		SERIAL_SUPER(RVComp)
	SERIAL_END
};

} // end of ctbn namespace

#endif
