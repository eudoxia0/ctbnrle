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
#include "ctbn.h"
#include "nullptr03.h"


namespace ctbn {

using namespace std;

SOBJCLASSDEF(CTBN)

CTBN::CTBN(RV *startdist, Dynamics *dyn) : Markov(startdist, dyn) {
}

CTBN::CTBN(std::istream &is) : Markov(is) {

}

CTBN::CTBN(const Markov &m) : Markov(m) {

}

CTBN *CTBN::Clone() const {
	return new CTBN(*this);
}

CTBN &CTBN::operator=(const CTBN &ctbn) {
	if(&ctbn==this) return *this;
	if(p0 != nullptr03) delete p0;
	if(d != nullptr03) delete d;
	p0 = ctbn.p0->Clone();
	d = ctbn.d->Clone();
	return *this;
}

const MarkovDyn* CTBN::Node(int i) const {
	return dynamic_cast<const MarkovDyn*>(CTBNDynamic()->Node(i));
}

const CTBNDyn* CTBN::CTBNDynamic() const {
	return dynamic_cast<const CTBNDyn*>(d);
}

} // end of ctbn namespace

