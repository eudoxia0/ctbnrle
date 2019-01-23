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
#include "markovdyn.h"



namespace ctbn {

SOBJCLASSDEF(MarkovDyn)

MarkovDyn::MarkovDyn(const Context &var, const Context &cvar, bool is_toggle)  : DynComp(var,cvar,(is_toggle? tmbase:mbase)) {
}

MarkovDyn::MarkovDyn(std::istream &is) : DynComp(is) {
}

MarkovDyn::~MarkovDyn() {
}

MarkovDyn *MarkovDyn::Clone() const {
	return new MarkovDyn(*this);
}

const MarkovSimple *MarkovDyn::mbase = new MarkovSimple();
const MarkovSimpleToggle *MarkovDyn::tmbase = new MarkovSimpleToggle();

} // end of ctbn namespace
