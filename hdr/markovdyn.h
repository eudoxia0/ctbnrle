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
#ifndef CTBNRLE_MARKOVDYN_H
#define CTBNRLE_MARKOVDYN_H

#include "markovsimple.h"
#include "dyncomp.h"



namespace ctbn {

// just a wrapper to make things simpler
// A Markov Dynamics is a composition of a set of MarkovSimple processes
// (either toggle variables, or -- by default -- regular Markov processes)
class MarkovDyn : public DynComp {
	SOBJCLASSDECL(MarkovDyn)
  public:
	MarkovDyn(const Context &var=Context(), const Context &cvar=Context(), 
			bool is_toggle = false);
	MarkovDyn(std::istream &is);

	virtual ~MarkovDyn();

	virtual MarkovDyn *Clone() const;

	inline const MarkovSimple *operator()(const Instantiation &i) const {
		return dynamic_cast<const MarkovSimple *>((*this)[i]);
	}
	inline MarkovSimple *operator()(const Instantiation &i) {
		return dynamic_cast<MarkovSimple *>((*this)[i]);
	}

	inline const MarkovSimple *operator()(int i) const {
		return dynamic_cast<const MarkovSimple *>((*this)[i]);
	}
	inline MarkovSimple *operator()(int i) {
		return dynamic_cast<MarkovSimple *>((*this)[i]);
	}
  private:
	const static MarkovSimple *mbase;
	const static MarkovSimpleToggle *tmbase;

	SERIAL_START_V(MarkovDyn)
		SERIAL_SUPER(DynComp)
	SERIAL_END
};


} // end of ctbn namespace

#endif

