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
#ifndef CTBNRLE_CTBN_H
#define CTBNRLE_CTBN_H

#include "nullptr03.h"
#include "markov.h"
#include "structure.h" 
#include "markovdyn.h"
#include "bn.h"



namespace ctbn {

// A CTBN is just a Markov process in which the starting
// distribution is a BN and the dynamics are described by a CTBNDyn
// This class exists only for naming convenience
class CTBN : public Markov
{
	SOBJCLASSDECL(CTBN)
public:
	// startdist and dyn are now owned by this class!
	CTBN(RV * startdist = nullptr03, Dynamics * dyn = nullptr03);
	CTBN(std::istream &is);
	CTBN(const Markov &m);
	CTBN &operator=(const CTBN &ctbn);
	virtual ~CTBN() {};

	CTBN *Clone() const;

	// These are new methods for CTBNs only
	// They are just short cuts to the CTBNDyn
	const MarkovDyn* Node(int i) const;
	const CTBNDyn *CTBNDynamic() const;

	SERIAL_START_V(CTBN)
		SERIAL_SUPER(Markov)
	SERIAL_END
};
} // end of ctbn namespace

#endif
