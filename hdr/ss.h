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
#ifndef CTBNRLE_SS_H
#define CTBNRLE_SS_H

#include "streamserial.h"


namespace ctbn {

// ABC for Sufficient Statistics
// Really just a blank class as what goes in here really depends
//  greatly and the methods are all part of the relevant RVs and Processes

class SS : public StreamObj {
public:
	SS();
	SS(std::istream &is);
	virtual ~SS();

	// virtual copy constructor
	virtual SS *Clone() const = 0;
	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
    virtual void Scale(double w) = 0;
    virtual void AddSS(const SS* nss, double w=1.0)=0;

	SERIAL_START_V(SS)
	SERIAL_END

};

} // end of ctbn namespace

#endif
