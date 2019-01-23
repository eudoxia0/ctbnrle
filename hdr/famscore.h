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
#ifndef CTBNRLE_FAMSCORE_H
#define CTBNRLE_FAMSCORE_H

#include "suffstatsquery.h"
#include "context.h"



namespace ctbn {

// Classes to handle the scoring done in structure search.
// Each class takes in the relevant prior parameters and a 
// SuffStatsQuery object that provides the sufficient statistics
// for scoring.

class FamScore {
  public:
  	// The pointer is owned by the caller of the constructor
     // but the pointer must remain valid for the life of the FamScore object
	FamScore(double nTrans, double aTime, SuffStatsQuery* ssQ);
	virtual ~FamScore();
	virtual double GetScore(const Context& v, const Context& cv) = 0;
  protected:
	double numTrans;
	double amtTime;

	// This class does not own this pointer
	SuffStatsQuery* ssQuery;
};

class BNFamScore : public FamScore {
  public:
  	// The pointer is owned by the caller of the constructor
     // but the pointer must remain valid for the life of the object
	BNFamScore(double nTrans, SuffStatsQuery* ssQ);
	virtual ~BNFamScore();
	virtual double GetScore(const Context& v, const Context& cv);
};

class CTBNFamScore : public FamScore {
  public:
  	// The pointer is owned by the caller of the constructor
     // but the pointer must remain valid for the life of the object
	CTBNFamScore(double nTrans, double aTime, SuffStatsQuery* ssQ);
	virtual ~CTBNFamScore();
	virtual double GetScore(const Context& v, const Context& cv);
};

} // end of ctbn namespace

#endif
