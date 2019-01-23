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
#ifndef CTBNRLE_MARKOV_H
#define CTBNRLE_MARKOV_H

#include "process.h"
#include "dynamics.h"
#include "rv.h"
#include "ss.h"
#include "structure.h"
#include "extramath.h"



namespace ctbn {

class FBInf;

// See process.h for a full description of the method's uses
class Markov : public Process {
	friend class FBInf;
	SOBJCLASSDECL(Markov)
public:
	// this object now "owns" these two pointers
	Markov(RV * startdist = NULL, Dynamics * dyn = NULL);
	Markov(std::istream &is);
	Markov(const Markov &m);
	Markov &operator=(const Markov &m);
	virtual ~Markov();

	virtual Markov *Clone() const;

	virtual void Mult(const Process *p);
	virtual SS *BlankSS() const;
	virtual void Maximize(const SS *ss);

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	virtual void Sample(Trajectory &tr, Random &rand=randomizer) const;

	const Dynamics *GetDynamics() const { return d; }
	const RV *GetStartDist() const {return p0;}

	double LLH(const Trajectory &tr, double w, 
			double &p0llh, double &dynllh) const;
	double LLH(const std::vector<Trajectory> &tr, 
			const std::vector<double> &w) const; 
	virtual double LLH(const SS *ss) const;
	SS* SuffStats(const std::vector<Trajectory> &tr, 
			const std::vector<double> &w) const;
	virtual void Scramble(double a=1.0, double b=1.0, 
				double alpha=1.0, double degree=1.0, 
				Random &rand=randomizer);
	
	virtual double GetScore(double numTrans, double amtTime, 
				const SS* ss) const;

protected:
	RV * p0;
	Dynamics * d;

	SERIAL_START_V(Markov)
		SERIAL_SUPER(Process)
		SERIAL_VAR(RV *,p0)
		SERIAL_VAR(Dynamics *,d)
	SERIAL_END
public:
	void serial_preload();
};

class MarkovSS : public SS {
	friend class Markov;
	friend class MeanFieldInf;
	friend class FBInf;
	SOBJCLASSDECL(MarkovSS)
public:
	MarkovSS(std::istream &is);
	MarkovSS(SS *rvss=NULL, SS *dynss=NULL);
	MarkovSS(const MarkovSS &mss);
	MarkovSS &operator=(const MarkovSS &mss);

	virtual ~MarkovSS();

	virtual MarkovSS *Clone() const;

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;
	virtual void Scale(double w);
	virtual void AddSS(const SS* nss, double w=1.0);
	virtual double NodeSS(int id, int val1, int val2, 
						  const Instantiation &cond) const;
private:
	SS *p0ss,*dss;
	
	SERIAL_START_V(MarkovSS)
		SERIAL_SUPER(SS)
		SERIAL_VAR(SS*,p0ss)
		SERIAL_VAR(SS*,dss)
	SERIAL_END
public:
	void serial_preload();
};

} // end of ctbn namespace

#endif

