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
#ifndef CTBNRLE_TRAJECTORY_H
#define CTBNRLE_TRAJECTORY_H

#include <map>
#include <vector>
#include <queue>
#include "context.h"
#include <iostream>

// This file describes a Trajectory... a mapping over the values
// of variables over time.  A Trajectory might be incomplete in that
// the values of certain variables are not known are certain times (or
// durations).  -1 (just like in context.h) represents an unknown
// value.  The variable IDs need not be consecutive, however, their
// values are assumed to be non-negative integers.
// Trajectories may start and end at any time (provided the former
// is before the latter).

// "point evidence" (evidence of the value of a variable at exactly
// one time instant but not before or after that) is not currently
// handled exactly.  Instead, a small window of time is used.
// See AddPointEvidence()

namespace ctbn {

// A single variable trajectory is a map from a time to
// the new value of the variable at that time
// (we might want to make this a sorted vector for speed 
//  purposes later)


typedef std::map<double, int> VarTrajectory;


enum EvidenceChangeEnum {
	unchanged_with_respect_to_context = 0,
	only_variable_observability_changed = 1,
	value_of_an_observed_variable_changed = 2
};


// A trajectory over multiple variables is...
class Trajectory {
public:
	Trajectory();
	virtual ~Trajectory();

	typedef std::map<int, VarTrajectory> ivmap;

	inline double TimeBegin() const { return ts; }
	inline double TimeEnd() const { return te; }

	inline void SetBeginTime(double t) { ts = t; }
	inline void SetEndTime(double t) { te = t; }

	int Value(int varid, double time, bool inclusive=true) const;
	Instantiation Values(const Context &c,
			double time, bool inclusive=true) const;

	// Adding a transition to -1 causes the variable to be unobserved
	void AddTransition(int varid, double time, int newval);
	void AddTransition(const Instantiation &i, double time);

	// adds point evidence.  This is currently handled by adding
	// a very short window of evidence (deltime below)
	void AddPointEvidence(int varid, double time, int value,
				double deltime = 1e-10);

	void SaveOld(std::ostream &os) const;
	void LoadOld(std::istream &is);

	// Removes all information about variable varid
	// noinit=false means that an initial value of -1 *will* be added
	// noinit=true means that no initial value (or any record) will remain
	void SetUnknown(int varid, bool noinit=false);
	// Same as SetUnknown, except that initial value is not set to -1
	// (i.e. variable is never even mentioned
	void RemoveNodeTraj(int varid);
	// replace varid's trajectory in *this with varid's trajectory in tr
	void ReplaceNodeTraj(int varid, const Trajectory &tr);
	Trajectory ExtractNodeTraj(int varid) const;

	// Index is a way of moving through a trajectory (similar to an
	// iterator).  It is constructed by calling Trajectory::Begin()
	class Index {
	public:
		~Index();

		// moves to the next transition of any variable
		Index & operator++() {
			TestInc(Context());
			return *this;
		}

		const Index operator++(int) {
			Trajectory::Index ret(*this);
			++(*this);
			return ret;
		}

		// moves to the next transition of the context
		void Inc(const Context &c) { while(!done && !TestInc(c)); }

		// returns the current "absolute" time
		double Time() const { return t; };
		// returns the time until the next event
		double DeltaT() const {
			return nextvar.empty() ? (endt-t) : (-nextvar.top().first-t);
		}
		// returned reference is invalid if Index is changed
		const Instantiation &Values() const { return v; }

		bool Done() const { return done; }
		// I really want these to be private... hmmm...
		Index(const Trajectory &tr, const Context &c, double t0, double tend);

		// increments to next transition and returns whether context
		// was changed
		// return == 0 if nothing in c is changed
		// return == 1 if only change is from known value to hidden 
		//    value (or the reverse)
		// return == 2 if change at least one variable in c changes
		//    from one known value to another
		EvidenceChangeEnum TestInc(const Context &c);
		int NextVar() const {return nextvar.empty() ? (-1) : (nextvar.top().second);}
		int NextTransitionVal() const;

	private:
		Instantiation v;
		double t,endt;
		// maps varid onto the pair of next trans and last trans
		std::map<int,std::pair<VarTrajectory::const_iterator,
		                  VarTrajectory::const_iterator> > vari;
		std::priority_queue<std::pair<double, int> > nextvar;
		bool done;
	};

	// a reverse index, build from Trajectory::End()
	class RIndex {
	public:
		~RIndex();

		// moves to the previous transition of any variable
		RIndex &operator--() { TestDec(Context()); return *this; }
		const RIndex operator--(int) {
			Trajectory::RIndex ret(*this);
			--(*this);
			return ret;
		}

		// moves to the next transition of the context
		void Dec(const Context &c) { while(!done && !TestDec(c)); }

		// returns the current "absolute" time
		double Time() const { return t; };
		// returns the time until the previous event
		double DeltaT() const {
			return nextvar.empty() ? (t-startt) : (t-nextvar.top().first);
		}
		// returned reference is invalid if Index is changed
		// Note that *despite* the trajectory being left continuous
		// this is the value just *before* the current transition
		const Instantiation &Values() const { return v; }

		bool Done() const { return done; }
		// I really want these to be private... hmmm...
		RIndex(const Trajectory &tr, const Context &c, double t0, double tstart);
		// increments to next transition and returns whether context
		// was changed
		// return == 0 if nothing in c is changed
		// return == 1 if only change is from known value to hidden 
		//    value (or the reverse)
		// return == 2 if change at least one variable in c changes
		//    from one known value to another
		int TestDec(const Context &c);

	private:
		Instantiation v;
		double t,startt;
		// maps varid onto the pair of next trans and last trans
		std::map<int,std::pair<VarTrajectory::const_reverse_iterator,
		                  VarTrajectory::const_reverse_iterator> > vari;
		std::priority_queue<std::pair<double, int> > nextvar;
		bool done;

	};

	// Return an Index (see above) over only the variables mentioned
	// in c
	Index Begin(const Context &c) const { return Index(*this,c,TimeBegin(),TimeEnd()); }
	// Same as above, but the Index only runs from startt to endt
	Index Begin(const Context &c, double startt, double endt) const { return Index(*this,c,startt,endt); }

	// Same as above, but run through the trajectory in reverse:
	RIndex End(const Context &c) const { return RIndex(*this,c,TimeEnd(),TimeBegin()); }
	// firstt < startt
	RIndex End(const Context &c, double startt, double firstt) const { return RIndex(*this,c,startt,firstt); }

	// makes a nicer looking (but unloadable) ASCII version
	void Draw(std::ostream &os) const;
	

private:
	double ts,te;

	ivmap traj; 

	int Value(const ivmap::const_iterator &vari,
					double time, bool inclusive=true) const;

	SERIAL_START(Trajectory)
		SERIAL_VAR(double,ts)
		SERIAL_VAR(double,te)
		SERIAL_VAR(ivmap,traj)
	SERIAL_END
};

std::ostream &operator<<(std::ostream &os, const Trajectory &tr);
std::istream &operator>>(std::istream &is, Trajectory &tr);

} // end of ctbn namespace

#endif
