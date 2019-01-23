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
#ifndef CTBNRLE_TRAJSIMPLE_H_
#define CTBNRLE_TRAJSIMPLE_H_

#include "trajectory.h"

namespace ctbn {

class TrajSimple {
public:
	virtual ~TrajSimple();

	virtual double TimeBegin() const = 0;
	virtual double TimeEnd() const = 0;

	virtual void SetBeginTime(double t) = 0;
	virtual void SetEndTime(double t) = 0;

	virtual int Value(double t, bool inclusive=true) const = 0;

	virtual void AddTransition(int val, double t) = 0;

	class Index {
	public:
		virtual ~Index();

		virtual void Next() = 0;
		virtual double Time() const = 0;
		virtual double DeltaT() const = 0;
		virtual int Value() const = 0;
		virtual bool Done() const = 0;
	};

	virtual Index *Begin() const = 0;
	virtual Index *Begin(double startt, double endt) const = 0;
	virtual Index *End() const = 0;
	virtual Index *End(double startt, double endt) const = 0;
};


class TrajSimpleMap : public TrajSimple {
public:
	TrajSimpleMap();
	virtual ~TrajSimpleMap();

	virtual double TimeBegin() const { return ts; }
	virtual double TimeEnd() const { return te; }
	virtual void SetBeginTime(double t) { ts = t; }
	virtual void SetEndTime(double t) { te = t; }

	virtual int Value(double t, bool inclusive=true) const;

	virtual void AddTransition(int val, double t);

	class Index : public TrajSimple::Index {
	public:
		Index(const TrajSimpleMap &tsm, double t0, double tend);
		virtual ~Index();
		virtual void Next();
		virtual double Time() const;
		virtual double DeltaT() const;
		virtual int Value() const;
		virtual bool Done() const;

	private:
		double endt,t;
		int val;
		bool done;
		VarTrajectory::const_iterator loc,end;
	};

	class RIndex : public TrajSimple::Index {
	public:
		RIndex(const TrajSimpleMap &tsm, double t0, double tstart);
		virtual ~RIndex();
		virtual void Next();
		virtual double Time() const;
		virtual double DeltaT() const;
		virtual int Value() const;
		virtual bool Done() const;

	private:
		double startt,t;
		int val;
		bool done;
		VarTrajectory::const_reverse_iterator loc,end;
	};

	virtual Index *Begin() const;
	virtual Index *Begin(double startt, double endt) const;
	virtual RIndex *End() const;
	virtual RIndex *End(double startt, double endt) const;

protected:
	double ts,te;
	VarTrajectory tr;
};

// one could also implement this with a non-const Trajectory &
// in which case the "AddTransition" and "Set.*Time" methods would work
// this is a wrapper around a trajectory to make it look like
// just a flat sequence of numbers (timestamped) constrained to the
// indexes of the Context vars
class ConstFlatTraj : public TrajSimple {
public:
	ConstFlatTraj(const Trajectory &traj, const Context &vars);
	virtual ~ConstFlatTraj();

	virtual double TimeBegin() const;
	virtual double TimeEnd() const;

	virtual void SetBeginTime(double t);
	virtual void SetEndTime(double t);

	virtual int Value(double t, bool inclusive=true) const;

	virtual void AddTransition(int val, double t);

	class Index : public TrajSimple::Index {
	public:
		Index(const ConstFlatTraj &ft, double tstart, double tend);
		virtual ~Index();
		virtual void Next();
		virtual double Time() const;
		virtual double DeltaT() const;
		virtual int Value() const;
		virtual bool Done() const;
	private:
		Trajectory::Index impl;
	};
	class RIndex : public TrajSimple::Index {
	public:
		RIndex(const ConstFlatTraj &ft, double tstart, double tend);
		virtual ~RIndex();
		virtual void Next();
		virtual double Time() const;
		virtual double DeltaT() const;
		virtual int Value() const;
		virtual bool Done() const;
	private:
		Trajectory::RIndex impl;
	};

	virtual Index *Begin() const;
	virtual Index *Begin(double startt, double endt) const;
	virtual RIndex *End() const;
	virtual RIndex *End(double startt, double endt) const;

protected:
	const Trajectory &tr;
	Context c;
};

// and here is such an implementation...
class FlatTraj : public TrajSimple {
public:
	FlatTraj(Trajectory &traj, const Context &vars);
	virtual ~FlatTraj();

	virtual double TimeBegin() const;
	virtual double TimeEnd() const;

	virtual void SetBeginTime(double t);
	virtual void SetEndTime(double t);

	virtual int Value(double t, bool inclusive=true) const;

	virtual void AddTransition(int val, double t);

	class Index : public TrajSimple::Index {
	public:
		Index(const FlatTraj &ft, double tstart, double tend);
		virtual ~Index();
		virtual void Next();
		virtual double Time() const;
		virtual double DeltaT() const;
		virtual int Value() const;
		virtual bool Done() const;
	private:
		Trajectory::Index impl;
	};
	class RIndex : public TrajSimple::Index {
	public:
		RIndex(const FlatTraj &ft, double tstart, double tend);
		virtual ~RIndex();
		virtual void Next();
		virtual double Time() const;
		virtual double DeltaT() const;
		virtual int Value() const;
		virtual bool Done() const;
	private:
		Trajectory::RIndex impl;
	};

	virtual Index *Begin() const;
	virtual Index *Begin(double startt, double endt) const;
	virtual RIndex *End() const;
	virtual RIndex *End(double startt, double endt) const;

protected:
	Trajectory &tr;
	Context c;
};

} // end of ctbn namespace

#endif

