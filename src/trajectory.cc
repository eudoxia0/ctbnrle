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
#include "trajectory.h"

namespace ctbn {

using namespace std;

Trajectory::Trajectory() {
	ts = te = 0;
}

Trajectory::~Trajectory() {
}

int Trajectory::Value(int varid, double time, bool inclusive) const {
	if (time>te || (time==te && inclusive)) return -1;
	if (time<ts || (time==ts && !inclusive)) return -1;
	ivmap::const_iterator vari = traj.find(varid);
	if (vari==traj.end()) return -1;
	return Value(vari,time,inclusive);
}

int Trajectory::Value(const ivmap::const_iterator &vari,
				double time, bool inclusive) const {
	VarTrajectory::const_reverse_iterator i(vari->second.upper_bound(time));
	while (i!=vari->second.rend()
			&& (i->first>time || (!inclusive && i->first==time))) {
		++i;
	}
	return i==vari->second.rend() ? -1 : i->second;
}

Instantiation Trajectory::Values(const Context &c,
			double time, bool inclusive) const {
	Instantiation ret(c);
	if (time>te || (time==te && inclusive)) return ret;
	if (time<ts || (time==ts && !inclusive)) return ret;
	for(ivmap::const_iterator vari=traj.begin();
			vari!=traj.end();vari++)
		ret.SetVal(vari->first,Value(vari,time,inclusive));
	return ret;
}

void Trajectory::AddTransition(int varid, double time, int newval) {
	if (ts==te && traj.empty()) ts=te=time;
	else if (ts>time) ts = time;
	else if (te<time) te = time;

	ivmap::iterator vari = traj.find(varid);
	if (vari==traj.end())
		vari = traj.insert(make_pair(varid,VarTrajectory())).first;
	vari->second.insert(make_pair(time,newval));
}

void Trajectory::AddTransition(const Instantiation &i, double time) {
	Instantiation oldi = Values(i,time);
	vector<int> varlist = i.VarList();
	for(unsigned int j=0;j<varlist.size();j++) {
		int newvalue = i.Value(varlist[j]);
		if (oldi.Value(varlist[j])!=newvalue) 
			AddTransition(varlist[j],time,newvalue);
	}
}

void Trajectory::AddPointEvidence(int varid, double time, int value,
			double deltime) {
	// note that if prevval!=-1 this is really a strange 
	// function to call!
	int prevval = Value(varid,time);
	AddTransition(varid,time,value);
	AddTransition(varid,time+deltime,prevval);
}

void Trajectory::SaveOld(ostream &os) const {
	os << ts << os.fill() << te << os.fill() << traj;
}

void Trajectory::LoadOld(istream &is) {
	is >> ts >> te >> traj;
}

ostream &operator<<(ostream &os, const Trajectory &tr) {
	tr.SaveOld(os);
	return os;
}

istream &operator>>(istream &is, Trajectory &tr) {
	tr.LoadOld(is);
	return is;
}

Trajectory::Index::~Index() {
}

Trajectory::Index::Index(const Trajectory &tr, const Context &c, double t0, double tend) {
	done = false;
	t = t0;
	endt = tend;
	v = tr.Values(c,t0);
	for(ivmap::const_iterator vi = tr.traj.begin();
					vi != tr.traj.end(); vi++) {
		if(!c.HasId(vi->first)) continue;

		VarTrajectory::const_iterator i = vi->second.upper_bound(t0);
		if (i!=vi->second.end()) {
			vari.insert(make_pair(vi->first,make_pair(i,vi->second.end())));
			if (i->first<endt)
				nextvar.push(make_pair(-i->first,vi->first));
		}
	}
}

EvidenceChangeEnum Trajectory::Index::TestInc(const Context & c) {

	EvidenceChangeEnum ret(unchanged_with_respect_to_context);

	if (nextvar.empty()) {
		done = true;
		t = endt;
		if (-1 != v.Index()) {
			ret = only_variable_observability_changed;
		}
		v.SetAllVal(-1);
		return ret; 
	}

	double newt = -nextvar.top().first;
	while(!nextvar.empty() && -nextvar.top().first == newt) {
		int vi = nextvar.top().second;
		nextvar.pop();

		std::map<int,std::pair<VarTrajectory::const_iterator,
		 VarTrajectory::const_iterator> >::iterator spot = vari.find(vi);
		VarTrajectory::const_iterator &vti = spot->second.first;
		VarTrajectory::const_iterator &vtiend = spot->second.second;

		if (c.HasId(vi)) {
			if (vti->second != -1 && v.Value(vi) != -1) {
				ret = value_of_an_observed_variable_changed;
			} else if (ret == unchanged_with_respect_to_context) {
				ret = only_variable_observability_changed;
			}
		}

		v.SetVal(vi, vti->second);
		++vti;

		if (vti!=vtiend && vti->first<endt) {
			nextvar.push(make_pair(-vti->first, vi));
		}
	}
	t = newt;
	return ret;
}

Trajectory::RIndex::~RIndex() {
}

Trajectory::RIndex::RIndex(const Trajectory &tr, const Context &c, double t0, double tstart) {
	done = false;
	t = t0;
	startt = tstart;
	v = tr.Values(c,t0,false);
	for(ivmap::const_iterator vi = tr.traj.begin();
					vi != tr.traj.end(); vi++) {
		if(!c.HasId(vi->first)) continue;
		VarTrajectory::const_reverse_iterator i(vi->second.lower_bound(t0));
		if (i!=vi->second.rend()) {
			vari.insert(make_pair(vi->first,make_pair(i,vi->second.rend())));
			if (i->first>startt)
				nextvar.push(make_pair(i->first,vi->first));
		}
	}
}

int Trajectory::RIndex::TestDec(const Context &c) {
	if (nextvar.empty()) {
		done = true;
		t = startt;
		int retval;
		if(v.Index()==-1) retval = 0;
		else retval = 1;
		v.SetAllVal(-1);
		return retval;
	}
	int ret = 0;
	double newt = nextvar.top().first;
	while(!nextvar.empty() && nextvar.top().first==newt) {
		int vi = nextvar.top().second;
		nextvar.pop();

		std::map<int,std::pair<VarTrajectory::const_reverse_iterator,
			VarTrajectory::const_reverse_iterator> >::iterator spot = vari.find(vi);
		VarTrajectory::const_reverse_iterator &vti = spot->second.first;
		VarTrajectory::const_reverse_iterator &vtiend = spot->second.second;
		++vti;

		if (c.HasId(vi)) {
			if (vti->second != -1 && v.Value(vi)!=-1) ret = 2;
			else if (ret==0) ret = 1;
		}

		v.SetVal(vi,vti->second);

		if (vti!=vtiend && vti->first>startt) {
			nextvar.push(make_pair(vti->first,vi));
		}
	}
	t = newt;
	return ret;
}

void Trajectory::Draw(ostream &os) const {
	ios::fmtflags oldflags = os.flags();
	os.setf(ios::fixed);
	Context all;
	os << "time \\ var   ";
	for(ivmap::const_iterator i=traj.begin();i!=traj.end();++i) {
		os.width(5);
		os << i->first << "   ";
		all.AddVar(i->first,2);
	}
	os << endl;
	Instantiation oldv(all,-1);
	for(Index ii = Begin(all);!ii.Done();++ii) {
		//os << "_____________________________________________________________" << endl;
		os << ii.Time() << "     ";
		for(ivmap::const_iterator i=traj.begin();i!=traj.end();++i) {
			if (oldv.Value(i->first)!=ii.Values().Value(i->first)) {
				os.width(5);
				os << ii.Values().Value(i->first) << "   ";
			} else {
				os << "        ";
			}
		}
		os << endl;
		oldv = ii.Values();
	}
	os.flags(oldflags);
}

void Trajectory::SetUnknown(int varid, bool noinit) {
	traj.erase(varid);
	if (!noinit) AddTransition(varid, ts, -1);
}

void Trajectory::RemoveNodeTraj(int varid) {
	traj.erase(varid);
}

void Trajectory::ReplaceNodeTraj(int varid, const Trajectory &tr)
{
	map<int, VarTrajectory>::const_iterator iter = tr.traj.find(varid);
	if (iter != tr.traj.end())
		traj[varid] = iter->second;
}

Trajectory Trajectory::ExtractNodeTraj(int varid) const {
	map<int, VarTrajectory>::const_iterator iter = traj.find(varid);
	Trajectory tr;
	tr.ts = this->ts;
	tr.te = this->te;
	if(iter != traj.end()) tr.traj[varid] = iter->second;
	return tr;
}
} // end of ctbn namespace
