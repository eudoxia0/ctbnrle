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
#include "trajsimple.h"

#include <cassert>

namespace ctbn {

using namespace std;

TrajSimple::~TrajSimple() {
}

TrajSimple::Index::~Index() {
}

TrajSimpleMap::TrajSimpleMap() {
	te = ts = 0.0;
}

TrajSimpleMap::~TrajSimpleMap() {
}

int TrajSimpleMap::Value(double t, bool inclusive) const {
	if (t>te || (t==te && inclusive)) return -1;
	if (t<ts || (t==ts && !inclusive)) return -1;
	VarTrajectory::const_reverse_iterator i(tr.upper_bound(t));
	while(i!=tr.rend() && (i->first>t || (!inclusive && i->first==t)))
		++i;
	return i==tr.rend() ? -1 : i->second;
}

void TrajSimpleMap::AddTransition(int val, double t) {
	if (ts==te && tr.empty()) ts=te=t;
	else if (ts>t) ts = t;
	else if (te<t) te = t;
	tr[t] = val;
}

TrajSimpleMap::Index::Index(const TrajSimpleMap &tsm, double t0, double tend) {
	endt = tend;
	t = t0;
	end = tsm.tr.end();
	loc = tsm.tr.upper_bound(t0);
	VarTrajectory::const_reverse_iterator i(loc);
	++i;
	val = i->second;
	done = false;
}

TrajSimpleMap::RIndex::RIndex(const TrajSimpleMap &tsm, double t0, double tstart) : loc(tsm.tr.lower_bound(t0)) {
	startt = tstart;
	t = t0;
	end = tsm.tr.rend();
	val = loc.base()->second;
	done = false;
}

TrajSimpleMap::Index::~Index() {
}

TrajSimpleMap::RIndex::~RIndex() {
}

void TrajSimpleMap::Index::Next() {
	if (loc!=end) {
		t = loc->first;
		val = loc->second;
		++loc;
		if (loc->first>=endt) loc=end;
	} else {
		done = true;
		t = endt;
	}
}

void TrajSimpleMap::RIndex::Next() {
	if (loc!=end) {
		t = loc->first;
		++loc;
		val = loc->second;
		if (loc->first<startt) loc=end;
	} else {
		done = true;
		t = startt;
	}
}

double TrajSimpleMap::Index::Time() const {
	return t;
}

double TrajSimpleMap::RIndex::Time() const {
	return t;
}

double TrajSimpleMap::Index::DeltaT() const {
	return loc==end ? endt-t : loc->first-t;
}

double TrajSimpleMap::RIndex::DeltaT() const {
	return loc==end ? t-startt : t-loc->first;
}

int TrajSimpleMap::Index::Value() const {
	return done ? -1 : val;
}

int TrajSimpleMap::RIndex::Value() const {
	return done ? -1 : val;
}

bool TrajSimpleMap::Index::Done() const {
	return done;
}

bool TrajSimpleMap::RIndex::Done() const {
	return done;
}

TrajSimpleMap::Index *TrajSimpleMap::Begin() const {
	return new TrajSimpleMap::Index(*this,ts,te);
}

TrajSimpleMap::Index *TrajSimpleMap::Begin(double startt, double endt) const {
	return new TrajSimpleMap::Index(*this,startt,endt);
}

TrajSimpleMap::RIndex *TrajSimpleMap::End() const {
	return new TrajSimpleMap::RIndex(*this,te,ts);
}

TrajSimpleMap::RIndex *TrajSimpleMap::End(double startt, double endt) const {
	return new TrajSimpleMap::RIndex(*this,endt,startt);
}

//----------

ConstFlatTraj::ConstFlatTraj(const Trajectory &traj, const Context &vars) :
	tr(traj), c(vars) {
}

ConstFlatTraj::~ConstFlatTraj() {
}

double ConstFlatTraj::TimeBegin() const {
	return tr.TimeBegin();
}

double ConstFlatTraj::TimeEnd() const {
	return tr.TimeEnd();
}

void ConstFlatTraj::SetBeginTime(double t) {
	assert(0); // not possible
}

void ConstFlatTraj::SetEndTime(double t) {
	assert(0); // not possible
}

int ConstFlatTraj::Value(double t, bool inclusive) const {
	return tr.Values(c,t,inclusive).Index();
}

void ConstFlatTraj::AddTransition(int val, double t) {
	assert(0); // not possible
}

ConstFlatTraj::Index::Index(const ConstFlatTraj &ft,
		double tstart, double tend) : impl(ft.tr,ft.c,tstart,tend) {
}

ConstFlatTraj::RIndex::RIndex(const ConstFlatTraj &ft,
		double tstart, double tend) : impl(ft.tr,ft.c,tend,tstart) {
}

ConstFlatTraj::Index::~Index() {
}

ConstFlatTraj::RIndex::~RIndex() {
}

void ConstFlatTraj::Index::Next() {
	++impl;
}

void ConstFlatTraj::RIndex::Next() {
	--impl;
}

double ConstFlatTraj::Index::Time() const {
	return impl.Time();
}

double ConstFlatTraj::RIndex::Time() const {
	return impl.Time();
}

double ConstFlatTraj::Index::DeltaT() const {
	return impl.DeltaT();
}

double ConstFlatTraj::RIndex::DeltaT() const {
	return impl.DeltaT();
}

int ConstFlatTraj::Index::Value() const {
	return impl.Values().Index();
}

int ConstFlatTraj::RIndex::Value() const {
	return impl.Values().Index();
}

bool ConstFlatTraj::Index::Done() const {
	return impl.Done();
}

bool ConstFlatTraj::RIndex::Done() const {
	return impl.Done();
}

ConstFlatTraj::Index *ConstFlatTraj::Begin() const {
	return new ConstFlatTraj::Index(*this,TimeBegin(),TimeEnd());
}

ConstFlatTraj::Index *ConstFlatTraj::Begin(double startt, double endt) const {
	return new ConstFlatTraj::Index(*this,startt,endt);
}

ConstFlatTraj::RIndex *ConstFlatTraj::End() const {
	return new ConstFlatTraj::RIndex(*this,TimeBegin(),TimeEnd());
}

ConstFlatTraj::RIndex *ConstFlatTraj::End(double startt, double endt) const {
	return new ConstFlatTraj::RIndex(*this,startt,endt);
}

//--------

FlatTraj::FlatTraj(Trajectory &traj, const Context &vars) :
	tr(traj), c(vars) {
}

FlatTraj::~FlatTraj() {
}

double FlatTraj::TimeBegin() const {
	return tr.TimeBegin();
}

double FlatTraj::TimeEnd() const {
	return tr.TimeEnd();
}

void FlatTraj::SetBeginTime(double t) {
	tr.SetBeginTime(t);
}

void FlatTraj::SetEndTime(double t) {
	tr.SetEndTime(t);
}

int FlatTraj::Value(double t, bool inclusive) const {
	return tr.Values(c,t,inclusive).Index();
}

void FlatTraj::AddTransition(int val, double t) {
	Instantiation i(c);
	i.SetIndex(val);
	tr.AddTransition(i,t);
}

FlatTraj::Index::Index(const FlatTraj &ft,
		double tstart, double tend) : impl(ft.tr,ft.c,tstart,tend) {
}

FlatTraj::RIndex::RIndex(const FlatTraj &ft,
		double tstart, double tend) : impl(ft.tr,ft.c,tend,tstart) {
}

FlatTraj::Index::~Index() {
}

FlatTraj::RIndex::~RIndex() {
}

void FlatTraj::Index::Next() {
	++impl;
}

void FlatTraj::RIndex::Next() {
	--impl;
}

double FlatTraj::Index::Time() const {
	return impl.Time();
}

double FlatTraj::RIndex::Time() const {
	return impl.Time();
}

double FlatTraj::Index::DeltaT() const {
	return impl.DeltaT();
}

double FlatTraj::RIndex::DeltaT() const {
	return impl.DeltaT();
}

int FlatTraj::Index::Value() const {
	return impl.Values().Index();
}

int FlatTraj::RIndex::Value() const {
	return impl.Values().Index();
}

bool FlatTraj::Index::Done() const {
	return impl.Done();
}

bool FlatTraj::RIndex::Done() const {
	return impl.Done();
}

FlatTraj::Index *FlatTraj::Begin() const {
	return new FlatTraj::Index(*this,TimeBegin(),TimeEnd());
}

FlatTraj::Index *FlatTraj::Begin(double startt, double endt) const {
	return new FlatTraj::Index(*this,startt,endt);
}

FlatTraj::RIndex *FlatTraj::End() const {
	return new FlatTraj::RIndex(*this,TimeBegin(),TimeEnd());
}

FlatTraj::RIndex *FlatTraj::End(double startt, double endt) const {
	return new FlatTraj::RIndex(*this,startt,endt);
}

} //end of ctbn namespace
