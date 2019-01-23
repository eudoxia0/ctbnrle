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
#include "utils.h"
#include "trajectory.h"
#include "context.h"
#include "random.h"
#if defined __linux__ || __CYGWIN__
#include <sys/time.h>
#include <sys/resource.h>
#endif
#include <cmath>
#include "params.h"



namespace ctbn {

using namespace std;

bool IsTrajSame(const Trajectory &tr1, const Trajectory &tr2, const Context &c) {
	Trajectory::Index i1 = tr1.Begin(c);
	Trajectory::Index i2 = tr2.Begin(c);
	Instantiation inst1, inst2;
	for (;;) {
		if (i1.Done() || i2.Done())
			break;
		inst1 = i1.Values();
		inst2 = i2.Values();
		Context diff(inst1, inst2);
		if (diff.VarList().size() > 0)
			return false;
		i1++;
		i2++;
	}
	return true;
}


bool IsTrajValid(const Trajectory &tr, const Context &c) {
	Trajectory::Index i = tr.Begin(c);
	Instantiation oldinst, newinst;
	oldinst = i.Values();
	for (;;) {
		++i;
		if (i.Done())
			break;
		newinst = i.Values();
		Context diff(oldinst, newinst);
		oldinst = newinst;
		if (diff.VarList().size() != 1)
			return false;
	}
	return true;
}

void RemoveInformation(Trajectory &tr, const Context &c, int var, double start, double end) {
	Trajectory newtr;
	newtr.SetBeginTime(tr.TimeBegin());
	newtr.SetEndTime(tr.TimeEnd());
	for(Trajectory::Index i = tr.Begin(c);!i.Done();++i) {
		Instantiation inst(i.Values());
		if (i.Time()>=start && i.Time()<end)
			inst.SetVal(var,-1);
		newtr.AddTransition(inst,i.Time());
		if (i.Time()<start && i.Time()+i.DeltaT()>start) {
			inst.SetVal(var,-1);
			newtr.AddTransition(inst,start);
		}
		if (i.Time()>start && i.Time()<end && i.Time()+i.DeltaT()>=end) {
			inst.SetVal(var,i.Values().Value(var));
			newtr.AddTransition(inst,end);
		}
	}
	tr = newtr;
}

void RemoveInformation(Trajectory &tr, const Context &c, int nvar, int nit, double frac) {
	double len = tr.TimeEnd()-tr.TimeBegin();
	for(int i=0;i<nit;i++) {
		int var = randomizer.RandInt(nvar);
		double st = randomizer.RandReal()*len*(1-frac)+tr.TimeBegin();
		RemoveInformation(tr, c, var, st, st+len*frac);
		//                cout << "var: " << var << ", st: " << st << ", end: " << st+len*frac << endl;
	}
}

void RemoveNodesInformation(Trajectory &tr, const Context &c, int nvar, int nit, double frac) {
  int excludenode = ParamInt("ExcludeNode", -1);
	double len = tr.TimeEnd()-tr.TimeBegin();
	for(int i=0;i<nit;i++) {
                int var = i % nvar;//randomizer.RandInt(nvar);
                if(var==excludenode) continue;
                //cout << "var: " << var << endl;
		double st = randomizer.RandReal()*len*(1-frac)+tr.TimeBegin();
		RemoveInformation(tr, c, var, st, st+len*frac);
		//                cout << "var: " << var << ", st: " << st << ", end: " << st+len*frac << endl;
	}
}

void SelectNodesInformation(Trajectory &tr, const Context &c, int nvar, double thresh, double frac) {
	int excludenode = ParamInt("ExcludeNode", -1);
	double begint = tr.TimeBegin();
	double endt = tr.TimeEnd();
	for(int var=0; var<nvar; var++) {
		priority_queue<double> select;
		double sum = SelectTime(select, frac, begint, endt, thresh);
		double res = sum - (endt-begint)*thresh;
		double seg = (endt-begint)*frac;
		while(!select.empty()) {
			double st = -select.top();
			double nextt = st + seg;
			select.pop();
			if(select.empty()) nextt -= res;
			RemoveInformation(tr, c, var, st, nextt);
		}
	}
}


double SelectTime(priority_queue<double> &select, double ratio, double begint, double endt, double thresh) {
	priority_queue<double> Q;
	double sum=0;
	double len = endt-begint;
	double seg = ratio*len;
	while(sum<thresh*len) {
		double st = randomizer.RandReal()*(len-seg) + begint;
		Q.push(-st);
		//cout << st << endl;
		sum = 0.0;
		priority_queue<double> tmpq(Q);
		double top = -tmpq.top();
		tmpq.pop();
		while(!tmpq.empty()) {
			double next = -tmpq.top();
			if(next-top<seg) sum += next-top;
			else sum += seg;
			top = next;
			tmpq.pop();
		}
		sum += seg;
		//cout << "sum=" << sum << endl;
	}
	select = Q;
	return sum;
}

double getcputime(void) {
	#if defined __linux__ || __CYGWIN__
	struct timeval tim;
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	tim=ru.ru_utime;
	double t=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
	tim=ru.ru_stime;
	t+=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
	return t;
	#elif defined _WIN32
	return 0.0;
	#endif
}

void Normalize(vector<double> &logw) {
	double maxweight = logw[0];
	for(unsigned int i = 1; i < logw.size(); i++)
		if(logw[i] > maxweight) maxweight = logw[i];
	double sumweight = 0;
	for(unsigned int i = 0; i < logw.size(); i++) {
		logw[i] = exp(logw[i] - maxweight);
		sumweight += logw[i];
	}
	for(unsigned int i = 0; i < logw.size(); i++)
		logw[i] /= sumweight;
	for (unsigned int i = 0; i < logw.size(); i++)
		logw[i] = log(logw[i]);
}

double mean(const vector<double> &v) {
	double ret = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		ret += v[i];
	return ret/v.size();
}

double relative_bias(const vector<double> &v, double vtrue) {
	double ret = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		ret += v[i];
	ret /= v.size();
	return fabs(ret-vtrue)/vtrue;
}

double relative_std(const vector<double> &v, double vtrue) {
	double avg = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		avg += v[i];
	avg /= v.size();
	double ret = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		ret += (v[i] - avg)*(v[i] - avg);
	return sqrt(ret/v.size())/vtrue;
}

} // end of ctbn namespace
