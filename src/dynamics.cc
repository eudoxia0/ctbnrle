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
#include "dynamics.h"

namespace ctbn {

using namespace std;

Dynamics::Dynamics(const Context &var, const Context &cvar) :
	StreamObj(), v(var), cv(cvar) {
}

Dynamics::Dynamics(istream &is) : StreamObj(), v(is), cv(is) {
}

Dynamics::~Dynamics() {
}

void Dynamics::LoadOld(istream &is) {
	v.LoadOld(is);
	cv.LoadOld(is);
}

void Dynamics::SaveOld(ostream &os) const {
	v.SaveOld(os);
	os << os.fill();
	cv.SaveOld(os);
}

void Dynamics::SampleTrajectory(Trajectory &tr, double t, Random &rand) const {
	double endt = tr.TimeEnd();
	Instantiation newi(Domain());
	Instantiation i(tr.Values(Domain(),t));
	while(1) {
		double newt;
		SampleNextEvent(i,t,newt,newi,rand);
		if (newt>=endt) break;
		tr.AddTransition(newi,newt);
		t = newt;
		i = newi;
	}
}

void Dynamics::Mult(const Dynamics *x) {
	v = v+x->v;
	cv = (cv+x->cv)-v;
}

SS* Dynamics::SuffStats(const vector<Trajectory> &tr,
                       const vector<double> &w) const {
	SS *ss = BlankSS();
	const Context &c = this->Domain() + this->CondDomain();
	const Context &owndomain = Domain();
	for(unsigned int i=0; i<tr.size(); i++) {
		Instantiation curri = tr[i].Values(c, tr[i].TimeBegin());
		Trajectory::Index index = tr[i].Begin(c);
		while(!index.Done()) {
			//curri.SaveOld(cout); cout << endl;
			double t = index.Time();
			double deltat = index.DeltaT();
			this->AddSS(curri, t, deltat, ss, w[i]);
			
			int sign = index.TestInc(owndomain);
			if(!index.Done()) {
				Instantiation nexti = index.Values();
				t = index.Time();
				if (sign==2) 
					this->AddTransSS(curri,nexti,
									 t,ss,w[i]);
				curri = nexti;
			}
		}
	}
	return ss;
}

} // end of ctbn namespace
