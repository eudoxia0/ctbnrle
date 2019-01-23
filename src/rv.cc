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
#include "rv.h"



namespace ctbn {

using namespace std;

RV::RV(const Context &var, const Context &cvar) :
			StreamObj(), v(var), cv(cvar) {
}

RV::RV(istream &is) : StreamObj(), v(is), cv(is) {
}

RV::~RV() {
}

void RV::LoadOld(istream &is) {
	v.LoadOld(is);
	cv.LoadOld(is);
}

void RV::SaveOld(ostream &os) const {
	v.SaveOld(os);
	os << os.fill();
	cv.SaveOld(os);
}

SS* RV::SuffStats(const vector<Trajectory> &tr,
                      const vector<double> &w) const {
        SS* ss = this->BlankSS();
	const Context &c = this->Domain() + this->CondDomain();
        for(unsigned int i=0; i<tr.size(); i++) {
                Instantiation curri = tr[i].Values(c, tr[i].TimeBegin());
                this->AddSS(curri, ss, w[i]);
        }
        return ss;
}

} // end of ctbn namespace
