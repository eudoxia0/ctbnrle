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

#include "bn.h"
#include "nullptr03.h"
#include "rvsimple.h"


namespace ctbn {


using namespace std;

RVSimple::RVSimple() : StreamObj() {
}

RVSimple::RVSimple(istream &is) : StreamObj() {
}

RVSimple::~RVSimple() {
}

void RVSimple::LoadOld(istream &is) {
}

void RVSimple::SaveOld(ostream &os) const {
}


std::ostream & RVSimple::PrintJoint(std::ostream & out, RV * the_initial_distribution_of_the_markov_process) const {
    vectr v(the_initial_distribution_of_the_markov_process->Domain().Size(), 0.0);
    double logf;
    this->GetDist(v, logf);

    //Printing scaled to when max = 1
    vectr scaled = v/v.max();
    out <<"\nJoint Distribution:\t";
    for(int i=0; i<scaled.length(); i++){
        out <<"["<<i<<"]="<<scaled[i]<<"\t";
    }
    out << endl;
    //

    //Open this up once you're done.
    out << "Exact Logf: " << logf << endl;
    v.niceprint(out);

    return out;
}


std::ostream & RVSimple::PrintMarginal(std::ostream & out, RV * the_initial_distribution_of_the_markov_process) const {
    vectr v(the_initial_distribution_of_the_markov_process->Domain().Size(), 0.0);
    double logf;
    this->GetDist(v, logf);
    out << "Exact Logf: " << logf << endl;
    BN const * const initial_distribution(dynamic_cast<BN *>(the_initial_distribution_of_the_markov_process));
    int nn = (nullptr03 == initial_distribution) ? 0 : initial_distribution->NumofNodes();
//    int nn = dynamic_cast<BN*>(the_initial_distribution_of_the_markov_process)->NumofNodes();
    int div = 1;
    int ds, ids;
    vector<vectr> mdist;
    mdist.resize(nn);
    for(int i=0; i<nn; i++){
        ds = dynamic_cast<BN*>(the_initial_distribution_of_the_markov_process)->NodeByVar(i)->Domain().Size();
        mdist.at(i) = vectr(ds, 0.0);
    }
    for (int ind=0; ind<v.length(); ind++){
        div = 1;
        for (int i=0; i<nn; i++){
            ds = mdist.at(i).length();
            ids = (ind/div)%ds;
            (mdist.at(i))[ids] = (mdist.at(i))[ids] + v[ind];
            div *= ds;
        }
    }
    for(int i=0; i<nn; i++){
        mdist.at(i).normalize();
        (mdist.at(i)).niceprint(out);
    }

    return out;
}


std::ostream & RVSimple::Print(std::ostream & out, RV * the_initial_distribution_of_the_markov_process) const {
    out << "0 ";
    PrintMarginal(out, the_initial_distribution_of_the_markov_process);
    return out;
}

//------

RVCondSimple::RVCondSimple() : StreamObj() {
}

RVCondSimple::RVCondSimple(std::istream &is) : StreamObj() {
}

RVCondSimple::~RVCondSimple() {
}


std::ostream & PrintCond(std::ostream & out) {
    throw "Not yet implemented for this subtype of RVCondSimple.";
}


void RVCondSimple::LoadOld(istream &is) {
}

void RVCondSimple::SaveOld(ostream &os) const {
}

//------

SOBJCLASSDEF(RVCSCompSS)

RVCSCompSS::RVCSCompSS(int nc) : ss(nc) {
	for(int i=0;i<nc;i++) ss[i] = nullptr03;
}

RVCSCompSS::RVCSCompSS(const RVCSCompSS &x) : SS(), ss(x.ss.size()) {
	int nc = x.ss.size();
	for(int i=0;i<nc;i++) if (x.ss[i] != nullptr03) ss[i] = x.ss[i]->Clone();
}

RVCSCompSS::RVCSCompSS(istream &is) {
	int nc;
	is >> nc;
	for(int i=0;i<nc;i++) {
		SS *sptr = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
		ss.push_back(sptr);
	}
}

RVCSCompSS::~RVCSCompSS() {
	int nc = ss.size(); 
	for(int i=0;i<nc;i++)
		if (ss[i] != nullptr03) delete ss[i];
}

void RVCSCompSS::LoadOld(istream &is) {
	int nc = ss.size(); 
	for(int i=0;i<nc;i++)
		if (ss[i] != nullptr03) delete ss[i];
	is >> nc;
	ss.resize(nc);
	for(int i=0;i<nc;i++)
		ss[i] = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
}

void RVCSCompSS::SaveOld(ostream &os) const {
	int nc = ss.size();
	os << nc;
	for(int i=0;i<nc;i++) {
		os << os.fill();
		ss[i]->SaveOldPtr(os);
	}
}

RVCSCompSS *RVCSCompSS::Clone() const {
	return new RVCSCompSS(*this);
}

RVCSCompSS &RVCSCompSS::operator=(const RVCSCompSS &x) {
	if (&x!=this) {
		int nc = ss.size();
		for(int i=0;i<nc;i++)
			if (ss[i] != nullptr03) delete ss[i];
		nc = x.ss.size();
		ss.resize(nc);
		for(int i=0;i<nc;i++)
			if (x.ss[i] == nullptr03) ss[i] = nullptr03;
			else ss[i] = x.ss[i]->Clone();
	}
	return *this;
}

void RVCSCompSS::Scale(double w) {
	int n = ss.size();
	for(int i=0; i<n; i++)
		ss[i]->Scale(w);
}

void RVCSCompSS::AddSS(const SS* nss, double w) {
	const RVCSCompSS* rss = dynamic_cast<const RVCSCompSS*>(nss);
	int n = ss.size();
	for(int i=0; i<n; i++)
		ss[i]->AddSS(rss->ss[i], w);
}

} // end of ctbn namespace
