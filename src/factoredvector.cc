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
#include "factoredvector.h"
#include "matrix.h"
#include "ctbn.h"
#include "ctbndyn.h"
#include "context.h"
#include "multirv.h"
#include "notyetimplementederror.h"
#include "nullptr03.h"

#include <string>


namespace ctbn {

using namespace std;


static const string the_class_name("FactoredMatrix");


FactoredVector::FactoredVector() : RVSimple()
{
	logf = 0.0;
}

FactoredVector::FactoredVector(const FactoredVector &v) : RVSimple() {
	unsigned int n = v.dists.size();
	dists.resize(n);
	for (unsigned int i = 0; i != n; ++i) {
		dists[i] = v.dists.at(i);
	}
	logf = v.logf;
}

FactoredVector::FactoredVector(RV *rv) : RVSimple() {
	BN *bn = dynamic_cast<BN*>(rv); 
	int n = bn->NumofNodes();
	dists.resize(n);
	int domainsize; 
	logf = 0.0;
	for (int i=0; i<n; i++){
		MultiRV *newrv = (MultiRV*) bn->NodeByVar(i);
		domainsize = newrv->Domain().Size();
		vectr d(domainsize, -1.0);
		double logfactor;
		newrv->operator[](0).GetDist(d, logfactor);
		dists[i] = d;
		logf += logfactor;

	}
}

FactoredVector::~FactoredVector() {

}

FactoredVector * FactoredVector::Clone() const {

	return new FactoredVector(*this);
}

void FactoredVector::Condition(const map<int, int> &ev) {
	map<int, int>::const_iterator it;
	int id, val;
	double tmp;
	for (it = (ev).begin(); it != ev.end(); it++){
		id = (*it).first;
		val = (*it).second;
		RestrictById(id, val);
	}
	
}

void FactoredVector::Restrict(const std::vector<int, std::allocator<int> >&) {
	cout << "FactoredVector::Restrict not implemented" << endl;
}

void FactoredVector::RestrictById(int id, int val) {
	double tmp = dists.at(id)[val];
	dists.at(id) = 0.0;
	dists.at(id)[val] = tmp;
}

void FactoredVector::RestrictById(int id, int val, double x) {
	dists.at(id) = 0.0;
	dists.at(id)[val] = x;
}

void FactoredVector::MultBy(const RVSimple *x) {
	const FactoredVector *v = dynamic_cast<const FactoredVector *>(x);
	for(unsigned int i=0; i<dists.size(); i++) {
		dists[i].multby(v->dists[i]);
	}
}

void FactoredVector::MakeUniform() {
	for (unsigned int i=0; i<dists.size(); i++){
		dists.at(i) = 1.0;
	}
	logf = 0.0;
	
}

double FactoredVector::Normalize() {

	double mofsum, tmp, diff;
	double eps = 0.00001;

	tmp = dists[0].normalize();
	mofsum = tmp;
	for (unsigned int i=1; i<dists.size(); i++){
		tmp = dists[i].normalize();

		diff = tmp - mofsum;
		if (!((diff < eps)  && (diff > - eps))){
			mofsum *= tmp;

		}
	
	}
	double ret= exp(logf)*mofsum;
	logf = 0.0;
	return ret;
	
	/*
	double tmp = 1.0;	
	for (unsigned int i=0; i<dists.size(); i++){
		tmp *= dists[i].normalize();	
	}
	
	double ret = exp(logf)*tmp;
	logf = 0.0;
	return ret;
	*/
}

double FactoredVector::absnormalize(vectr &v) {
	double tmpsum = 0.0;
	for (int i=0; i<v.length(); i++){
		tmpsum += fabs(v[i]);
	}
	for (int i=0; i<v.length(); i++){
		v[i] = fabs(v[i]) / tmpsum;
	}
	return tmpsum;
}

double FactoredVector::maxnormalize(vectr &v) {
	double max = v[0];
	double min = v[0];
	for (int i=0; i<v.length(); i++){
		if (max < v[i])
			max = v[i];
		if (min > v[i])
			min = v[i];
	}
	for (int i=0; i<v.length(); i++){
		v[i] = (v[i] - min)/ (max-min);
	}
	return max-min;
}

double FactoredVector::NormalizeNeg(){
	double mofsum, diff;
//  double tmp;
	double eps = 0.00001;
	mofsum = 1.0;

	for (unsigned int i=0; i<dists.size(); i++) {
		mofsum *= dists.at(i).sum();
		for (int j=0; j<(dists[i]).length(); j++) {
			if ((dists[i])[j] < 0.0){

				dists[i][j] = 0.0;
			}
		}
		
//		tmp = dists[i].normalize();	
	}			


	double ret= exp(logf)*mofsum;

	logf = 0.0;
	return ret;
	
}


double FactoredVector::GetLogf() const {
	return logf;
}


void FactoredVector::SetLogf(double logfactor) {
	logf = logfactor;
}
	

double FactoredVector::Sum(bool islog) const {
	double nc = 0.0;
	for (unsigned int i=0; i<dists.size(); i++) {
		nc += log(dists[i].sum());
	}
	if (islog) {
	    return (logf+nc);
	} else {
		return (exp(logf+nc));
	}
	
}

//Assumes both vectors have the same length
void FactoredVector::AddMult(const FactoredVector & v, double x) {
	
/*
	const FactoredVector * a = dynamic_cast<const FactoredVector *>(&v);
	
	const FactoredVector *minfv, *maxfv; 
	double difflogf, maxlogf;
	double alogf = a->logf + log(x);
	if(logf < alogf){
		minfv = this;
		maxfv = a;
		difflogf = alogf - logf;
		maxlogf = alogf;
	}
	else{
		minfv = a;
		maxfv = this;
		difflogf = logf - alogf;
		maxlogf = logf;
	}

	double expf = exp(difflogf);
	for (unsigned int i=0; i<dists.size(); i++) {
		dists.at(i) = (minfv->dists.at(i))/expf + maxfv->dists.at(i);
		dists.at(i) = dists.at(i)/(1+1/expf);
	}
	logf = maxlogf + log(1+1/expf);
*/
/*
	if (logf!=0.0) {
		double m = exp(logf);
		for(unsigned int i=0;i<dists.size();i++)
			dists.at(i) *= m;
	}
	double m = exp(v.logf)*x;
	vector<double> sums(dists.size());
	double sprod = 1.0;
	for(unsigned int i=0;i<v.dists.size();i++)
		sprod *= (sums[i] = v.dists.at(i).sum());
	for(unsigned int i=0;i<dists.size();i++)
		dists.at(i) += v.dists.at(i)*m*sprod/sums[i];
	logf = 0.0;
*/
	if (logf!=0.0) {
		double m = exp(logf);
		for(unsigned int i=0;i<dists.size();i++)
			dists.at(i) *= m;
	}
	double m = exp(v.logf)*x;
	vector<double> sums(dists.size());
	for(unsigned int i=0;i<v.dists.size();i++)
		sums[i] = v.dists.at(i).sum();
	for(unsigned int i=0;i<dists.size();i++) {
		double sprod = m;
		for(unsigned int j=0;j<dists.size();j++)
			if (j!=i) sprod *= sums[i];
		dists.at(i) += v.dists.at(i)*sprod;
	}
	logf = 0.0;
}
	
void FactoredVector::Add(const RVSimple *x) {
	const FactoredVector * a = dynamic_cast<const FactoredVector *>(x);
	
	const FactoredVector *minfv, *maxfv; 
	double difflogf, maxlogf;
	if(logf < a->logf){
		minfv = this;
		maxfv = a;
		difflogf = a->logf - logf;
		maxlogf = a->logf;
		
	}
	else{
		minfv = a;
		maxfv = this;
		difflogf = logf - a->logf;
		maxlogf = logf;
	}

	double expf = exp(difflogf);
	for (unsigned int i=0; i<dists.size(); i++) {
		dists.at(i) = (minfv->dists.at(i))/expf + maxfv->dists.at(i);
		dists.at(i) = dists.at(i)/(1+1/expf);
	}
	logf = maxlogf + log(1+1/expf);
	
}

void FactoredVector::Mult(double x) {

   	if(x > 0){
	    logf += log(x);
	}
	else if( x == 0){
		logf = 0;
	 	for (unsigned int i=0; i<dists.size(); i++) {
		    dists[i] *= 0;
	    }
	}	
	else {
		cout << "Multiplying with < 0 !!!" << endl;
		cout << __LINE__ << " logf " << logf << endl;
		/*
	    for (unsigned int i=0; i<dists.size(); i++) {
		    dists[i] *= x;
	    }
	    */
	}
}

void FactoredVector::GetDist(vectr &v, int varid, double &logfactor) {
    unsigned varid_as_uint(varid);
	if (/*varid_as_uint >= 0 &&*/ varid_as_uint < dists.size()) {
		v = dists[varid_as_uint];
		logfactor = logf;
	}
}

	
void FactoredVector::SetDist(const vectr &v, int varid, double logfactor) {
    unsigned varid_as_uint(varid);
    if ( varid_as_uint < dists.size()){
		dists[varid_as_uint] = v;	
		logf = logfactor;
	}
//	cout << __LINE__ << " logf " << logf << endl;
    /*
	v.niceprint(cout);
	double vsum = v.sum();
	bool allzero = false;
	if (vsum == 0){
		allzero = true;
		cout << "sum zero" << endl;
		for(int j=0; j<v.getm(); j++){
			if(v[j] > 0)	allzero = false;
		}
	}
		
		
	if ( varid_as_uint < dists.size()) {//varid_as_uint >= 0 &&
		if (allzero){	
			dists[varid_as_uint] = v+1e-10;	
			vsum = 2*(1e-10);
		}
		else
			dists[varid_as_uint] = v;
			
		logf = logfactor;
		
		if(vsum > 0){
			for(uint i=0; i<dists.size(); i++){
				if(i != varid_as_uint)	dists[i] = dists[i]*vsum; 
			}
		}
	}
	*/
		
}

double FactoredVector::GetDistInstNorm(int varid, int index) {
    double dsum = dists[varid].sum();
    if(dsum > 0)
        return (dists[varid])[index]/dsum;
    else
	    return (dists[varid])[index];
}

double FactoredVector::GetDistInst(int varid, int index) {
	return (dists[varid])[index];
}

void FactoredVector::SetDistInst(int varid, int index, double value) {

	(dists[varid])[index] = value;
}

double FactoredVector::GetMargMin() {
	double minval = GetMin(0);
	double tmp = 0.0;
	for(unsigned int i=1; i<dists.size(); i++) {
		tmp = GetMin(i);
		if (minval > tmp) minval = tmp;
	}
	return minval;
} 

double FactoredVector::GetJointMin() {
	double minval = 1.0;
	for (unsigned int i=0; i<dists.size(); i++) {
		minval *= GetMin(i);
	}
	return minval;
}

double FactoredVector::GetMin(int varid) {
	double minval = fabs((dists[varid])[0]);
	for (int i=0; i<dists[varid].length(); i++) {
		if (minval < fabs((dists[varid])[i]))
			minval = (dists[varid])[i];
	}
	return minval;
}

double FactoredVector::GetMargMax() {
	double maxval = GetMax(0);
	double tmp = 0.0;
	for(unsigned int i=1; i<dists.size(); i++) {
		tmp = GetMax(i);
		if (maxval < tmp) maxval = tmp;
	}
	return maxval;
} 

double FactoredVector::GetJointMax() {
	double maxval = 1.0;
	for(unsigned int i=0; i<dists.size(); i++) {
		maxval *= GetMax(i);
	}
	return maxval;
}

double FactoredVector::GetMax(int varid) {
	double maxval = fabs((dists[varid])[0]);
	for (int i=0; i<dists[varid].length(); i++) {
		if (maxval > fabs((dists[varid])[i]))
			maxval = (dists[varid])[i];
	}
	return maxval;
}

int FactoredVector::GetDistSize(int varid) {
	return (dists[varid]).length();
}

int FactoredVector::Size() { return dists.size(); }


std::ostream & FactoredVector::Print(std::ostream & out, RV * ignored = nullptr03) const {
	for (unsigned int i=0; i<dists.size(); i++){
		
		vectr tmp = dists.at(i);
		tmp.normalize();
		tmp.niceprint(out);	
	}
	out << endl;
	return out;
}


std::ostream & FactoredVector::PrintSimple(std::ostream & out) const {
    double logf = this->GetLogf();
    out << "FactoredVector logf: " << logf << endl;
    this->Print(out);
    return out;
}


//From RVSimple
void FactoredVector::GetDist(vectr &d, double &logfactor) const {
    throw not_yet_implemented_error(the_class_name, "GetDist");
}


void FactoredVector::SetDist(const vectr &d, double logfactor) {
    throw not_yet_implemented_error(the_class_name, "SetDist");
}


double FactoredVector::Prob(int ind, bool log) const{
    throw not_yet_implemented_error(the_class_name, "Prob");
}


void FactoredVector::Load(std::istream &is) {
    throw not_yet_implemented_error(the_class_name, "Load");
}

void FactoredVector::Save(std::ostream &os) const {
    throw not_yet_implemented_error(the_class_name, "Save");
}

void FactoredVector::Reindex(const std::vector<std::vector<int> > &ind) {
    throw not_yet_implemented_error(the_class_name, "Reindex");
}

void FactoredVector::Add(const RVSimple *x, double w) {
    throw not_yet_implemented_error(the_class_name, "Add");
}

SS *FactoredVector::BlankSS() const {
    throw not_yet_implemented_error(the_class_name, "BlankSS");
}

void FactoredVector::AddSS(int x, SS *ss, double w) const {
    throw not_yet_implemented_error(the_class_name, "AddSS");
}

void FactoredVector::AddExpSS(SS *ss, double w) const {
    throw not_yet_implemented_error(the_class_name, "AddExpSS");
}

void FactoredVector::AddSS(const SS *toadd, const RVSimple* rvs,
			const std::vector<std::vector<int> > &mapping, 
			SS *ss, double w) const {
    throw not_yet_implemented_error(the_class_name, "AddSS");
}

int FactoredVector::Sample(Random &rand) const {
    throw not_yet_implemented_error(the_class_name, "Sample");
}

void FactoredVector::Maximize(const SS *ss) {
    throw not_yet_implemented_error(the_class_name, "Maximize");
}

void FactoredVector::Scramble(double alpha, double degree, Random &rand) {
    throw not_yet_implemented_error(the_class_name, "Scramble");
} 

double FactoredVector::LLH(const SS *ss) const {
    throw not_yet_implemented_error(the_class_name, "LLH");
}

double FactoredVector::GetScore(double numTrans, const SS* ss) const {
    throw not_yet_implemented_error(the_class_name, "GetScore");
}

}

