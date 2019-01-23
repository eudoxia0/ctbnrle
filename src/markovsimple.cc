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
#include "markovsimple.h"
#include "params.h"
#include "multisimple.h"
#include "rk.h"
#include "extramath.h"
#include "rvsimple.h"

namespace ctbn {

using namespace std;

SOBJCLASSDEF(MarkovSimple)

MarkovSimple::MarkovSimple(int n) : DynSimple(), Q(n,n,0.0) {
}

MarkovSimple::MarkovSimple(istream &is) : DynSimple(is) {
	is >> Q;
}

MarkovSimple::~MarkovSimple() {
}

void MarkovSimple::LoadOld(istream &is) {
	DynSimple::LoadOld(is);
	is >> Q;
}

void MarkovSimple::SaveOld(ostream &os) const {
	DynSimple::SaveOld(os);
	os << os.fill() << Q;
}

MarkovSimple *MarkovSimple::Clone() const {
	return new MarkovSimple(*this);
}

MarkovSimple *MarkovSimple::MakeNew(int nstates) const {
	return new MarkovSimple(nstates);
}

void MarkovSimple::Normalize() {
	int n = Q.getm();
	for(int i=0;i<n;i++) {
		double sum = 0.0;
		for(int j=0;j<n;j++)
			if (i!=j) sum += Q[i][j];
		Q[i][i] = -sum;
	}
}

void MarkovSimple::Restrict(const vector<int> &ind) {
	int n = ind.size();
	int m = Q.getm();
	int ii = 0;
	for(int i=0;i<m;i++) {
		if (ind[ii]==i) ii++;
		else {
			for(int j=0;j<n;j++) 
				Q[ind[j]][i] = 0.0;
		}
	}
}

void MarkovSimple::Expand(const vector<vector<int> > &ind, int newn) {
    matrix newQ(newn,newn,0.0);
    const int n = Q.getm();
    for (int i=0; i<n; i++) {
        const int ns = ind[i].size();
        for (int j=0; j<n; j++)
            if (Q[i][j]!=0)
                for (int s=0; s<ns; s++)
                    newQ[ind[i][s]][ind[j][s]] = Q[i][j];
        /*
            for (int s=0; s<ns; s++) {
                    newQ[ind[i][s]][ind[i][s]] = Q[i][i];
                    for(int j=0;j<n;j++)
                            newQ[ind[i][s]][ind[j][s]] = Q[i][j];
            }
        */
    }

    Q.swap(newQ);
    //      Q = newQ;
}


CondTransQ<const matrix &> *MarkovSimple::Cond(double t0, double t1) const {
	return new CondTransQ<const matrix &>(Q,t1-t0);
}


SparseCondTransQ<matrix> *MarkovSimple::CondRestrict(double t0, double t1,
		const vector<int> &ind) const {
	return new SparseCondTransQ<matrix>(matrix(Q,ind),t1-t0,ind);
}


CondTransQ2<matrix> *MarkovSimple::Cond(double t) const {
	matrix newQ(Q);
	int n = Q.getm();
	for(int i=0;i<n;i++) newQ[i][i] = 0.0;
	return new CondTransQ2<matrix>(newQ);
}

SparseCondTransQ2<matrix> *MarkovSimple::CondRestrict(double t,
		const vector<int> &fromind, const vector<int> &toind,
		bool transition) const {
	if (transition)
		return new SparseCondTransQ2<matrix>(matrix(Q,fromind,toind),
										fromind,toind);
	else {
		matrix m(fromind.size(),toind.size(),0.0);
		unsigned int i=0, j=0;
		while(i<fromind.size() && j<toind.size()) {
			if (fromind[i]==toind[j]) {
				m[i][j] = 1.0;
				i++;
				j++;
			} else if (fromind[i]<toind[j]) {
				i++;
			} else {
				j++;
			}
		}
		return new SparseCondTransQ2<matrix>(m,fromind,toind);
	}
}

MarkovSimpleSS *MarkovSimple::BlankSS() const {
	return new MarkovSimpleSS(Q.getm());
}

void MarkovSimple::Mult(const DynSimple *x) {
	const MarkovSimple *xx = dynamic_cast<const MarkovSimple *>(x);
	assert(xx!=NULL); // otherwise, not sure what to do...
	if (Q.getm()!=xx->Q.getm()) Q = xx->Q;
	else Q += xx->Q;
}

void MarkovSimple::AddSS(int x, double t0, double deltat, 
					SS *ss, double w) const {
	MarkovSimpleSS *msss = dynamic_cast<MarkovSimpleSS *>(ss);
	msss->c[x][x] += deltat*w;
}

void MarkovSimple::AddTransSS(int x1,int x2, double t, 
				SS *ss,double w) const {
	MarkovSimpleSS *msss = dynamic_cast<MarkovSimpleSS *>(ss);
	if (x1 != x2)
		msss->c[x1][x2] += w;
}

void MarkovSimple::AddExpSS(const RVSimple *alpha, const RVSimple *beta,
		double t0, double deltat, SS *ss, double w) const {
	MarkovSimpleSS *msss = dynamic_cast<MarkovSimpleSS *>(ss);
	const SparseMultiZSimple *a = 
			dynamic_cast<const SparseMultiZSimple *>(alpha);
	assert(a!=NULL);  // otherwise, need to implement more code
	const SparseMultiZSimple *b = 
			dynamic_cast<const SparseMultiZSimple *>(beta);
	assert(b!=NULL); // ditto
	// a->Indexes() should be the same as b->Indexes()!
	const vector<int> &ind = a->Indexes();
	int n = ind.size();
	matrix c;
	double z = suffstat(a->Dist(),b->Dist(),c,matrix(Q,ind),deltat);

	double sum = 0.0;
	for(int i=0;i<n;i++) sum += c[i][i];
	for(int i=0;i<n;i++) for(int j=0;j<n;j++)
		if (i!=j) c[i][j] *= Q[ind[j]][ind[i]]; // c is transposed!
	// normalize c to condition on remaining within the subsystem
	// (and multiply by w)
	c *= w*deltat/sum;

	// add in (with transpose)
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++) 
			msss->c[ind[j]][ind[i]] += c[i][j];
}

void MarkovSimple::AddExpTransSS(const RVSimple *x1, const RVSimple *x2,
		double t, SS *ss, double w) const {
	MarkovSimpleSS *msss = dynamic_cast<MarkovSimpleSS *>(ss);
	const SparseMultiZSimple *a = 
				dynamic_cast<const SparseMultiZSimple *>(x1);
	assert(a!=NULL);  // otherwise, need to implement more code
	const SparseMultiZSimple *b = 
				dynamic_cast<const SparseMultiZSimple *>(x2);
	assert(b!=NULL); // ditto

	const vector<int> &aind = a->Indexes();
	const vector<int> &bind = b->Indexes();
	matrix c(Q,aind,bind);
	c.multbycol(a->Dist());
	c.multbyrow(b->Dist());
	int mm = c.getm();
	int nn = c.getn();
	double csum = 0.0;
	for(int i=0;i<mm;i++)
		for (int j=0;j<nn;j++)
			if (aind[i]==bind[j])
				c[i][j] = 0.0;
			else csum += c[i][j];
	/*
	for(int i=0;i<nn;i++)
		c[i][i] = 0.0;
	// condition on exactly 1 transition happening
	// (and multiply by w)
	double csum = c.sum();
	*/
	if (csum<=0.0) return;
	c *= w/csum;
	int m = aind.size(), n = bind.size();
	for(int i=0;i<m;i++) for(int j=0;j<n;j++)
		msss->c[aind[i]][bind[j]] += c[i][j];
}

void MarkovSimple::AddSS(const SS *toadd, const DynSimple *dyn,
		const std::vector<std::vector<int> > &mapping,
		SS *ss, double w) const {
	const MarkovSimple *m = dynamic_cast<const MarkovSimple *>(dyn);
	assert(m!=NULL); // otherwise, need double-dispatch and more code
	MarkovSimpleSS *mss = dynamic_cast<MarkovSimpleSS *>(ss);
	const MarkovSimpleSS *addmss =
		dynamic_cast<const MarkovSimpleSS *>(toadd);

	int n = Q.getm();
	for(int i=0;i<n;i++) {
		int li = mapping[i].size();
		for(int ki=0;ki<li;ki++) {
			mss->c[i][i] +=
				w*addmss->c[mapping[i][ki]][mapping[i][ki]];

// new version that assumes this comes from a process in which two
// variables do not change at exactly the same time:
			for(int j=0;j<n;j++) {
				if (j==i) continue;
				int lj = mapping[j].size();
				mss->c[i][j] += w*addmss->c[mapping[i][ki]][mapping[j][ki]];
/* old version:
            for(int j=0;j<n;j++) {
               if (j==i) continue;
               int lj = mapping[j].size();
               for(int kj=0;kj<lj;kj++)
                   mss->c[i][j] += w*addmss->c[mapping[i][ki]][mapping[j][kj]];
*/
			}
		}
	}
}


void MarkovSimple::SampleNextEvent(int ind, double t,
		double &newt, int &newind, Random &rand) const {
	newt = t + rand.SampleExp(-Q[ind][ind]);
	double sum = 0.0;
	int n = Q.getn();
	for(int i=0;i<n;i++) if (i!=ind) sum += Q[ind][i];
	newind = rand.SampleMultinomial(&Q[ind][0],n,sum);
}

void MarkovSimple::Maximize(const SS *ss) {
	static const double divzeroq = ParamDouble("DivZeroQ",1e-6);

	const MarkovSimpleSS *msss = dynamic_cast<const MarkovSimpleSS *>(ss);
	int m = Q.getm(), n = Q.getn();
	for(int i=0;i<m;i++) {
		double numtr = 0.0;
		for(int j=0;j<n;j++)
			if (j!=i) numtr += msss->c[i][j];
		double tottime = msss->c[i][i];
		if (tottime<=0) {
			for(int j=0;j<n;j++)
				if (j==i) Q[i][j] = -divzeroq*(n-1);
				else Q[i][j] = divzeroq;
		} else {
			for(int j=0;j<n;j++)
				if (j==i) Q[i][j] = -numtr/tottime;
				else Q[i][j] = msss->c[i][j]/tottime;
		}
	}
	for (int i = 0; i < m; ++i) {
		if (fabs(Q[i][i]) < divzeroq) {
			for (int j = 0; j < n; ++j)
				Q[i][j] = divzeroq;
			Q[i][i] = -divzeroq*(n-1);
		}
		else{
			for (int j = 0; j < n; ++j)
				if(i!=j && fabs(Q[i][j]) < divzeroq)
					{
                        Q[i][j] = divzeroq;
                        Q[i][i] -= divzeroq;
					}
		}
	}
/*
	//DEBUG
	for(int i=0; i<m;++i) {
		if(Q[i][i] > 0) {
			cout << "stopped" << endl;
			cin.get();
		}
	}
*/
}

//by Yu
void MarkovSimple::Scramble(double a, double b, double alpha, 
				double degree, Random &rand) {
	int fix = ParamInt("FixScramble", -1);
	if(fix!=-1) { 
		int n = Q.getn();
		int m = Q.getm();
		for (int i = 0; i < m; ++i)
			for (int j = 0; j < n; ++j)
				if(i==j) Q[i][j] = -a;
				else Q[i][j] = a/(n-1);
		return;
	}

	int n = Q.getn();
	int m = Q.getm();
	double *in = new double[n-1];
	matrix original(Q);

	for(int i=0;i<n-1;i++)
		in[i] = alpha;
	for(int i=0;i<m;i++) {
		double q = rand.SampleGamma(a,b);
		rand.SampleDirichlet(in,n-1,&(Q[i][0]));
		Q[i][m-1] = Q[i][i];
		Q[i][i] = -q;
		for(int j=0;j<n;j++)
			if (j!=i) Q[i][j] *= q;
	}
	delete []in;
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			Q[i][j] = Q[i][j] 
				  * degree 
				  + original[i][j] 
				  * (1.0 - degree);
}

double MarkovSimple::LLH(const SS *ss) const {
	const MarkovSimpleSS *mss = dynamic_cast<const MarkovSimpleSS *>(ss);
	double llh = 0.0;
	matrix c = mss->c;
	unsigned int n = c.getn();
	unsigned int m = c.getm();
	for (unsigned int i=0; i<m; i++)
		for (unsigned int j=0; j<n; j++)
			if (c[i][j]>0.0) {
				if (i!=j) {
                    llh += c[i][j] * log(Q[i][j]);
				} else {
                    llh += c[i][j] * Q[i][j];
                }
            }
	return llh;
}

// Score over a single matrix. (sum over states)
double MarkovSimple::GetScore(int parentCard, double numTrans, double amtTime,
				const SS* ss) const {

	double localScore(0);
	const MarkovSimpleSS *mss = dynamic_cast<const MarkovSimpleSS *>(ss);
	const matrix& c = mss->c;

	int size = c.getm();

	double alpha_xx_u(numTrans / ((size-1)*size*parentCard));
	double alpha_x_u((size-1)*alpha_xx_u);
	double tau_x_u(amtTime / (size*parentCard));
/*
    cout << "size=" << size << endl;
    cout << "amtTime=" << amtTime << endl;
    cout << "numTrans=" << numTrans << endl;
    cout << "parentCard=" << parentCard << endl;
    cout << "alpha_xx_u=" << alpha_xx_u << endl;
    cout << "alpha_x_u=" << alpha_x_u << endl;
    cout << "tau_x_u=" << tau_x_u << endl;
    cout << "c=" << c << endl;
*/

	for(int i = 0; i < size; i++) {
		// Sum over each state of the variable
		double m_x_u(0);
		double t_x_u = c[i][i];

		// Theta score (transitions)
		for(int j = 0; j < size; j++) {
			if(i == j) continue;
			m_x_u += c[i][j];
			localScore += lngamma(alpha_xx_u + c[i][j]);
		}
		localScore += lngamma(alpha_x_u) - lngamma(alpha_x_u + m_x_u)
				- (size-1)*lngamma(alpha_xx_u);

		// q score (time)
		localScore += lngamma(alpha_x_u + m_x_u + 1)
				+ (alpha_x_u + 1)*log(tau_x_u)
				- lngamma(alpha_x_u + 1)
				- (alpha_x_u + m_x_u + 1)
				* log(tau_x_u + t_x_u);
	}
    //cout << "\nlocalScore=" << localScore << endl;
	return localScore;
}
		
//-----

SOBJCLASSDEF(MarkovSimpleSS)

MarkovSimpleSS::MarkovSimpleSS(int n) : c(n,n,0.0) {
}

MarkovSimpleSS::MarkovSimpleSS(istream &is) {
	LoadOld(is);
}

MarkovSimpleSS::~MarkovSimpleSS() {
}

MarkovSimpleSS *MarkovSimpleSS::Clone() const {
	return new MarkovSimpleSS(*this);
}

void MarkovSimpleSS::LoadOld(istream &is) {
	is >> c;
}

void MarkovSimpleSS::SaveOld(ostream &os) const {
	os << c;
}

void MarkovSimpleSS::Scale(double w) { 
	for(int i=0; i<c.getm(); i++)
		for(int j=0; j<c.getn(); j++)
			c[i][j] /= w;

}

void MarkovSimpleSS::AddSS(const SS* nss, double w) { 
	const MarkovSimpleSS *mss = dynamic_cast<const MarkovSimpleSS*>(nss);
	for(int i=0; i<c.getm(); i++)
		for(int j=0; j<c.getn(); j++)
			c[i][j] += mss->c[i][j] * w;

}

MarkovSimpleToggle::MarkovSimpleToggle(int n) : MarkovSimple(n) {
}

MarkovSimpleToggle::MarkovSimpleToggle(istream &is) : MarkovSimple(is) {
}

MarkovSimpleToggle::~MarkovSimpleToggle() {
}

MarkovSimpleToggle *MarkovSimpleToggle::Clone() const {
	return new MarkovSimpleToggle(*this);
}

MarkovSimpleToggle *MarkovSimpleToggle::MakeNew(int nstates) const {
	return new MarkovSimpleToggle(nstates);
}

void MarkovSimpleToggle::Maximize(const SS *ss) {
	static const double divzeroq = ParamDouble("DivZeroQ",1e-6);

	const MarkovSimpleSS *msss = dynamic_cast<const MarkovSimpleSS *>(ss);
	matrix &Q = Intensity();

	int m = Q.getm(), n = Q.getn();
	double tottr = 0.0, tottime = 0.0;
	for (int i = 0; i < m; ++i) 
		for (int j = 0; j < n; ++j)
			if (j == i)
				tottime += msss->c[i][i];
			else
				tottr += msss->c[i][j];
	double toggleq = divzeroq;
	if (tottime > 0) toggleq = tottr/tottime;
	for (int i = 0; i < m; ++i) 
		for (int j = 0; j < n; ++j) {
			if (j == i) Q[i][i] = -(n-1)*toggleq;
			else Q[i][j] = toggleq;
		}
	for (int i = 0; i < m; ++i) {
		if (fabs(Q[i][i]) < divzeroq) {
			for (int j = 0; j < n; ++j)
				Q[i][j] = divzeroq;
			Q[i][i] = -divzeroq*(n-1);
		}
		else{
			for (int j = 0; j < n; ++j)
				if(i!=j && fabs(Q[i][j]) < divzeroq) {
					Q[i][j] = divzeroq;
					Q[i][i] -= divzeroq;
				}
		}

	}
}

void MarkovSimpleToggle::Scramble(double a, double b, double alpha, 
				double degree, Random &rand) {
	int fix = ParamInt("FixScramble", -1);
	if(fix!=-1) { 
		int n = Q.getn();
		int m = Q.getm();
		for (int i = 0; i < m; ++i)
			for (int j = 0; j < n; ++j)
				if(i==j) Q[i][j] = -a;
				else Q[i][j] = a/(n-1);
		return;
	}

	int n = Q.getn();
	int m = Q.getm();
	matrix original(Q);
	double q = rand.SampleGamma(a,b);

	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			Q[i][j] = i==j ? -(n-1)*q : q;
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			Q[i][j] = Q[i][j] 
				  * degree 
				  + original[i][j] 
				  * (1.0 - degree);
}

} // end of ctbn namespace
