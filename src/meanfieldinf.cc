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

#include "meanfieldinf.h"
#include "markovdyn.h"
#include "rk.h"
#include "params.h"
#include "nullptr03.h"
#include <iostream>
#include <deque>
#include <algorithm>

namespace ctbn {

using namespace std;

// The precision for the energy convergence accuracy
const double PRECISION = 1e-5;


MeanFieldInf::MeanFieldInf() : Inference(), factored(nullptr03),
				initialized(false), optimized(false), 
				eps(1e-5) {
}

MeanFieldInf::~MeanFieldInf() {
	mus.clear();
	gammas.clear();
	delete factored;
}

MeanFieldInf* MeanFieldInf::Clone() const {
	MeanFieldInf* cloneInf = new MeanFieldInf(*this);
	if(this->factored) cloneInf->factored = this->factored->Clone();
	else cloneInf->factored = nullptr03;
	return cloneInf;
}


void MeanFieldInf::SetProcess(const Process *p) {
	Reset();
	this->p = nullptr03;
	this->p = dynamic_cast<const Markov *>(p);
}

void MeanFieldInf::SetTrajectory(const Trajectory *tr) {
	Reset();
	this->tr = nullptr03;
	this->tr = tr;
}

double MeanFieldInf::Filter(const Instantiation &x, double t, 
				bool log) {
	return 0.0;
}

double MeanFieldInf::Smooth(const Instantiation &x, double t,
				bool log) {
	if(!optimized) Optimize();
	return GetProb(x,t,log);
}

double MeanFieldInf::Prob(double t, bool log) {
	cout << "not implemented" << endl;
	return 0.0;
}

void MeanFieldInf::AddExpSuffStats(SS *ss, double w) {
	if(!optimized) Optimize();
	MarkovSS *mss = dynamic_cast<MarkovSS *>(ss);
	CTBNDynSS *css = dynamic_cast<CTBNDynSS *>(mss->dss);
	const BN* bn = dynamic_cast<const BN *>(p->GetStartDist());
	Instantiation b(bn->Domain()+bn->CondDomain());
	b.SetAllVal(0);
	do {
		double startprob = Smooth(b, tr->TimeBegin());
		if(startprob>0) bn->AddSS(b, mss->p0ss, startprob*w);
	} while(b.Inc());
	
	const CTBNDyn* cdyn = 
		dynamic_cast<const CTBNDyn *>(p->GetDynamics());
	Context c = cdyn->Domain();
	vector<int> varL = c.VarList();

	vector<double> cutPoints;
	for(unsigned int i=0; i<varL.size(); i++) {
		vector<double> varCP = 
			mus.find(varL[i])->second.GetCutPoints();
		for(unsigned int j=0; j<varCP.size()-1; j++)
			cutPoints.push_back(varCP[j]);
	}
	sort(cutPoints.begin(), cutPoints.end());
	vector<double>::iterator itv;
	itv = unique(cutPoints.begin(),cutPoints.end());
	cutPoints.resize(itv - cutPoints.begin());

	for(unsigned int i=0; i<varL.size(); i++) {
		const MarkovDyn* node = dynamic_cast<const MarkovDyn* >(
						cdyn->NodeByVar(varL[i]));
		const Context& v = node->Domain();
		const Context& cv = node->CondDomain();
		Instantiation idx(v+cv);
		idx.SetAllVal(0);
		do {
			double timeSS = 0;
			for(unsigned int j=0; j<cutPoints.size()-1; j++) {
				timeSS +=
					mfSuffStat(this,idx,
						cutPoints[j],
						cutPoints[j+1],eps);
			}
			node->AddSS(idx, 0, w, css->dynss[i],timeSS);
		} while(idx.Inc());
		Instantiation v1(v);
		Instantiation v2(v);
		Instantiation parents(cv);
		v1.SetAllVal(0);
		v2.SetAllVal(0);
		parents.SetAllVal(0);
		Instantiation idx1(v+cv);
		Instantiation idx2(v+cv);

		Instantiation curri = tr->Values(v, tr->TimeBegin());
		Instantiation nexti;
		
		Trajectory::Index index = tr->Begin(v);
		while(!index.Done()) {
			double t = index.Time();
			int sign = index.TestInc(c);
			if(!index.Done()) {
				nexti = index.Values();
				t = index.Time();
				if(sign==2) {
					do {
						v1.SetVal(curri);
						v2.SetVal(nexti);
						idx1.SetVal(v1);
						idx1.SetVal(parents);
						idx2.SetVal(v2);
						idx2.SetVal(parents);
						double transSS = 
							Smooth(parents, t);
						node->AddTransSS(idx1,
							idx2,t,
							css->dynss[i],
							transSS*w);
					} while(parents.Inc());
				}
				curri = nexti;
			}
		}

		v1.SetAllVal(0);
		v2.SetAllVal(0);
		parents.SetAllVal(0);

		do {
			idx1.SetVal(parents);
			idx2.SetVal(parents);
			do {
				idx1.SetVal(v1);
				do {
					idx2.SetVal(v2);
					double transSS = 0;
					for(unsigned int j=0; j<cutPoints.size()-1; j++) {
						transSS +=
						mfTransSuffStat(this,varL[i],
								idx1,idx2,
								cutPoints[j],
								cutPoints
								[j+1],eps);
					}
				if(transSS < 0) {
					PrintMu(); cout << endl;
					PrintGamma(); cout << endl;
					cin.get();
				}
					node->AddTransSS(idx1,idx2,
						0, css->dynss[i],transSS*w);
				} while(v2.Inc());
			} while(v1.Inc());
		} while(parents.Inc());
	}
}


void MeanFieldInf::AddExpSuffStats(const Dynamics *dyn, SS *ss,
					double w) {
	if(!optimized) Optimize();
	
	const Context &v = dyn->Domain();
	int vid = v.MinVar();
	const Context &cv = dyn->CondDomain();
	Context allv(v+cv);

	vector<int> varL = allv.VarList();

	vector<double> cutPoints;
	for(unsigned int i=0; i<varL.size(); i++) {
		vector<double> varCP = 
			mus.find(varL[i])->second.GetCutPoints();
		for(unsigned int j=0; j<varCP.size()-1; j++)
			cutPoints.push_back(varCP[j]);
	}
	sort(cutPoints.begin(), cutPoints.end());
	vector<double>::iterator itv;
	itv = unique(cutPoints.begin(),cutPoints.end());
	cutPoints.resize(itv - cutPoints.begin());

	Instantiation idx(v+cv);
	idx.SetAllVal(0);

	do {
		double timeSS = 0;
		for(unsigned int j=0; j<cutPoints.size()-1; j++) {
			timeSS +=
				mfSuffStat(this,idx,
					cutPoints[j],
					cutPoints[j+1],eps);
		}
		dyn->AddSS(idx, 0, timeSS, ss, w);
	} while(idx.Inc());
	Instantiation v1(v);
	Instantiation v2(v);
	Instantiation parents(cv);
	v1.SetAllVal(0);
	v2.SetAllVal(0);
	parents.SetAllVal(0);
	Instantiation idx1(v+cv);
	Instantiation idx2(v+cv);

	Instantiation curri = tr->Values(v, tr->TimeBegin());
	Instantiation nexti;
	
	Trajectory::Index index = tr->Begin(v);
	while(!index.Done()) {
		double t = index.Time();
		int sign = index.TestInc(allv);
		if(!index.Done()) {
			nexti = index.Values();
			t = index.Time();
			if(sign==2) {
				do {
					v1.SetVal(curri);
					v2.SetVal(nexti);
					idx1.SetVal(v1);
					idx1.SetVal(parents);
					idx2.SetVal(v2);
					idx2.SetVal(parents);
					double transSS = 
						Smooth(parents, t);
					dyn->AddTransSS(idx1,
						idx2,t,
						ss,
						transSS*w);
				} while(parents.Inc());
			}
			curri = nexti;
		}
	}

	v1.SetAllVal(0);
	v2.SetAllVal(0);
	parents.SetAllVal(0);

	do {
		idx1.SetVal(parents);
		idx2.SetVal(parents);
		do {
			idx1.SetVal(v1);
			do {
				idx2.SetVal(v2);
				double transSS = 0;
				for(unsigned int j=0; j<cutPoints.size()-1; j++) {
					transSS +=
					mfTransSuffStat(this,vid,
							idx1,idx2,
							cutPoints[j],
							cutPoints
							[j+1],eps);
				}
				dyn->AddTransSS(idx1,idx2,
					0, ss,transSS*w);
			} while(v2.Inc());
		} while(v1.Inc());
	} while(parents.Inc());
}

void MeanFieldInf::AddExpSuffStats(const RV *p0, SS *ss,
					double w) {
	Instantiation i(p0->Domain()+p0->CondDomain());
	i.SetAllVal(0);
	do {
		double startprob = Smooth(i, tr->TimeBegin());
		if(startprob>0) p0->AddSS(i, ss, startprob*w);
	} while(i.Inc());
}

void MeanFieldInf::TransToFactored() {
	if(!p || factored) return;
	const CTBNDyn* cdyn = dynamic_cast <const CTBNDyn *>(p->GetDynamics());
	Context c = cdyn->Domain();
	vector<int> varL = c.VarList();

	//Create a factored version of the dynamics
	CTBNDyn* factoredCDyn = new CTBNDyn();
	Context nullc;
	for(unsigned int i=0; i<varL.size(); i++) {
		Context cur;
		cur.AddVar(varL[i], c.Cardinality(varL[i]));
		factoredCDyn->AddNode(new MarkovDyn(cur, nullc));
	}
	factoredCDyn->FillParams(*cdyn,true);
	factored = new Markov(p->GetStartDist()->Clone(),factoredCDyn);
}

void MeanFieldInf::Initialize() {
	if(!p || !tr || initialized) return;
	TransToFactored();
	const Context &c = factored->GetDynamics()->Domain();
	vector<int> varL = c.VarList();

	for(unsigned int ii=0; ii<varL.size(); ii++) {
		int curVarId = varL[ii];
		const Dynamics *curNode = 
			dynamic_cast<const CTBNDyn *>(
				factored->GetDynamics())->
				NodeByVar(curVarId);
		const Context &nodeCon = curNode->Domain();
		const MarkovDyn *curmdyn =
			dynamic_cast<const MarkovDyn *>(curNode);

		const matrix &Q = (*curmdyn)(0)->Intensity();
		int states = nodeCon.Size();

		Trajectory trNode = tr->ExtractNodeTraj(curVarId);
		double ts = trNode.TimeBegin();
		double te = trNode.TimeEnd();
		Instantiation curri =
			trNode.Values(nodeCon, trNode.TimeBegin());
		Trajectory::Index i = trNode.Begin(nodeCon);

		vector<double> times;
		vector<vectr> vals;
		vector<int> transtype;

		times.push_back(trNode.TimeBegin());
		vectr e0(states,0.0);
		int val = i.Values().Value(curVarId);
		bool startHidden = false;
		if(val == -1) {
			startHidden = true;
			e0 = vectr(states,1.0/states);
		}
		else {
			e0[val] = 1;
		}
		transtype.push_back(0);
		vals.push_back(e0);

		while(!i.Done()) {
			int trans = i.TestInc(nodeCon);
			times.push_back(i.Time());
			int state = i.Values().Value(curVarId);

			vectr evid(states,0.0);
			if(state >= 0) evid[state] = 1.0;
			if(i.Done() && state == -1) {
				evid = vectr(states,1.0);
			}
			transtype.push_back(trans);
			vals.push_back(evid);
		}
		AddVar(curVarId);
		transTimes.insert(make_pair<int,vector<double> >(
			curVarId,times));

		bool inHidden = false;
		for(unsigned int i=0; i<vals.size(); i++) {
			if(i == vals.size()-1) {
				if(transtype[i] == 1) {
					AddMuVal(curVarId,times[i],
						vals[i-1]);
					matrix gammaVal(states,states,0.0);
					AddGammaVal(curVarId,times[i],
							gammaVal);
				}
				continue;
			}
			else if(!(startHidden && i==0) &&
				(transtype[i]!=1 || inHidden)) {
				AddMuVal(curVarId,times[i],vals[i]);
				AddMuVal(curVarId,times[i+1]-1e-11,vals[i]);
				matrix gammaVal(states,states,0.0);
				AddGammaVal(curVarId,times[i],gammaVal);
				AddGammaVal(curVarId,times[i+1]-1e-11,gammaVal);
				inHidden = false;
				continue;
			}
			inHidden = true;
			vectr alphastart;
			if(startHidden && i==0)
				alphastart = vals[i];
			else
				alphastart = vals[i-1];

			vectr betastart = vals[i+1];
			double starttime = times[i];
			double endtime = times[i+1];
			AddHiddenInterval(curVarId,starttime,endtime);
			int steps = 20;
			double deltat = (endtime-starttime)/steps;
			vector<vectr> alpha;
			deque<vectr> beta;
			vectr curAlpha = alphastart;
			vectr curBeta = betastart;
			alpha.push_back(curAlpha);
			beta.push_front(curBeta);

			for(int j=0; j<steps; j++) {
				vexpmt(curAlpha,Q,deltat);
				alpha.push_back(curAlpha);

				expmtv(curBeta,Q,deltat);
				beta.push_front(curBeta);
			}

			double t = starttime;
			AddMuVal(curVarId,t,alphastart);

			for(unsigned int j=0; j<alpha.size(); j++, t+=deltat) {
				vectr dist = alpha[j].dotstar(beta[j]);
				double probEvid = dist.normalize();

				matrix outer(alpha[j],beta[j]);
				matrix gammaVal = outer.dotstar(Q);
				gammaVal = gammaVal / probEvid;
				for(int k=0; k<gammaVal.getn(); k++) {
					double sumrow = 0;
					for(int l=0;l<gammaVal.getm();l++) {
						if(k==l) continue;
						sumrow += gammaVal[k][l];
					}
					gammaVal[k][k] = -sumrow;
				}
				if(j==alpha.size()-1 && abs(te-t) > 1e-11) 
					t -= 1e-11;
				AddMuVal(curVarId,t,dist);
				AddGammaVal(curVarId,t,gammaVal);
			}
		}
	}
	initialized = true;
}

void MeanFieldInf::Optimize() {
	// A bound on the maximum number of iterations for the optimization of the 
	// ODEs
	const int KMAX = ParamInt("KMAX",30);

	if(!initialized) {
		PrintMu();
		Initialize();
	}
	if(optimized) return;

	const Context &c = factored->GetDynamics()->Domain();
	vector<int> varL = c.VarList();
	bool converged = false;
	int runNumber = 1;
	double prevEnergy = -INT_MAX;
	double curEnergy = -INT_MAX;
	vectr deltaEnergies(3,INT_MAX);

	while(!converged) {
		prevEnergy = curEnergy;
		for(unsigned int i=0; i<varL.size(); i++) {
			int curVarId = varL[i];
			ContFunction<vectr> &curMu = 
				mus.find(curVarId)->second;
			ContFunction<matrix> &curGamma = 
				gammas.find(curVarId)->second;

			const vector<double_pair> &hiddenPairs =
				hiddenIntervals.find(curVarId)->second;
			int states = GetNumStates(curVarId);

				
			for(unsigned int j=0; j<hiddenPairs.size(); j++) {
				double t0 = hiddenPairs[j].first;
				double t1 = hiddenPairs[j].second;

				//Get children transitions
				map<double,int> childTrans = 
					GetChildTransTimes(curVarId,t0,t1);
				map<double,int>::reverse_iterator rit = 
					childTrans.rbegin();

				rhoTemp.Clear();
				if(j==hiddenPairs.size()-1 && 
					tr->Value(varL[i],tr->TimeEnd())) {
					rhoTemp.AddVal(t1,vectr(states,1));
				}
				else {
					rhoTemp.AddVal(t1,curMu.GetVal(t1));
				}


				// Solve rho backward
				// Integrate from t1 to t0
				double t1temp = t1;

				// Take care of observed child transitions
				while(rit != childTrans.rend() && 
					rit->first > t0) {
					double tPlus = rit->first+1e-11;
					int chVarId = rit->second;
					int chstates = GetNumStates(chVarId);
					matrix ceQ(chstates,chstates,0.0);
					mfBackward(this,curVarId,
						t1,tPlus,eps);

					pair<int,int> trans = 
						GetTrans(chVarId,
							rit->first);
					int y0 = trans.first;
					int y1 = trans.second;

					vectr rhotPlus = 
						rhoTemp.GetVal(tPlus);
					for(int k=0; k<states; k++) {
						CalcExpectQ(chVarId,
							rit->first,
							ceQ,curVarId,k,
							false);
						rhotPlus[k] *= ceQ[y0][y1];
					}
					t1=rit->first;
					rhoTemp.AddVal(t1,rhotPlus);
					++rit;
				}

				// Finish up the backward integration
				mfBackward(this,curVarId,t1,t0,eps);

				t1=hiddenPairs[j].second;
				muTemp.Clear();

				// Run mean field for BNs on t_0 if unobserved
				if(t0 == tr->TimeBegin() &&
					tr->Value(curVarId,
						tr->TimeBegin()) == -1) {
					vectr mut0(states,0.0);
					vectr rhot0 = rhoTemp.GetVal(t0);
					CalcMuStart(curVarId,rhot0,mut0);
					muTemp.AddVal(t0,mut0);
					
				}
				else {
					muTemp.AddVal(t0,curMu.GetVal(t0));
				}

				gammaTemp.Clear();
				matrix gammaVal(states,states,0.0);
				CalcGamma(curVarId,t0,
					muTemp.GetVal(t0),
					rhoTemp.GetVal(t0),gammaVal);
				gammaTemp.AddVal(t0,gammaVal);

				// Solve mu forward
				// Integrate from t0 to t1
				mfForward(this,curVarId,t0,t1,eps);
			
				if(j==hiddenPairs.size()-1 && 
					tr->Value(curVarId,
						tr->TimeEnd()) == -1) {
					muTemp.AddVal(t1,
						muTemp.GetVal(t1-1e-11));
					muTemp.EraseVal(t1-1e-11);
					gammaTemp.AddVal(t1,
						gammaTemp.GetVal(t1-1e-11));
					gammaTemp.EraseVal(t1-1e-11);
				}
				else {
					muTemp.AddVal(t1,curMu.GetVal(t1));
					gammaVal *= 0;
					gammaTemp.AddVal(t1,gammaVal);
				}

				curMu.Replace(muTemp);
				curGamma.Replace(gammaTemp);
			}
		}
		curEnergy = CalcEnergy();
		deltaEnergies[runNumber%3] = curEnergy - prevEnergy;
		converged = runNumber >= KMAX
				|| abs(deltaEnergies.sum() /
				(3 * prevEnergy)) < 0.01 
				|| curEnergy == 0
				|| abs(curEnergy - prevEnergy) < PRECISION;
		runNumber++;
	}
	optimized = true;
}

void MeanFieldInf::CalcRhoGrad(int varid, double t, const vectr &rhoVal,
				vectr &rhoGrad) {
	int states = GetNumStates(varid);
	assert(rhoGrad.getm() == states);

	matrix qBack(states,states,0.0);
	CalcExpectQ(varid,t,qBack);

	vectr psi(states,0.0);
	CalcPsi(varid,t,psi);
	
	for(int xi=0; xi<states; xi++) {
		rhoGrad[xi] = -rhoVal[xi]*(qBack[xi][xi] + psi[xi]);
		for(int yi=0; yi<states; yi++) {
			if(xi==yi) continue;
			rhoGrad[xi] -= rhoVal[yi]*qBack[xi][yi];
		}
	}
}

void MeanFieldInf::CalcMuGrad(int varid, double t, const matrix &gammaVal,
				vectr &muGrad) {
	int states = GetNumStates(varid);
	assert(muGrad.getm() == states);
	muGrad *= 0;

	for(int xi=0; xi<states; xi++) {
//		muGrad[xi] = 0;
		for(int yi=0; yi<states; yi++) {
			muGrad[xi] += gammaVal[yi][xi];
		}
	}
}

void MeanFieldInf::CalcGamma(int varid, double t, 
				const vectr &muVal, const vectr &rhoVal,
				matrix &gamma) {
	int states = GetNumStates(varid);
	assert(gamma.getm() == states);

	matrix eQ(states,states,0.0);
	CalcExpectQ(varid,t,eQ);

	for(int xi=0; xi<states; xi++) {
		double sum = 0;
		for(int yi=0; yi<states; yi++) {
			if(xi==yi) continue;
			double numerator = eQ[xi][yi]*rhoVal[yi]*muVal[xi];
			double res = 
				rhoVal[xi] < 1e-5 ? 0 : numerator/rhoVal[xi];
			gamma[xi][yi] = res;
			sum -= res;
		}
		gamma[xi][xi] = sum;
	}
}

void MeanFieldInf::CalcExpectQ(int varid, double t, matrix &expectQ, 
	int jvarid, int xj, bool log) {
	const CTBNDyn* cdyn = 
		dynamic_cast<const CTBNDyn* >(p->GetDynamics());	
	const MarkovDyn* node = dynamic_cast<const MarkovDyn* >(
					cdyn->NodeByVar(varid));
	const Context& v = node->Domain();
	const Context& cv = node->CondDomain();

	int states = v.Size();
//	expectQ = matrix(states,states,0.0);
//	expectQ *= 0;
	
	Context jv;
	bool cond = false;
	if(jvarid != -1 && xj != -1) {
		jv = cdyn->NodeByVar(jvarid)->Domain();
		cond = true;
	}

	Instantiation i(cv-jv);
	vector<int> paList = i.VarList();
	i.SetAllVal(0);

	Instantiation jvi(jv);
	if(cond) jvi.SetVal(jvarid, xj);

	Instantiation idx(cv);

	for(int xi=0;xi<states;xi++) {
		for(int yi=0;yi<states;yi++) {
			double val = 0;
			do {
				double muprod = 1.0;
				for(unsigned int k=0; k<paList.size(); k++) {
					double thismu = 
						mus.find(paList[k])->
						second.GetVal(t)
						[i.Value(paList[k])];
					muprod *= thismu;
				}
				if (muprod==0.0) continue;
				idx.SetVal(i);
				idx.SetVal(jvi);
				double rate = 
					(*node)(idx)->Intensity()[xi][yi];
				if(xi!=yi && log) rate = ::log(rate);
				val += muprod*rate;
			} while(i.Inc());
			expectQ[xi][yi] = (xi!=yi && log ? exp(val) : val);
		}
	}
}

// Needs mu, gamma, eQ updated
void MeanFieldInf::CalcPsi(int varid, double t, vectr &psi) {
	int states = GetNumStates(varid);
	psi *= 0;

	const CTBNDyn* cdyn = 
		dynamic_cast<const CTBNDyn* >(p->GetDynamics());	
	vector<int> children = cdyn->GetChildrenByVar(varid);

	for(unsigned int j=0; j<children.size(); j++) {
		int chvarid = children[j];
		int chstates = cdyn->NodeByVar(chvarid)->Domain().Size();

		vectr muVal = mus.find(chvarid)->second.GetVal(t);
		matrix gammaVal = gammas.find(chvarid)->second.GetVal(t);

		for(int xi=0; xi<states; xi++) {
			matrix eQ(chstates,chstates,0.0);

			// Conditioned on varid=xi
			CalcExpectQ(chvarid,t,eQ,varid,xi);
			for(int xj=0;xj<chstates;xj++) {
				psi[xi] += muVal[xj]*eQ[xj][xj];
				for(int yj=0;yj<chstates;yj++) {
					if(xj==yj || gammaVal[xj][yj]==0) continue;
					psi[xi] += gammaVal[xj][yj]
							* log(eQ[xj][yj]);
				}
			}
		}
	}
}

double MeanFieldInf::CalcPointCompEnergy(int varid, double t) {
	int states = GetNumStates(varid);

	matrix eQ(states,states,0.0);
	CalcExpectQ(varid,t,eQ);

	vectr muVal = mus.find(varid)->second.GetVal(t);
	matrix gammaVal = gammas.find(varid)->second.GetVal(t);

	double res = 0;

	for(int x=0; x<states; x++) {
		if(muVal[x] <= 0) continue;
		res += ((muVal[x]*eQ[x][x]) - 
			gammaVal[x][x]*(1+log(muVal[x])));
		for(int y=0; y<states; y++) {
			if(x==y) continue;
			res += gammaVal[x][y]*
				(log(eQ[x][y]) -
				log(gammaVal[x][y]));
		}
	}
	double ret = isnan(res) ? 0 : res;
	assert(!isnan(ret) && !isinfinite(ret));
	return ret;
}

double MeanFieldInf::PointTransSuffStat(int varid, const Instantiation &x1,
					const Instantiation &x2,
					double t) {
	const CTBNDyn* cdyn = dynamic_cast<const CTBNDyn *>(
		p->GetDynamics());
	const MarkovDyn* node = dynamic_cast<const MarkovDyn *>(
		cdyn->NodeByVar(varid));
	const Context& v = node->Domain();
	const Context& cv = node->CondDomain();
	int i1 = v.Index(x1);
	int i2 = v.Index(x2);
	Instantiation parents(x1-v);
	parents.SetVal(x1);

	double ret = 1;
	if(i1!=i2) {
		ret *= gammas.find(varid)->second.GetVal(t)[i1][i2];
		ret *= Smooth(parents,t);
	}
	else 
		ret = 0;
	return ret;
}

double MeanFieldInf::CalcCompEnergy(int varid) {
	ContFunction<vectr> &curMu = 
		mus.find(varid)->second;
	map<int,vector<double> >::iterator it = transTimes.find(varid);
	double res = 0;
	double finalres = 0;
	vector<double> cutPoints = curMu.GetCutPoints();
	for(unsigned int i=0; i<cutPoints.size()-1; i++) {
		res = mfEnergy(this,varid,cutPoints[i],cutPoints[i+1],eps);
		finalres += res;
	}
	return finalres;
}

double MeanFieldInf::CalcEnergy() {
	double res = 0;
	const Context &c = factored->GetDynamics()->Domain();
	vector<int> varL = c.VarList();
	for(unsigned int i=0; i<varL.size(); i++) 
		res += CalcCompEnergy(varL[i]);
	return res;
}

void MeanFieldInf::GetRhoTempVal(double t, vectr &rhoVal) {
	rhoVal = rhoTemp.GetVal(t);
}

void MeanFieldInf::GetMuTempVal(double t, vectr &muVal) {
	muVal = muTemp.GetVal(t);
}

void MeanFieldInf::GetGammaTempVal(double t, matrix &gammaVal) {
	gammaVal = gammaTemp.GetVal(t);
}

void MeanFieldInf::AddRhoTempVal(double t, const vectr &rhoVal) {
	rhoTemp.AddVal(t,rhoVal);
}

void MeanFieldInf::AddMuTempVal(double t, const vectr &muVal) {
	muTemp.AddVal(t,muVal);
}

void MeanFieldInf::AddGammaTempVal(double t, const matrix &gammaVal) {
	gammaTemp.AddVal(t,gammaVal);
}

void MeanFieldInf::PrintMu() const {
	map<int,ContFunction<vectr> >::const_iterator it = mus.begin();
	for(; it!=mus.end(); ++it) {
		cout << "var: " << it->first << endl;
		it->second.Print();
		cout << endl << "====================" << endl << endl;
	}
}

void MeanFieldInf::PrintGamma() const {
	map<int,ContFunction<matrix> >::const_iterator it = 
		gammas.begin();
	for(; it!=gammas.end(); ++it) {
		cout << "var: " << it->first << endl;
		it->second.Print();
		cout << endl << "====================" << endl << endl;
	}
}

void MeanFieldInf::PrintHidden() const {
	map<int,vector<double_pair> >::const_iterator it = 
		hiddenIntervals.begin();
	for(; it!=hiddenIntervals.end(); ++it) {
		cout << "var: " << it->first << endl;
		cout << it->second;
		cout << endl << "====================" << endl << endl;
	}
}


double MeanFieldInf::CalcQuery(QueryCalculator &calc) {
	SS *ss = p->BlankSS();
	AddExpSuffStats(ss);
	double ret = calc.Calculate(ss);
	delete ss;
	return ret;
}

pair<int,int> MeanFieldInf::GetTrans(int varid, double t) const {
	int from = tr->Value(varid,t-1e-11);
	int to = tr->Value(varid,t+1e-11);
	return make_pair<int,int>(from,to);
}

map<double,int> MeanFieldInf::GetChildTransTimes(int varid,
					double start, double end) const {
	const CTBNDyn* cdyn = dynamic_cast<const CTBNDyn *>(
		p->GetDynamics());
	vector<int> children = cdyn->GetChildrenByVar(varid);

	map<double,int> ret;
	for(unsigned int i=0; i<children.size(); i++) {
		vector<double> transTimes = 
			GetObservedTransTimes(children[i],start,end);
		for(unsigned int j=0; j<transTimes.size(); j++) {
			ret.insert(make_pair<double,int>(
				transTimes[j],children[i]));
		}
	}
	return ret;
}

vector<double> MeanFieldInf::GetObservedTransTimes(int varid, 
					double start, double end) const {
	const Context &nodeCon  = 
		dynamic_cast<const CTBNDyn *>(
				p->GetDynamics())->
				NodeByVar(varid)->Domain();
	Trajectory trNode = tr->ExtractNodeTraj(varid);
	Instantiation curri =
		trNode.Values(nodeCon, trNode.TimeBegin());
	Trajectory::Index i = trNode.Begin(nodeCon);

	vector<double> times;

	while(!i.Done()) {
		if(i.TestInc(nodeCon)!=2) continue;
		double t = i.Time();
		if(t>start && t<end) times.push_back(t);
	}
	return times;
}

void MeanFieldInf::CalcMuStart(int varid, const vectr &rho0, vectr &mu0) {
	const BN *bn = dynamic_cast<const BN *>(p->GetStartDist());
	const RV *node = bn->NodeByVar(varid); 
	const Context &v = bn->NodeByVar(varid)->Domain();

	vector<int> varL = bn->GetChildrenByVar(varid);
	varL.push_back(varid);
	int states = v.Size();

	vectr factor(states,0.0);

	for(unsigned int i=0; i<varL.size(); i++) {
		const RV *curNode = bn->NodeByVar(varL[i]);
		Instantiation x(curNode->Domain()+curNode->CondDomain());
		x.SetAllVal(0);
		Instantiation nx(x-v);
		do {
			nx.SetVal(x);
			double prob = GetProb(nx,tr->TimeBegin());
			if (prob>0)
				factor[x.Value(varid)] += prob*curNode->Prob(x,true);
		} while(x.Inc());
	}
	for(int i=0; i<factor.getm(); i++) {
		mu0[i] = rho0[i]*exp(factor[i]);
	}
	double musum = mu0.sum();
	if (musum==0.0) mu0 = 1.0/mu0.getm();
	else mu0 /= musum;
	//mu0.normalize();
}

double MeanFieldInf::GetProb(const Instantiation &x, double t,
				bool log) {
	double ret = 1.0;
	vector<int> vList = x.VarList();
	for(unsigned int i=0; i<vList.size();i++) {
		int vid = vList[i];
		int val = x.Value(vid);
		if (val==-1) continue;
		ret *= mus.find(vid)->second.GetVal(t)[val];
		if(ret==0.0) break;
	}
	return !log ? ret : ::log(ret);
}

void MeanFieldInf::AddVar(int varid) {
	mus.insert(make_pair<int,ContFunction<vectr> >(
			varid, ContFunction<vectr>()));

	gammas.insert(make_pair<int,ContFunction<matrix> >(
			varid, ContFunction<matrix>()));

	hiddenIntervals.insert(make_pair<int,vector<double_pair> >(
			varid, vector<double_pair>()));
}

void MeanFieldInf::AddMuVal(int varid, double t, const vectr &vals) {
	map<int,ContFunction<vectr> >::iterator it = mus.find(varid);
	if(it == mus.end()) return;
	it->second.AddVal(t, vals);
}

void MeanFieldInf::AddGammaVal(int varid, double t, const matrix &vals) {
	map<int,ContFunction<matrix> >::iterator it = gammas.find(varid);
	if(it == gammas.end()) return;
	it->second.AddVal(t, vals);
}

void MeanFieldInf::AddHiddenInterval(int varid, double start, double end) {
	map<int,vector<double_pair> >::iterator it = 
		hiddenIntervals.find(varid);
	if(it == hiddenIntervals.end()) return;
	it->second.push_back(double_pair(start,end));
}

void MeanFieldInf::Reset() {
	initialized = false;
	optimized = false;
	if(factored) delete factored;
	factored = nullptr03;

	mus.clear();
	gammas.clear();

	hiddenIntervals.clear();
	transTimes.clear();
}

void MeanFieldInf::SetEPS(double neweps) {
	Reset();
	eps = neweps;
}

} // end of ctbn namespace
