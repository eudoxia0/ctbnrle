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
#include "em.h"
#include "params.h"
#include "utils.h"
#include "markov.h"
#include "structure.h"
#include "ctbndyn.h"
#include "dynamics.h"
#include "expsuffstatsquery.h"
#include "nonexpsuffstatsquery.h"
#include "brutestructuresearch.h"
#include "grapheditsearch.h"
#include "famscore.h"
#include "nullptr03.h"

#ifdef USE_PTHREADS
#include <pthread.h>
#endif

namespace ctbn {

using namespace std;

SS *EStep(const vector<Trajectory> &data, int si, int ei,
		Process *dist, Inference *inf) {
	SS *ss = dist->BlankSS();
	inf->SetProcess(dist);
	for(int j=si;j<ei;j++) {
		inf->SetTrajectory(&data[j]);
		inf->AddExpSuffStats(ss);
	}
	return ss;
}

#ifdef USE_PTHREADS
struct EThreadInfo {
	int si,ei;
	Process *dist;
	Inference *inf;
	const vector<Trajectory> *data;
	pthread_t thread;
};

void *RunEThread(void *iptr) {
	EThreadInfo *tinfo = (EThreadInfo *)iptr;
	SS *ret = EStep(*(tinfo->data),tinfo->si,tinfo->ei,tinfo->dist,tinfo->inf);
	return ret;
}
#endif

void EM(const vector<Trajectory> &data, Process *dist, Inference *inf) {
	// max # of iterations
	int nitt = ParamInt("EMItt", 5); 
	// LLH improvement (per example) below which loop will terminate early
	// (negative value means never terminate early)
	double llhimprov = ParamDouble("EMLLHImprov",-1);
	int verbose = ParamInt("VerboseEM",0);
	int ndata = data.size();
#ifdef USE_PTHREADS
	int nthreads = ParamInt("EMNumThreads",1);
	EThreadInfo *tinfo = nthreads<=1 ? nullptr03 : new EThreadInfo[nthreads-1];
	pthread_attr_t tattr;
	pthread_attr_init(&tattr);
	pthread_attr_setdetachstate(&tattr,PTHREAD_CREATE_JOINABLE);
	for(int i=0;i<nthreads-1;i++) {
		tinfo[i].si = ndata*(i+1)/nthreads;
		tinfo[i].ei = ndata*(i+2)/nthreads;
		tinfo[i].dist = nullptr03;
		tinfo[i].data = &data;
		tinfo[i].inf = inf->Clone();
	}
#endif

	double oldllh = -INFINITY;
	for(int i = 0; i < nitt; i++) {
		if (verbose) cout << "EMItt=" << i+1 << endl;

		SS *ss;
#ifdef USE_PTHREADS
		if (nthreads<=1) { // single thread...
#endif
			ss = EStep(data,0,ndata,dist,inf);
#ifdef USE_PTHREADS
		} else { // multiple threads...
			for(int i=0;i<nthreads-1;i++) { // fire up
				if (tinfo[i].dist != nullptr03) delete tinfo[i].dist;
				tinfo[i].dist = dist->Clone();
				int rc = pthread_create(&(tinfo[i].thread),&tattr,
					RunEThread,(void *)(&(tinfo[i])));
				if (rc!=0) {
					cerr << "Error # " << rc << " in starting thread in EM" << endl;
					exit(1);
				}
			}
			inf->SetProcess(dist); // run first batch in main thread
			ss = EStep(data,0,ndata/nthreads,dist,inf);
			for(int i=0;i<nthreads-1;i++) { // collect results
				SS *nss;
				int rc = pthread_join(tinfo[i].thread,(void **)(&nss));
				if (rc!=0) {
					cerr << "Error # " << rc << " in ending thread in EM" << endl;
					exit(1);
				}
				ss->AddSS(nss);
				delete nss;
			}
		}
#endif
		if (verbose>2) ss->SaveV(cout);
		if (verbose>1) cout << "\tLLH/example=" << dist->LLH(ss)/ndata << endl << "...maximize..." << endl;
		dist->Maximize(ss);
		if (llhimprov>=0) {
			double newllh = dist->LLH(ss)/ndata;
			if (verbose) cout << "\tLLH/example=" << newllh << endl;
			if (newllh-oldllh<llhimprov) {
				delete ss;
				break;
			}
			oldllh = newllh;
		} else if (verbose) {
			double newllh = dist->LLH(ss)/ndata;
			if (verbose) cout << "\tLLH/example=" << newllh << endl;
		}
		if (verbose>1) dist->SaveV(cout);
		delete ss;
	}
#ifdef USE_PTHREADS
	for(int i=0;i<nthreads-1;i++) {
		if (tinfo[i].dist != nullptr03) delete tinfo[i].dist;
		delete tinfo[i].inf;
	}
	pthread_attr_destroy(&tattr);
	if (tinfo != nullptr03) delete[] tinfo;
#endif
}

void SEM(const vector<Trajectory> &data, 
		Process *&dist, 
		Inference *inf) {
	int nitt = ParamInt("SEMItt", 10);
	double nTrans = ParamDouble("alpha", 1.0);
	double nTime = ParamDouble("tau", 1.0);
	int verbose = ParamInt("VerboseSEM",ParamInt("VerboseEM",0));
	Structure s;

	// Initial parameter update
	EM(data,dist,inf);

	int strmod = 1;
	for(int i = 0; i < nitt && strmod>0; i++) {
		if (verbose) cout << "SEMItt=" << i+1 << endl;
		ExpSuffStatsQuery* essQuery = 
			new ExpSuffStatsQuery(inf->Clone(),true);
		essQuery->SetData(&data);
		essQuery->SetProcess(dist);

		// Structure modification step
		Markov* m = dynamic_cast<Markov *>(dist);
		const BN* bnet = dynamic_cast<const BN *>(m->GetStartDist());
		
		const CTBNDyn* cdyn = 
			dynamic_cast<const CTBNDyn *>(m->GetDynamics());
		
		BNFamScore* bfs = new BNFamScore(nTrans,essQuery);
		GraphEditSearch gesBN(bnet->Domain(),bfs,true);
		if (verbose) cout << "Learning BN..." << endl;
		s = gesBN.LearnStructure();
		delete bfs;
		
		RV* newrv = new BN(s, bnet->Domain());
		dynamic_cast<BN *>(newrv)->FillParams(*bnet);

		Structure olds;
		bnet->GetStructure(olds);
		strmod = olds.HammingDist(s);
		
		CTBNFamScore* cfs = new CTBNFamScore(nTrans,nTime,essQuery);
		GraphEditSearch ges(cdyn->Domain(),cfs);
		if (verbose) cout << "Learning CTBN..." << endl;
		s = ges.LearnStructure();
		delete cfs;

		delete essQuery;

		Dynamics* newdyn = new CTBNDyn(s, cdyn->Domain());
		dynamic_cast<CTBNDyn *>(newdyn)->FillParams(*cdyn);

		if (strmod==0) {
			cdyn->GetStructure(olds);
			strmod = olds.HammingDist(s);
		}

		Process* newdist = new Markov(newrv, newdyn);
		delete dist;
		dist = newdist;
		
		// Inference and parameter update step
		EM(data, dist, inf);
	}
}

} // end of ctbn namespace
