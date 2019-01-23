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
#include "ctbndyn.h"
#include "markovdyn.h"
#include "markovsimple.h"
#include <map>
#include "extramath.h"
#include "samplequeue.h"
#include "searchqueue.h"
#include <algorithm>
#include "params.h"

namespace ctbn {

using namespace std;

SOBJCLASSDEF(CTBNDyn)

CTBNDyn::CTBNDyn() : Dynamics(Context(),Context()) {
	joint = NULL;
}

CTBNDyn::CTBNDyn(std::istream &is) : Dynamics(is) {
	int n;
	is >> n;
	nodes.resize(n);
	for(int i=0;i<n;i++)
		nodes[i] = dynamic_cast<Dynamics *>(StreamObj::LoadOldPtr(is));
	joint = NULL;
}

CTBNDyn::CTBNDyn(const CTBNDyn &ctbn) : Dynamics(ctbn) {
	int n = ctbn.nodes.size();
	nodes.resize(n);
	for(int i=0;i<n;i++)
		nodes[i] = ctbn.nodes[i]->Clone();
	joint = NULL;
}

CTBNDyn::CTBNDyn(const Structure &s, const Context &vars) : 
				Dynamics(Context(), Context()) {
	joint = NULL;
	vector<vector<int> > parentlist = s.parentlist;
	int n = parentlist.size();

	for(int i=0;i<n;i++) {
		vector<int> currentParentlist = parentlist[i];
		Context v;
		Context cv;
		if(vars.HasId(i)) 
			v.AddVar(i,vars.Cardinality(i));
		else continue;

		for(unsigned int j=0;j<currentParentlist.size();j++) {
			Context temp;
			int varid = currentParentlist[j];
			if(vars.HasId(varid)) {
				temp.AddVar(varid,vars.Cardinality(varid));
				cv = Context(cv,temp);
			}
		}
		AddNode(new MarkovDyn(v,cv));
	}
}
	

CTBNDyn &CTBNDyn::operator=(const CTBNDyn &ctbn) {
	if (&ctbn!=this) {
		if (joint) delete joint;
		joint = NULL;
		int n = nodes.size();
		for(int i=0;i<n;i++)
			delete nodes[i];
		n = assocss.size();
		for(int i=0;i<n;i++)
			assocss[i]->dyn = NULL;

		n = ctbn.nodes.size();
		nodes.resize(n);
		for(int i=0;i<n;i++)
			nodes[i] = ctbn.nodes[i]->Clone();
		ClearNode2IdMaps();
	}
	return *this;
}

CTBNDyn *CTBNDyn::Clone() const {
	return new CTBNDyn(*this);
}

CTBNDyn::~CTBNDyn() {
	if (joint) delete joint;
	int n = nodes.size();
	for(int i=0;i<n;i++)
		delete nodes[i];
	n = assocss.size();
	for(int i=0;i<n;i++)
		assocss[i]->dyn = NULL;
}

void CTBNDyn::LoadOld(istream &is) {
	if (joint) delete joint;
	joint = NULL;
	Dynamics::LoadOld(is);
	int n = nodes.size();
	for(int i=0;i<n;i++)
		delete nodes[i];
	n = assocss.size();
	for(int i=0;i<n;i++)
		assocss[i]->dyn = NULL;
	assocss.resize(0);

	is >> n;
	for(int i=0;i<n;i++)
		nodes[i] = dynamic_cast<Dynamics *>(StreamObj::LoadOldPtr(is));

	ClearNode2IdMaps();
}

void CTBNDyn::serial_preload() {
	if (joint) delete joint;
	joint = NULL;
	int n = nodes.size();
	for(int i=0;i<n;i++) delete nodes[i];
	n = assocss.size();
	for(int i=0;i<n;i++)
		assocss[i]->dyn = NULL;
	assocss.resize(0);
	ClearNode2IdMaps();
}

void CTBNDyn::SaveOld(ostream &os) const {
	Dynamics::SaveOld(os);
	int n = nodes.size();
	os << os.fill() << n;
	for(int i=0;i<n;i++) {
		os << os.fill();
		nodes[i]->SaveOldPtr(os);
	}
}

void CTBNDyn::Normalize() {
	int n = nodes.size();
	if (joint) delete joint;
	joint = NULL;
	for(int i=0;i<n;i++)
		nodes[i]->Normalize();
}

void CTBNDyn::Restrict(const Instantiation &x) {
	int n = nodes.size();
	if (joint) delete joint;
	joint = NULL;
	for(int i=0;i<n;i++)
		nodes[i]->Restrict(x);
	ClearNode2IdMaps(); // could change...
}

void CTBNDyn::Mult(const Dynamics *x) {
	// not really sure how to do this...
	assert(false);
}

RVCondSimple *CTBNDyn::Cond(double t0, double t1, 
				const Instantiation &x) const {
	static bool usesparse = ParamInt("AvoidJoint",1);
	if (!usesparse) return Joint()->Cond(t0,t1,x);

	Instantiation knownx(x.KnownVars(),x);
	// assumes that only one variable can change at a time
	Context missing = Domain()+CondDomain()-knownx;
	int nmiss = missing.Size();
	matrix subQ(nmiss,nmiss,0.0);

	vector<int> ix(nmiss);
	map<int,int> rev;
	int ii = 0;
	forallassign_subset(Instantiation &i1,missing,knownx) {
		ix[ii] = i1.Index();
		rev[i1.Index()] = ii;
		ii++;
	}

	ii = 0;
	forallassign_subset(Instantiation &i1,missing,knownx) {
		for(vector<Dynamics *>::const_iterator n=nodes.begin();
						n!=nodes.end();++n) {
			const DynComp *dyn = dynamic_cast<const DynComp *>(*n);
			if (dyn==NULL) return Joint()->Cond(t0,t1,x);
			const MarkovSimple *simp =
				dynamic_cast<const MarkovSimple *>((*dyn)[i1]);
			if (simp==NULL) return Joint()->Cond(t0,t1,x);
			const matrix &localQ = simp->Intensity();
			const Context &locald = dyn->Domain();
			int lii1 = locald.Index(i1);
			subQ[ii][ii] += localQ[lii1][lii1];
			Context ch = locald / missing;
			foralldiffone_subset(Instantiation &i2,i1,ch) {
				int lii2 = locald.Index(i2);
				subQ[ii][rev[i2.Index()]] += localQ[lii1][lii2];
			}
		}	
		ii++;
	}

	return new SparseCondTransQ<matrix>(subQ,t1-t0,ix);
}

RVCondSimple *CTBNDyn::Cond(double t, const Instantiation &from,
		const Instantiation &to, bool transition) const {
	static bool usesparse = ParamInt("AvoidJoint",1);
	if (!usesparse) return Joint()->Cond(t,from,to,transition);


	vector<int> fromind,toind;
	Context all = Domain()+CondDomain();
	all.ConsistentIndexes(fromind,from);
	all.ConsistentIndexes(toind,to);
	matrix subQ(fromind.size(),toind.size(),0.0);
	// assumes that only one variable can change at a time
	if (transition) {
		Instantiation knownfrom(from.KnownVars(),from);
		Instantiation knownto(to.KnownVars(),to);
		vector<int> inbothvar = (knownfrom / knownto).VarList();
		int chvar = -1;
		for(vector<int>::iterator v=inbothvar.begin();v!=inbothvar.end();++v)
			if (knownfrom.Value(*v)!=-1 && knownto.Value(*v)!=-1 &&
						knownfrom.Value(*v) != knownto.Value(*v)) {
				chvar = *v;
				break;
			}
		assert(chvar != -1);
		PrepNode2IdMaps();
		const DynComp *dyn =
				dynamic_cast<const DynComp *>(nodes[var2node[chvar]]);
		if (dyn==NULL) return Joint()->Cond(t,from,to,transition);
		const Context &locald = dyn->Domain();
			
		map<int,int> torev;
		for(unsigned int i=0;i<toind.size();i++)
			torev[toind[i]] = i;

		int newv = to.Value(chvar);
		
		int ii = 0;
		forallassign_except(Instantiation &ifrom,
							Domain()+CondDomain(),knownfrom) {
			const MarkovSimple *simp =
					dynamic_cast<const MarkovSimple *>((*dyn)[ifrom]);
			if (simp==NULL) return Joint()->Cond(t,from,to,transition);
			const matrix &localQ = simp->Intensity();

			Instantiation ito(ifrom);
			ito.SetVal(chvar,newv);
			int liito = locald.Index(ito);
			int liifrom = locald.Index(ifrom);
			subQ[ii][torev[ito.Index()]] += localQ[liifrom][liito];
			ii++;
		}
		
	} else {
		// ripped straight from markovsimple.cc (should perhaps merge):
          unsigned int i=0,j=0;
          while(i<fromind.size() && j<toind.size()) {
               if (fromind[i]==toind[j]) {
                    subQ[i][j] = 1.0;
                    i++;
                    j++;
               } else if (fromind[i]<toind[j]) {
                    i++;
               } else {
                    j++;
               }
          }
	}
	return new SparseCondTransQ2<matrix>(subQ,fromind,toind);
}


CTBNDynSS *CTBNDyn::BlankSS() const {
	CTBNDynSS *ret = new CTBNDynSS(this);
	assocss.push_back(ret);
	int n = nodes.size();
	ret->dynss.resize(n);
	for(int i=0;i<n;i++)
		ret->dynss[i] = nodes[i]->BlankSS();
	ret->jss = NULL;
	return ret;
}

void CTBNDyn::AddSS(const Instantiation &x, double t0, double deltat,
		SS *ss, double w) const {
	CTBNDynSS *css = dynamic_cast<CTBNDynSS *>(ss);
	int n = nodes.size();
	for(int i=0;i<n;i++)
		nodes[i]->AddSS(x,t0,deltat,css->dynss[i],w);
}

void CTBNDyn::AddTransSS(const Instantiation &x1, const Instantiation &x2,
		double t, SS *ss, double w) const {
	CTBNDynSS *css = dynamic_cast<CTBNDynSS *>(ss);
	int n = nodes.size();
	for(int i=0;i<n;i++)
		nodes[i]->AddTransSS(x1,x2,t,css->dynss[i],w);
}

void CTBNDyn::AddExpSS(const RV *alpha, const RV *beta,
		double t0, double deltat, SS *ss, double w) const {
	CTBNDynSS *css = dynamic_cast<CTBNDynSS *>(ss);
	int n = nodes.size();
	for(int i=0;i<n;i++)
		nodes[i]->AddExpSS(alpha,beta,t0,deltat,css->dynss[i],w);
}

void CTBNDyn::AddExpTransSS(const RV *x1, const RV *x2,
		const Context &changevar,
		double t, SS *ss, double w) const {
	if (!changevar.IsOverlap(Domain())) return;
	CTBNDynSS *css = dynamic_cast<CTBNDynSS *>(ss);
	int n = nodes.size();
	for(int i=0;i<n;i++)
		nodes[i]->AddExpTransSS(x1,x2,changevar,t,css->dynss[i],w);
}

void CTBNDyn::AddExpSS(const Instantiation &condi, const RVSimple *alpha,
		const RVSimple *beta, double t0, double deltat, SS *ss,
		double w) const {
    CTBNDynSS *css = dynamic_cast<CTBNDynSS *>(ss);
    if (NULL == css->jss) css->jss = Joint()->BlankSS();
    Joint()->AddExpSS(condi, alpha, beta, t0, deltat, css->jss, w);
/*
	SS *allss = Joint()->BlankSS();
	Joint()->AddExpSS(condi,alpha,beta,t0,deltat,allss,w);
	CTBNDynSS *css = dynamic_cast<CTBNDynSS *>(ss);
	int n = nodes.size();
	for(int i=0;i<n;i++)
		nodes[i]->AddSS(allss,Joint(),css->dynss[i]);
	delete allss;
*/
}

void CTBNDyn::AddSS(const SS *toadd, const Dynamics *dyn,
		SS *ss, double w) const {
	int n = nodes.size();
	CTBNDynSS *ctbndynss = dynamic_cast<CTBNDynSS *>(ss);
	for(int i=0;i<n;i++)
		nodes[i]->AddSS(toadd,dyn,ctbndynss->dynss[i],w);
}

void CTBNDyn::AddExpSS(int condi, const RVSimple *alpha, const RVSimple *beta,
		double t0, double deltat, SS *ss, double w) const {
	Instantiation i(CondDomain());
	i.SetIndex(condi);
	AddExpSS(i,alpha,beta,t0,deltat,ss,w);
}

void CTBNDyn::AddExpTransSS(const Instantiation &condi, const RVSimple *x1,
		const RVSimple *x2, const Context &changevar,
		double t, SS *ss, double w) const {
	RVSimple *myx1 = x1->Clone();
	RVSimple *myx2 = x2->Clone();
	myx1->Normalize();
	myx2->Normalize();
	CTBNDynSS *css = dynamic_cast<CTBNDynSS *>(ss);
	int n = nodes.size();
	for(int i=0;i<n;i++) {
		if (!nodes[i]->Domain().IsOverlap(changevar)) continue;
		Context cond = nodes[i]->CondDomain()-CondDomain();
		vectr condprob(nodes[i]->CondDomain().Size(), 0.0);
		forallassign_subset(Instantiation &transii,cond,condi) {
			vector<vector<int> > mapping(nodes[i]->Domain().Size());
			forallassign_subset(Instantiation &ii,nodes[i]->Domain(),transii)
				Domain().ConsistentIndexes(
					mapping[nodes[i]->Domain().Index(ii)],ii);

			RVSimple *t1 = myx1->Clone();
			t1->Reindex(mapping);
			RVCondSimple *p1 = nodes[i]->Cond(t, transii, transii, true);
			p1->Mult(t1);
			RVSimple *t2 = myx2->Clone();
			t2->Reindex(mapping);
			t1->MultBy(t2);
			condprob[cond.Index(transii)] = t1->Sum();
			delete t2;
			delete p1;
			delete t1;
		}
		condprob.normalize();
		forallassign_subset(Instantiation &transii,cond,condi) {
			vector<vector<int> > mapping(nodes[i]->Domain().Size());
			forallassign_subset(Instantiation &ii,nodes[i]->Domain(),transii)
				Domain().ConsistentIndexes(
					mapping[nodes[i]->Domain().Index(ii)],ii);
			RVSimple *t1 = myx1->Clone();
			t1->Reindex(mapping);
			t1->Normalize();
			RVSimple *t2 = myx2->Clone();
			t2->Reindex(mapping);
			t2->Normalize();
			double condw = w*condprob[cond.Index(transii)];
			if (condw>0.0)
				nodes[i]->AddExpTransSS(transii,t1,t2,changevar,t,
									css->dynss[i],condw);
			delete t2;
			delete t1;
		}
	}

	delete myx1;
	delete myx2;
}

void CTBNDyn::AddExpTransSS(int condi, const RVSimple *x1, const RVSimple *x2,
		const Context &changevar,
		double t, SS *ss, double w) const {
	Instantiation i(CondDomain());
	i.SetIndex(condi);
	AddExpTransSS(i,x1,x2,changevar,t,ss,w);
}

void CTBNDyn::SampleNextEvent(const Instantiation &i, double t, double &nextt,
		Instantiation &nexti, Random &rand) const {
	int n = nodes.size();
	nextt = INFINITY;
	for(int j=0;j<n;j++) {
		double tt = t;
		Instantiation ni(i);
		nodes[j]->SampleNextEvent(i,t,tt,ni,rand);
		if (finite(tt) && tt<nextt) {
			nextt = tt;
			nexti = ni;
		}
	}
}

void CTBNDyn::Maximize(const SS *ss) {
	const CTBNDynSS *css = dynamic_cast<const CTBNDynSS *>(ss);
	int n = nodes.size();
	css->Compact();
	for(int i=0;i<n;i++)
		nodes[i]->Maximize(css->dynss[i]);
	if (joint) delete joint;
	joint = NULL;
}

void CTBNDyn::PrepNode2IdMaps() const {
	if (node2var.size() == nodes.size()) return;

	ClearNode2IdMaps();
	int m = nodes.size();
	node2var.resize(m);
	for(int j=0;j<m;j++) {
		node2var[j] = nodes[j]->Domain().VarList()[0];
		var2node[node2var[j]] = j;

		vector<int> par = nodes[j]->CondDomain().VarList();
		for(unsigned int k=0;k<par.size();k++)
			children[par[k]].push_back(j);
	}
	for(int j=0;j<m;j++) {
		if(children.find(node2var[j]) == children.end())
			children.insert(make_pair(node2var[j],
				vector<int>()));
	}
}

void CTBNDyn::ClearNode2IdMaps() const {
	node2var.clear();
	var2node.clear();
	children.clear();
}

// We could just keep the default SampleTrajectory.  However that has 
// a running time O(nm) where n is the number of samples and m is the number
// of variables.  This method is O(m+nkln(m)) where k is the maximal number of
// parents in the graph.
// This currently makes the assumption that each node only corresponds to
// one variable
void CTBNDyn::SampleTrajectory(Trajectory &tr, double t, Random &rand) const {
	int m = nodes.size();
	SampleQueue events(m);

	double endt = tr.TimeEnd();
	Instantiation i(tr.Values(Domain(),t));

	PrepNode2IdMaps();

	for(int j=0;j<m;j++) {
		Instantiation newv;
		double newt;
		nodes[j]->SampleNextEvent(i,t,newt,newv,rand);
		events.Add(SampleQueue::Event(j,newv.Value(node2var[j]),
						newt));
	}

	while(1) {
		SampleQueue::Event e;
		if (!events.Head(e) || e.time>endt) return;
		int varid = node2var[e.var];
		tr.AddTransition(varid, e.time, e.value);
		i.SetVal(varid,e.value);
		t = e.time;

		events.Remove(e.var);
		map<int,vector<int> >::iterator ci = 
			children.find(node2var[e.var]);
		vector<int> c;
		if (ci != children.end() && ci->first == node2var[e.var])
			c = ci->second;
		for(unsigned int j=0;j<c.size();j++)
			events.Remove(c[j]);

		Instantiation newv;
		double newt;
		nodes[e.var]->SampleNextEvent(i,t,newt,newv,rand);
		events.Add(SampleQueue::Event(e.var,
					newv.Value(node2var[e.var]),newt));
		for(unsigned int j=0;j<c.size();j++) {
			Instantiation newv;
			nodes[c[j]]->SampleNextEvent(i,t,newt,newv,rand);
			events.Add(SampleQueue::Event(c[j],
					newv.Value(node2var[c[j]]),newt));
		}
	}
}

void CTBNDyn::AddNode(Dynamics *d) {
	nodes.push_back(d);
	v = v+d->Domain();
	cv = (cv+d->CondDomain())-v;
	if (joint!=NULL) {
		delete joint;
		joint = NULL;
	}
	ClearNode2IdMaps();
}

const Dynamics *CTBNDyn::Joint() const {
	if (!joint) {
		int n = nodes.size();
		if (n==0) return NULL;
		joint = nodes[0]->Clone();
		for(int i=1;i<n;i++) {
			joint->Mult(nodes[i]);
		}
	}
	return joint;
}

void CTBNDyn::GetStructure(Structure &s) const {
	PrepNode2IdMaps();
	s.var2node = var2node;
	s.node2var = node2var;
	s.children = children;
	s.parentlist.clear();
	int n = nodes.size();
	for (int i=0; i<n; i++) {
		Context c = nodes[i]->CondDomain();
		s.parentlist.push_back(c.VarList());
	}
}

void CTBNDyn::Scramble(double a, double b, double alpha, double degree, 
			Random &rand) {
	//Added by Yu
	int nc = nodes.size();
	for(int i=0;i<nc;i++)
		nodes[i]->Scramble(a, b, alpha, degree, rand);
}

double CTBNDyn::LLH(const SS *ss) const {
	const CTBNDynSS *css = dynamic_cast<const CTBNDynSS *>(ss);
	css->Compact();
	int nc = nodes.size();
	double llh = 0.0;
	for(int i=0; i<nc; i++)
		llh += nodes[i]->LLH(css->dynss[i]);
	return llh;
}

matrix CTBNDyn::JointMatrix() const {
	const MarkovDyn *jointp = dynamic_cast<const MarkovDyn *>(Joint());
	matrix m = (*jointp)(0)->Intensity();
	return m;
}


double CTBNDyn::GetScore(double numTrans, 
			 double amtTime, const SS *ss) const {
	const CTBNDynSS *css = dynamic_cast<const CTBNDynSS* >(ss);
	css->Compact();
	double score(0);
	int vars = css->dynss.size();
	double localscore;
	// Sum over each variable
	for(int i=0;i<vars;i++) {
		localscore = nodes[i]->
			GetScore(numTrans, amtTime, css->dynss[i]);
		score += localscore;
	}
	return score;
}

// Assumes the source CTBNDyn has the same organization of nodes.
void CTBNDyn::FillParams(const CTBNDyn &src, bool useRandom) {
	for(int i = 0; i < NumofNodes(); i++) {
		const MarkovDyn *csdyn = dynamic_cast<MarkovDyn *>(src.nodes[i]);
		MarkovDyn *cdyn = dynamic_cast<MarkovDyn *>(nodes[i]);
		const Context &csrc = csdyn->CondDomain();
		const Context &cur = cdyn->CondDomain();

		if (useRandom) {
			forallassign(Instantiation &curins,cur) {
				int idx = randomizer.RandInt(csrc.Size());
				(*cdyn)(curins)->Intensity() = (*csdyn)(idx)->Intensity();
			}
		} else {
			int nstates = cdyn->Domain().Size();
			forallassign(Instantiation &curins,cur) {
				matrix AIM(nstates,nstates,0.0);
				// Compute average IM over extra vars
				int count = 0;
				forallassign_except(Instantiation &ins,csrc,curins) {
					AIM += (*csdyn)(ins)->Intensity();
					count++;
				}
				AIM /= count;
				(*cdyn)(curins)->Intensity() = AIM;
			}
		}
	}
}


const Dynamics *CTBNDyn::NodeByVar(int i) const {
	PrepNode2IdMaps();
	int idx = var2node[i];
	return nodes[idx];
}

const SS *CTBNDyn::NodeSSByVar(int i, SS* ss) const {
	PrepNode2IdMaps();
	int idx = var2node[i];
	const CTBNDynSS *css = dynamic_cast<const CTBNDynSS* >(ss);
	css->Compact();
	return css->dynss[idx];
}

vector<int> CTBNDyn::GetChildrenByVar(int varid) const {
	PrepNode2IdMaps();
	vector<int> ret;
	map<int,vector<int> >::iterator it = children.find(varid);
	if(it == children.end()) return ret;
	ret = it->second;
	for(unsigned int i = 0; i < ret.size(); i++) 
		ret[i] = node2var[ret[i]];
	return ret;
}

SOBJCLASSDEF(CTBNDynSS)

CTBNDynSS::CTBNDynSS(istream &is) {
	int n;
	is >> n;
	dynss.resize(n);
	for(int i=0;i<n;i++)
		dynss[i] = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
	bool hasjss;
	is >> hasjss;
	if (hasjss) jss = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
	else jss = NULL;
	dyn = NULL;
}

CTBNDynSS::CTBNDynSS(const CTBNDyn *d) {
	jss = NULL;
	dyn = d;
}

CTBNDynSS::CTBNDynSS(const CTBNDynSS &ctbnss) : SS() {
	int n = ctbnss.dynss.size();
	dynss.resize(n);
	for(int i=0;i<n;i++)
		dynss[i] = ctbnss.dynss[i]->Clone();
    if (NULL == ctbnss.jss) {
        jss = NULL;
    } else {
        jss = ctbnss.jss->Clone();
    }
	dyn = ctbnss.dyn;
	if (dyn) dyn->assocss.push_back(this);
}

CTBNDynSS::~CTBNDynSS() {
	int n = dynss.size();
	for(int i=0;i<n;i++)
		delete dynss[i];
	if (jss) delete jss;
	if (dyn) dyn->assocss.erase(remove(dyn->assocss.begin(),dyn->assocss.end(),this),dyn->assocss.end());
}

CTBNDynSS &CTBNDynSS::operator=(const CTBNDynSS &ctbnss) {
	if (&ctbnss!=this) {
		int n = dynss.size();
		for(int i=0;i<n;i++)
			delete dynss[i];
		if (jss) delete jss;
		n = ctbnss.dynss.size();
		dynss.resize(n);
		for(int i=0;i<n;i++)
			dynss[i] = ctbnss.dynss[i]->Clone();
		if (ctbnss.jss == NULL) jss = NULL;
		else jss = ctbnss.jss->Clone();
		if (dyn) dyn->assocss.erase(remove(dyn->assocss.begin(),dyn->assocss.end(),this),dyn->assocss.end());
		dyn = ctbnss.dyn;
		if (dyn) dyn->assocss.push_back(this);
	}
	return *this;
}

CTBNDynSS *CTBNDynSS::Clone() const {
	return new CTBNDynSS(*this);
}

void CTBNDynSS::Compact() const {
	if (jss != NULL) {
		assert(dyn!=NULL);
//		cerr << "compacting jss: " << endl;
//		jss->SaveV(cerr);
		int n = dyn->nodes.size();
		for(int i=0;i<n;i++) {
			dyn->nodes[i]->AddSS(jss,dyn->Joint(),dynss[i]);
//			cerr << i << ": " << endl;
//			dynss[i]->SaveV(cerr);
		}
		delete jss;
		jss = NULL;
	}
}

void CTBNDynSS::LoadOld(istream &is) {
	int n = dynss.size();
	for(int i=0;i<n;i++)
		delete dynss[i];
	if (dyn) dyn->assocss.erase(remove(dyn->assocss.begin(),dyn->assocss.end(),this),dyn->assocss.end());

	is >> n;
	dynss.resize(n);
	for(int i=0;i<n;i++)
		dynss[i] = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
	bool hasjss;
	is >> hasjss;
	if (hasjss) jss = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
	else jss = NULL;
	dyn = NULL;
}

void CTBNDynSS::serial_preload() {
	int n = dynss.size();
	for(int i=0;i<n;i++)
		delete dynss[i];
	if (dyn) dyn->assocss.erase(remove(dyn->assocss.begin(),
				dyn->assocss.end(),this),dyn->assocss.end());
	dyn = NULL;
}

void CTBNDynSS::serial_presave() const {
	Compact();
}

void CTBNDynSS::SaveOld(ostream &os) const {
	int n = dynss.size();
	os << n;
	for(int i=0;i<n;i++) {
		os << os.fill();
		dynss[i]->SaveOldPtr(os);
	}
	bool hasjss = (jss!=NULL);
	os << os.fill() << hasjss;
	if (hasjss) {
		os << os.fill();
		jss->SaveOldPtr(os);
	}
}

void CTBNDynSS::Scale(double w) {
	int n = dynss.size();
	for (int i = 0; i<n; i++)
		dynss[i]->Scale(w);
    if (jss) {
        jss->Scale(w);
    }
}

double CTBNDynSS::NodeSS(int id, int val1, int val2, 
				const Instantiation &cond) const {
	vector<int> index;
	Compact();
	//cond.ConsistentIndexes(index);
	assert(dyn!=NULL);
	dyn->PrepNode2IdMaps();
	return dynamic_cast<const DynCompSS *>(dynss[dyn->var2node[id]])->
		Element(index, val1, val2);
}

void CTBNDynSS::AddSS(const SS* nss, double w) {
	int n = dynss.size();
	const CTBNDynSS *css = dynamic_cast<const CTBNDynSS*>(nss);
	for(int i=0; i<n; i++)
		dynss[i]->AddSS(css->dynss[i], w);
	if (css->jss!=NULL) {
		if (!jss) jss = dyn->Joint()->BlankSS();
		jss->AddSS(css->jss,w);
	}
}

} // end of ctbn namespace
