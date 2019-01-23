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
#include "multirv.h"
#include "params.h"
#include "multisimple.h"

namespace ctbn {

using namespace std;

SOBJCLASSDEF(BN)

BN::BN() : RV(Context(),Context()) {
}

BN::BN(istream &is) : RV(is) {
	int n;
	is >> n;
	nodes.resize(n);
	for(int i=0;i<n;i++)
		nodes[i] = dynamic_cast<RV *>(StreamObj::LoadOldPtr(is));
}

BN::BN(const BN &bn) : RV(bn), nodes(bn.nodes.size()) {
	int n = nodes.size();
	for(int i=0;i<n;i++)
		nodes[i] = bn.nodes[i]->Clone();
}

BN::BN(const Structure &s, const Context &vars) : 
				RV(Context(), Context()) {
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
		AddNode(new MultiRV(v,cv));
	}
}
	
BN &BN::operator=(const BN &bn) {
	if (&bn!=this) {
		RV::operator=(bn);
		int n = nodes.size();
		for(int i=0;i<n;i++)
			delete nodes[i];
		n = bn.nodes.size();
		nodes.resize(n);
		for(int i=0;i<n;i++)
			nodes[i] = bn.nodes[i]->Clone();
		sampleorder = bn.sampleorder;
	}
	return *this;
}

BN *BN::Clone() const {
	return new BN(*this);
}

BN::~BN() {
	int n = nodes.size();
	for(int i=0;i<n;i++)
		delete nodes[i];
}

void BN::LoadOld(istream &is) {
	RV::LoadOld(is);
	int n = nodes.size();
	for(int i=0;i<n;i++)
		delete nodes[i];
	is >> n;
	nodes.resize(n);
	for(int i=0;i<n;i++)
		nodes[i] = dynamic_cast<RV *>(StreamObj::LoadOldPtr(is));
	sampleorder.clear();
}

void BN::serial_preload() {
	int n = nodes.size();
	for(int i=0;i<n;i++)
		delete nodes[i];
	sampleorder.clear();
	node2var.clear();
	var2node.clear();
	children.clear();
}
	

void BN::SaveOld(ostream &os) const {
	RV::SaveOld(os);
	int n = nodes.size();
	os << os.fill() << n;
	for(int i=0;i<n;i++) {
		os << os.fill();
		nodes[i]->SaveOldPtr(os);
	}
}

void BN::Normalize() {
	assert(0); // not implemented (requires changing structure of network)
	// below would work, but would create one large node:
	RV *v = Joint(Context()); // don't condition on anything
	v->Normalize();
	sampleorder.clear();
	int n = nodes.size();
	for(int i=0;i<n;i++)
		delete nodes[i];
	nodes.clear();
	nodes.push_back(v);
}

RV *BN::Joint(const Instantiation &x) const {
	//Context rest = (Domain()+CondDomain())-x.KnownVars();

	RV *ans = nodes[0]->Clone();
	ans->Restrict(x);
	//ans->Project(ans->Domain()/rest, ans->CondDomain()/rest);

	int n = nodes.size();
	for(int i=1;i<n;i++) {
		RV *v = nodes[i]->Clone();
		v->Restrict(x);
		//v->Project(v->Domain()/rest, v->CondDomain()/rest);
		ans->MultBy(v);
		delete v;
	}
	return ans;
}

// This is the variable elimination that does
// not sum out until the very end (not efficient)
double BN::Prob(const Instantiation &x, bool log) const {
	RVSimple *v = MakeSimple(x,false);
	double ret = v->Sum(log);
	delete v;
	return ret;
}

void BN::Restrict(const Instantiation &x) {
	int n= nodes.size();
	for(int i=0;i<n;i++)
		nodes[i]->Restrict(x);
}

void BN::Project(const Context &c, const Context &cc) {
	assert(0);
	// not yet implemented
	RV::Project(c,cc);
}

void BN::MultBy(const RV *x) {
	assert(0);
	// again, not yet implemented
}

BNSS *BN::BlankSS() const {
	BNSS *ret = new BNSS();
	int n = nodes.size();
	ret->rvss.resize(n);
	for(int i=0;i<n;i++)
		ret->rvss[i] = nodes[i]->BlankSS();
	return ret;
}

void BN::AddSS(const Instantiation &x, SS *ss, double w) const {
	BNSS *bnss = dynamic_cast<BNSS *>(ss);
	// x should be a complete instantiation, so this should work
	// fine (and it not complex at all)
	int n = nodes.size();
	for(int i=0;i<n;i++)
		nodes[i]->AddSS(x,bnss->rvss[i],w);
}


void BN::AddExpSS(const Instantiation &x, SS *ss, double w) const {
	RVSimple *rvs = MakeSimple(x);
	AddExpSS(x,rvs,ss,w);
	delete rvs;
}

// Again, this one is correct (hopefully) but in no way efficient!
void BN::AddExpSS(const Instantiation &x, const RVSimple *rvs, 
		SS *ss, double w) const {
	BNSS *bnss = dynamic_cast<BNSS *>(ss);
	int n = nodes.size();
	Instantiation knowncond(CondDomain(),x);
	for(int i=0;i<n;i++) {
		Context var = nodes[i]->Domain()+nodes[i]->CondDomain();
		Instantiation inst(var,-1);
		Context conditt = nodes[i]->CondDomain() - CondDomain();
		inst.SetVal(knowncond);
		forallassign_subset(Instantiation &instcond,
					nodes[i]->CondDomain()-CondDomain(),knowncond) {
		
			vector<vector<int> > reind;
			forallassign_subset(Instantiation &inst,
						nodes[i]->Domain(),instcond) {
				vector<int> vv;
				Domain().ConsistentIndexes(vv,inst);
				reind.push_back(vv);
			}
			RVSimple *subrv = rvs->Clone();
			subrv->Reindex(reind);
			nodes[i]->AddExpSS(instcond,subrv,bnss->rvss[i],w);
			delete subrv;
		}
	}
}

void BN::Sample(Instantiation &x, Random &rand) const {
	// Added for sampling more efficiently with no evidence:
	if (x.NumVars()==x.NumMissingVars()) ImportanceSample(x,false,rand);
	else {
		RV *v = Joint(x);
		v->Normalize();
		Instantiation sample(v->Domain());
		v->Sample(sample,rand);
		x.SetVal(sample,true);
		delete v;
	}
}

void BN::Maximize(const SS *ss) {
	const BNSS *bnss = dynamic_cast<const BNSS *>(ss);
	int n = nodes.size();
	for(int i=0;i<n;i++)
		nodes[i]->Maximize(bnss->rvss[i]);
}

RVSimple *BN::MakeSimple(const Instantiation &x, bool normalize) const {
	const bool usesparse = ParamInt("AvoidJoint",1);
	if (usesparse) return MakeSimpleSparse(x,normalize);
	RV *v = Joint(x);
	RVSimple *ret = v->MakeSimple(x,normalize);
	delete v;
	return ret;
}

RVSimple *BN::MakeSimpleSparse(const Instantiation &x, bool normalize) const {
	Instantiation knownx(x.KnownVars(),x);
	Context missing = Domain() - knownx;
	vectr p(missing.Size());
	vector<int> ix(missing.Size());
	Instantiation i(Domain()+CondDomain(),0);
	i.SetVal(knownx);
	int ii = 0;
	forallassign_subset(Instantiation &i, missing, knownx) {
		p[ii] = 1.0;
		for(vector<RV *>::const_iterator n=nodes.begin();n!=nodes.end();++n)
			p[ii] *= (*n)->Prob(i);
		ix[ii] = i.Index();
		ii++;
	}
	SparseMultiZSimple *ret = new SparseMultiZSimple(p,ix);
	if (normalize) ret->Normalize();
	return ret;
}

void BN::AddNode(RV *rv) {
	nodes.push_back(rv);
	v = v+rv->Domain();
	cv = (cv+rv->CondDomain())-v;
	sampleorder.clear();
}


void BN::SetSampleOrder() const {
	sampleorder.clear();
	unsigned int numnode = nodes.size();
	bool *mat = new bool[numnode*numnode];
	bool *added = new bool[numnode];
	int *numpar = new int[numnode];
	for(unsigned int i=0; i<numnode; i++) {
		added[i] = false;
		numpar[i] = 0;
		for(unsigned int j=0; j<numnode; j++)
			mat[i*numnode+j] = 0;
	}

	for(unsigned int i=0; i<numnode; i++) {
		Context conddomain = nodes[i]->CondDomain();
		if(conddomain.VarList().size()==0) continue;
		for(unsigned int j=0; j<numnode; j++) {
			if(i==j) continue;
			Context nodedomain = nodes[j]->Domain();
			if(conddomain.IsSubset(nodedomain)) {
				mat[j*numnode+i] = 1; numpar[i]++;
			}
		}
	}

	bool change = true;
	while (change) {
		change = false;
		for(unsigned int i=0; i<nodes.size(); i++)
			if(!added[i] && numpar[i]==0) {
				added[i]=true;
				sampleorder.push_back(i);
				for(unsigned int j=0; j<numnode; j++)
					if(i!=j && mat[i*numnode+j]) numpar[j]--;
				change = true;
			}
	}
	delete []mat;
	delete []numpar;
	delete []added;
}


double BN::ImportanceSample(Instantiation &x, bool log, Random &rand) const {
	double ret = log?0.0:1.0;

	if (sampleorder.size() < nodes.size())
		SetSampleOrder();

	for(vector<int>::const_iterator i=sampleorder.begin();
			i!=sampleorder.end();++i) {
		// If there is only one variable per node, this is faster;
		if(nodes[*i]->Domain().NumVars()==1) {
			if (x.Value(nodes[*i]->Domain().VarList()[0])!=-1) {
				if (log) ret += nodes[*i]->Prob(x,true);
				else ret *= nodes[*i]->Prob(x,false);
			} else {
				nodes[*i]->Sample(x,rand);
			}
		} else {
			// If not, this is necessary:
			if (log) ret += nodes[*i]->Prob(x,true);
			else ret *= nodes[*i]->Prob(x,false);
			nodes[*i]->Sample(x,rand);
		}
	}
	return ret;
}

void BN::Scramble(double alpha, double degree, Random &rand) {
	int n = nodes.size();
	for(int i=0; i<n; i++)
		nodes[i]->Scramble(alpha, degree, rand);
}

double BN::LLH(const SS *ss) const {
	double llh = 0.0;
	const BNSS *bnss = dynamic_cast<const BNSS *>(ss);
	int n = nodes.size();
	for(int i=0; i<n; i++) {
		llh += nodes[i]->LLH(bnss->rvss[i]);
	}
	return llh;
}

double BN::GetScore(double numTrans, const SS* ss) const {
	const BNSS *bnss = dynamic_cast<const BNSS *>(ss);
	double score(0);
	int numNodes = NumofNodes();
	for(int i = 0; i < numNodes; i++)
		score += nodes[i]->GetScore(numTrans, bnss->rvss[i]);
	return score;
}

// Assumes the source BN has the same organization of nodes.
void BN::FillParams(const BN &src) {
	for(int i = 0; i < NumofNodes(); i++) {
		const MultiRV *csrv = 
			dynamic_cast<MultiRV *>(src.nodes[i]);
		MultiRV *crv = 
			dynamic_cast<MultiRV *>(nodes[i]);
	
		const Context &csrc = csrv->CondDomain();
		const Context &cur = crv->CondDomain();

		Instantiation ins(csrc);
		ins.SetAllVal(0);
		
		int nstates = crv->Domain().Size();

		forallassign(Instantiation &curins, cur) {
			vectr avgDist(nstates,0.0);
			int count = 0;
			// Compute average over extra conditioning vars
			forallassign_subset(Instantiation &ins, csrc-cur ,curins) {
				vectr thisDist(nstates,0.0);
				double logfactor = 0.0;
				(*csrv)[ins].GetDist(thisDist,logfactor);
				avgDist += thisDist;
				count++;			
			} 
			avgDist /= count;
			(*crv)[curins].SetDist(avgDist);
		}
	}
}


void BN::PrepNode2IdMaps() const {
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

void BN::ClearNode2IdMaps() const {
	node2var.clear();
	var2node.clear();
	children.clear();
}

const RV* BN::NodeByVar(int i) const {
	PrepNode2IdMaps();
	int idx = var2node[i];
	return nodes[idx];
}

vector<int> BN::GetChildrenByVar(int varid) const {
	PrepNode2IdMaps();
	vector<int> ret;
	map<int,vector<int> >::iterator it = children.find(varid);
	if(it == children.end()) return ret;
	ret = it->second;
	for(unsigned int i=0; i<ret.size(); i++) ret[i] = node2var[ret[i]];
	return ret;
}

void BN::GetStructure(Structure &s) const {
	PrepNode2IdMaps();
	s.var2node = var2node;
	s.node2var = node2var;
	s.children = children;
	
	int n = nodes.size();
	s.parentlist.clear();
	for (int i=0; i<n; i++) {
		const Context &c = nodes[i]->CondDomain();
		s.parentlist.push_back(c.VarList());
	}	
}



SOBJCLASSDEF(BNSS)

BNSS::BNSS() {
}

BNSS::BNSS(istream &is) {
	int n;
	is >> n;
	rvss.resize(n);
	for(int i=0;i<n;i++)
		rvss[i] = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
}

BNSS::BNSS(const BNSS &bnss) : SS() {
	int n = rvss.size();
	for(int i=0;i<n;i++)
		delete rvss[i];
	n = bnss.rvss.size();
	rvss.resize(n);
	for(int i=0;i<n;i++)
		rvss[i] = bnss.rvss[i]->Clone();
}

BNSS &BNSS::operator=(const BNSS &bnss) {
	if (&bnss!=this) {
		int n = rvss.size();
		for(int i=0;i<n;i++)
			delete rvss[i];
		n = bnss.rvss.size();
		rvss.resize(n);
		for(int i=0;i<n;i++)
			rvss[i] = bnss.rvss[i]->Clone();
	}
	return *this;
}

BNSS *BNSS::Clone() const {
	return new BNSS(*this);
}

BNSS::~BNSS() {
	int n = rvss.size();
	for(int i=0;i<n;i++)
		delete rvss[i];
}

void BNSS::LoadOld(istream &is) {
	int n = rvss.size();
	for(int i=0;i<n;i++)
		delete rvss[i];
	is >> n;
	rvss.resize(n);
	for(int i=0;i<n;i++)
		rvss[i] = dynamic_cast<SS *>(StreamObj::LoadOldPtr(is));
}

void BNSS::serial_preload() {
	int n = rvss.size();
	for(int i=0;i<n;i++)
		delete rvss[i];
}

void BNSS::SaveOld(ostream &os) const {
	int n = rvss.size();
	os << n;
	for(int i=0;i<n;i++) {
		os << os.fill();
		rvss[i]->SaveOldPtr(os);
	}
}


void BNSS::Scale(double w) {
	int n = rvss.size();
	for(int i=0; i<n; i++)
		rvss[i]->Scale(w);
}

void BNSS::AddSS(const SS* nss, double w) {
	int n = rvss.size();
	const BNSS *bss = dynamic_cast<const BNSS*>(nss);
	for(int i=0; i<n; i++)
		rvss[i]->AddSS(bss->rvss[i], w);
}

} // end of ctbn namespace

