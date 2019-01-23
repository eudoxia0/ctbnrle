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
#include "graphic.h"
#include "params.h"

namespace ctbn {

using namespace std;

// If the "dot" command (from graphviz) is not in your
// path, or if you wish to define an explicit path to it,
// change the string below accordingly
#define DOTCOMMAND "/usr/bin/dot"

/*
 * GraphicEvent is class for generating trajectory graph
 */

//return the difference between two Instantiation
vector<int> diffvalues(const Instantiation &i1, const Instantiation &i2) {
     Context both(i1,i2,Context::UNION);
     vector<int> vars = both.VarList();

     vector<int> ret;
     for(vector<int>::iterator i=vars.begin();i!=vars.end();++i)
          if (i1.Value(*i) != i2.Value(*i))
               ret.push_back(*i);

     return ret;
}


Graphic::GraphicEvent::GraphicEvent() {
	// Auto-generated constructor stub
	change_node = 0;
	time = 0;
	state.clear();
}

void Graphic::GraphicEvent::SetChanges(int node,
				       double change_time,
				       const vector<int> &change_state) {
	change_node = node;
	time = change_time;
	state = change_state;
}

double Graphic::GraphicEvent::GetTime() const {
	return time;
}

int Graphic::GraphicEvent::GetChange() const {
	return change_node;
}

vector <int> Graphic::GraphicEvent::GetState() const {
	return state;
}

void Graphic::GraphicEvent::SetName(string insertname) {
	name = insertname;
}

string Graphic::GraphicEvent::GetName() const{
	return name;
}


/*
 * GraphicNode class is for drawing the CTBN graph.
 * this class will save all the data need to be generate the CTBN graph
 */

Graphic::GraphicNode::GraphicNode() {
	a = "";
	b = 0;
	check_string = false;
}

string Graphic::GraphicNode :: GetName() {
	return a;
}

int Graphic::GraphicNode :: GetSize() {
	return matching.size();
}

void Graphic::GraphicNode::Add(const string & name) {
	a = name;
}

void Graphic::GraphicNode::Add(const int & name) {
	exist = false;
	check_string = true;
	stringstream ss;
	ss << name;
	ss >> a;
	ss.clear(stringstream::goodbit);
}

void Graphic::GraphicNode::SetMatrix(const matrix & matrix, int & conn, int & locate) {
	dyn.push_back(matrix);
	matching.push_back(conn);
	place.push_back(locate);
}

vector <int> Graphic::GraphicNode::GetPlace() {
	return place;
}

vector <int> Graphic::GraphicNode::GetMatching() {
	return matching;
}

void Graphic::GraphicNode::SetMatching(vector <int> matching) {
	exist = true;
	this->matching = matching;
}

void Graphic::GraphicNode::SetPlace(int location) {
	place.push_back(location);
}

void Graphic::GraphicNode::SetMatrix(const matrix & matrix) {
	dyn.push_back(matrix);
}

void Graphic::GraphicNode::SetBool() {
	exist = true;
}

bool Graphic::GraphicNode::GetExist() {
	return exist;
}

matrix Graphic::GraphicNode::GetMatrix(const int & i) {
	return dyn[i];
}

bool Graphic::GraphicNode::CheckStrings() {
	return check_string;
}

void Graphic::GraphicNode::SetName(const string & name) {
	a = name;
	check_string = true;
}


//initialize graphic class
Graphic :: Graphic() {
	d = NULL;
	size = 0;
	//size of x coordinate the trajectory visualization
	size_x = ParamInt("traj_win_size", 800);
	//size of the font on trajectory
	change_point = ParamInt("traj_font", 15);
	//size of changed time mark in trajectory graph
	size_distance = ParamInt("traj_ruler", 5);
	//size of the font on ctbn graph
	size_font = ParamInt("ctbn_font", 16);
	//start time of the trajectory
	start = 0;
	//ending time of the trajectory
	finish = 2;
	largest_tr = 0;
	change_name = false;
}

//outputing both graph representation of Markov and trajectory

void Graphic::DrawBoth(const Markov & A, Trajectory & traj,
		       ofstream & ls, string & name, ofstream & lb) {
	this->start_point = traj.TimeBegin();
	this->ending_point = traj.TimeEnd();
	const CTBNDyn * B = dynamic_cast<const CTBNDyn *>(A.GetDynamics());
	size = B->NumofNodes();
	Context context = Context();
	for (int i = 0; i < B->NumofNodes(); ++i)
		context = context + B->Node(i)->Domain();
	Trajectory::Index i = traj.Begin(context);
	ExtractInfo(A);
	Changes();
	BuildRank(ls);
	Create(i);
	large_point = 0;
	getlocation.clear();
	string command2 = string(DOTCOMMAND) + " -Tdot " + name + " > temp.dot";
	const char * commend;
	commend = command2.c_str();
	if (system(commend)<0) return; // report error? (assume >0 not error code)
	GetLocation("temp.dot");
	BuildStateRank(lb);
}

//extracting location of each placement in graph
void Graphic::GetLocation(string name) {
	ifstream obtainlocation;
	double get_distance;
	obtainlocation.open(name.c_str());
	if(!obtainlocation.is_open()) return;
	int count = size;
	vector <double> location;
	string textmatch = "rects=";
	string::iterator it = textmatch.begin();

	while(!obtainlocation.eof()) {
		if( it == textmatch.end()) {
			obtainlocation.get();
			obtainlocation.get();
			obtainlocation >> get_distance;
			obtainlocation.get();
			obtainlocation >> get_distance;

			if(large_point < get_distance) {
				large_point = get_distance;
			}
			location.push_back(get_distance);
			obtainlocation.get();
			obtainlocation >> get_distance;
			obtainlocation.get();
			obtainlocation >> get_distance;
			if(large_point < get_distance) {
				large_point = get_distance;
			}
			location.push_back(get_distance);
			count--;
			getlocation.push_back(location);
			location.clear();
			if(count == 0) {
				break;
			}
			it = textmatch.begin();
		}
		else if(* it == obtainlocation.get())
			it ++;
		else
			it = textmatch.begin();
	}
	obtainlocation.close();
}

//output Markov graph in one arrange line
void Graphic::DrawCTBNRank(const Markov & A, ofstream & ls,
			   string & name, ofstream & lb) {
	ExtractInfo(A);
	Changes();
	BuildRank(ls);
	large_point = 0;
	getlocation.clear();
	string command2 = string(DOTCOMMAND) + " -Tdot " + name + " > temp.dot";
	const char * commend;
	commend = command2.c_str();
	if (system(commend)<0) return; // report error? (assume >0 not error code)
	GetLocation("temp.dot");
	Create(A, start, finish);
	BuildStateRank(lb);
}


void Graphic::DrawCTBN(const Markov & A, ofstream & ls) {
	ExtractInfo(A);
	Changes();
	Build(ls);
}

void Graphic::DrawCTBNMatrix(const Markov & A, ofstream & ls) {
	ExtractInfo(A);
	Changes();
	BuildM(ls);
}

void Graphic ::SetSizeTrajectory(const int & size_x) {
	this -> size_x = size_x;
}

void Graphic :: SetFontTrajectory(const int & change_font, const int & change_font1) {
	if(change_font > 0 && change_font1 > 0) {
		change_point = change_font;
		size_distance = change_font1;
	}
}

void Graphic ::SetFontCTBN(const int & size_font) {
	this -> size_font = size_font;
}

//extract the data from Markov
void Graphic :: ExtractInfo(const Markov & A) {
	graph.clear();
	nodes.clear();
	d = A.GetDynamics();
	const CTBNDyn * B = dynamic_cast<const CTBNDyn *>(d);

	int n = B->NumofNodes();
	size = B->NumofNodes();
	for(int i = 0; i < n; i++) {
		nodes.push_back(dynamic_cast
				<const MarkovDyn *>(B->Node(i)));
	}
	Context x, y;
	vector <int> domain;
	vector <int> condomain;
	std::map<int, GraphicNode>::iterator it;
	for(int i = 0; i < n; i ++) {
		x = nodes[i]->Domain();
		y = nodes[i]->CondDomain();
		domain = x.VarList();
		if(graph.find(domain[0]) != graph.end())
		{
			it = graph.find(domain[0]);
		}
		else
		{
			GraphicNode newnode;
			newnode.Add(domain[0]);
			graph[domain[0]] = newnode;
			it = graph.find(domain[0]);
		}
		if(y.NumVars() > 0) {
			condomain = y.VarList();
			it->second.SetMatching(condomain);
			Instantiation inst(y, 0);
			bool good = true;
			while(good) {
				matrix cim = (*nodes[i])(inst)->Intensity();
				for(unsigned int j = 0; j < condomain.size(); j ++)
				{
					graph[domain[0]].SetPlace(inst.Value(condomain[j]));
				}
				it->second.SetMatrix(cim);
				good = inst.Inc();
			}
		} else {
			it->second.SetBool();
			Instantiation inst(y, 0);
			matrix cim = (*nodes[i])(inst)->Intensity();
			it->second.SetMatrix(cim);
			it->second.SetMatching(domain);
			it->second.SetPlace(0);
		}
	}
}

void Graphic::Changes (vector <string> new_name) {
	change_name = true;
	change_names = new_name;
}

void Graphic::ClearName()
{
	change_names.clear();
	change_name = false;
}

void Graphic::Changes() {
	if(change_name) {
		std::map<int, GraphicNode>::iterator it = graph.begin();
		for(unsigned int i = 0; i < change_names.size()
			&& it != graph.end(); i++, it++) {
			if(it->second.GetExist()) {
				it->second.SetName(change_names[i]);
			}
		}
	}
}

void Graphic ::DrawTrajectory(const Markov & A, ofstream & ls) {
	ExtractInfo(A);
	Create(A, start, finish);
	Changes();
	BuildState(ls);
}

void Graphic::Largest(int big) {
	if (big > largest_tr) {
		largest_tr = big;
	}
}

void Graphic ::DrawTrajectory(Trajectory traj,const Context & context,
							ofstream & ls) {
	this->start_point = traj.TimeBegin();
	this->ending_point = traj.TimeEnd();

	Trajectory::Index i = traj.Begin(context);
	vector <int> domain = context.VarList();
	this ->size = domain.size();
	graph.clear();
	for(int k = 0; k < size; k ++)
	{
		GraphicNode newItem;
		newItem.Add(domain[k]);
		newItem.SetBool();
		graph[domain[k]] = newItem;
	}
	Create(i);
	Changes();
	BuildState(ls);
}

void Graphic ::Create(Trajectory::Index & i) {
	graph_change.clear();
	stringstream ss;
	string name;
	Instantiation oldinst(i.Values());
	vector <int> state;
	for(int n = 0; n < size; n++) {
		state.push_back(oldinst.Value(n));
		Largest(oldinst.Value(n));
	}

	GraphicEvent newItem;
	newItem.SetChanges(size+1, i.Time(), state);
	newItem.SetName("0");
	graph_change.push_back(newItem);
	int count = 1;
	while(!i.Done()) {
		++i;
		GraphicEvent newItem;
		Instantiation newinst(i.Values());
		vector<int> change_id = diffvalues(newinst, oldinst);
		if(change_id.size() > 0 && ending_point != i.Time()) {
			for(unsigned int j = 0; j < change_id.size(); j ++ )
			{
				state[change_id[j]] = newinst.Value(change_id[j]);
				newItem.SetChanges(change_id[j], i.Time(),  state);
				Largest(newinst.Value(change_id[j]));
				ss << count;
				ss >> name;
				ss.clear(stringstream::goodbit);
				newItem.SetName(name);
				graph_change.push_back(newItem);
				count ++;
			}
		}

		oldinst = newinst;
	}


}


void Graphic ::Create(const Markov & A, double & start, double & finish) {
	largest_tr = 0;
	graph_change.clear();
	ending_point = finish;
	stringstream ss;
	d = A.GetDynamics();
	const CTBNDyn * B = dynamic_cast<const CTBNDyn *>(d);
	size = B->NumofNodes();
	Context context = Context();
	for (int i = 0; i < B->NumofNodes(); ++i)
		context = context + B->Node(i)->Domain();
	string name;
	Trajectory tr;
	tr.SetBeginTime(start);
	tr.SetEndTime(finish);
	A.Sample(tr);
	Trajectory::Index i = tr.Begin(context);
	Instantiation oldinst(i.Values());
	vector <int> state;
	for(int n = 0; n < size; n++) {
		state.push_back(oldinst.Value(n));
		Largest(oldinst.Value(n));
	}
	GraphicEvent newItem;
	newItem.SetChanges(size+1, i.Time(), state);
	newItem.SetName("0");
	graph_change.push_back(newItem);
	int count = 1;
	while(!i.Done()) {
		++i;
		GraphicEvent newItem;
		Instantiation newinst(i.Values());
		vector<int> change_id = diffvalues(newinst, oldinst);
		if(change_id.size() > 0 && ending_point != i.Time()) {
			for(unsigned int j = 0; j < change_id.size(); j ++ )
			{
				state[change_id[j]] = newinst.Value(change_id[j]);
				newItem.SetChanges(change_id[j], i.Time(),  state);
				Largest(newinst.Value(change_id[j]));
				ss << count;
				ss >> name;
				ss.clear(stringstream::goodbit);
				newItem.SetName(name);
				graph_change.push_back(newItem);
				count ++;
			}
		}

		oldinst = newinst;
	}

}

void Graphic::BuildMatrix(matrix cim, ofstream& ls) {
	int matrix_size = 0;
	matrix_size = cim.getm();
	for(int m = 0; m < matrix_size; m++) {
		if(m != 0) {
			ls << "{";
		}
		for (int n = 0; n < matrix_size; n++) {
			ls << cim[n][m];
			ArrangeFind(cim);
			Arrange(ls,(int_size[n]- arrangement[m][n]));
			if(n +1 != matrix_size) ls<< " | ";
		}
		ls <<  "}";
		if ((m+1) != matrix_size) ls << "|";
	}
}

int Graphic::Largest(vector<int> length) {
	int large = length[0];
	for(unsigned int i = 1; i < length.size(); i++)
	{
		if(length[i] > large) large = length[i];
	}
	return large;
}

void Graphic::Arrange(ofstream & ls, int space) {
	for(int i = 0; i < space; i++) ls << "\\ ";
}

void Graphic::ArrangeFind(matrix cim) {
	arrangement.clear();
	int_size.clear();
	vector<int> length;
	int matrix_size = cim.getm();
	string size;
	int size_count;
	for(int m = 0; m < matrix_size; m++) {
		length.clear();
		for (int n = 0; n < matrix_size; n++) {
			stringstream ss;
			ss << cim[n][m];
			ss >> size;
			size_count = 0;
			for(unsigned int i = 0; i < size.size(); i ++) {
				if(size[i] == '.') {
					size_count = i +1;
				}
				if(size_count == 0 && i +1 == size.size()) {
					size_count = i +1;
				}
			}
			length.push_back((size.size()-size_count));
			ss.clear(stringstream::goodbit);
		}
		int_size.push_back(Largest(length));
		arrangement.push_back(length);
	}

}

void Graphic::CTBNGraphicIntro(std::ofstream & ls)
{
	ls << "digraph g {" << endl;
	ls << "ordering=out;" << endl;

	ls << "graph [ " << endl << "rankdir = "
	   << "\"" << "LR" << "\"" << endl
	   << "];" << endl;
	ls << "node [" << endl << "fontsize = "
	   << "\"" << size_font << "\"" << "\n" << "shape = "
	   << "\"" << "eclipse" << "\"" << "];" << endl;
	ls << "edge [" << endl << "];" << endl;
}

void Graphic::BuildM(ofstream& ls) {
	CTBNGraphicIntro(ls);
	string name = "";
	string node = "node";
	stringstream ss;
	vector <int> temp;
	vector <int> places;
	int count;
	std::map<int, GraphicNode>::iterator it;
	for(it = graph.begin(); it != graph.end(); it ++) {
		if(it->second.GetExist()) {
		node = node + it->second.GetName();
		ls << "\"" << node << "\"" << "[" << endl;
		ls << "label = " << "\"" << it->second.GetName();
		temp = it->second.GetMatching();
		places = it->second.GetPlace();
		count = 0;

		int matrix_count = 0;
		for(unsigned int x = 0; x < places.size(); x ++) {
			if (x ==0) {
				ls << "| " << "{";
			}
			if(places.size() != 1) {
				if(count == 0) {
					ls << graph.find(temp[count])->second.GetName()
					   << "= " << places[x] << " ";
				}
				else {
					ls << ", "
					   << graph.find(temp[count])->second.GetName()
					   << "= " << places[x] << " ";
				}
			}
			count ++;
			if((unsigned) count >= temp.size()) {
				count = 0;
				if((x+1) != places.size()) {

					ls << "| {";
					matrix cim = it->second.
						     GetMatrix(matrix_count);
					BuildMatrix(cim, ls);
					ls <<"} || {";
					matrix_count++;
				}
			}
		}
		if(places.size() != 1) {
			ls << "| ";
		}
		ls << "{";
		matrix cim = it->second.GetMatrix(matrix_count);
		BuildMatrix(cim, ls);
		ls <<"}";
		ls << "\"" << endl;
		ls << "shape = " << "\"" << "record" << "\"" << endl;
		ls << "];" << endl;
		}
		node = "node";
	}
	string nod1 ="node";
	count = 0;
	temp.clear();


	for(it = graph.begin(); it != graph.end(); it ++) {
		node = node + it->second.GetName();
		temp = it->second.GetMatching();
		for(unsigned int j = 0; j < temp.size(); j++) {
			if(it->second.CheckStrings()) {
				nod1 = nod1 + graph.find(temp[j])->second.GetName();
			}
			else {
				nod1 = nod1 + name;
			}
			ss << count;
			ss >> name;
			if (nod1 != node) {
			ls << "\"" << nod1 << "\"" << "-> "
			   << "\"" << node << "\"" << " [ id = "
			   << name << "];" << endl;
			}
			nod1 = "node";
			ss.clear(stringstream::goodbit);
			count ++;
		}
		node = "node";
	}
	ls << "overlap=false" << endl <<"}" << endl;
}

void Graphic::Build(ofstream & ls) {
	CTBNGraphicIntro(ls);
	string name = "";
	string node = "node";
	stringstream ss;
	vector <int> temp;
	vector <int> places;
	std::map<int, GraphicNode>::iterator it;
	int count;
	for(it = graph.begin(); it != graph.end(); it ++) {
		if(it->second.GetExist()) {
		node = node + it->second.GetName();
		ls << "\"" << node << "\"" << "[" << endl;
		ls << "label = " << "\"" << it->second.GetName();
		temp = it->second.GetMatching();
		places = it->second.GetPlace();
		count = 0;
		ls << "\"" << endl;
		ls << "shape = " << "\"" << "record"
		   << "\"" << endl;
		ls << "];" << endl;
		}
		node = "node";
	}
	string nod1 ="node";
	count = 0;
	temp.clear();
	for(it = graph.begin(); it != graph.end(); it ++) {
		node = node + it->second.GetName();
		temp = it->second.GetMatching();
		for(unsigned int j = 0; j < temp.size(); j++) {
			if(it->second.CheckStrings()) {
				nod1 = nod1 + graph.find(temp[j])->second.GetName();
			}
			else {
				nod1 = nod1 + name;
			}
			ss << count;
			ss >> name;
			if (nod1 != node) {
			ls << "\"" << node << "\"" << " -> "
			   << "\"" << nod1 << "\""
			   << " [ dir=back, id = " << name << "];" << endl;
			}
			nod1 = "node";
			ss.clear(stringstream::goodbit);
			count ++;
		}
		node = "node";
	}
	ls << "overlap=false" << endl <<"}" << endl;
}

void Graphic::BuildRank(ofstream & ls) {

	CTBNGraphicIntro(ls);

	string name = "";
	string node = "node";
	stringstream ss;
	vector <int> temp;
	vector <int> places;
	int count;
	std::map<int, GraphicNode>::iterator it;
	for(it = graph.begin(); it != graph.end(); it ++) {
		if(it->second.GetExist()) {
		node = node + it->second.GetName();
		ls << "\"" << node << "\"" << "[" << endl;
		ls << "label = " << "\"" << it->second.GetName();
		temp = it->second.GetMatching();
		places = it->second.GetPlace();
		count = 0;
		ls << "\"" << endl;
		ls << "shape = " << "\"" << "record"
		   << "\"" << endl;
		ls << "];" << endl;
		}
		node = "node";
	}
	string nod1 ="node";
	count = 0;
	temp.clear();
	ls << "subgraph ctbn { " << endl;
	ls << "rank = same; " << endl;
	for(it = graph.begin(); it != graph.end(); it ++) {
		node = node + it->second.GetName();
		temp = it->second.GetMatching();
		for(unsigned int j = 0; j < temp.size(); j++) {
			if(it->second.CheckStrings()) {
				nod1 = nod1 + graph.find(temp[j])->second.GetName();
			}
			else {
				nod1 = nod1 + name;
			}
			ss << count;
			ss >> name;
			if (nod1 != node) {
			ls << "\"" << node << "\"" << " -> "
			   << "\"" << nod1 << "\""
			   << " [dir = back id = " << name << "];" << endl;
			}
			nod1 = "node";
			ss.clear(stringstream::goodbit);
			count ++;
		}
		node = "node";
	}
	ls << endl << "}";
	ls << endl <<"}" << endl;
}

vector< vector<double> > Graphic::setcolor() {
	int change;
	int count;
	vector< vector<double> > color_set;
	vector<double> color_sequence;

	for(int i = 0; (largest_tr - i) > -1; i ++) {
		change = largest_tr -i;
		count = 0;
		color_sequence.clear();
		while(change > 7) {
			change = change -7;
			count++;
		}
		if(change >= 4) {
			change = change - 4;
			if(count >= 3 && count <= 5) {
				color_sequence.push_back(
							9.00 - (.15 * count));
			}
			else {
				color_sequence.push_back(1.00);
			}
		}
		else {
			if(count >= 3 && count <= 5) {
				color_sequence.push_back(.18 * count);
			}
			else {
				color_sequence.push_back(0.00);
			}
		}
		if(change >= 2) {
			change= change - 2;
			if(count >= 1 && count <= 3) {
				color_sequence.push_back(
							1.00 - (.15 * count));
			}
			else {
				color_sequence.push_back(1.00);
			}
		}
		else {
			if(count >= 1 && count <= 3) {
				color_sequence.push_back(.18 * count);
			}
			else {
				color_sequence.push_back(0.10);
			}
		}
		if(change >= 1) {
			change= change - 1;
			color_sequence.push_back(1.00 - (.15 * count));
		}
		else {
			if(count >= 1) {
				color_sequence.push_back(.18 * count);
			}
			else {
				color_sequence.push_back(0.00);
			}
		}
		color_set.push_back(color_sequence);
	}
	return color_set;
}

int Graphic ::BiggestNumber() {
	int size_text = 0;
	int compare = size;
	while(compare >= 0) {
		compare = compare - 10;
		size_text ++;
	}
	return size_text;
}

int Graphic::BiggestText() {
	int size_text = graph.begin()->second.GetName().length();
	std::map<int, GraphicNode>::iterator it;
	for(it = graph.begin()++; it != graph.end(); it ++) {
		if(it->second.GetName().length() > (unsigned) size_text) {
			size_text = it->second.GetName().length();
		}
	}
	return size_text;
}


void Graphic::BuildState(ofstream & ls) {
	vector <vector <double> > space_place;
	vector <double> space_inn;
	int size_text;
	if(graph.begin()->second.CheckStrings()) {
		size_text = BiggestText();
		size_text = (int) ((size_text * 0.55) + 0.5);
	}
	else {
		size_text= BiggestNumber();
	}

	space_inn.push_back(((change_point * size_text)) + change_point);
	space_inn.push_back(((size+size_distance)
			  * change_point - ( change_point * 2.5 )));


	int state_size = graph_change.size();
	ls << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
	ls << "%%Title: CTBN Graphic" << endl;
	ls << "%%Creator: UCR Research project" << endl;
	time_t timenow;
	struct tm * timeinfo;

	time ( &timenow );
	timeinfo = localtime ( &timenow );
	ls << "%%CreationDate: "<< asctime (timeinfo);
	ls << "%%BoundingBox: 0 0 "
	   << (size_x + change_point)
	   <<" " <<(size+size_distance)* change_point << endl;
	ls << "%%EndComments" << endl;
	ls << "%%Document-Fonts: Times-Roman" << endl;

	ls << "/Times-Roman findfont" << endl;
	ls << change_point << " scalefont" << endl;
	ls << "setfont \n" << endl;
	ls << "1 setlinewidth " << endl;

	int init_box = (change_point * 2);

	int count_line = 0;
	vector< vector<double> > color_set = setcolor();
	for(unsigned int i = 0; i < color_set.size(); i++) {
		if(init_box > size_x) {
			init_box = (change_point * 2);
			count_line++;
		}
		ls << init_box - change_point << " "
		   << ((size+size_distance)* change_point)
		      -(count_line * change_point ) << " "
		   << "newpath moveto" << endl;
		ls << init_box << " "
		   << ((size+size_distance)* change_point)
		      - (count_line * change_point ) << " "
		   << "lineto" << endl;
		ls << init_box  << " "
		   << ((size+size_distance)* change_point)
		      - (count_line * change_point ) - change_point
		   << " " << "lineto" << endl;
		ls << init_box - change_point<< " "
		   << ((size+size_distance)* change_point)
		      - (count_line * change_point ) - change_point
		   << " " << "lineto" << endl;
		ls << "closepath" << endl;
		ls << "gsave" << endl;
		ls << color_set[color_set.size() - i -1][0]
		   << " " << color_set[color_set.size() - i -1][1]
		   << " " << color_set[color_set.size() - i -1][2]
		   << " setrgbcolor fill" << endl;
		ls << "grestore" << endl;
		ls << "stroke" << endl;
		ls << init_box << " "
		   << ((size+size_distance)* change_point)
		      - (count_line * change_point )- change_point
		   << " " << "newpath moveto" << endl;
		ls << "( =" << i << " )" << " show" << endl << endl;
		init_box = init_box + 4 * change_point;
	}


	for(int i = 0; i < size; i ++) {
		space_inn[1] = space_inn[1] - (change_point * i);
		space_place.push_back(space_inn);
		ls << space_inn[0] - (change_point * size_text)
		   << " " << space_inn[1] << " "
		   << "newpath moveto" << endl;
		ls << space_inn[0] << " "
		   << space_inn[1] << " " << "lineto" << endl;
		ls << space_inn[0]  << " "
		   << space_inn[1] - change_point
		   << " " << "lineto" << endl;
		ls << space_inn[0] - (change_point * size_text)
		   << " " << space_inn[1] - change_point
		   << " " << "lineto" << endl;
		ls << "closepath" << endl;
		ls << "stroke" << endl;
		ls << space_inn[0]- change_point* size_text +(change_point / 4)
		   << " " << space_inn[1]
		      - change_point+(change_point / 4)
		   << " " << "newpath moveto" << endl;
		ls << "(" << graph[i].GetName()
		   << ")" << " show" << endl << endl;
		space_inn[1] = ((size+size_distance)* change_point
				       - ( change_point * 2.5 ));
	}
	int measure;
	vector <int> check_value;
	for(int i = 0; i < size; i ++) {
		for(int j = 0; j < state_size; j ++) {
			if (graph_change[j].GetChange() == i) {
				measure = (int) (((size_x -
						space_inn[0])*(graph_change[j].
                                                GetTime() / (ending_point
                                                - start_point))) +
                                                space_inn[0]);
				ls << space_place[i][0] << " "
				   << space_place[i][1] << " "
				   << "newpath moveto" << endl;
				ls << measure << " "
				   << space_place[i][1]
				   << " " << "lineto" << endl;
				ls << measure << " "
				   << space_place[i][1] - change_point
				   << " " << "lineto" << endl;
				ls << space_place[i][0]<< " "
				   << space_place[i][1] - change_point
				   << " " << "lineto" << endl;
				ls << "closepath" << endl;
				ls << "gsave" << endl;
				check_value = graph_change[j-1].GetState();
				if(check_value[i] != -1) {
					ls << color_set[color_set.size() - 1
					   - check_value[i]][0] << " "
					   << color_set[color_set.size() - 1
					   - check_value[i]][1] << " "
					   << color_set[color_set.size() - 1
					   - check_value[i]][2]
					   << " setrgbcolor fill" << endl;
				}
				ls << "grestore" << endl;
				ls << "stroke" << endl << endl;
				ls << "1 setlinewidth " << endl;
				space_place[i][0] = measure;
			}

		}
		ls << space_place[i][0]
		   << " " << space_place[i][1]
		   << " " << "newpath moveto" << endl;
		ls << size_x << " " << space_place[i][1]
		   << " " << "lineto" << endl;
		ls << size_x << " " << space_place[i][1] - change_point
		   << " " << "lineto" << endl;
		ls << space_place[i][0]<< " "
		   << space_place[i][1] - change_point
		   << " " << "lineto" << endl;
		ls << "closepath" << endl;
		ls << "gsave" << endl;
		check_value =graph_change[state_size-1].GetState();

		if(check_value[i] != -1) {
			ls << color_set[color_set.size()
			      - 1 - check_value[i]][0] << " "
			   << color_set[color_set.size()
			      - 1 - check_value[i]][1] << " "
			   << color_set[color_set.size()
			      - 1 - check_value[i]][2]
			   << " setrgbcolor fill" << endl;
		}
		ls << "grestore" << endl;
		ls << "stroke" << endl << endl;
	}
	ls << "/Times-Roman findfont" << endl;
	ls << size_distance << " scalefont" << endl;
	ls << "setfont \n" << endl;
	int first_row = (int) (space_place[size-1][1]
	                - (change_point + size_distance));
	int next_row = size_distance;
	vector <double> firstline_row;
	firstline_row.push_back(0);
	unsigned int line = 0;
	int present_column;
	int rejust;
	double justed;
	for(int i = 0; i < state_size; i ++) {
		if(i > 0 && graph_change[i].GetTime() != graph_change[i-1].GetTime())
		{
			measure = (int) ((size_x - space_inn[0])
				  * (graph_change[i].GetTime()
				  / (ending_point - start_point))
				  + space_inn[0]);

			line = 0;
			while(firstline_row[line]  > measure - (size_distance * 2.5)) {
				line ++;
				if(line >= firstline_row.size())
					firstline_row.push_back(0);
			}

			present_column = first_row -(next_row * line);
			ls << measure << " "
			<< present_column << " "
			<< "newpath moveto" << endl;
			firstline_row[line] = measure;
			
			rejust = (int) (graph_change[i].GetTime() * 1000);
			justed = rejust;
			justed = justed / 1000;
			ls << "(" << justed << ")" << " show" << endl << endl;
		}
	}


	ls << "showpage"<< endl;
	ls << "end";
}

void Graphic::BuildStateRank(ofstream & ls) {
	vector <double> space_inn;

	int size_text;
	if(graph[0].CheckStrings()) {
		size_text = BiggestText();
		size_text = (int) ((size_text * 0.55) + 0.5);
	}
	else {
		size_text= BiggestNumber();
	}

	int state_size = graph_change.size();
	ls << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
	ls << "%%Title: CTBN Graphic" << endl;
	ls << "%%Creator: UCR Research project" << endl;
	time_t timenow;
	struct tm * timeinfo;

	time ( &timenow );
	timeinfo = localtime ( &timenow );
	ls << "%%CreationDate: "<< asctime (timeinfo);
	ls << "%%BoundingBox: 0 0 "
	   << (size_x + change_point) <<" "
	   << large_point + (9 * change_point) << endl;
	ls << "%%EndComments" << endl;
	ls << "%%Document-Fonts: Times-Roman" << endl;

	ls << "/Times-Roman findfont" << endl;
	ls << change_point << " scalefont" << endl;
	ls << "setfont \n" << endl;
	ls << "1 setlinewidth " << endl;

	int init_box = (change_point* 2);
	double position_y = large_point + (9 * change_point);
	int count_line = 0;
	vector< vector<double> > color_set = setcolor();
	for(unsigned int i = 0; i < color_set.size(); i++) {
		if(init_box > size_x) {
			init_box = (change_point * 2);
			count_line++;
		}
		ls << init_box - change_point<< " "
		   << position_y -(count_line * change_point )
		   << " " << "newpath moveto" << endl;
		ls << init_box << " "
		   << position_y -(count_line * change_point )
		   << " " << "lineto" << endl;
		ls << init_box  << " "
		   << position_y -
		      (count_line * change_point ) - change_point
		   << " " << "lineto" << endl;
		ls << init_box - change_point << " "
		   << position_y -(count_line * change_point )
		      - change_point
		   << " " << "lineto" << endl;
		ls << "closepath" << endl;
		ls << "gsave" << endl;
		ls << color_set[color_set.size() - i -1][0] << " "
		   << color_set[color_set.size() - i -1][1] << " "
		   << color_set[color_set.size() - i -1][2]
		   << " setrgbcolor fill" << endl;
		ls << "grestore" << endl;
		ls << "stroke" << endl;
		ls << init_box << " "
		   << position_y -(count_line * change_point ) - change_point
		   << " " << "newpath moveto" << endl;
		ls << "( =" << i << " )" << " show" << endl << endl;
		init_box = init_box + 4 * change_point;
	}

	init_box = ((change_point * size_text)) + change_point;
	for(int i = 0; i < size; i ++) {
		ls << init_box - (change_point * size_text) << " "
		   << getlocation[i][1] + (5 * change_point) << " "
		   << "newpath moveto" << endl;
		ls << init_box  << " "
		   << getlocation[i][1] + (5 * change_point)
		   << " " << "lineto" << endl;
		ls << init_box  << " "
		   << getlocation[i][0] + (5 * change_point) << " "
		   << "lineto" << endl;
		ls << init_box - (change_point * size_text) << " "
		   << getlocation[i][0] + (5 * change_point) << " "
		   << "lineto" << endl;
		ls << "closepath" << endl;
		ls << "stroke" << endl;
		ls << init_box - change_point * size_text +(change_point / 4)
		   << " " << getlocation[i][0] + (5 * change_point)
		      + ((getlocation[i][1] - getlocation[i][0]) / 4)
		   << " " << "newpath moveto" << endl;
		ls << "(" << graph[i].GetName() << ")"
		   << " show" << endl << endl;
		space_inn.push_back(init_box);
	}
	int measure;
	vector <int> check_value;
	for(int i = 0; i < size; i ++) {
		for(int j = 0; j < state_size; j ++) {
			if (graph_change[j].GetChange() == i) {
				measure = (int) ((size_x -
						init_box)*(graph_change[j].
                                                GetTime() / (ending_point
                                                - start_point)) +
                                                init_box);
				ls << space_inn[i] << " "
				   << getlocation[i][1] + (5 * change_point)
				   << " " << "newpath moveto" << endl;
				ls << measure << " "
				   << getlocation[i][1] + (5 * change_point)
				   << " "
				   << "lineto" << endl;
				ls << measure << " "
				   << getlocation[i][0] + (5 * change_point)
				   << " "
				   << "lineto" << endl;
				ls << space_inn[i] << " "
				   << getlocation[i][0] + (5 * change_point)
				   << " "
				   << "lineto" << endl;
				ls << "closepath" << endl;
				ls << "gsave" << endl;
				check_value =graph_change[j-1].GetState();
				if(check_value[i] != -1) {
					ls << color_set[color_set.size()
					      - 1 - check_value[i]][0] << " "
					   << color_set[color_set.size()
					      - 1 - check_value[i]][1] << " "
					   << color_set[color_set.size()
					      - 1 - check_value[i]][2]
					   << " setrgbcolor fill" << endl;
				}
				else
				{
					ls << 1 << " "
					   << 1 << " "
					   << 1
					   << " setrgbcolor fill" << endl;
				}
				ls << "grestore" << endl;
				ls << "stroke" << endl << endl;
				ls << "1 setlinewidth " << endl;
				space_inn[i] = measure;
			}

		}
		ls << space_inn[i]<< " "
		   << getlocation[i][1] + (5 * change_point) << " "
		   << "newpath moveto" << endl;
		ls << size_x << " "
		   << getlocation[i][1] + (5 * change_point) << " "
		   << "lineto" << endl;
		ls << size_x << " "
		   << getlocation[i][0] + (5 * change_point) << " "
		   << "lineto" << endl;
		ls << space_inn[i]<< " "
		   << getlocation[i][0] + (5 * change_point) << " "
		   << "lineto" << endl;
		ls << "closepath" << endl;
		ls << "gsave" << endl;
		check_value =graph_change[state_size-1].GetState();

		if(check_value[i] != -1) {
			ls << color_set[color_set.size()
			      - 1 - check_value[i]][0] << " "
			   << color_set[color_set.size()
			      - 1 - check_value[i]][1] << " "
			   << color_set[color_set.size()
			      - 1 - check_value[i]][2]
			   << " setrgbcolor fill" << endl;
		}
		ls << "grestore" << endl;
		ls << "stroke" << endl << endl;
	}
	ls << "/Times-Roman findfont" << endl;
	ls << size_distance << " scalefont" << endl;
	ls << "setfont \n" << endl;
	int first_row = (int) getlocation[0][0]
					 + size_distance;
	int next_row = size_distance;
	vector <double> firstline_row;
	firstline_row.push_back(0);
	unsigned int line = 0;
	int present_column;
	int rejust;
	double justed;
	for(int i = 0; i < state_size; i ++) {
		if(i > 0 && graph_change[i].GetTime() != graph_change[i-1].GetTime())
		{
			measure = (int) ((size_x - init_box)
				  * (graph_change[i].GetTime()
				  / (ending_point - start_point))
				  + init_box);

			line = 0;
			while(firstline_row[line]  > measure - (size_distance * 2.5)) {
				line ++;
				if(line >= firstline_row.size())
					firstline_row.push_back(0);
			}

			present_column = first_row -(next_row * line);
			ls << measure << " "
			<< present_column << " "
			<< "newpath moveto" << endl;
			firstline_row[line] = measure;

			rejust = (int) (graph_change[i].GetTime() * 1000);
			justed = rejust;
			justed = justed / 1000;
			ls << "(" << justed << ")" << " show" << endl << endl;
		}
	}


	ls << "showpage"<< endl;
	ls << "end";
}

} // end of ctbn namespace
