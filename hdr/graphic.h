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
#ifndef CTBNRLE_GRAPHIC_H
#define CTBNRLE_GRAPHIC_H

#include "markov.h"
#include "dynamics.h"
#include "context.h"
#include "ctbndyn.h"
#include "matrix.h"
#include "markovdyn.h"
#include "trajectory.h"

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace ctbn {

/*
 * this class is devoted to drawing a ctbn graph and a trajectory
 * given the start time and ending time it will generate the trajectory graph
 * or output the Trajectory::Index graph.
 *
 *
 * In int main you need to use InitParams(argc, argv); before anything else
 * params.h
 * to change the window x coordinate size use the "traj_win_size"
 * to change the font on trajectory graph use the "traj_font"
 * to change the font on ctbn graph use the "ctbn_font"
 * to change the font on trajectory graph ruler use the "traj_ruler"
 */

class Graphic {
  public:
	Graphic();
	//output CTBN graph without the matrices
	void DrawCTBN(const Markov &, std::ofstream &);
	//output trajectory graph
	void DrawTrajectory(Trajectory,const Context &,
			            std::ofstream &);
	//output trajectory graph by generate new trajectory
	void DrawTrajectory(const Markov &, std::ofstream &);
	//output CTBN graph with the matrices
	void DrawCTBNMatrix (const Markov &, std::ofstream &);
	//output both graph (i.e. it require name of ctbn graph as string)
	void DrawCTBNRank(const Markov &,
			           std::ofstream &, std::string &,
			           std::ofstream &);
	//output both graph (i.e. it require name of ctbn graph as string)
	void DrawBoth(const Markov &, Trajectory &, 
		      std::ofstream &,std::string &, std::ofstream &);
	//change size for name, and time of change
	void SetFontTrajectory(const int & change_font = 15 , 
			       const int & change_font1 = 5 );
	//set size of page of trajectory
	void SetSizeTrajectory(const int &);
	//set font size on CTBN graph
	void SetFontCTBN(const int &);
	//change name of the variable.
	void Changes(std::vector <std::string>);
	//clearing the changed name
	void ClearName();

  private:
	//class use for save up the trajectory graphic
	class GraphicEvent {
	  public:
		GraphicEvent();
		void SetChanges(int, double, const std::vector<int> &);
		void SetName(std::string);
		std::string GetName() const;
		int GetChange() const;
		double GetTime() const;
		std::vector <int> GetState() const;

	  private:
		std::vector <int> state;
		std::string name;
		int change_node;
		double time;
	};

	//class use for save up the ctbn graphic
	class GraphicNode {
	  public:
		GraphicNode();
		void Add(const int &);
		void Add(const std::string &);
		std::string GetName();
		int GetSize();
		void SetName(const std::string &);
		void SetMatrix(const matrix &, int &, int &);
		void SetPlace(int);
		void SetMatrix(const matrix &);
		void SetBool();
		void SetMatching(std::vector <int>);
		void SetMatching(std::vector <std::string>);
		std::vector <int> GetPlace();
		std::vector <int> GetMatching();
		bool GetExist();
		matrix GetMatrix(const int &);
		bool CheckStrings();

	  private:
		bool exist;
		bool check_string;
		std::string a;
		int b;
		std::vector <int> connect;
		std::vector <matrix> dyn;
		std::vector <int> matching;
		std::vector <int> place;
	};

	//generate CTBN process with time limit for trajectory
	void Create(const Markov &, double &, double &);
	//create the vector to generate CTBN graph
	void ExtractInfo(const Markov &);
	void GetLocation(std::string);
	void Create(Trajectory::Index &);
	//function created for attempt align each number together
	void Arrange(std::ofstream& ls, int);
	void ArrangeFind(matrix);
	//get largest value
	int Largest(std::vector<int>);
	//get larger value
	void Largest(int);
	//get the different color for each number
	std::vector< std::vector<double> > setcolor();
	//get size of text
	int BiggestText();
	//get size of number
	int BiggestNumber();
	//change name of the variable.
	void Changes();
	//build graph of CTBN with matrixes
	void BuildM(std::ofstream & ls);
	//build graph of CTBN without matrixes
	void Build(std::ofstream & ls);
	//build graph of CTBN without matrixes with rank
	void BuildRank(std::ofstream & ls);
	//build trajectory
	void BuildState(std::ofstream & ls);
	//part of the buildm which generate matrix on the graph
	void BuildMatrix(matrix, std::ofstream& ls);
	void BuildStateRank(std::ofstream & ls);

	void CTBNGraphicIntro(std::ofstream & ls);

	std::map <int, GraphicNode> graph;
	std::vector <GraphicEvent> graph_change;
	std::vector <std::string> change_names;
	int size;
	//size of x coordinate the trajectory visualization
	int size_x;
	//size of the font on trajectory
	int change_point;
	//size in trajectory change mark
	int size_distance;
	//start time of the trajectory
	double start_point;
	//ending time of the trajectory
	double ending_point;
	double large_point;
	int largest_tr;
	//size of the font on ctbn graph
	int size_font;
	double start;
	bool change_name;
	double finish;
	std::vector< std::vector<double> > getlocation;
	std::vector < std::vector <int> > arrangement;
	std::vector < std::vector <int> > troject;
	std::vector <int> int_size;
	//They draw the pointers D and nodes form the CTBN which are own by this class
	const Dynamics *d;
	std::vector <const MarkovDyn *> nodes;
};


} // end of ctbn namespace

#endif /*GRAPHIC_H_*/

