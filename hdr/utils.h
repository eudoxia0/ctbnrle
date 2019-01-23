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
#ifndef CTBNRLE_UTILS_H
#define CTBNRLE_UTILS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <queue>


// A set of useful utility functions

namespace ctbn {

class Trajectory;
class Context;

bool IsTrajSame(const Trajectory &tr1, const Trajectory &tr2, const Context &c);
bool IsTrajValid(const Trajectory &tr, const Context &c);
void RemoveInformation(Trajectory &tr, const Context &c, int var, double start, double end);
void RemoveInformation(Trajectory &tr, const Context &c, int nvar, int nit, double frac);
void RemoveNodesInformation(Trajectory &tr, const Context &c, int nvar, int nit, double frac);

void SelectNodesInformation(Trajectory &tr, const Context &c, int nvar, double thresh, double frac);
double SelectTime(std::priority_queue<double> &select, double ratio, double begint, double endt, double thresh );

double getcputime(void);
void Normalize(std::vector<double> &logw);
double mean(const std::vector<double> &v);
double relative_bias(const std::vector<double> &v, double vtrue);
double relative_std(const std::vector<double> &v, double vtrue);
bool LoadNetworkTraj(std::istream &in, std::vector<Trajectory> &results);
bool SplitPortTraj(const std::vector<Trajectory> &whole, std::vector<std::vector<Trajectory> > &splitted);
std::vector<Trajectory> PartitionTraj(const Trajectory &whole);

Trajectory LoadTraj(std::ifstream &in);

} // end of ctbn namespace

#endif
