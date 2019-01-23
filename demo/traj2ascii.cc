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
// this file reads in a trajectory from stdin and "draws" an
// ASCII representation of it on stdout
#include "trajectory.h"
#include <iostream>
#include "params.h"
#include "streamextra.h"

using namespace std;
using namespace ctbn;

int main(int argc, char **argv)
{
	InitParams(argc, argv);

	addnan(cout);
	addnan(cin);
 
	Trajectory tr;
	tr.Load(cin);
	tr.Draw(cout);
	cout << endl;
	return 0;
}
