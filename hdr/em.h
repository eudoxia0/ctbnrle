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
#ifndef CTBNRLE_EM_H
#define CTBNRLE_EM_H

#include "inference.h"
#include "process.h"
#include "trajectory.h"
#include "structure.h"

// All of the pointers in these functions are owned by the caller

namespace ctbn {

//  Perform Expectation-Maximization with the data data
//  on the Process dist with inference method inf
void EM(const std::vector<Trajectory> &data, Process *dist, Inference *inf);

// Perform Structural Expectation-Maximization
// In this case, dist may need to be changed in which case the original
// one will be "deleted" (freed) and replaced with a new one
void SEM(const std::vector<Trajectory> &data, 
		Process *&dist, 
		Inference *inf);

} // end of ctbn namespace

#endif
