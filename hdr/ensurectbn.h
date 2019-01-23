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
#ifndef ENSURECTBN_H
#define ENSURECTBN_H

#include "ctbn.h"
#include "bn.h"
#include "markov.h"
#include "multirv.h"
#include "ctbn.h"

namespace ctbn {
ENSURECLASS(CTBN)
ENSURECLASS(BN)
ENSURECLASS(Markov)
ENSURECLASS(CTBNDyn)
ENSURECLASS(MultiRV)
ENSURECLASS(MarkovDyn)
ENSURECLASS(MarkovSimple)
ENSURECLASS(MarkovSimpleToggle)
ENSURECLASS1(RVCondSimpleComp,MultiZSimple)
}

#endif
