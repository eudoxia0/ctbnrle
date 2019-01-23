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
#include "streamextra.h"

#include <cstring>


namespace {
	char infstring[] = "inf";
	char nanstring[] = "NaN";

	// the "C" means classic, and it should really be based
	// on the current locale, but currently the "C" is ignored anyway!
	nannum_put<char,infstring,nanstring> *putfacet = 
			new nannum_put<char,infstring,nanstring>("C");
	nannum_get<char,infstring,nanstring> *getfacet = 
			new nannum_get<char,infstring,nanstring>("C");
}

void addnan(std::basic_ios<char> &str) {
	std::locale withread(str.getloc(),putfacet);
	std::locale withboth(withread,getfacet);
	str.imbue(withboth);
}

