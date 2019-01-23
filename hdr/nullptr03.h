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
#ifndef CTBNRLE_NULLPTR03_H
#define CTBNRLE_NULLPTR03_H

// A la Scott Meyers...
namespace ctbn {

	const class nullptr_t {
	public:
		template <typename T> inline operator T*() const { return 0; }
		template <class C, typename T> inline operator T C::*() const { return 0; }
	private:
		void operator&() const;
	} nullptr03 = {};
}

#endif
