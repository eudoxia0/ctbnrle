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
#ifndef CTBNRLE_STREAMEXTRA_H
#define CTBNRLE_STREAMEXTRA_H

#include <math.h>
#include <string.h>
#include <algorithm>
#include "defines.h"

/* This header file does two separated (but related) things:
   1. It defines a way to save and load via streams (ie serialize)
      maps and vectors from STL.  Because of scoping issues it
	 has to be defined in the std namespace (which is unfortunate)

   2. It defines a function addnan which can be called on any
      input or output character stream.  It makes sure that stream
	 can read or write "NaN" and infinite values properly.
*/

#include <map>
#include <vector>
#include <iostream>
#include <iterator>
#include <string>
#include <typeinfo>

// first declaration of some helper methods to push maps and vectors
// automatically through operator<< and operator>>
// It would be *great* if this could be taken out of the std namespace!
// Unfortunately, with current C++ compilers (March 2010), that is not
// possible.
namespace std {
	template<typename T, typename S>
	ostream &operator<<(ostream &os, const pair<T,S> &p) {
		return os << p.first << os.fill() << p.second;
	}
	template<typename T, typename S>
	istream &operator>>(istream &is, pair<T,S> &p) {
		return is >> p.first >> p.second;
	}
	
	// basically the same as istream_iterator, but with a limit
	// on the number of items that can be read.
	template<typename T, typename CharT = char,
		typename Traits = char_traits<CharT>,
		typename Dist = ptrdiff_t>
	class istream_lim_iterator
	: public iterator<input_iterator_tag,T,Dist,const T*, const T&> {
		private:
			int togo;
			basic_istream<CharT,Traits> *s;
			T t;
		public:
			istream_lim_iterator() {
				togo = -1;
				s = NULL;
			}
			istream_lim_iterator(basic_istream<CharT,Traits> &is,
						 int limit)
					: togo(limit), s(&is) {
				read();
			}

			istream_lim_iterator& operator++() {
				read();
				return *this;
			}
			const istream_lim_iterator operator++(int) {
				istream_lim_iterator tmp = *this;
				(*this)++;
				return tmp;
			}
			const T operator*() const {
				return t;
			}
			const T* operator->() const {
				return &(operator*());
			}
			bool isequal(const istream_lim_iterator &ili) const {
				return (togo==-1 && ili.togo==-1) || 
						(togo!=-1 && ili.togo!=-1);
			}
		private:
			void read() {
				if (togo>=0) {
					if (togo>0 && s && s->good()) {
						*s >> t;
						togo--;
						if (!s || !s->good()) 
							togo = -1;
					} else togo=-1;
				}
			}
	};

	template<typename T, typename CharT, typename Traits, typename Dist>
	inline bool operator==(
			const istream_lim_iterator<T,CharT,Traits,Dist> &x,
			const istream_lim_iterator<T,CharT,Traits,Dist> &y) {
		return x.isequal(y);
	}

	template<typename T, typename CharT, typename Traits, typename Dist>
	inline bool operator!=(
			const istream_lim_iterator<T,CharT,Traits,Dist> &x,
			const istream_lim_iterator<T,CharT,Traits,Dist> &y) {
		return !x.isequal(y);
	}

	// And now we have the definitions of the shift operators:
	template<typename K, typename T, typename C, typename A>
	ostream &operator<<(ostream &os, const map<K,T,C,A> &m) {
		os << m.size() << os.fill();
		copy(m.begin(),m.end(),ostream_iterator<pair<K,T> >(os," "));
		return os;
	}

	template<typename K, typename T, typename C, typename A>
	istream &operator>>(istream &is, map<K,T,C,A> &m) {
		m.clear();
		int s;
		is >> s;
		m.insert(istream_lim_iterator<pair<K,T> >(is,s),
				istream_lim_iterator<pair<K,T> >());
		return is;
	}

	template<typename T, typename A>
	ostream &operator<<(ostream &os, const vector<T,A> &v) {
		os << v.size() << os.fill();
		copy(v.begin(),v.end(),ostream_iterator<T>(os," "));
		return os;
	}

	template<typename T, typename A>
	istream &operator>>(istream &is, vector<T,A> &v) {
		int s;
		is >> s;
		v.assign(istream_lim_iterator<T>(is,s),
				istream_lim_iterator<T>());
		return is;
	}

}

// Here is some code to define a new locale so that floating point
// numbers like -inf and Nan can be read and written.

// You can forget about all of it.  What you need to do is call
// this function (addnan) on any input or output stream before
// using it.  It will make sure that stream reads or writes "inf"
// and "NaN" values properly
void addnan(std::basic_ios<char> &str);


// This should work, but someone didn't add num_put_byname
// (there is time_put_byname however!)
//template <class charT, class OutputIterator>
//class nannum_put : public std::num_put_byname<charT> {
// so instead we do...
template <class charT, charT *infstr, charT *nanstr, 
		class OutputIterator = std::ostreambuf_iterator<charT> >
class nannum_put : public std::num_put<charT,OutputIterator> {
	public:
		virtual ~nannum_put() throw() {};
		explicit nannum_put(const char *locname)
			//: std::num_put_byname<charT>(locname) {}
			: std::num_put<charT,OutputIterator>() {}
	protected:
		virtual OutputIterator do_put(OutputIterator itt, 
						std::ios_base &iob, 
						charT fill, 
						double d) const {
			if (!isfinite(d)) { 
				itt = writeit<double>(itt,d);
			} else {
				//itt = std::num_put_byname<charT>::
				//	do_put(itt,iob,fill,d);
				itt = std::num_put<charT,OutputIterator>::
						do_put(itt,iob,fill,d);
			}
			return itt;
		}
		virtual OutputIterator do_put(OutputIterator itt, 
				std::ios_base &iob, charT fill, 
				long double d) const {
			if (!isfinite(d)) { 
				itt = writeit<long double>(itt,d);
			} else {
				//itt = std::num_put_byname<charT>::
				//		do_put(itt,iob,fill,d);
				itt = std::num_put<charT,OutputIterator>::
					do_put(itt,iob,fill,d);
			}
			return itt;
		}
		virtual OutputIterator do_put(OutputIterator itt, 
						std::ios_base &iob, 
						charT fill, 
						float f) const {
			if (!isfinite(f)) { 
				itt = writeit<float>(itt,f);
			} else {
				//itt = std::num_put_byname<charT>::
				//		do_put(itt,iob,fill,f);
				itt = std::num_put<charT,OutputIterator>::
					do_put(itt,iob,fill,f);
			}
			return itt;
		}
	private:
		template<typename FP>
			static OutputIterator writeit(OutputIterator itt,
							FP f) {
				static const char *minus = "-";
				// commented-out code works for g++'s STL, but 
				// not for others.  Uncommented-out code is fine
				// more generally
				if (isinf(f)) {
					if (f<0) 
						//itt = std::__write(itt,"-",1);
						std::copy(minus,minus+1,itt);
					std::copy(infstr,
						infstr
						 +std::char_traits<charT>::length(infstr),
						itt);
					/*
					itt = std::__write(itt,
							infstr,
							std::
							char_traits<charT>::
							length(infstr));
							*/
				} else
					std::copy(nanstr,
						nanstr
						 +std::char_traits<charT>::length(nanstr),
						itt);
					/*
					itt = std::__write(itt,
							nanstr,
							std::
							char_traits<charT>::
							length(nanstr));
							*/
				return itt;
			}
};

// This should work, but someone didn't add num_get_byname
// (there is time_get_byname however!)
//template <class charT, class InputIterator>
//class nannum_get : public std::num_put_byname<charT> {
// so instead we do...
template <class charT, charT *infstr, charT *nanstr, 
		class InputIterator = std::istreambuf_iterator<charT> >
class nannum_get : public std::num_get<charT,InputIterator> {
	public:
		virtual ~nannum_get() throw() {};
		typedef std::basic_string<charT> strtype;

		explicit nannum_get(const char *locname)
			//: std::num_get_byname<charT>(locname) {}
			: std::num_get<charT,InputIterator>() {}
	protected:
		virtual InputIterator do_get(InputIterator beg, 
						InputIterator end, 
						std::ios_base &iob, 
						std::ios_base::iostate &iost,
						 double &d) const {
			bool negate = false;
			if (*beg == '-') { negate=true; ++beg; }
			if (*beg == '+') { negate=false; ++beg; }
			if (*beg==*infstr) {
				if (checkstr(beg,end,infstr,iost))
					d = negate ? -INFINITY : INFINITY;
				return beg;
			} else if (*beg==*nanstr) {
				if (checkstr(beg,end,nanstr,iost))
					d = nan("");
				return beg;
			} else {
				InputIterator ret = 
					std::num_get<charT,InputIterator>::
						do_get(beg,end,iob,iost,d);
				if (negate) d = -d;
				return ret;
			}
		}
		virtual InputIterator do_get(InputIterator beg, 
						InputIterator end, 
						std::ios_base &iob, 
						std::ios_base::iostate &iost,
						 long double &d) const {
			bool negate = false;
			if (*beg == '-') { negate = true; ++beg; }
			if (*beg == '+') { negate = false; ++beg; }
			if (*beg==*infstr) {
				if (checkstr(beg,end,infstr,iost))
					d = negate ? -INFINITY : INFINITY;
				return beg;
			} else if (*beg==*nanstr) {
				if (checkstr(beg,end,nanstr,iost))
					d = nanl("");
				return beg;
			} else {
				InputIterator ret = 
					std::num_get<charT,InputIterator>::
						do_get(beg,end,iob,iost,d);
				if (negate) d = -d;
				return ret;
			}
		}
		virtual InputIterator do_get(InputIterator beg, 
						InputIterator end, 
						std::ios_base &iob, 
						std::ios_base::iostate &iost,
						float &f) const {
			bool negate = false;
			if (*beg == '-') { negate = true; ++beg; }
			if (*beg == '+') { negate = false; ++beg; }
			if (*beg==*infstr) {
				if (checkstr(beg,end,infstr,iost))
					f = negate ? -INFINITY : INFINITY;
				return beg; } else if (*beg==*nanstr) {
				if (checkstr(beg,end,nanstr,iost))
					f = nanf("");
				return beg;
			} else {
				InputIterator ret = 
					std::num_get<charT,InputIterator>::
						do_get(beg,end,iob,iost,f);
				if (negate) f = -f;
				return ret;
			}
		}
	private:
		static bool checkstr(InputIterator &beg, 
					InputIterator end, 
					charT *str, 
					std::ios_base::iostate &iost) {
			while(*str && beg!=end) {
				if (*beg==*str) { ++str; ++beg; }
				else {
					iost |= std::ios_base::failbit;
					return false;
				}
			}
			return *str==0;
		}
};
// End of stuff for addnan

#endif

