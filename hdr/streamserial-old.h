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

#ifndef CTBNRLE_STREAMSERIAL_OLD_H
#define CTBNRLE_STREAMSERIAL_OLD_H

// NOTE NOTE NOTE NOTE NOTE
/* This header file is old.  It will be removed in subsequent releases
 * it is here to provide a parallel method for loading and saving so
 * that files can be converted from the old format to the new
 */

/* This header file does two things:
   1. It defines a SteamObj class.  This class can save and load
      any of its descendents to or from a stream without knowing the
	 type at compile time.  That is, you can load and save pointers
	 without having the object "sliced" to the type of the pointer.

	 Using it is slightly complex.  Basically, you need to define
	 a class that inherits (directly or indirectly) from StreamObj.
	 This class needs two things:
          i. a constructor with parameter "std::istream &", and
	     ii. a virtual method SaveOld that returns void and takes in
		      a "std::ostream &os"
      The former serves as a "load" method and the later as a save
	 method.

	 In the class declaration (on the first line) you need to
	 add the macro "SOBJCLASSDECL(xxx)" where xxx is the name of
	 the class.  Similarly, in the definition (.cc) file, you
	 need to add the macro "SOBJCLASSDEF(xxx)"

	 Below is a simple example.  The first constructor is not
	 required by the class.

	-------------in the .h file:

	class A : public StreamObj {
		SOBJCLASSDECL(A)
		
		public:
			A(int xx, int yy);
			A(std::istream &is);
			virtual void SaveOld(std::ostream &os) const;

		private:
			int x,y;
	};

	--------------in the .cc file:

	SOBJCLASSDEF(A)

	A::A(int xx, int yy) {
		x = xx; y = yy;
	}

	A::A(std::istream &is) {
		is >> x >> y;
	}

	void A::SaveOld(std::ostream &os) const {
		os << x << os.fill() << y;
	}

	---------------

	 To use the capabilities, simply save pointers to objects
	 by calling "SaveOldPtr" and load them using StreamObj::LoadOldPtr.

	 For example:

	 A *aptr = new A(2,3);
	 aptr->SaveOldPtr(cout);  // save object pointed to by aptr to cout
	 	// Note that if aptr pointed to a derived class of A,
		// provided it also had the form above, it would be saved
		// in its entirety (not just the "class A" portion).

	 StreamObj *objptr = StreamObj::LoadOldPtr(cin); // load obj from cin
	 	// Note that the resulting object will be of the 
		// "correct" type and not necessarily a StreamObj object.

	 2. it defines the marco ENSURECLASSOLD(...) that makes sure
	    that the argument (a class name) is loaded in.  This is necessary
	    if your code never makes explicit reference to the class, but
	    expects it to be able to be loaded dynamically from a file
*/

#include "defines.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <typeinfo>
#include <vector>

namespace ctbn {

// This is the StreamObj.  The macros below take care of all of the
// internals for inherited classes.  Do NOT overload any of these
// methods except SaveOld!  Just define SaveOld and a constructor (from 
// an input stream) for your new class and ignore the rest.
//
// The top of this file has instructions and an example
class StreamObj {
	public:
		StreamObj();
		virtual ~StreamObj();
		
		// This is the one that saves the info and should
		// be called on all ptrs
		void SaveOldPtr(std::ostream &os) const;
		// This is standard (and the one you should overwrite)
		virtual void SaveOld(std::ostream &os) const;

		// this one does all of the "RTTI-ish" stuff
		// this is the one you want to call to "load a ptr"
		static StreamObj *LoadOldPtr(std::istream &is);

		// this is the straight-forward one
		static StreamObj *Create(std::istream &is);

		typedef StreamObj* (*CreatorFn)(std::istream &);
		static const char *addID(const char *classname, CreatorFn fn);
		static const char *addID(const char *classname, 
					 	const char *classname2, 
						CreatorFn fn);
		//static char *killSpaces(const char *classname);
		static const char *IDname();
	protected:
		virtual const char *myID() const;

	private:
		static std::map<std::string, CreatorFn> &loadtable();
		static const char *_forcename;
};

#define SOBJCLASSDECL(classname) \
	private: \
		static const char *_forcename; \
	protected: \
		virtual const char *myID() const; \
	public: \
		static const char *IDname(); \
		static StreamObj *Create(std::istream &is);

#define SOBJCLASSDEF(classname) \
	const char *classname::IDname() { \
		static const char *ret = StreamObj::addID(#classname,classname::Create); \
		return ret; \
	} \
	StreamObj *classname::Create(std::istream &is) { \
		return new classname(is); \
	} \
	const char *classname::myID() const { return IDname(); } \
	const char *classname::_forcename = classname::IDname();

//StreamObj::addID(#classname,tclass::IDname,classname<tclass>::Create); 
// This isn't that portable as typeid( ).name() isn't gauranteed to
// return anything useful and it certainly doesn't return the same thing
// across compilers.
// Another fix doesn't really seem possible just now.
// (asking class tclass for its name gets into a circular mess
//  where that name might not be set up yet -- web search for
//  "static fiasco" if you are interested)

/* Old Method:
#define SOBJCLASSDEFTEMPL1(classname,tclass) \
	template <class tclass> \
	const char *classname<tclass>::IDname = \
			StreamObj::addID(typeid(classname<tclass>).name(),classname<tclass>::Create); \
	template <class tclass> \
	StreamObj *classname<tclass>::Create(std::istream &is) { \
		return new classname<tclass>(is); \
	} \
	template <class tclass> \
	const char *classname<tclass>::myID() const { return IDname; } \
*/

char *striptypename(const char *name);

/* New solution: */
#define SOBJCLASSDEFTEMPL1(classname,tclass)\
    template <class RVS>\
    const char * RVCondSimpleComp<RVS>::IDname() {\
        static const char * ret = StreamObj::addID(\
         #classname,\
         RVS::IDname(),\
         RVCondSimpleComp<RVS>::Create);\
        /* line below so that even older objects still load */\
        static const char * ret2 = StreamObj::addID(\
         striptypename(typeid(RVCondSimpleComp<RVS>).name()),\
         RVCondSimpleComp<RVS>::Create);\
        /* line below so that old objects still load */\
        static const char * ret3 = StreamObj::addID(\
         typeid(RVCondSimpleComp<RVS>).name(),\
         RVCondSimpleComp<RVS>::Create);\
        return ret;\
    }\
    template <class RVS>\
    StreamObj * RVCondSimpleComp<RVS>::Create(std::istream &is) {\
        return new RVCondSimpleComp<RVS>(is);\
    }\
    template <class RVS>\
    const char * RVCondSimpleComp<RVS>::myID() const { return IDname(); }\
    template <class RVS>\
    const char * RVCondSimpleComp<RVS>::_forcename = RVCondSimpleComp<RVS>::IDname();


// can add SOBJCLASSDEFTEMPL2 for classes with two template params if needed.

#define ENSURECLASSOLD(classname) \
	namespace ctbn_ensureclass { \
		classname *_ ## classname ## _dummy_construct(std::istream &is) { \
			return new classname(is); \
		} \
	}

} // end of ctbn namespace

#endif

