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
#ifndef CTBNRLE_SERIAL_H
#define CTBNRLE_SERIAL_H

#include "nullptr03.h"

/* I think the following remark is based upon a misunderstanding
 * of "exec" within the template metaprogramming below.
// need prototype for exec()
#ifdef _MSC_VER
#   include <process.h>
#else
#   include <unistd.h>
#endif
*/

#include <exception>
#include <string>

#ifndef SERIALNAMESPACE
#define SERIALNAMESPACE serial
#endif

namespace SERIALNAMESPACE {

    using ctbn::nullptr03;

	// can be thrown by the "Load" method
	class streamexception : public std::exception {
	public:
		streamexception(std::string txt) throw() : expl(txt) { }
		streamexception(const streamexception &e) throw()
			{ expl = e.expl; }
		streamexception &operator=(const streamexception &e) throw() 
			{ if (this != &e) expl = e.expl; return *this; }
		virtual ~streamexception() throw() { }
		virtual const char *what() const throw() { return expl.c_str(); }
	private:
		std::string expl;
	};

}

/* Here starts "helper" macros: */

#define SERIAL_GETTERCOMMON(expr,cexpr,vtype,vname) \
			typedef vtype valtype; \
			static inline vtype &getvalue(SERIAL_BASETYPE *o) \
				{ return expr; } \
			static inline const char *getname(SERIAL_BASETYPE *o) \
				{ return vname; } \
			\
			static inline vtype const &getvalue(const SERIAL_BASETYPE *o) \
				{ return cexpr; } \
			static inline const char *getname(const SERIAL_BASETYPE *o) \
				{ return vname; } \

#define SERIAL_GETTERSTRUCT(expr,cexpr,vtype,vname) \
		struct { SERIAL_GETTERCOMMON(expr,cexpr,vtype,vname) }

#define SERIAL_GETTERSTRUCT_DEFAULT(expr,cexpr,vtype,vname,dval) \
		struct { \
			SERIAL_GETTERCOMMON(expr,cexpr,vtype,vname) \
			static inline void setdefault(vtype &v) { v = dval; } \
		} 

#define SERIAL_NOTV_COMMON \
	friend struct SERIALNAMESPACE::IsEmpty<SERIAL_BASETYPE,void>; \
	private: \
		template<typename SERIAL__T,bool SERIAL__B> \
		struct serial__vallocT { \
			const static bool valid = !SERIAL__B; \
			typedef SERIAL_BASETYPE rettype; \
			typedef SERIAL_BASETYPE *(*createfntype)(void); \
			typedef std::map<std::string,createfntype> allocmaptype; \
			inline static SERIAL_BASETYPE *\
					allocbyname(const std::string &name) { \
				return nullptr03; } \
			inline static allocmaptype &alloctable() { \
				static allocmaptype table; \
				return table; \
			} \
			inline static const char * \
					addalloc(const char *name, createfntype fn) \
				{ return name; } \
		}; \

#define SERIAL_V_COMMON(myname) \
	friend struct SERIALNAMESPACE::IsEmpty<SERIAL_BASETYPE,void>; \
	private: \
		template<typename SERIAL__T, bool SERIAL__B> \
		struct serial__vallocT { \
			const static bool valid = SERIAL__B; \
			typedef SERIAL_BASETYPE rettype; \
			typedef SERIAL_BASETYPE *(*createfntype)(void); \
			typedef std::map<std::string,createfntype> allocmaptype; \
			inline static SERIAL_BASETYPE * \
					allocbyname(const std::string &name) { \
				typename allocmaptype::iterator i = alloctable().find(name); \
				if (i==alloctable().end()) return nullptr03; \
				else return i->second(); \
			} \
			inline static allocmaptype &alloctable() { \
				static allocmaptype table; \
				return table; \
			} \
			inline static const char * \
					addalloc(const char *name, createfntype fn) { \
				alloctable()[std::string(name)] = fn; \
				return name; \
			} \
		}; \
	public: \
		virtual void SaveV(std::ostream &os, const char *name="", int indent=0) const { \
			Save(os,name,indent); \
		} \
		inline static SERIAL_BASETYPE *LoadV(std::istream &is) { \
			SERIAL_BASETYPE *ret; SERIALNAMESPACE::LoadWrapper(ret,is); \
			return ret; \
		} \
		virtual void serial__loadwrapv(std::istream &is, \
					const SERIALNAMESPACE::XMLTagInfo &info) { \
			SERIALNAMESPACE::LoadWrapper(*this,info,is); \
		}

#define SERIAL_SUPER_COMMON(cname) \
		serial__class##cname##lst; \
		template<bool SERIAL__B> \
		struct serial__vallocT<typename SERIALNAMESPACE::Type_If<SERIAL__B && cname::serial__valloc::valid,SERIAL_BASETYPE>::type,SERIAL__B> { \
			const static bool valid = SERIAL__B; \
			typedef cname::serial__valloc::rettype rettype; \
			typedef cname::serial__valloc::createfntype createfntype; \
			inline static SERIAL_BASETYPE *allocbyname(const std::string &name) { \
				cname::serial__valloc::rettype *cret = cname::serial__valloc::allocbyname(name); \
				if (cret == nullptr03) return nullptr03; \
				SERIAL_BASETYPE *ret = dynamic_cast<SERIAL_BASETYPE *>(cret); \
				if (ret == nullptr03) return nullptr03; \
				return ret; \
			} \
			inline static const char *addalloc(const char *name, createfntype fn) \
				{ return cname::serial__valloc::addalloc(name,fn); } \
		}; \

#define SERIAL_NAME0(cname) \
		#cname

#define SERIAL_NAME1(cname,tname1) \
		SERIALNAMESPACE::TName(#cname,1, \
					SERIALNAMESPACE::TypeInfo<tname1>::namestr())

#define SERIAL_NAME2(cname,tname1,tname2) \
		SERIALNAMESPACE::TName(#cname,2, \
					SERIALNAMESPACE::TypeInfo<tname1>::namestr(), \
					SERIALNAMESPACE::TypeInfo<tname2>::namestr())

#define SERIAL_NAME3(cname,tname1,tname2,tname3) \
		SERIALNAMESPACE::TName(#cname,3, \
					SERIALNAMESPACE::TypeInfo<tname1>::namestr(), \
					SERIALNAMESPACE::TypeInfo<tname2>::namestr(), \
					SERIALNAMESPACE::TypeInfo<tname3>::namestr())

#define SERIAL_NAME4(cname,tname1,tname2,tname3,tname4) \
		SERIALNAMESPACE::TName(#cname,4, \
					SERIALNAMESPACE::TypeInfo<tname1>::namestr(), \
					SERIALNAMESPACE::TypeInfo<tname2>::namestr(), \
					SERIALNAMESPACE::TypeInfo<tname3>::namestr(), \
					SERIALNAMESPACE::TypeInfo<tname4>::namestr())

#define SERIAL_ID0(cname) \
	static const char *serial_IDname() { \
		static const char *ret = SERIAL_NAME0(cname); \
		return ret; \
	} \

#define SERIAL_ID1(cname,tname1) \
	static const char *serial_IDname() { \
		static const char *ret = SERIAL_NAME1(cname,tname1); \
		return ret; \
	}

#define SERIAL_ID2(cname,tname1,tname2) \
	static const char *serial_IDname() { \
		static const char *ret = SERIAL_NAME2(cname,tname1,tname2); \
		return ret; \
	}

#define SERIAL_ID3(cname,tname1,tname2,tname3) \
	static const char *serial_IDname() { \
		static const char *ret = SERIAL_NAME3(cname,tname1,tname2,tname3); \
		return ret; \
	}

#define SERIAL_ID4(cname,tname1,tname2,tname3,tname4) \
	static const char *serial_IDname() { \
		static const char *ret = \
			SERIAL_NAME4(cname,tname1,tname2,tname3,tname4); \
		return ret; \
	}

#define SERIAL_SNAME0(cname) \
	static const char *serial__shiftname() { \
		static const char *ret = SERIAL_NAME0(cname); \
		return ret; \
	} \

#define SERIAL_SNAME1(cname,tname1) \
	static const char *serial__shiftname() { \
		static const char *ret = SERIAL_NAME1(cname,tname1); \
		return ret; \
	}

#define SERIAL_SNAME2(cname,tname1,tname2) \
	static const char *serial__shiftname() { \
		static const char *ret = SERIAL_NAME2(cname,tname1,tname2); \
		return ret; \
	}

#define SERIAL_SNAME3(cname,tname1,tname2,tname3) \
	static const char *serial__shiftname() { \
		static const char *ret = SERIAL_NAME3(cname,tname1,tname2,tname3); \
		return ret; \
	}

#define SERIAL_SNAME4(cname,tname1,tname2,tname3,tname4) \
	static const char *serial__shiftname() { \
		static const char *ret = \
			SERIAL_NAME4(cname,tname1,tname2,tname3,tname4); \
		return ret; \
	}


#define SERIAL_ID0_V(cname) \
	static const char *serial_IDname() { \
		static const char *ret = serial__valloc::addalloc( \
			SERIAL_NAME0(cname), \
			SERIALNAMESPACE::creator<SERIAL_BASETYPE>::create<serial__valloc::rettype>); \
		return ret; \
	} \

#define SERIAL_ID1_V(cname,tname1) \
	static const char *serial_IDname() { \
		static const char *ret = serial__valloc::addalloc( \
			SERIAL_NAME1(cname,tname1), \
			SERIALNAMESPACE::template creator<SERIAL_BASETYPE>::template create<typename serial__valloc::rettype>); \
		return ret; \
	}

#define SERIAL_ID2_V(cname,tname1,tname2) \
	static const char *serial_IDname() { \
		static const char *ret = serial__valloc::addalloc( \
			SERIAL_NAME2(cname,tname1,tname2), \
			SERIALNAMESPACE::template creator<SERIAL_BASETYPE>::template create<typename serial__valloc::rettype>); \
		return ret; \
	}

#define SERIAL_ID3_V(cname,tname1,tname2,tname3) \
	static const char *serial_IDname() { \
		static const char *ret = serial__valloc::addalloc( \
			SERIAL_NAME3(cname,tname1,tname2,tname3), \
			SERIALNAMESPACE::template creator<SERIAL_BASETYPE>::template create<typename serial__valloc::rettype>); \
		return ret; \
	}

#define SERIAL_ID4_V(cname,tname1,tname2,tname3,tname4) \
	static const char *serial_IDname() { \
		static const char *ret = serial__valloc::addalloc( \
			SERIAL_NAME4(cname,tname1,tname2,tname3,tname4), \
			SERIALNAMESPACE::creator<SERIAL_BASETYPE>::create<typename serial__valloc::rettype>); \
		return ret; \
	}

/* Here starts "USESHIFT" macros: */

#define SERIAL_USESHIFT(cname) \
	private: \
		typedef cname SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME0(cname) \
		enum { SERIAL_ISSHORT = false };
		//inline static bool serial__isshort() { return false; }

#define SERIAL_USESHIFT_SHORT(cname) \
	private: \
		typedef cname SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME0(cname) \
		enum { SERIAL_ISSHORT = true };
		//inline static bool serial__isshort() { return true; }

#define SERIAL_USESHIFT1(cname,tname1) \
	private: \
		typedef cname<tname1> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME1(cname,tname1) \
		inline static bool serial__isshort() { return false; }

#define SERIAL_USESHIFT1_SHORT(cname) \
	private: \
		typedef cname<tname1> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME1(cname,tname1) \
		inline static bool serial__isshort() { return true; }

#define SERIAL_USESHIFT2(cname,tname1,tname2) \
	private: \
		typedef cname<tname1,tname2> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME2(cname,tname1,tname2) \
		inline static bool serial__isshort() { return false; }

#define SERIAL_USESHIFT2_SHORT(cname,tname1,tname2) \
	private: \
		typedef cname<tname1,tname2> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME2(cname,tname1,tname2) \
		inline static bool serial__isshort() { return true; }

#define SERIAL_USESHIFT3(cname) \
	private: \
		typedef cname<tname1,tname2,tname3> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME3(cname,tname1,tname2,tname3) \
		inline static bool serial__isshort() { return false; }

#define SERIAL_USESHIFT3_SHORT(cname,tname1,tname2,tname3) \
	private: \
		typedef cname<tname1,tname2,tname3> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME3(cname,tname1,tname2,tname3) \
		inline static bool serial__isshort() { return true; }

#define SERIAL_USESHIFT4(cname,tname1,tname2,tname3,tname4) \
	private: \
		typedef cname<tname1,tname2,tname3,tname4> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME4(cname,tname1,tname2,tname3,tname4) \
		inline static bool serial__isshort() { return false; }

#define SERIAL_USESHIFT4_SHORT(cname) \
	private: \
		typedef cname<tname1,tname2,tname3,tname4> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_SNAME4(cname,tname1,tname2,tname3,tname4) \
		inline static bool serial__isshort() { return true; }

/* Here starts "START" macros: */

#define SERIAL_START(cname) \
	private: \
		typedef cname SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_ID0(cname) \
	private: \
		typedef SERIALNAMESPACE::ListEnd 

#define SERIAL_START1(cname,tname) \
	private: \
		typedef cname<tname> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_ID1(cname,tname) \
	private: \
		typedef SERIALNAMESPACE::ListEnd

#define SERIAL_START2(cname,tname1,tname2) \
	private: \
		typedef cname<tname1,tname2> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_ID2(cname,tname1,tname2) \
	private: \
		typedef SERIALNAMESPACE::ListEnd

#define SERIAL_START3(cname,tname1,tname2,tname3) \
	private: \
		typedef cname<tname1,tname2,tname3> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_ID3(cname,tname1,tname2,tname3) \
	private: \
		typedef SERIALNAMESPACE::ListEnd

#define SERIAL_START4(cname,tname1,tname2,tname3,tname4) \
	private: \
		typedef cname<tname1,tname2,tname3,tname4> SERIAL_BASETYPE; \
	SERIAL_NOTV_COMMON \
	public: \
		SERIAL_ID4(cname,tname1,tname2,tname3,tname4) \
	private: \
		typedef SERIALNAMESPACE::ListEnd

/* Here starts "START_V" macros: */

#define SERIAL_START_V(cname) \
	private: \
		typedef cname SERIAL_BASETYPE; \
		SERIAL_V_COMMON(cname) \
	public: \
		SERIAL_ID0_V(cname); \
	private: \
		typedef SERIALNAMESPACE::ListEnd

#define SERIAL_START1_V(cname,tname) \
	private: \
		typedef cname<tname> SERIAL_BASETYPE; \
		SERIAL_V_COMMON(cname) \
	public: \
		SERIAL_ID1_V(cname,tname); \
	private: \
		typedef SERIALNAMESPACE::ListEnd

#define SERIAL_START2_V(cname,tname1,tname2) \
	private: \
		typedef cname<tname> SERIAL_BASETYPE; \
		SERIAL_V_COMMON(cname) \
	public: \
		SERIAL_ID2_V(cname,tname1,tname2); \
	private: \
		typedef SERIALNAMESPACE::ListEnd

#define SERIAL_START3_V(cname,tname1,tname2,tname3) \
	private: \
		typedef cname<tname> SERIAL_BASETYPE; \
		SERIAL_V_COMMON(cname) \
	public: \
		SERIAL_ID3_V(cname,tname1,tname2,tname3); \
	private: \
		typedef SERIALNAMESPACE::ListEnd

#define SERIAL_START4_V(cname,tname1,tname2,tname3,tname4) \
	private: \
		typedef cname<tname> SERIAL_BASETYPE; \
		SERIAL_V_COMMON(cname) \
	public: \
		SERIAL_ID4_V(cname,tname1,tname2,tname3,tname4); \
	private: \
		typedef SERIALNAMESPACE::ListEnd

/* Here starts "VAR" macros */

#define SERIAL_VAR(vtype,vname) \
		serial__##vname##lst; \
		typedef SERIAL_GETTERSTRUCT((o->vname),(o->vname),vtype,#vname) \
			serial__get##vname; \
		typedef SERIALNAMESPACE::List<serial__get##vname,serial__##vname##lst>

#define SERIAL_VAR_DEFAULT(vtype,vname,val) \
		serial__##vname##lst; \
		typedef SERIAL_GETTERSTRUCT_DEFAULT((o->vname),(o->vname),vtype,#vname,val) \
			serial__get##vname; \
		typedef SERIALNAMESPACE::List<serial__get##vname,serial__##vname##lst>

#define SERIAL_VAR_N(vtype,vname,sname) \
		serial__##vname##lst; \
		typedef SERIAL_GETTERSTRUCT((o->vname),(o->vname),vtype,sname) \
			serial__get##vname; \
		typedef SERIALNAMESPACE::List<serial__get##vname,serial__##vname##lst>

#define SERIAL_VAR_DEFAULT_N(vtype,vname,val,sname) \
		serial__##vname##lst; \
		typedef SERIAL_GETTERSTRUCT_DEFAULT((o->vname),(o->vname),vtype,sname,val) \
			serial__get##vname; \
		typedef SERIALNAMESPACE::List<serial__get##vname,serial__##vname##lst>

/* Here starts "SUPER" macros */

#define SERIAL_SUPER(cname) \
		SERIAL_SUPER_COMMON(cname) \
		typedef SERIAL_GETTERSTRUCT((*(static_cast<cname *>(o))),*(static_cast<const cname *>(o)),cname,#cname) \
			serial__getclass##cname; \
		typedef SERIALNAMESPACE::List<serial__getclass##cname,serial__class##cname##lst>
			
			
#define SERIAL_SUPER_DEFAULT(cname,val) \
	SERIAL_SUPER_COMMON(cname) \
	typedef SERIAL_GETTERSTRUCT_DEFAULT((*(static_cast<cname *>(o))), \
		*(static_cast<const cname *>(o)),cname,#cname,val) \
		serial__getclass##cname; \
	typedef SERIALNAMESPACE::List<serial__getclass##cname, \
		serial__class##cname##lst>
			
#define SERIAL_SUPER_N(cname,sname) \
		SERIAL_SUPER_COMMON(cname) \
		typedef SERIAL_GETTERSTRUCT((*(static_cast<cname *>(o))),*(static_cast<const cname *>(o)),cname,sname) \
			serial__getclass##cname; \
		typedef SERIALNAMESPACE::List<serial__getclass##cname,serial__class##cname##lst>
			
			
#define SERIAL_SUPER_DEFAULT_N(cname,val,sname) \
	SERIAL_SUPER_COMMON(cname) \
	typedef SERIAL_GETTERSTRUCT_DEFAULT((*(static_cast<cname *>(o))), \
		*(static_cast<const cname *>(o)),cname,sname,val) \
		serial__getclass##cname; \
	typedef SERIALNAMESPACE::List<serial__getclass##cname, \
		serial__class##cname##lst>
			

/* Here starts "ENSURECLASS" macros */

#define ENSURECLASS(cname) \
	namespace { static const char *namefor__##cname=cname::serial_IDname(); }

#define ENSURECLASS1(cname,tname1) \
	namespace { static const char *namefor__##cname##tname1 = \
		cname<tname1>::serial_IDname(); }

#define ENSURECLASS2(cname,tname1,tname2) \
	namespace { static const char *namefor__##cname##tname1##tname2 = \
		cname<tname1,tname2>::serial_IDname(); }

#define ENSURECLASS3(cname,tname1,tname2,tname3) \
	namespace { \
		static const char *namefor__##cname##tname1##tname2##tname3 = \
			cname<tname1,tname2,tname3>::serial_IDname(); \
	}

#define ENSURECLASS4(cname,tname1,tname2,tname3,tname4) \
	namespace { \
		static const char *namefor__##cname##tname1##tname2##tname3##tname4 \
			= cname<tname1,tname2,tname3,tname4>::serial_IDname(); \
	}

//#include "serial.tcc"

#include <map>
#include <set>
#include <vector>
#include <set>
#include <utility>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include <typeinfo>
#include <cstdarg>
#include <iostream>

namespace SERIALNAMESPACE {

    using ctbn::nullptr03;

    char *TName(const char *cname, int n, ...);
    void Indent(std::ostream &os, int indent);
    void *AllocByName(const std::string &name);
    typedef void *(*createfntype)(void);
    const char *AddAlloc(const char *name, createfntype fn);

    // a simple list:
    struct ListEnd {};

    template<typename H, typename T>
    struct List {};

    // pretty standard enable_if template (see Boost, for example)
    template<bool B, typename T>
    struct Type_If {
        typedef T type;
    };
    template<typename T>
    struct Type_If<false,T> {
    };

#define SERIAL_DECVAL(vname,expr) \
    enum { vname = (expr) };
    //static int const vname = (expr);

    // This one has appeared in one form or another all over the
    // internet

    template<typename T>
    struct TypeProp {
        template<typename C, C> struct type_check;

        template<typename S> static char
            (& chksave(type_check<void (S::*)(std::ostream&,int),
                           &S::serial__Save>*))[1];
        template<typename> static char
            (& chksave(...))[2];
        SERIAL_DECVAL(HasSave,sizeof(chksave<T>(0)) == 1)

        template<typename S> static char
            (& chkid(type_check<const char *(*)(), &S::serial_IDname>*))[1];
        template<typename> static char
            (& chkid(...))[2];
        SERIAL_DECVAL(HasIDname,sizeof(chkid<T>(0)) == 1)

        template<typename S> static char
            (& chkshift(type_check<const char *(*)(),
                    &S::serial__shiftname>*))[1];
        template<typename> static char
            (& chkshift(...))[2];
        SERIAL_DECVAL(HasShift,sizeof(chkshift<T>(0)) == 1)

        template<typename S> static char
            (& chkdef(type_check<void (*)(typename S::valtype &),
                &S::setdefault>*))[1];
        template<typename> static char
            (& chkdef(...))[2];
        SERIAL_DECVAL(HasDefault,sizeof(chkdef<T>(0)) == 1)

        template<typename S> static char
            (& chkv(type_check<T *(*)(std::istream &),
                &S::LoadV>*))[1];
        template<typename> static char
            (& chkv(...))[2];
        SERIAL_DECVAL(HasV,sizeof(chkv<T>(0)) == 1)

        template<typename S> static char
            (& chkpreload(type_check<void (S::*)(void),
                &S::serial_preload>*))[1];
        template<typename> static char
            (& chkpreload(...))[2];
        SERIAL_DECVAL(HasPreLoad,sizeof(chkpreload<T>(0)) == 1)

        template<typename S> static char
            (& chkpostload(type_check<void (S::*)(void),
                &S::serial_postload>*))[1];
        template<typename> static char
            (& chkpostload(...))[2];
        SERIAL_DECVAL(HasPostLoad,sizeof(chkpostload<T>(0)) == 1)

        template<typename S> static char
            (& chkpresave(type_check<void (S::*)(void) const,
                &S::serial_presave>*))[1];
        template<typename> static char
            (& chkpresave(...))[2];
        SERIAL_DECVAL(HasPreSave,sizeof(chkpresave<T>(0)) == 1)

        template<typename S> static char
            (& chkpostsave(type_check<void (S::*)(void) const,
                &S::serial_postsave>*))[1];
        template<typename> static char
            (& chkpostsave(...))[2];
        SERIAL_DECVAL(HasPostSave,sizeof(chkpostsave<T>(0)) == 1)

        //enum { BlankList = SameType<ListEnd,typename T::serial__alllist>::value };
    };

    template<typename T>
    struct HasDefCon {
        struct small { char x[1]; };
        struct large { char x[2]; };
        template<typename U>
            static small chk(U (*)[1]);
        template<typename U>
            static large chk(...);
        //SERIAL_DECVAL(value,sizeof(HasDefCon<T>::template chk<T>(0))==sizeof(small))
        SERIAL_DECVAL(value,sizeof(chk<T>(0))==sizeof(small))
    };

    template<typename T,typename Condition=void>
    struct IsEmpty {
        SERIAL_DECVAL(value,false)
    };

    template<typename L>
    struct ListEmpty {
        SERIAL_DECVAL(value,false)
    };

    template<>
    struct ListEmpty<ListEnd> {
        SERIAL_DECVAL(value,true)
    };

    template<typename H, typename T>
    struct ListEmpty<List<H,T> > {
        SERIAL_DECVAL(value,IsEmpty<typename H::valtype>::value && ListEmpty<T>::value)
    };

    template<typename T>
    struct IsEmpty<T,typename Type_If<ListEmpty<typename T::serial__alllist>::value,void>::type> {
        SERIAL_DECVAL(value,true)
    };


    template<typename AT>
    struct creator {
        template<typename RT>
        inline static typename Type_If<HasDefCon<AT>::value,RT *>::type
        create() { return new AT(); }

        template<typename RT>
        inline static typename Type_If<!HasDefCon<AT>::value,RT *>::type
        create() { return nullptr03; }
    };


    struct XMLTagInfo {
        std::string name;
        std::map<std::string,std::string> attr;
        bool isstart,isend;
    };

    template<typename G> struct SaveItem;
    template<typename L> struct LoadList;
    void ReadTag(std::istream &is,XMLTagInfo &info);
    void ReadEndTag(std::istream &is,const char *ename);
    void ReadStr(std::istream &is, std::string &ret,const char *endchar);
    void WriteStr(std::ostream &os, const std::string &s,bool escape=true);

    // a reverse iterator on a list for saving
    template<typename L,typename Condition=void>
    struct SaveItt {
        template<typename O>
        inline static void exec(O o,std::ostream &os,int indent) { }
    };

    template<typename H, typename T>
    struct SaveItt<List<H,T>,
        typename Type_If<!IsEmpty<typename H::valtype>::value,void>::type> {
        template<typename O>
        inline static void exec(O o,std::ostream &os,int indent) {
            SaveItt<T>::exec(o,os,indent);
            SaveItem<H>::apply(o,os,indent);
        }
    };

    template<typename H, typename T>
    struct SaveItt<List<H,T>,
        typename Type_If<IsEmpty<typename H::valtype>::value,void>::type> {
        template<typename O>
        inline static void exec(O o,std::ostream &os,int indent) {
            SaveItt<T>::exec(o,os,indent);
        }
    };

    template<>
    struct SaveItt<ListEnd,void> {
        template<typename O>
        inline static void exec(O o,std::ostream &os, int indent) { }
    };


    template<typename L>
    struct LoadOne {
        template<typename O>
        inline static bool exec(O o,std::istream &, const XMLTagInfo &) {
            return false;
        }
    };

    template<>
    struct LoadOne<ListEnd> {
        template<typename O>
        inline static bool exec(O o,std::istream &, const XMLTagInfo &) {
            return false;
        }
    };

    template<typename H, typename T>
    struct LoadOne<List<H,T> > {
        template<typename O>
        inline static bool exec(O o,std::istream &is,
                const XMLTagInfo &info) {
            std::map<std::string,std::string>::const_iterator ni
                =info.attr.find("name");
            if (ni!=info.attr.end() && ni->second == H::getname(o)) {
                LoadWrapper(H::getvalue(o),info,is);
                return true;
            } else
                return LoadOne<T>::exec(o,is,info);
        }
    };

    template<typename L,typename Condition=void>
    struct AllLoaded {
        template<typename O>
        inline static bool exec(O o,const std::set<std::string> &nset) {
            return true;
        }
    };

    template<>
    struct AllLoaded<ListEnd,void> {
        template<typename O>
        inline static bool exec(O o,const std::set<std::string> &nset) {
            return true;
        }
    };

    template<typename H, typename T>
    struct AllLoaded<List<H,T>,
      typename Type_If<TypeProp<H>::HasDefault,void>::type > {
        template<typename O>
        inline static bool exec(O o,const std::set<std::string> &nset) {
            if (nset.find(H::getname(o)) == nset.end())
                H::setdefault(H::getvalue(o));
            return AllLoaded<T>::exec(o,nset);
        }
    };

    template<typename H, typename T>
    struct AllLoaded<List<H,T>,
      typename Type_If<!TypeProp<H>::HasDefault
                      && IsEmpty<typename H::valtype>::value,void>::type > {
        template<typename O>
        inline static bool exec(O o,const std::set<std::string> &nset) {
            return AllLoaded<T>::exec(o,nset);
        }
    };

    template<typename H, typename T>
    struct AllLoaded<List<H,T>,
            typename Type_If<!TypeProp<H>::HasDefault
                && !IsEmpty<typename H::valtype>::value,void>::type > {
        template<typename O>
        inline static bool exec(O o,const std::set<std::string> &nset) {
            if (nset.find(H::getname(o)) == nset.end()) return false;
            else return AllLoaded<T>::exec(o,nset);
        }
    };

    template<typename T,typename Condition=void>
    struct IsShiftable {
        SERIAL_DECVAL(atall,false)
        SERIAL_DECVAL(quickly,false)
    };

    template<typename T>
    struct IsShiftable<const T,void> {
        SERIAL_DECVAL(atall,IsShiftable<T>::atall)
        SERIAL_DECVAL(quickly,IsShiftable<T>::quickly)
    };

    template<typename T>
    struct IsShiftable<T,typename Type_If<TypeProp<T>::HasShift,void>::type> {
        SERIAL_DECVAL(atall,true)
        SERIAL_DECVAL(quickly,T::SERIAL_ISSHORT)
    };

#define SERIAL__BASETYPE(tname,strname) \
    template<> \
    struct IsShiftable<tname,void> { \
        SERIAL_DECVAL(atall,true) \
        SERIAL_DECVAL(quickly,true) \
    }; \
    template<> \
    struct TypeInfo<tname,void> { \
        inline static const char *namestr() { return #strname; } \
        inline static void writeotherattr(std::ostream &, const tname &) { } \
        inline static bool isshort(const tname &) { return true; } \
        inline static bool isinline(const tname &) { return false; } \
        inline static void save(const tname &t,std::ostream &os, int indent) { \
            os << t; \
        } \
        inline static void load(tname &t, const XMLTagInfo &info, std::istream &is) {\
            std::map<std::string,std::string>::const_iterator vi  \
                =info.attr.find("value"); \
            if (vi!=info.attr.end()) { \
                std::istringstream ss(vi->second); \
                ss.copyfmt(is); \
                ss >> t; \
                if (info.isend) return; \
            } else { \
                is >> t; \
            } \
            ReadEndTag(is,namestr()); \
        } \
    };


    template<typename T,typename Condition=void>
    struct TypeInfo {
        enum { x = sizeof(T::__Attempting_to_serialize_type_without_serialization_information) };

        inline static const char *namestr() {
            throw streamexception("Streaming Error: type info asked about unknown class");
        }
        inline static void writeotherattr(std::ostream &, const T &) {
            throw streamexception("Streaming Error: type info asked about unknown class");
        }
        inline static bool isshort(const T &)  {
            throw streamexception("Streaming Error: type info asked about unknown class");
        }
        inline static bool isinline(const T &)  {
            throw streamexception("Streaming Error: type info asked about unknown class");
        }
        inline static void save(const T &t,std::ostream &os, int indent) {
            throw streamexception("Streaming Error: save called for unknown class");
        }
        inline static void load(T &t, const XMLTagInfo &info,
                    std::istream &is) {
            throw streamexception("Streaming Error: load called for unknown class");
        }
    };

    // These were moved upward to avoid ambiguity.
    SERIAL__BASETYPE(bool,bool)
    SERIAL__BASETYPE(char,char)
    SERIAL__BASETYPE(unsigned char,u_char)
    SERIAL__BASETYPE(short,short)
    SERIAL__BASETYPE(unsigned short,u_short)
    SERIAL__BASETYPE(int,int)
    SERIAL__BASETYPE(unsigned int,u_int)
    SERIAL__BASETYPE(long,long)
    SERIAL__BASETYPE(unsigned long,u_long)
    SERIAL__BASETYPE(long long,l_long)
    SERIAL__BASETYPE(float,float)
    SERIAL__BASETYPE(double,double)
    SERIAL__BASETYPE(long double,l_double)


    template<typename T>
    struct TypeInfo<const T,void> {
        inline static const char *namestr() {
            return TypeInfo<T>::namestr();
        }
        inline static void writeotherattr(std::ostream &os, const T &t) {
            TypeInfo<T>::writeotherattr(os,t);
        }
        inline static void writeotherattr(std::ostream &os, const T &t,
                const char *s1, const char *s2) {
            TypeInfo<T>::otherattr(os,t,s1,s2);
        }
        inline static bool isshort(const T &t)  {
            return TypeInfo<T>::isshort(t);
        }
        inline static bool isinline(const T &t)  {
            return TypeInfo<T>::isinline(t);
        }
        inline static void save(const T &t,std::ostream &os, int indent) {
            TypeInfo<T>::save(t,os,indent);
        }
        inline static void load(const T &t, const XMLTagInfo &info,
                    std::istream &is) {
            throw streamexception("Streaming Error: load called for constant type");
        }
    };


    template<typename T>
    inline void SaveWrapper(const T &v,const char *vname,
            std::ostream &os,int indent) {
        Indent(os,indent);
        os << "<" << TypeInfo<T>::namestr();
        TypeInfo<T>::writeotherattr(os,v);
        if (vname && vname[0]!=0) os << " name=\"" << vname << "\"";
        if (TypeInfo<T>::isshort(v)) {
            os << " value=\"";
            TypeInfo<T>::save(v,os,indent);
            os << "\" \\>" << std::endl;
        } else if (TypeInfo<T>::isinline(v)) {
            os << " \\>" << std::endl;
        } else {
            os << ">";
            TypeInfo<T>::save(v,os,indent);
            os << "<\\" << TypeInfo<T>::namestr() << ">" << std::endl;
        }
    }

    template<typename T>
    inline typename Type_If<TypeProp<T>::HasV,void>::type
    SaveWrapper(T * const &v, const char *vname,
            std::ostream &os, int indent) {
        if (v == nullptr03) {
            Indent(os,indent);
            os << "<" << TypeInfo<T>::namestr();
            if (vname && vname[0]!=0) os << " name=\"" << vname << "\"";
            os << " isnull=\"1\" \\>" << std::endl;
        } else v->SaveV(os,vname,indent);
    }

    template<typename T>
    inline typename Type_If<!TypeProp<T>::HasV,void>::type
    SaveWrapper(T * const &v, const char *vname,
            std::ostream &os, int indent) {
        if (v == nullptr03) {
            Indent(os,indent);
            os << "<" << TypeInfo<T>::namestr();
            if (vname && vname[0]!=0) os << " name=\"" << vname << "\"";
            os << " isnull=\"1\" \\>" << std::endl;
        } else SaveWrapper(*v,vname,os,indent);
    }

    template<typename T>
    inline void SaveWrapper(const T &v,const char *vname,
            std::ostream &os,int indent, const char *e1, const char *e2) {
        Indent(os,indent);
        os << "<" << TypeInfo<T>::namestr();
        TypeInfo<T>::writeotherattr(os,v,e1,e2);
        if (vname && vname[0]!=0) os << " name=\"" << vname << "\"";
        if (TypeInfo<T>::isinline(v)) {
            os << " \\>" << std::endl;
        } else {
            os << ">";
            TypeInfo<T>::save(v,os,indent,e1,e2);
            os << "<\\" << TypeInfo<T>::namestr() << ">" << std::endl;
        }
    }


    template<typename T>
    struct CheckTag {
        inline static void check(const XMLTagInfo &info) {
            if (!info.isstart)
                throw streamexception(std::string("Stream Input Format Error: expected start tag for ")+TypeInfo<T>::namestr()+", received end tag for "+info.name);
            if (info.name != TypeInfo<T>::namestr())
                throw streamexception(std::string("Stream Input Format Error: expected start tag for ")+TypeInfo<T>::namestr()+", received start tag for "+info.name);
        }
    };

    template<typename T>
    inline void LoadWrapper(T &v, const XMLTagInfo &info, std::istream &is) {
        CheckTag<T>::check(info);
        TypeInfo<T>::load(v,info,is);
    }

    template<typename T>
    inline typename Type_If<TypeProp<T>::HasV,void>::type
    LoadWrapper(T *&v, const XMLTagInfo &info, std::istream &is) {
        if (!info.isstart)
            throw streamexception(std::string("Stream Input Format Error: expected start tag, received end tag for ")+info.name);
        std::map<std::string,std::string>::const_iterator vi
            =info.attr.find("isnull");
        if (vi!=info.attr.end() && vi->second=="1") {
            v = nullptr03;
            if (!info.isend) {
                XMLTagInfo einfo;
                ReadTag(is,einfo);
                if (einfo.name!=info.name || !info.isend || info.isstart)
                    throw streamexception(std::string("Stream Input Format Error: null pointer object for type ")+TypeInfo<T>::namestr()+" and name "+info.name+" has non-empty contents");
            }
            return;
        }
        T *vv = T::serial__valloc::allocbyname(info.name);
        if (vv == nullptr03)
            throw streamexception(std::string("Stream Input Format Error: expected start tag for subtype of ")+TypeInfo<T>::namestr()+", received start tag for type "+info.name+" which is either unknown or not a subtype");
        v = vv;
        v->serial__loadwrapv(is,info);
    }

    template<typename T>
    inline typename Type_If<!TypeProp<T>::HasV,void>::type
    LoadWrapper(T *&v, const XMLTagInfo &info, std::istream &is) {
        std::map<std::string,std::string>::const_iterator vi
            =info.attr.find("isnull");
        if (vi!=info.attr.end() && vi->second=="1") {
            v = nullptr03;
            if (!info.isend) {
                XMLTagInfo einfo;
                ReadTag(is,einfo);
                if (einfo.name!=info.name || !info.isend || info.isstart)
                    throw streamexception(std::string("Stream Input Format Error: null pointer object for type ")+TypeInfo<T>::namestr()+" and name "+info.name+" has non-empty contents");
            }
            return;
        }
        v = new T();
        LoadWrapper(*v,info,is);
    }

    template<typename T>
    inline void LoadWrapper(T &v, std::istream &is) {
        XMLTagInfo info;
        ReadTag(is,info);
        LoadWrapper(v,info,is);
    }

    template<typename T>
    inline typename Type_If<!TypeProp<T>::HasPreLoad
                         && !TypeProp<T>::HasPostLoad,void>::type
    SerialLoadWrap(T &t, std::istream &is) {
        t.serial__Load(is);
    }

    template<typename T>
    inline typename Type_If<TypeProp<T>::HasPreLoad
                         && !TypeProp<T>::HasPostLoad,void>::type
    SerialLoadWrap(T &t, std::istream &is) {
        t.serial_preload();
        t.serial__Load(is);
    }

    template<typename T>
    inline typename Type_If<!TypeProp<T>::HasPreLoad
                         && TypeProp<T>::HasPostLoad,void>::type
    SerialLoadWrap(T &t, std::istream &is) {
        t.serial__Load(is);
        t.serial_postload();
    }

    template<typename T>
    inline typename Type_If<TypeProp<T>::HasPreLoad
                         && TypeProp<T>::HasPostLoad,void>::type
    SerialLoadWrap(T &t, std::istream &is) {
        t.serial_preload();
        t.serial__Load(is);
        t.serial_postload();
    }

    template<typename T>
    inline typename Type_If<!TypeProp<T>::HasPreSave
                    && !TypeProp<T>::HasPostSave,void>::type
    SerialSaveWrap(const T &t, std::ostream &os, int indent) {
        t.serial__Save(os,indent);
    }

    template<typename T>
    inline typename Type_If<TypeProp<T>::HasPreSave
                    && !TypeProp<T>::HasPostSave,void>::type
    SerialSaveWrap(const T &t, std::ostream &os, int indent) {
        t.serial_presave();
        t.serial__Save(os,indent);
    }

    template<typename T>
    inline typename Type_If<!TypeProp<T>::HasPreSave
                    && TypeProp<T>::HasPostSave,void>::type
    SerialSaveWrap(const T &t, std::ostream &os, int indent) {
        t.serial__Save(os,indent);
        t.serial_postsave();
    }

    template<typename T>
    inline typename Type_If<TypeProp<T>::HasPreSave
                    && TypeProp<T>::HasPostSave,void>::type
    SerialSaveWrap(const T &t, std::ostream &os, int indent) {
        t.serial_presave();
        t.serial__Save(os,indent);
        t.serial_postsave();
    }

    template<>
    struct TypeInfo<std::string,void> {
        inline static const char *namestr() { return "string"; }
        inline static void writeotherattr(std::ostream &, const std::string &) { }
        inline static bool isshort(const std::string &s) { return s.length()<20; }
        inline static bool isinline(const std::string &s) { return false; }
        inline static void save(const std::string &s, std::ostream &os,
                int indent) {
            WriteStr(os,s,isshort(s));
        }
        inline static void load(std::string &s, const XMLTagInfo &info,
                std::istream &is) {
            std::map<std::string,std::string>::const_iterator vi
                =info.attr.find("value");
            if (vi!=info.attr.end()) {
                std::istringstream ss(vi->second);
                ss.copyfmt(is);
                ReadStr(ss,s,"");
                if (info.isend) return;
            } else ReadStr(is,s,"<");
            ReadEndTag(is,namestr());
        }
    };

    template<typename T>
    struct TypeInfo<T, typename Type_If<TypeProp<T>::HasIDname,void>::type> {
        inline static const char *namestr() { return T::serial_IDname(); }
        inline static void writeotherattr(std::ostream &,const T &) { }
        inline static bool isshort(const T &) { return false; }
        inline static bool isinline(const T &) { return false; }
        inline static void save(const T &t,std::ostream &os,int indent) {
            os << std::endl;
            SerialSaveWrap(t,os,indent+1);
            Indent(os,indent);
        }
        inline static void load(T &t,const XMLTagInfo &info,
                std::istream &is) {
            SerialLoadWrap(t,is);
        }
    };

    template<typename T>
    struct TypeInfo<T, typename Type_If<TypeProp<T>::HasShift,void>::type> {
        inline static const char *namestr() { return T::serial__shiftname(); }
        inline static void writeotherattr(std::ostream &, const T &) {}
        inline static bool isshort(const T &) { return T::SERIAL_ISSHORT; }
        inline static bool isinline(const T &) { return false; }
        inline static void save(const T &t,std::ostream &os,int indent) {
            os << t;
        }
        inline static void load(T &t, const XMLTagInfo &info,
                std::istream &is) {
            std::map<std::string,std::string>::const_iterator vi
                =info.attr.find("value");
            if (vi!=info.attr.end()) {
                std::istringstream ss(vi->second);
                ss.copyfmt(is);
                ss >> t;
                if (info.isend) return;
            } else {
                is >> t;
            }
            ReadEndTag(is,namestr());
        }
    };

    template<typename T1,typename T2>
    struct TypeInfo<std::pair<T1,T2>,
        typename Type_If<!IsShiftable<T1>::quickly
                      || !IsShiftable<T2>::quickly,void>::type> {
        inline static const char *namestr() { return "pair"; }
        inline static void writeotherattr(std::ostream &, const std::pair<T1,T2> &,
                const char * fn = nullptr03, const char * sn = nullptr03) {
        }
        inline static bool isshort(const std::pair<T1,T2> &) { return false; }
        inline static bool isinline(const std::pair<T1,T2> &) { return false; }
        inline static void save(const std::pair<T1,T2> &p,
                std::ostream &os,int indent,
                    const char *firstname= "first",
                    const char *secondname="second") {
            os << std::endl;
            SaveWrapper(p.first,firstname,os,indent+1);
            SaveWrapper(p.second,secondname,os,indent+1);
            Indent(os,indent);
        }

        class PName {
        public:
            PName(std::pair<T1,T2> &pa,const char *n1, const char *n2)
                : p(pa), s1(n1), s2(n2) {}
            std::pair<T1,T2> &p;
            const char *s1,*s2;
        };

        struct Get1 {
            inline static T1 &getvalue(PName &o)
                { return o.p.first; }
            inline static const char *getname(PName &o)
                { return o.s1; }
        };
        struct Get2 {
            inline static T2 &getvalue(PName &o)
                { return o.p.second; }
            inline static const char *getname(PName &o)
                { return o.s2; }
        };

        typedef List<Get1,List<Get2,ListEnd> > PList;

        inline static void load(std::pair<T1,T2> &p,
                const XMLTagInfo &info, std::istream &is,
                    const char *firstname = "first",
                    const char *secondname = "second") {
            PName pn(p,firstname,secondname);
            LoadList<PList>::exec(pn,is,"pair"); \
        }
    };

    template<typename T1,typename T2>
    struct TypeInfo<std::pair<T1,T2>,
        typename Type_If<IsShiftable<T1>::quickly
                      && IsShiftable<T2>::quickly,void>::type> {
        inline static const char *namestr() {
            static std::string ret = std::string("pair.")
                +TypeInfo<T1>::namestr()
                +"."+TypeInfo<T2>::namestr();
            return ret.c_str();
        }
        inline static void writeotherattr(std::ostream &os, const std::pair<T1,T2> &p,
                const char *firstname="first",
                const char *secondname="second") {
            os << " " << firstname << "=\"" << p.first << "\" "
                << secondname << "=\"" << p.second << "\"";
        }
        inline static bool isshort(const std::pair<T1,T2> &) { return false; }
        inline static bool isinline(const std::pair<T1,T2> &) { return true; }
        inline static void save(const std::pair<T1,T2> &p,
                std::ostream &os,int indent,
                const char *firstname= "first",
                const char *secondname="second") {
        }

        inline static void load(std::pair<T1,T2> &p,
                const XMLTagInfo &info, std::istream &is,
                    const char *firstname = "first",
                    const char *secondname = "second") {
            std::map<std::string,std::string>::const_iterator vi
                =info.attr.find(firstname);
            if (vi!=info.attr.end()) {
                std::istringstream ss(vi->second);
                ss.copyfmt(is);
                ss >> p.first;
            } else
                throw streamexception(std::string("Stream Input Format Error: missing attribute ")+firstname+" in tag "+info.name);
            vi = info.attr.find(secondname);
            if (vi!=info.attr.end()) {
                std::istringstream ss(vi->second);
                ss.copyfmt(is);
                ss >> p.second;
            } else
                throw streamexception(std::string("Stream Input Format Error: missing attribute ")+secondname+" in tag "+info.name);
            if (!info.isend) {
                XMLTagInfo einfo;
                ReadTag(is,einfo);
                if (einfo.name!=info.name || !info.isend || info.isstart)
                    throw streamexception(std::string("Stream Input Format Error: tag ")+info.name+" has non-empty contents");
            }
        }
    };


    template<typename T,typename A>
    struct TypeInfo<std::vector<T,A>,
            typename Type_If<!IsShiftable<T>::atall,void>::type> {
        inline static const char *namestr() { return "vector"; }
        inline static void writeotherattr(std::ostream &os, const std::vector<T,A> &v) {
            os << " nelem=" << v.size();
        }
        inline static bool isshort(const std::vector<T,A> &) { return false; }
        inline static bool isinline(const std::vector<T,A> &) { return false; }
        inline static void save(const std::vector<T,A> &v,
                std::ostream &os,int indent) {
            os << std::endl;
            int c=0;
            for(typename std::vector<T,A>::const_iterator i=v.begin();
                i!=v.end();++i,++c) {
                SaveWrapper(*i,"",os,indent+1);
            }
            Indent(os,indent);
        }
        inline static void load(std::vector<T,A> &v, const XMLTagInfo &info,
                    std::istream &is) {
            std::map<std::string,std::string>::const_iterator ni
                    = info.attr.find("nelem");
            int n = 0;
            if (ni != info.attr.end())
                v.resize(n = atoi(ni->second.c_str()));
            else v.resize(0);
            XMLTagInfo eleminfo;
            int i=0;
            while(1) {
                ReadTag(is,eleminfo);
                if (eleminfo.isend && !eleminfo.isstart) {
                    if (eleminfo.name == namestr()) {
                        if (i!=n) v.resize(i);
                        return;
                    }
                    throw streamexception(std::string("Stream Input Format Error: expected end tag for ")+namestr()+", received end tag for "+eleminfo.name);
                }
                if (i>=n) v.resize(n=i+1);
                LoadWrapper(v[i++],eleminfo,is);
            }
        }
    };

    template<typename T,typename A>
    struct TypeInfo<std::vector<T,A>,
            typename Type_If<IsShiftable<T>::atall,void>::type> {
        inline static const char *namestr() {
            static char *ret = TName("vector",1,TypeInfo<T>::namestr());
            return ret;
        }
        inline static void writeotherattr(std::ostream &os, const std::vector<T,A> &v) {
            os << " nelem=" << v.size();
        }
        inline static bool isshort(const std::vector<T,A> &) { return false; }
        inline static bool isinline(const std::vector<T,A> &) { return false; }
        inline static void save(const std::vector<T,A> &v,
                std::ostream &os, int indent) {
            for(typename std::vector<T,A>::const_iterator i=v.begin();i!=v.end();++i)
                os << *i << ' ';
        }
        inline static void load(std::vector<T,A> &v, const XMLTagInfo &info,
                std::istream &is) {
            std::map<std::string,std::string>::const_iterator ni
                    = info.attr.find("nelem");
            if (ni == info.attr.end())
                throw streamexception("Stream Input Format Error: ivector needs nelem attribute");
            int n = atoi(ni->second.c_str());
            v.resize(n);
            for(int i=0;i<n;i++)
                is >> v[i];
            ReadEndTag(is,namestr());
        }
    };

    template<typename T,typename C,typename A>
    struct TypeInfo<std::set<T,C,A>,
            typename Type_If<!IsShiftable<T>::atall,void>::type> {
        inline static const char *namestr() { return "set"; }
        inline static void writeotherattr(std::ostream &,const std::set<T,C,A> &) { }
        inline static bool isshort(const std::set<T,C,A> &) { return false; }
        inline static bool isinline(const std::set<T,C,A> &) { return false; }
        inline static void save(const std::set<T,C,A> &s,
                std::ostream &os,int indent) {
            os << std::endl;
            int c=0;
            for(typename std::set<T,C,A>::const_iterator i=s.begin();
                i!=s.end();++i,++c) {
                SaveWrapper(*i,"",os,indent+1);
            }
            Indent(os,indent);
        }
        inline static void load(std::set<T,C,A> &s, const XMLTagInfo &info,
                    std::istream &is) {
            s.clear();
            XMLTagInfo eleminfo;
            int i=0;
            while(1) {
                ReadTag(is,eleminfo);
                if (eleminfo.isend && !eleminfo.isstart) {
                    if (eleminfo.name == namestr()) return;
                    throw streamexception(std::string("Stream Input Format Error: expected end tag for ")+namestr()+", received end tag for "+eleminfo.name);
                }
                T temp;
                LoadWrapper(temp,eleminfo,is);
                s.insert(temp);
            }
        }
    };

    template<typename T,typename C,typename A>
    struct TypeInfo<std::set<T,C,A>,
            typename Type_If<IsShiftable<T>::atall,void>::type> {
        inline static const char *namestr() {
            static char *ret = TName("set",1,TypeInfo<T>::namestr());
            return ret;
        }
        inline static void writeotherattr(std::ostream &os, const std::set<T,C,A> &v) {
            os << " nelem=" << v.size();
        }
        inline static bool isshort(const std::set<T,C,A> &) { return false; }
        inline static bool isinline(const std::set<T,C,A> &) { return false; }
        inline static void save(const std::set<T,C,A> &s,
                std::ostream &os, int indent) {
            for(typename std::set<T,C,A>::const_iterator i=s.begin();i!=s.end();++i)
                os << *i << ' ';
        }
        inline static void load(std::set<T,C,A> &s, const XMLTagInfo &info,
                std::istream &is) {
            std::map<std::string,std::string>::const_iterator ni
                    = info.attr.find("nelem");
            if (ni == info.attr.end())
                throw streamexception("Stream Input Format Error: ivector needs nelem attribute");
            int n = atoi(ni->second.c_str());
            s.clear();
            for(int i=0;i<n;i++) {
                T temp;
                is >> temp;
                s.insert(temp);
            }
            ReadEndTag(is,namestr());
        }
    };

    template<typename K,typename T,typename C, typename A>
    struct TypeInfo<std::map<K,T,C,A>, void> {
        inline static const char *namestr() { return "map"; }
        inline static void writeotherattr(std::ostream &,const std::map<K,T,C,A> &) { }
        inline static bool isshort(const std::map<K,T,C,A> &) { return false; }
        inline static bool isinline(const std::map<K,T,C,A> &) { return false; }
        inline static void save(const std::map<K,T,C,A> m,
                std::ostream &os,int indent) {
            os << std::endl;
            int c=0;
            for(typename std::map<K,T,C,A>::const_iterator i=m.begin();
                i!=m.end();++i,++c) {
                SaveWrapper(*i,"",os,indent+1,"key","value");
            }
            Indent(os,indent);
        }
        inline static void load(std::map<K,T,C,A> &m, const XMLTagInfo &info,
                            std::istream &is) {
            m.clear();
            XMLTagInfo eleminfo;
            int i=0;
            while(1) {
                ReadTag(is,eleminfo);
                if (eleminfo.isend && !eleminfo.isstart) {
                    if (eleminfo.name == namestr()) return;
                    throw streamexception(std::string("Stream Input Format Error: expected end tag for ")+namestr()+", received end tag for "+eleminfo.name);
                }
                std::pair<K,T> elem; // not as efficient as I would like
                TypeInfo<std::pair<K,T> >::load(elem,eleminfo,is,"key","value");
                m.insert(elem); // copying elements here (C++0x might help)
            }
        }
    };


    // how to save an item where G is a "getter" for value,type,name
    template<typename G>
    struct SaveItem {
        template<typename T>
        inline static void apply(T o,std::ostream &os,int indent) {
            SaveWrapper(G::getvalue(o),G::getname(o),os,indent);
        }
    };

    template<typename G>
    struct LoadItem {
        template<typename T>
        inline static bool apply(T o,std::istream &is,
                const XMLTagInfo &info) {
            std::map<std::string,std::string>::const_iterator ni
                =info.attr.find("name");
            if (ni!=info.attr.end() || ni->second != G::getname(o))
                return false;
            LoadWrapper(G::getvalue(o),info,is);
            return true;
        }
    };

    template<typename L>
    struct LoadList {
        template<typename O>
        inline static void exec(O o, std::istream &is, const char *cname) {
            XMLTagInfo info;
            std::set<std::string> loadedmem;
            while(1) {
                ReadTag(is,info);
                if (info.isend && !info.isstart) {
                    if (info.name != cname)
                        throw streamexception(std::string("Stream Input Format Error: expected end tag for ")+cname+", received end tag for "+info.name);
                    if (!AllLoaded<L>::exec(o,loadedmem))
                        throw streamexception(std::string("Stream Input Format Error: not all fields present for ")+cname);
                    return;
                }
                if (LoadOne<L>::exec(o,is,info))
                    loadedmem.insert(info.attr["name"]);
                else
                    throw streamexception(std::string("Extra field ")+info.attr["name"]+" of type "+info.name+" in object "+cname);
            }
        }
    };


/* Here starts "END" macros (only 1, actually) */

#define SERIAL_END \
        serial__tmp__alllist; \
    public: \
        typedef serial__tmp__alllist serial__alllist; \
        typedef serial__vallocT<SERIAL_BASETYPE,true> serial__valloc; \
        inline void serial__Save(std::ostream &os, int indent=0) const { \
            SERIALNAMESPACE::SaveItt<serial__alllist>::exec(this,os,indent); } \
        inline void Save(std::ostream &os, const char *name="", \
                        int indent=0) const { \
            SERIALNAMESPACE::SaveWrapper(*this,name,os,indent); \
        } \
        inline void serial__Load(std::istream &is) { \
            SERIALNAMESPACE::LoadList<serial__alllist>::exec \
                (this,is,serial_IDname()); \
        } \
        inline void Load(std::istream &is) { \
            SERIALNAMESPACE::LoadWrapper(*this,is); \
        } \
        inline static SERIAL_BASETYPE *LoadPtr(std::istream &is) { \
            SERIAL_BASETYPE *ret; \
            SERIALNAMESPACE::LoadWrapper(ret,is); \
            return ret; \
        }

}

#endif
