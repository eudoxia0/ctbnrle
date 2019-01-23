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
#ifndef CTBNRLE_DEFINES_H
#define CTBNRLE_DEFINES_H


/* This file describes changes that need to be made for compilation
 * on different platforms.  
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <climits>
#include <cfloat>

typedef long long int longlong;

#if defined __linux__ || __CYGWIN__
// Linux/Cygwin section...
#define TIME_STRUCT longlong
#define GET_TIME    GetLinuxTime()
TIME_STRUCT GET_TIME;

#elif defined _WIN32
// Windows section...
#include <limits>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN 
#include <windows.h>

#ifdef UNION
#undef UNION
#endif

#ifdef INTERSECTION
#undef INTERSECTION
#endif

#ifdef DIFFERENCE
#undef DIFFERENCE
#endif 


#define TIME_STRUCT longlong
#define GET_TIME GetWindowsTime()
/*
#define TIME_STRUCT longlong
#define GET_TIME	GetTime()
*/
TIME_STRUCT GET_TIME;
#else
// Mac OS section...

//   #define TIME_STRUCT hrtime_t
//   #define GET_TIME    gethrtime()
// for mac os x
#include <mach/mach_time.h>
#define GET_TIME    mach_absolute_time()
#endif


#if defined(_WIN32)
// more Windows section...
#   ifndef finite
#       define finite(x) (_finite (x))
#   endif

#   ifndef isfinite 
#       define isfinite(x) (_finite (x))
#   endif

#   ifndef isnan
#       define isnan(d) (_isnan(d))
#   endif

#   ifndef __func__
#       define __func__ __FUNCTION__
#   endif

#   ifndef isinf 
//#     define isinf(x) (!_finite(x) && !_isnan(x)) 
#       inline bool isinf(double x) { return !_finite(x) && !_isnan(x); } 
#   endif

#   ifndef isinfinite
#       inline bool isinfinite(double x) { return isinf(x); }
#   endif

#   ifndef INFINITY
#       define INFINITY (std::numeric_limits<double>::infinity())
#   endif

    extern const unsigned long win32nan[2];
    double nan(const char *);
    float nanf(const char *);
    long double nanl(const char *);
#else
// non-Windows solution to isfinite being C99 and not C++
// should be changed when C++0x is adopted
#   include <cmath>
#   ifdef __INTEL_COMPILER
#       include <mathinf.h>
        inline bool isinfinite(double x) { return isinf(x); }
//        inline bool isfinite(double x) { return !isinf(x) && !isnan(x); }
        inline bool finite(double x) { return !isinf(x); }
#   else
        inline bool isinfinite(double x) { return std::isinf(x); }
 //       inline bool isfinite(double x) { return !std::isinf(x) && !std::isnan(x); }
#   endif
#endif

#endif // _DEFINES_H
