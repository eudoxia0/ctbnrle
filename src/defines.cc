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
#include "defines.h"

#if defined __linux__ || __CYGWIN__
#include <sys/time.h>

TIME_STRUCT GetLinuxTime() {

	struct timeval t;

	struct timezone z;
	gettimeofday(&t, &z);
	return (TIME_STRUCT)t.tv_sec * (TIME_STRUCT)1000000000 +
	       (TIME_STRUCT)t.tv_usec * (TIME_STRUCT)1000;
}


#endif

#if defined(_WIN32)

longlong GetWindowsTime() {
     SYSTEMTIME syst;
     GetSystemTime(&syst);
     FILETIME filet;
     SystemTimeToFileTime(&syst,&filet);
     ULARGE_INTEGER larget;
     larget.LowPart = filet.dwLowDateTime;
     larget.HighPart = filet.dwHighDateTime;

     return larget.QuadPart/10;
}

const unsigned long win32nan[2]={0xffffff, 0x7ffffff};
double nan(const char *) { return *((double *)win32nan);}
float nanf(const char *) { return (float)(*((double *)win32nan));}
long double nanl(const char *) { return (long double)(*((double *)win32nan));}

#endif
