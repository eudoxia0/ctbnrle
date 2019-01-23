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
#ifndef CTBNRLE_PARAMS_H
#define CTBNRLE_PARAMS_H

#include <string>



namespace ctbn {

// This file creates a method for setting global parameters in a variety
// of ways and preventing the values from needing to be passed all around
// if they are more-or-less "assumed."

// To make this work, the function "InitParams" (below) should be the
// first thing called by the main function.  This will initialize the
// parameters (argc and argv should be the corresponding arguments to 
// the main function).
void InitParams(int argc, char **argv);

// Then, whenever the program needs a parameter, it can call one of the
// three functions immediately following.  The first argument is the
// name of the parameter.  The second value is the value that should be
// returned in case the parameter was not specified.
std::string ParamStr(const std::string &name, const std::string &dflt="");
int ParamInt(const std::string &name, int dflt=0);
double ParamDouble(const std::string &name, double dflt=0.0);

// How this works:
//  InitParams sets up a global variable (see top of .cc file) who scope
//  is hidden from all except the functions above.
//  It reads parameters from two places: parameter files and the command
//  line.  Parameter files take the form of one parameter per line.  
//  Lines beginning with "#" or ";" are ignored.  Whitespace at the
//  beginning of a line is ignored.  The first token on the line
//  (delimited by whitespace) is the name of the parameter.  The remainder
//  of the line (after the separating whitespace) is the parameter value.
//  
//  So, the file
//  #-----------
//  parameter1  3.0
//    par two words
//  #do not use
//  q 4
//  #-----------
//  defines three parameters:
//     "parameter1" => 3.0
//     "par" => "two words"
//     "q" => 4
//
//  InitParams reads the parameters in in the following order with more
//  recent values overriding older values:
//  1.  Parameter file "params" (in the current directory)
//  2.  Parameter file "<exe>.params" (where <exe> is the executed name of
//       the program)
//  3.  Command line
//
//  On the command line, -D<paramname> <value> defines paramname to have
//  value value.  The argument "-param <filename>" reads in file filename.
//  And the argument "-listparams" will cause the program to read in all
//  of the parameters, list them all to stdout, and then quit 
//  with error code 0.

// This sets a parameter.  Use with caution as it slightly destroys
// the value of *user* set parameters.  Best probably if set before
// the call to InitParams so that user values can override this value
// (but this value overrides the default given by the associated "Param*"
//  call)
void SetParam(const std::string &name, const std::string &val);

} // end of ctbn namespace

#endif

