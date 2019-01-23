#! /usr/bin/python

# Continuous Time Bayesian Network Reasoning and Learning Engine
# Copyright (C) 2009 The Regents of the University of California
#
# Authored by Yu Fan, William Lam, Joon Lee, Christian Shelton, and Jing Xu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os;


# execute the Makefile
os.system("make")

# Run the EXES the file from EXES= in Makefile

#open the Makefile file to read
filename = "Makefile"
FILE = open(filename, "r")

#start the loop
while 1:
	#read one line from the file
	sLine = FILE.readline()
	#check if the EXES= is in the line just being passed in
	if sLine.startswith("EXES="):
		while 1:
			#seperate the every word by space
			thisline = sLine.split(" ")
			#go through the array of string
			for element in thisline:
				# if it isn't EXES= run those file 
				if "EXES=" not in element:
					if '\\' not in element:
						if " " not in element:
							#take out new line
							oPipe = os.popen("./" + element.strip('\n'))
							#start other loop
							while 1:
								#read the out put of the file
								osLine = oPipe.read()
								#check if pass is in the output first line 
								if "PASS" in osLine: 
									print "\n" + element.strip('\n') + " has PASSED."
									break
								else:
									print "\n" + element.strip('\n') + " has FAILED."
									print osLine
									#output everything form the program
									for osLine in oPipe.readlines():
										print osLine
									break
			#read one more time if the sLine end with \\
			if '\\' not in sLine:
				sline = ""
				break
			else:
				sLine = FILE.readline()
	#same as eof()
   	if len(sLine) is 0:
		break

#new line
print "\n"
#os.system("make clean")
