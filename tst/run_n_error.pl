#!/usr/bin/perl -w

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


sub main {
	#system("rm -fr bhps_query_results");
	#system("mkdir bhps_query_results");
	@nsample = (10,20,50,100,200,500,1000); #,2000,5000,10000,20000,50000,100000,200000,500000);
	for (my $i = 0; $i <= $#nsample; $i++) {
		system("./approxinf drug.ctbn drug.traj_q.txt $nsample[$i] >drug_query_results/gibbs.$nsample[$i].txt &");
	}
}

&main;
