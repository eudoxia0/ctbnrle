How to cite:
============

If you wish to cite this code, please cite the JMLR paper

Christian R. Shelton, Yu Fan, William Lam, Joon Lee, and Jing Xu.
Continuous Time Bayesian Network Reasoning and Learning Engine.
Journal of Machine Learning Research. 11(Mar), 1137-1140, 2010.

a suitable bibtex entry would be

@article{CTBNRLE,
	author = "Christian R. Shelton and Yu Fan and William Lam and Joon Lee
and Jing Xu",
	title = "Continuous Time {B}ayesian Network Reasoning and Learning
Engine",
	journal = "Journal of Machine Learning Research",
	volume = 11,
	number = "Mar',
	pages = "1137--1140",
	year = 2010
}

How to compile:
===============
see the file "INSTALL" in the root directory

Structure:
==========

The base directory is split into multiple directorys:

demo/  This contains example programs on how to use the libraries.  They
	are often useful on their own (e.g. sample from a CTBN, learn a CTBN
	from samples, etc).  The README.txt file also has information on 
	the file extensions used for data files.

doc/   This directory contains information on how to read the code and 
	use it, literature to read, and funding sources that made this code
	base possible.

src/   This directory contains all of the non-header files.

hdr/   This directory contains all of the header files.  This, plus the ctbn.a
	created in /src are all that is necessary to link a new program to
	the libraries here

tst/   This directory contains regression tests that are used to verify that
	the code is working

eigen/ and glpk/ These directories contain the source trees for Eigen
	(a matrix package -- headers only) and GLPK (a linear program solver)

