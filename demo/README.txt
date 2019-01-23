The excutables in this folder give examples of how to use the 
CTBN-RLE libraries. CTBN models are usually loaded from files
in the excutables. The following is a quick description of the 
file formats and their extensions  

File Ext        Description
.ctbn           A CTBN model saved by method Save()
.ctbnptr        A CTBN model saved by method SavePtr()
.ctbndyn        A CTBNDyn model saved by method Save()
.ctbndynptr     A pointer of CTBNDyn model saved by method SavePtr()
.traj           A single trajectory saved by method Save()

The example files drug.ctbn, drug.ctnbptr, and queryinput.data are
created by the Makefile.  This is because, for this version of the code,
the file formats depend on the RTTI information supplied by the compiler,
which differs across compilers.


Demo executables
================

Constructing a network
----------------------
twonode: use as "twonode"
	Creates a simple two-node network:  The initial distribution BN has
	the structure 0->1 and the dynamics have the structure 0<->1.
	The parameters are arbitrarily chosen.  Good to read to understand
	how to build a CTBN.

makedrug:  use as "makedrug [options] > <filename>"
	Creates the drug network as used in Nodelman et al. 2002 and
	saves it to stdout. See makedrug.cc for details on optional 
	parameters.

makechain:  use as "makechain [option] > <filename>"
	Creates a chain structured network and saves it to stdout. See 
	makechain.cc for details on optional parameters.


Generating trajectories
-------------------------
sampleprocess:   use as "sampleprocess <endtime> < <ctbnfile>"
	Loads a process from stdin and outputs (to stdout) a sample from
	0 until time endtime.

gen_fulldata:   use as "gen_fulldata <.ctbn file> <output file> [option]"
	Generates complete trajectories from a given CTBN and saves the 
	trajectories in a file. The trajectories can be used as fully 
	observed data for learning CTBN parameters. See gen_fulldata.cc for 
	details on optional paramters.

gen_partialdata: use as "gen_partialdata [option] < .ctbnfile > outfile"
	Generates incomplete trajectories from a given CTBN and save the 
	trajectories in a file. The trajectories can be used as partially 
	observed data for learning CTBN parameters using EM algorithm. See 
	gen_partialdata.cc for details on optional patameters.



Queries
-------
exactquery: use as "exactquery"
	Load in model, evidence, true query answers from the file 
	"queryinput.data" and compute queries using exact Markov inference 
	and compare to the true answers. See exactquery.cc for details.

importancequery: use as "importancequery"
	Load in model, evidence, true query answers from the file 
	"queryinput.data" Compute queries using approximate inference via 
	importance sampler and compare to the true answers. See 
	importancequery.cc for details.

gibbsquery: use as "gibbsquery"
	Load in model, evidence, true query answers from the file 
	"queryinput.data" and compute queries using approximate inference via 
	gibbs sampler and compare to the true answers. See gibbsquery.cc for 
	details.

Note there are other inference methods also available.  See the tst 
directory for other examples.


Parameter learning
------------------
learnparams:    use as "learnparams <.ctbn file> <datafile>"
	Demonstates parameter learning of a given CTBN with complete data. 
	The algorithm was described in Nodelman et al., 2003. The data can be
	generated using gen_fulldata.

emlearn:    use as "emlearn <.ctbn file> <datafile>"
	Demonstate parameter learning of a given CTBN with incomplete data.
	The algorithm was described in Nodelman et al., 2005. The data can be
	generated using gen_partialdata.


Structure learning
------------------
structurelearn:  use as "structurelearn"
	Demonstrates structure learning on the drug network used in Nodelman
	et al. 2002. Reports on the difference in number of parameters for 
	the learned BN and Hamming distance for the learned CTBN. See 
	structurelearn.cc for details on optional parameters.

sem:  use as "sem"
	Demonstrates structural EM on the drug network used in Nodelman
	et al. 2002. Reports on the difference in number of parameters for 
	the learned BN and Hamming distance for the learned CTBN. See sem.cc
	for details on optional parameters.


Constructing a clique tree
--------------------------
clique_example: use as "clique_example <.ctbn file>"
	Generates a clique tree from a CTBN.


Visualization
-------------
traj2ascii: use as "traj2ascii < .trajfile"
	Takes a trajectory (say, produced from sampleprocess) and tries to
	make a nice "ASCII" picture of the trajectory.

using_graphic:  use as "using_graphic"
	Produces a .dot file for graphviz from a CTBN
	Draws a trajectory as a postscript file.
	see using_graphic.cc for details on optional parameters.
