17 September 2014
Document written by Samuel Flaxman, samuel.flaxman@colorado.edu
revised 14 May 2015

README: A summary of the nLocusSim program, written initially for version 1.0.2 (most current C source code now in file nLocusSim_v110.c).  

ATTRIBUTION:  
+ The source code (nLocusSim_vXXX.c) was written by Samuel Flaxman.  
+ The .R wrappers were written by Aaron Comeault.  
+ The random number generation code in the "MT" directory is from: Saito, M., and Matsumoto, M. (2006). SIMD-oriented Fast Mersenne Twister: a 128-bit pseudorandom number generator. In Monte Carlo and Quasi-Monte Carlo Methods, A. Keller, S. Heinrich, and H. Niederreiter, eds. (Heidelberg, Berlin: Springer-Berlin), pp. 607–622.  The version used here is dSMFT v2.1.  See: http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/ for versions, downloads, and descriptions of the random number generator code.

WHAT THE SOURCE CODE IS FOR:
This github repo provides the source code that was used for simulations that appear as part of a scientific publication on the population genetics of adaptation walking stick insects.  The bibliographic citation is:
A. A. Comeault et al.  2015.  Selection on a genetic polymorphism counteracts ecological speciation in a stick insect.  Current Biology [volume, pages, and DOI forthcoming]

The full dataset used in that publication is archived at datadryad.org [DOI forthcoming].


KEY UPDATES FROM PREVIOUS VERSIONS:
+ MULTIRUN mode (see V. below)
+ Some command line options have been added and others changed!  Read them carefully after you compile by running ./nLocusSim_v110 -?
+ Some new files are produced; see below
+ names of some parameters have been changed for clarity 


I.  Description
	A.  Model: 
This is an individual-based model of population genetic dynamics in Timema.  Multiple loci are considered, each with two alleles.  Alleles are coded as zeros and ones in a large "genotypes" array.  Recombination, migration, and differential reproduction are all stochastic.    Migration is NOT conditional in any way except that melanistic and green individuals may have different migration rates.  The program is initialized in a state of maximum divergence at pattern, and neutral loci.  At these loci, the individuals on Adenostoma hosts have "0" alleles, and the individuals on Ceanothus hosts have "1" alleles.  For the color locus, the starting frequency is a parameter, and for the sex locus, we assume a 50:50 sex ratio.

	B.  A note on indexing in data arrays in C:  
Numbered indexes (e.g., individuals, loci, deme #, etc) all start with 0 in C, not with 1.  So the first element of an array is element #0 (not #1).  This should help with understanding some of the outputs in data files and comments below.  This may mean that indexes output from the program have to be incremented by 1 to be usable in common analysis programs such as R.

	C.  Loci and alleles:
		1.  The "genome" setup is hard coded.  There are 6 loci:
	LocusIDnumber	LocusTypeCode	LocusType	RecombProbabilityWithPrevLocus
	0		0		NEUTRAL_LOCUS	NA
	1		0		NEUTRAL_LOCUS	0.5
	2		1		COLOR_LOCUS	0.001
	3		2		PATTERN_LOCUS	COLOR_PATTERN_RECOMB_PROBABILITY
	4		0		NEUTRAL_LOCUS	0.001
	5		3		SEX_LOCUS	0.5

		2.  Note: "LocusTypeCode" is an integer coding used internally in the program and spit out some places.  The column "LocusType" in the above table is how these codes should be decoded.
		3.  In this scheme, locus #0 is an unlinked neutral locus, locus #1 is a neutral locus tightly linked to the color locus (locus #2), locus #4 is a neutral locus tightly linked to the pattern locus (locus #3), and locus #5 is like a sex chromosome for an XX-XO system.  We assume here that none of the other loci are sex-linked.  COLOR_PATTERN_RECOMB_PROBABILITY is a parameter that can be altered on the command line with the option "-r" .

	D.  Random number generation:
Pseudorandom numbers are generated using the Mersenne Twister; the source code for this is in the directory called MT.  The citation is: Saito, M. and M. Matsumoto. 2006. SIMD-oriented Fast Mersenne Twister: a 128-bit pseudorandom number generator. Pp. 607-622 in A. Keller, S. Heinrich, and H. Niederreiter, eds. Monte Carlo and Quasi-Monte Carlo Methods. Springer Berlin Heidelberg, Berlin.

II. Compiling the program.
	A. Requirements:  a UNIX-style terminal and gcc (the GNU C compiler) and standard C libraries
	B. Instructions:  [ NOTE: instructions beginning with "$" denote UNIX commands to issue on the command line in your terminal ]
		1.  download the source (.c), Makefile, and MT directory all to the same directory
		2.  Open a UNIX terminal and cd to the directory
		3.  $ make
		4.  That should make an executable.  If it doesn't, contact Sam.

III.  Running the program
	A.  Requirements:  UNIX terminal, successfully compiled executable, RnumSeed.txt file present in same directory as executable
	B.  Instructions:
		1.  cd to the directory with the executable
		2.  COMMAND LINE OPTIONS: to see a list of all available command line options for the program, at your command prompt, type:
			$ ./nLocusSim_v110 -? 
		3.  read the options, their usage, and the defaults.  For example, to run the simulation with a population size of 10000, which is NOT the default, you would do:  
			$ ./nLocusSim_v110 -N 10000
		4.  Run the program once with defaults by typing:
			$ ./nLocusSim_v110
		5.  With the defaults, it should take ~10 seconds to complete.
		6.  Your directory should now be populated with 10 new data files (.csv) as well as one .R file with useful parameter values and codes.

IV.  Data outputs from the program in regular mode.
	A.  File types:  
		1.  most data are plain text, comma-delimited (.csv) files.  By default, these have text headers.  Headers can be suppressed with the -h 0 option (see III.B.2.)
		2.  Files with info that can be sourced in R are plain text R files with R-script-style formatting.  These files include comments about what they have; you should be able to look at those files in the RStudio program OR in any text editor.  
	B.  Contents of files: the data files are time series of various divergence metrics
		1.  "AlleleFreqTimeSeries.csv" :  First column is the time step; subsequent columns are allele frequencies.  The allele frequencies for a given locus (see I.C. above) are the frequencies of the "1" allele at that locus.  So, the frequency of the "0" allele is one minus this.  These are global allele frequencies (over all demes combined).
		2.  "AlleleFreqByPatch.csv" :  Same as the previous one, but these are frequencies in each deme (patch).  Headers tell which deme and which locus.
		3.  "FSTtimeSeries.csv" :  Like "AlleleFreqTimeSeries.csv", but FST values instead of allele frequencies.
		4 & 5.  "InitialPopulation.csv" and "FinalPopulation.csv" :  The genotypes, phenotypes, patch (deme) locations, and coordinates of all individuals at the beginning and end (respectively) of simulations.  Coordinates are in (continuous) Euclidian space, where 0,0 corresponds to the upper left of the rectangular environment.
		6.  "LDtimeSeries.csv" : Like "AlleleFreqTimeSeries.csv", with LD values between the two loci.
		7.  "MigrationRecord.csv" : Overall proportion of individuals changing patches (demes) each time step.
		8.  "PatchAbundances.csv" : Numbers of individuals in each patch at each time
		9.  "GenomeScheme.csv" : gives the table shown above in I.C.
		10.  "InterHostMatingTS.csv" :  gives a time series of the frequency of matings between males and females that originated on different host types
		11.  "parameters.R" : an R-readable script that contains information about states of all parameters that affected the run and the data recording.

V.  Running program in MULTIRUN mode.
	A.  Description: 
		New to this version, the program is set up to be iterated multiple times (in sequence) directly from the command line or with a wrapper.  The basic idea is that final states from multiple independent simulations, along with parameter values of each, are stored in a single .csv file. This should facilitate analysis and parameter studies.  If the master data file gets too big, it should be easy to filter it with standard UNIX command line tools (awk, grep, etc.).
	B.  Running in MULTIRUN mode:  run the program with the option "-w 1", e.g., 
		$ ./nLocusSim_v110 -w 1 
		(note that all other usual options work as well)
	C.  Files produced from MULTIRUN mode.  Only 3 files are produced:
		1.  "LastRunParameters.R" : a record of only the parameters from the last run that was completed. This is probably redundant and unnecessary
		2.  "MultirunEndpointData.csv" : this is the master data file.  It contains columns detailing parameters followed by columns detailing results, i.e., info about the state of the system at the end of a given simulation run.
		3.  "MultirunCount.txt" : this is a file used by the program simply to index individual runs.  If this is deleted from the directory, the program will start counting at "1" again but will not erase any previous data.
		4.  RnumSeed.txt is automatically modified to ensure that independent runs have different random number seeds.
	D.  Multirun mode never over-writes the data in "MultirunEndpointData.csv", but rather simply appends new data to it.
	E.  Wrapping: 
		1.  An very simple example wrapper script is provided ("exampleWrapper.sh").  This wrapper is a very simple Bash shell script.  After the code has been compiled (above), this script or anything similar should work.  To use this wrapper:
			$ sh exampleWrapper.sh
		2.  Wrappers we actually used to generate datasets are provided as two R scripts: "model_nLocusSim_v110_simulations.R" and "model_run_nLocusSim_v108_simulations.R".  These are the wrappers that generated the simulations for the Current Biology publication cited above.  The full data are archived at Dryad (also cited above).
	F.  Monitoring progress of a large number of runs: in MULTIRUN mode, the program prints messages about when it starts and finishes a run to stderr (i.e., prints messages to screen unless redirected).

