#!/bin/bash

# example wrapper script for Multirun use of nLocusSim_v110
# suppose we want to iterate simulations over two variables:
# -S controls selection against melanistic morph in green patch
# -d controls proportion of host patches that are "dark niche"

# to loop in a certain increment, in BASH we can use the seq command:
# suppose we want to vary -S in increments of 0.05 from 0.05 to 0.25:
Svals=`seq 0.05 0.05 0.25`
# suppose we want to vary -d in increments of 0.1 from 0.1 to 0.4:
dvals=`seq 0.1 0.1 0.4`
# these increments are too large and ranges too small for a real
# parameter study, but it keeps this example from going too long:
# this should give us 5 * 4 = 20 parameter combos = 20 runs.

# if we want to use all other parameters as defaults, the only
# additional thing we need to do is remember to use the -w 1 option to run in "MULTIRUN" mode

for i in $Svals
do
	for j in $dvals
	do
		./nLocusSim_v110 -w 1 -S $i -d $j
	done
done

