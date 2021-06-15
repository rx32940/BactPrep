#!/bin/bash

# The following inputs are used:
#inputFile = './data3_many/alignment-rec.fa';
#outputFile = './data3_many/data3_many_res.mat';
#inputSpecsFile = './fG_input_specs.txt';
#fileWithOrdering = './data3_many/data3_given_ordering.txt'
#
# RUNNING: sh ./example_script.sh

##############################################################################
# Run fastGEAR
# 
# USAGE:
# ./run_fastGEAR.sh <path to matlab/mcr> <inputFile> <outputFile> <inputSpecsFile>
#
./run_fastGEAR.sh /m/fs/software/matlab/r2016a ./data3_many/alignment-rec.fa ./data3_many/data3_many_res.mat ./fG_input_specs.txt

#
#
# RESULTS as text can now be found in files "lineage_information.txt", "recombinations_recent.txt", and "recombinations_ancestral.txt" in the output directory which can be found in the same directory as the output file specified.



##############################################################################
# Plot recombinations in strains
#
# USAGE:
#./run_plotRecombinations.sh <path to matlab/mcr> <outputFile> <type (1=recent, 2=ancestral)> <strainOrder (1=original, 2=by cluster, or fileWithOrdering)>
#
# Plot recent recombinations using the original ordering (Note that the script will continue only after you close the figure window)
./run_plotRecombinations.sh /m/fs/software/matlab/r2016a ./data3_many/data3_many_res.mat 1 1
#
#
# Plot ancestral recombinations, order strains according to a user-specified ordering
./run_plotRecombinations.sh /m/fs/software/matlab/r2016a ./data3_many/data3_many_res.mat 2 ./data3_many/data3_given_ordering.txt



##############################################################################
#For reference, show colors used in the plots
#
# USAGE:
# ./run_plotColors.sh <path to matlab/mcr> <number of lineages, use - for the number detected)> <outputFile>
#
./run_plotColors.sh /m/fs/software/matlab/r2016a - ./data3_many/data3_many_res.mat


####################
#Plot marginal lineage probabilities for a strain. Note that this works only if in the input specifications "produce complete output" has been selected.
#
# USAGE:
# ./run_plotMarginalsForStrain.sh <path to matlab/mcr> <outputFile> <strainIndex> <donorLineage (if =0, then shows all lineages)>
#
#./run_plotMarginalsForStrain.sh /m/fs/software/matlab/r2016a ./data3_many/data3_many_res.mat 1 0
