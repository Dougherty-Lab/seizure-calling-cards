#!/usr/bin/env bash

# USAGE:
# Make sure you first make all changes to versions/paths
# noted in the script
#
# run_1234.csv is the input samplesheet
# sbatch run_nf_callingcards.sh run_1234.csv

#SBATCH --job-name=nf-callingcards
#SBATCH --mem-per-cpu=5G
#SBATCH --output=nf-callingcards_%j.out
#SBATCH --error=nf-callingcards_%j.err
#SBATCH --mail-type=ALL

# load system dependencies -- on HTCF, we use spack
eval $(spack load --sh singularityce@3.8.0)
eval $(spack load --sh nextflow@23.04.4)

tmp=./tmp$(mktemp -d /tmp/$USER-singularity-XXXXXX)

mkdir local_tmp
mkdir tmp

export NXF_HOME=/path/to/your/nextflow/directory #change this to your nextflow directory
export NXF_SINGULARITY_CACHEDIR=singularity_images
export SINGULARITY_TMPDIR=$tmp
export SINGULARITY_CACHEDIR=$tmp

config=wustl_htcf.config # Change this to the correct path of your config file, which may be different for your computing environment
params=params.json # Change this to the path to your params.json

# The first cmd line input to this script is expected to be
# the samplesheet csv. The samplesheet's basename is used to
# create the outdir and workdir directory names.
# For example, if your samplesheet is called
# run_1234.csv, then the results directory will be
# run_1234_results, and the work directory will be
# run_1234_work. Once the workflow completes, you
# can and should delete the work directory
input=$1
input_basename=$(basename $input)
outdir=${input_basename%.csv}_results
workdir=${input_basename%.csv}_work

# The following is the nextflow run
# command. `-resume` may be left on,
# even if this hasn't been run before.
# if the workflow exits unexpectedly,
# as long as you don't delete the _work
# and .nextflow and .nextflow.log directories/files
# in your launch directory, you can resubmit
# this same script and the pipeline will
# resume from where it left off
nextflow pull nf-core/callingcards -r dev
nextflow run nf-core/callingcards \
    -revision dev \
    -profile default_mammals,singularity \
    -c $config \
    -params-file $params \
    --input $input \
    --outdir $outdir \
    -work-dir $workdir \
    -resume
