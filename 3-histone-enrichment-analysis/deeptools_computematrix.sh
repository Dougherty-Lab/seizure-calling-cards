#!/bin/bash
#SBATCH --job-name=deeptools_peaks
#SBATCH --output=deeptools_compute%j.out
#SBATCH --error=deeptools_compute%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Load py-deeptools
eval $( spack load --sh py-deeptools )

# Define number of processors to use (adjust as needed)
NUM_CORES=8

# Compute values for the matrix using the original and shifted bed files and parallel processing
computeMatrix reference-point --referencePoint center \
    -S ENCFF422JMO.bigWig \
    -R peaks.bed \
    -b 10000 -a 10000 \
    -o matrix_center_fc_H3K27ac.out \
    -p $NUM_CORES

computeMatrix reference-point --referencePoint center \
    -S ENCFF742UNH.bigWig \
    -R peaks.bed \
    -b 10000 -a 10000 \
    -o matrix_center_fc_H3K27me3.out \
    -p $NUM_CORES

computeMatrix reference-point --referencePoint center \
    -S ENCFF326LFP.bigWig \
    -R peaks.bed \
    -b 10000 -a 10000 \
    -o matrix_center_fc_H3K4me1.out \
    -p $NUM_CORES

computeMatrix reference-point --referencePoint center \
    -S ENCFF422JMO.bigWig \
    -R peaks_shifted_20kb.bed \
    -b 10000 -a 10000 \
    -o matrix_center_fc_H3K27ac_20kb.out \
    -p $NUM_CORES

computeMatrix reference-point --referencePoint center \
    -S ENCFF742UNH.bigWig \
    -R peaks_shifted_20kb.bed \
    -b 10000 -a 10000 \
    -o matrix_center_fc_H3K27me3_20kb.out \
    -p $NUM_CORES

computeMatrix reference-point --referencePoint center \
    -S ENCFF326LFP.bigWig \
    -R peaks_shifted_20kb.bed \
    -b 10000 -a 10000 \
    -o matrix_center_fc_H3K4me1_20kb.out \
    -p $NUM_CORES

