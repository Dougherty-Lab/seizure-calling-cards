#!/bin/bash
#SBATCH --job-name=deeptools_plot
#SBATCH --output=deeptools_plot%j.out
#SBATCH --error=deeptools_plot%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load py-deeptools
eval $( spack load --sh py-deeptools )

# Create the heatmaps
plotHeatmap -m matrix_center_fc_H3K27ac.out \
    -out center_fc_H3K27ac.pdf \
    --heatmapHeight 5.8 \
    --regionsLabel "" \
    --plotTitle "" \
    --xAxisLabel "" \
    --yAxisLabel "" \
    --yMin 0.7 \
    --yMax 2.7 \
    --legendLocation none \
    --missingDataColor "gray" \
    --colorMap RdBu_r \
    --zMin 0 \
    --zMax 7

plotHeatmap -m matrix_center_fc_H3K4me1.out \
    -out center_fc_H3K4me1.pdf \
    --heatmapHeight 5.8 \
    --regionsLabel "" \
    --plotTitle "" \
    --xAxisLabel "" \
    --yAxisLabel "" \
    --yMin 0.7 \
    --yMax 2.7 \
    --legendLocation none \
    --missingDataColor "gray" \
    --colorMap RdBu_r \
    --zMin 0 \
    --zMax 7

plotHeatmap -m matrix_center_fc_H3K27me3.out \
    -out center_fc_H3K27me3.pdf \
    --heatmapHeight 5.8 \
    --regionsLabel "" \
    --plotTitle "" \
    --xAxisLabel "" \
    --yAxisLabel "" \
    --yMin 0.7 \
    --yMax 2.7 \
    --legendLocation none \
    --missingDataColor "gray" \
    --colorMap RdBu_r \
    --zMin 0 \
    --zMax 7

plotHeatmap -m matrix_center_fc_H3K27ac_20kb.out \
    -out center_fc_H3K27ac_20kb.pdf \
    --heatmapHeight 5.8 \
    --regionsLabel "" \
    --plotTitle "" \
    --xAxisLabel "" \
    --yAxisLabel "" \
    --yMin 0.7 \
    --yMax 2.7 \
    --legendLocation none \
    --missingDataColor "gray" \
    --colorMap RdBu_r \
    --zMin 0 \
    --zMax 7

plotHeatmap -m matrix_center_fc_H3K4me1_20kb.out \
    -out center_fc_H3K4me1_20kb.pdf \
    --heatmapHeight 5.8 \
    --regionsLabel "" \
    --plotTitle "" \
    --xAxisLabel "" \
    --yAxisLabel "" \
    --yMin 0.7 \
    --yMax 2.7 \
    --legendLocation none \
    --missingDataColor "gray" \
    --colorMap RdBu_r \
    --zMin 0 \
    --zMax 7

plotHeatmap -m matrix_center_fc_H3K27me3_20kb.out \
    -out center_fc_H3K27me3_20kb.pdf \
    --heatmapHeight 5.8 \
    --regionsLabel "" \
    --plotTitle "" \
    --xAxisLabel "" \
    --yAxisLabel "" \
    --yMin 0.7 \
    --yMax 2.7 \
    --legendLocation none \
    --missingDataColor "gray" \
    --colorMap RdBu_r \
    --zMin 0 \
    --zMax 7
