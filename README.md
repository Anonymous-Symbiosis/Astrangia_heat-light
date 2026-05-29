# Astrangia_heat-light
This repository contains the data and code associated with the manuscript: "Light environment shapes divergent coral thermal responses across symbiotic states."

## Organization
This repository is organized into three folders: scripts, data, and figures.

The scripts folder contains all code used for data processing, analysis, and figure generation. This includes scripts for:
1. Experimental design visualization (Apoculata_ExpDesign.R)
2. Physiological analyses (Apoculata_Physiology.R, Apoculata_Phys_PCA.R)
3. Tag-seq read mapping (read_mapping.sh)
4. Gene expression analyses (Apoculata_DESeq.R, Apoculata_DEGs.R, Apoculata_GE_PCA.R, Apoculata_Baseline.R, Apoculata_VennDiagrams.R, Apoculata_ComparativeAnalysis.R, Apoculata_GOBubblePlots.R, Apoculata_WGCNA.R, Apoculata_Heatmaps.R).

Each script includes comments at the beginning describing its role in the workflow, and all required R functions are included within the scripts folder.

The data folder contains all input files required to run the analyses. All raw sequencing data are deposited in NCBI's SRA (BioProject PRJNA1472096).

The figures folder includes all final figures presented in the manuscript, along with supplementary figures.
