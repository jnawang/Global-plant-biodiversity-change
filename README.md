# Description of the code function for modeling global plant biodiversity change
1. species distribution model: Run_SDMs_seq.R and MadMax.R
2. predict species-specific migration velocity: predict_shift_rate.R
3. predict realized species range: overlay_colonizable_suitable2023.R
4. calculate global species richness: cal_richness_more.R
5. calculate global beta diversity: beta_diversity.R
6. calculate global phylogenetic diversity: phylo_diversity.R
7. calculate non-analog community: map2bin.R readbin.cpp  species_pair_matrix.cpp species_pair_overlap.cpp species_pair_overlap_chunk.cpp novel_species_pair.cpp novel_species_pair_chunk.cpp
# System requirements
1. Most analyses were performed on R version 4.2.1 except the analyses of the global distribution of novel communities and the calculation of species co-occurrence strength that are performed using cpp in order to improve computational efficiency.
2. The code requires an operation system that can run R version 4.2.1 or later, and a cpp compiler such as gcc.
# Installation guide
1. to run the code, the user needs to install a series of R papckages such as terra, raster, SDMtune, MuMin, lme4, randomForest, V.PhyloMaker, PVR, Rphylopars and so on.
2. Installation of these R packages is typically fast, and can be finished within a few minutes.
# Demo
1. the code can be run easily with a R compiler or cpp compiler.
# Instructions for use
1. The code has to be run in the order of the code function description or following the description in the Materials and Methods section of this paper.
