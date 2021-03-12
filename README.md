# The Natural History of Clonal Haematopoiesis

## Motivation

This is the code required to run 

## Code map

### Requirements

#### Software

* `clonex` - a multi-purpose and efficient implementation of Fisher-Wright simulations. Available [here](https://github.com/josegcpa/clonex)
* `R` 3.6.3 - other versions will probably work
* `Rstudio`

#### R library

`reticulate`,`greta`,`ggplot2`,`ggpubr`,`tidyverse`,`bayesplot`,`openxlsx`,`gtools`,`cowplot`,`ggsci`,`ggrepel`,`extraDistr`,`default`,`dendextend`,`ape`,`ggtree`,`grid`,`scatterpie`,`reghelper`,`phylodyn`

### Notebooks

*Notebooks with complete runs are already present in [LINK TO UPDATE](), but if one wants to run this locally the following steps can be taken:*

*Please note that this entails installing your own version of `greta`, which is the package used for MCMC sampling, and the adequate alteration of paths in `vaf_dynamics_functions.R`*

*Please note this assumes that you are running this script from the root directory of this project*

1. Technical overdispersion estimation
    a. Run `Notebook_Overdispersion.Rmd` - R notebook containing the overdispersion estimation from technical replicates.

2. Method validation
    a. Run `simulate_range.sh` (`sh simulate_range.sh`) - this will run a set of Fisher-Wright simulations with different driver fitness advantages and mutation rates. As it is, the script will submit jobs to a LSF job scheduler - if no such job scheduler is available, one should adjust accordingly by removing the line containing `bsub`.
    b. Run `Notebook_Simulations.Rmd` (in `Rstudio`) - R notebook containing the method validation using Fisher-Wright simulations

3. Data analysis - growth rate coefficient and age at onset inference, possible associations with phenotype
    a. Run `Scripts/run_model.R` (`Rscript Scripts/run_model.R`) - this runs the model that infers all growth rate coefficients
    b. Run `Notebook_GrowthCoefficients_AgeAtOnset_PossibleAssociations.Rmd` - this is the notebook containing the bulk of the analysis

4. Analysis of the historical growth effect and poor fits
    a. Run `Notebook_HistoricalGrowth_PoorFits.Rmd` - this runs analyses which factors - technical and biological - may be determinant of the difference between historical and inferred growth and a fit being poor (having one or more outlier)

5. Phylogenetic trees from single cell colines, the determination of growth per year and age at onset from these trees and comparison with estimates from longitudinal data
    a. Run `Scripts/plot_tree.R` (`Rscript Scripts/plot_tree.R`) - this runs all of the aforementioned analysis and plots it as displayed in the manuscript

*Optional: run `Scripts/plots_for_initial_panel.R`*
