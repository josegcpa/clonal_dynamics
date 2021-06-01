# The Natural History of Clonal Haematopoiesis

## Motivation

This is the repository for [The Natural History of Clonal Haematopoiesis](). In this work, we investigate the fitness differences between genes, sites and individual clones using longitudinal sequencing data from >300 individuals between 54 and 103 years old.

## Code map

### Requirements

#### Software

* `R` 3.6.3 - other versions may work, but this was the version it was tested on
* `Rstudio` - to run the R notebooks
* `clonex` - a multi-purpose and efficient implementation of Fisher-Wright simulations. Available [here](https://github.com/josegcpa/clonex)

#### R library

`reticulate`,`greta`,`ggplot2`,`ggpubr`,`tidyverse`,`bayesplot`,`openxlsx`,`gtools`,`cowplot`,`ggsci`,`ggrepel`,`extraDistr`,`default`,`dendextend`,`ape`,`ggtree`,`grid`,`scatterpie`,`reghelper`,`phylodyn`,`survival`,`survminer`

### Accessing the analysis

*Being served by `htmlpreview.github.io`*

1. [Estimating the technical overdispersion](https://htmlpreview.github.io/?ttps://github.com/josegcpa/vaf_dynamics/blob/master/Notebooks/Notebook_Overdispersion.html)
2. [Validating the approach to estimate annual growth from longitudinal data with Wright-Fisher simulations](https://htmlpreview.github.io/?https://github.com/josegcpa/vaf_dynamics/blob/master/Notebooks/Notebook_Simulations.html)
3. [Validating the approach to estimate annual growth from single-cell phylogenies with Wright-Fisher simulations](https://htmlpreview.github.io/?https://github.com/josegcpa/vaf_dynamics/blob/master/Notebooks/Notebook_BNPRFit.html)
4. [Growth rate coefficients and age at onset inference, possible associations with phenotype](https://htmlpreview.github.io/?https://github.com/josegcpa/vaf_dynamics/blob/master/Notebooks/Notebook_GrowthCoefficients_AgeAtOnset_PossibleAssociations.html)
5. [Analysing the trees obtained from single cell colonies in old individuals in Mitchell et al. (2021)](https://htmlpreview.github.io/?https://github.com/josegcpa/vaf_dynamics/blob/master/Notebooks/Notebook_Mitchell.html)
6. [Investigating the historical growth effect and poor fits](https://htmlpreview.github.io/?https://github.com/josegcpa/vaf_dynamics/blob/master/Notebooks/Notebook_HistoricalGrowth_PoorFits.html)

### Running the analysis

*Please note that this entails installing your own version of `greta`, which is the package used for MCMC sampling, and the adequate alteration of paths in `vaf_dynamics_functions.R`*

*Please note this assumes that you are running this script from the root directory of this project and that steps where the instructions mention running "in `Rstudio`" can be replaced by `./knit NOTEBOOK.Rmd` where `NOTEBOOK.Rmd` is the relevant Rmarkdown notebook. This generates `NOTEBOOK.html` and output files which are stored in `data_output`*

*Notebooks with complete runs are already present above, but if one wants to run this locally the following steps can be taken:*

1. Technical overdispersion estimation
    1. Run `Notebook_Overdispersion.Rmd` (in `Rstudio`) - R notebook containing the overdispersion estimation from technical replicates.

2. Longitudinal modelling validation (please note that, due to the stochastic nature of the simulations, results may differ slightly)
    1. Run `simulate_range.sh` (`sh simulate_range.sh`) - this will run a set of Fisher-Wright simulations with different driver fitness advantages and mutation rates. The `CLONEX_PATH` should be updated. As it is, the script will submit jobs to a LSF job scheduler - if no such job scheduler is available, one should adjust accordingly by removing the line containing `bsub`.
    2. Run simulations - this will run the model that uses the Fisher-Wright simulations to validate our approach
        * Run `Scripts/run_simulation_50k_1.R` (`Rscript Scripts/run_simulation_50k_1.R`)
        * Run `Scripts/run_simulation_100k_5.R` (`Rscript Scripts/run_simulation_100k_5.R`)
        * Run `Scripts/run_simulation_200k_200.R` (`Rscript Scripts/run_simulation_200k_13.R`)
    3. Run `Notebook_Simulations.Rmd` (in `Rstudio`) - R notebook containing the method validation using Fisher-Wright simulations
    4. Additional validation (regarding estimation using early and late parts of the trajectory and the effect of clonal competition on inference can also be done)
        1. Run `simulate_range_2.sh` (`sh simulate_range_2.sh`)
        2. Run `Scripts/investigate_simulations_early_late.R` (`Rscript Scripts/investigate_simulations_early_late.R`)
        3. Run `Scripts/investigate_simulations_competition.R` (`Rscript Scripts/investigate_simulations_competition.R`)

3. Phylogenetic population size and annual growth rates validation
    1. Run `simulate_few_complete.sh` (`sh simulate_few_complete.sh`)
    2. Run `Notebook_BNPRFit.Rmd` (in `Rstudio`)

4. Data analysis - growth rate coefficient and age at onset inference, possible associations with phenotype
    1. Run `Scripts/run_model.R` (`Rscript Scripts/run_model.R`) - this runs the model that infers all growth rate coefficients
    2. Run `Notebook_GrowthCoefficients_AgeAtOnset_PossibleAssociations.Rmd` (in `Rstudio`) - this is the notebook containing the bulk of the analysis

5. Phylogenetic trees from single cell colonies, the determination of growth per year and age at onset from these trees and comparison with estimates from longitudinal data
    1. Run `Scripts/plot_tree.R` (`Rscript Scripts/plot_tree.R`) - this runs all of the aforementioned analysis and plots it as displayed in the manuscript
    2. Run `Notebook_Mitchell.Rmd` (in `Rstudio`) - analyses the phylogenetic data from [Michell et al. (2021)]()

6. Analysis of the historical growth effect and poor fits
    1. Run `Notebook_HistoricalGrowth_PoorFits.Rmd` (in `Rstudio`) - this runs analyses which factors - technical and biological - may be determinant of the difference between historical and inferred growth and a fit being poor (having one or more outlier)

*Optional: run `Scripts/plots_for_initial_panel.R`*


