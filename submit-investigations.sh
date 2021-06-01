bsub -n 4 -M 16000 Rscript Scripts/investigate_simulations_competition.R

bsub -n 4 -M 16000 Rscript Scripts/investigate_simulations_developmental.R

bsub -n 4 -M 16000 Rscript Scripts/investigate_simulations_early_late.R
