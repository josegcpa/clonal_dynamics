bsub -n 4 -M 8000 -o /dev/null -e /dev/null Rscript Scripts/run_simulation_50k_1.R
bsub -n 4 -M 8000 -o /dev/null -e /dev/null Rscript Scripts/run_simulation_100k_5.R
bsub -n 4 -M 8000 -o /dev/null -e /dev/null Rscript Scripts/run_simulation_200k_13.R
