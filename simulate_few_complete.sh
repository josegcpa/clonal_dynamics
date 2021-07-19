CLONEX_PATH=/homes/josegcpa/clonex/clonex
N_DRIVERS=50
R=20

DIR=hsc_output_bnpr_complete # output directory
mkdir -p $DIR
#rm -rf hsc_output_bnpr_complete/*

N=200000
g=800 # ~ 80 years
G=20 # every 2 years
NP=7500
MR="20e-6"

for fmr in 0.005_200 0.01_50 0.015_20 0.02_15 0.025_10 0.03_6
do
    fitness=$(echo $fmr | cut -d '_' -f 1)
    mutation_rate=$(echo $fmr | cut -d '_' -f 2)
    for i in $(seq 1 $R);
    do
        dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"_$i
        mkdir -p $dir_name
        job_name=SIM_"$fitness"_$i
        bsub \
            -J $job_name \
            -M 2000 -n 1 \
            -o /dev/null -e /dev/null \
            "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p $NP -d $N_DRIVERS -u "$mutation_rate"e-9 -v $MR -G $G -s $fitness -R 1 -r $RANDOM"
        bsub \
            -J TRAJ_$job_name \
            -M 2000 -n 1 \
            -o /dev/null -e /dev/null \
            -w "ended($job_name)" \
            "./clonex_trajectory $dir_name/r001.csv $N_DRIVERS > $dir_name/driver_trajectory"
        bsub \
            -J LAST_$job_name \
            -M 2000 -n 1 \
            -o /dev/null -e /dev/null \
            -w "ended($job_name)" \
            "cat $dir_name/r001.csv | awk -v g=$g -F'\t' '{ if (\$1 == g) print \$0 }' > $dir_name/last_generation"
    done
done
