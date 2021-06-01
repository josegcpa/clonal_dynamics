CLONEX_PATH=/homes/josegcpa/clonex/clonex
N_DRIVERS=25
R=20

DIR=hsc_output_bnpr # output directory
mkdir -p $DIR

N=200000
g=800 # ~ 80 years
G=400 # every ~ 40 years
NP=7500
MR="20e-6"

fitness=0.005
mutation_rate=400
for i in $(seq 1 $R);
do
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"_$i
    mkdir -p $dir_name
    bsub -J SIM_$fitness -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p $NP -d $N_DRIVERS -u "$mutation_rate"e-9 -v $MR -G $G -s $fitness -R 1 -r $RANDOM"
done

fitness=0.01
mutation_rate=50
for i in $(seq 1 $R);
do
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"_$i
    mkdir -p $dir_name
    bsub -J SIM_$fitness -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p $NP -d $N_DRIVERS -u "$mutation_rate"e-9 -v $MR -G $G -s $fitness -R 1 -r $RANDOM"
done

fitness=0.015
mutation_rate=20
for i in $(seq 1 $R);
do
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"_$i
    mkdir -p $dir_name
    bsub -J SIM_$fitness -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p $NP -d $N_DRIVERS -u "$mutation_rate"e-9 -v $MR -G $G -s $fitness -R 1 -r $RANDOM"
done

fitness=0.02
mutation_rate=15
for i in $(seq 1 $R);
do
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"_$i
    mkdir -p $dir_name
    bsub -J SIM_$fitness -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p $NP -d $N_DRIVERS -u "$mutation_rate"e-9 -v $MR -G $G -s $fitness -R 1 -r $RANDOM"
done

fitness=0.025
mutation_rate=10
for i in $(seq 1 $R);
do
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"_$i
    mkdir -p $dir_name
    bsub -J SIM_$fitness -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p $NP -d $N_DRIVERS -u "$mutation_rate"e-9 -v $MR -G $G -s $fitness -R 1 -r $RANDOM"
done

fitness=0.03
mutation_rate=5
for i in $(seq 1 $R);
do
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"_$i
    mkdir -p $dir_name
    bsub -J SIM_$fitness -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p $NP -d $N_DRIVERS -u "$mutation_rate"e-9 -v $MR -G $G -s $fitness -R 1 -r $RANDOM"
done
