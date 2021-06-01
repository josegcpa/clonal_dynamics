DIR=hsc_output_developmental # output directory
CLONEX_PATH=/homes/josegcpa/clonex/clonex
mkdir -p $DIR

mutation_rate=300
mutation_rate=0$(echo "scale=3; $mutation_rate/100" | bc)
for fitness in $(seq 1 10) 12 15 20
do
    fitness=0$(echo "scale=4; $fitness/1000" | bc)
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"
    mkdir -p $dir_name
    bsub -M 1000 -n 1 -o /dev/null -e /dev/null\
        $CLONEX_PATH -a 64 -n 1 -N 200000 -g 1300 -f $dir_name -t 0 -p 1000 -d 100 -u "$mutation_rate"e-9 -v 2e-6 -G 20 -s $fitness -R 100
done
