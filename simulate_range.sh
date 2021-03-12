mut_rate_range="10 5 350"
fitness_range="1 1 25"

DIR=hsc_output # output directory
CLONEX_PATH=/homes/josegcpa/clonex/clonex
mkdir -p $DIR

for mutation_rate in $(seq $mut_rate_range)
  do
  mutation_rate=0$(echo "scale=3; $mutation_rate/100" | bc)
  for fitness in $(seq $fitness_range)
    do
    fitness=0$(echo "scale=4; $fitness/1000" | bc)
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"
    if [[ ! -e $dir_name  ]]; then mkdir $dir_name; fi
    bsub -M 1000 -n 1 -o /dev/null -e /dev/null\
        $CLONEX_PATH -n 200000 -N 200000 -g 1300 -f $dir_name -t 0 -p 200 -d 20 -u "$mutation_rate"e-9 -v 2e-6 -G 20 -s $fitness -R 3
  done
done
