mrl=$(seq 1 30 | tr ' ' '\n' | awk '{print $1 / 2}' | xargs)
mrh=$(seq 10 10 300)
fl="1 2 3 4 5"
fh=$(seq 6 25)
CLONEX_PATH=/homes/josegcpa/clonex/clonex
N_DRIVERS=50
R=3

DIR=hsc_output_200k_13 # output directory
P_DIR=hsc_output_200k_13_processed
mkdir -p $DIR
mkdir -p $P_DIR

## N = 200K; g = 13
N=200000
g=1300
G=5

for mutation_rate in $mrl
  do
  for fitness in $fl
  do
    fitness=0$(echo "scale=4; $fitness/1000" | bc)
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"
    p_dir_name=$P_DIR/$(basename $dir_name).csv
    if [[ ! -e $dir_name  ]]; then mkdir $dir_name; fi
    echo bsub -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p 1000 -d $N_DRIVERS -u "$mutation_rate"e-9 -v 2e-6 -G $G -s $fitness -R $R; Rscript Scripts/process_simulation.R $dir_name 13 200000 5 40 $p_dir_name"
  done
done

for mutation_rate in $mrh
  do
  mutation_rate=0$(echo "scale=3; $mutation_rate/100" | bc)
  for fitness in $fh
  do
    fitness=0$(echo "scale=4; $fitness/1000" | bc)
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"
    p_dir_name=$P_DIR/$(basename $dir_name).csv
    if [[ ! -e $dir_name  ]]; then mkdir $dir_name; fi
    echo bsub -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p 1000 -d $N_DRIVERS -u "$mutation_rate"e-9 -v 2e-6 -G $G -s $fitness -R $R; Rscript Scripts/process_simulation.R $dir_name 13 200000 5 40 $p_dir_name"
  done
done

## N 100K; g = 5
DIR=hsc_output_100k_5 # output directory
P_DIR=hsc_output_100k_5_processed
mkdir -p $DIR
mkdir -p $P_DIR
N=100000
g=500
G=5

for mutation_rate in $mrl
  do
  mutation_rate=$(echo "scale=2; $mutation_rate*13/5" | bc)
  for fitness in $fl
  do
    fitness=0$(echo "scale=4; $fitness/1000*13/5" | bc)
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"
    p_dir_name=$P_DIR/$(basename $dir_name).csv
    if [[ ! -e $dir_name  ]]; then mkdir $dir_name; fi
    echo bsub -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p 1000 -d $N_DRIVERS -u "$mutation_rate"e-9 -v 5.2e-6 -G $G -s $fitness -R $R; Rscript Scripts/process_simulation.R $dir_name 5 100000 5 15 $p_dir_name"
  done
done

for mutation_rate in $mrh
  do
  mutation_rate=$(echo "scale=2; $mutation_rate*13/5" | bc)
  mutation_rate=0$(echo "scale=3; $mutation_rate/100" | bc)
  for fitness in $fh
  do
    fitness=0$(echo "scale=4; $fitness/1000*13/5" | bc)
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"
    p_dir_name=$P_DIR/$(basename $dir_name).csv
    if [[ ! -e $dir_name  ]]; then mkdir $dir_name; fi
    echo bsub -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p 1000 -d $N_DRIVERS -u "$mutation_rate"e-9 -v 5.2e-6 -G $G -s $fitness -R $R; Rscript Scripts/process_simulation.R $dir_name 5 100000 5 15 $p_dir_name"
  done
done

## N = 50K; g = 1
DIR=hsc_output_50k_1 # output directory
P_DIR=hsc_output_50k_1_processed
mkdir -p $DIR
mkdir -p $P_DIR
N=50000
g=100
G=1

for mutation_rate in $mrl
  do
  mutation_rate=$(echo "scale=2; $mutation_rate*13" | bc)
  for fitness in $fl
  do
    fitness=0$(echo "scale=4; $fitness/1000*13" | bc)
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"
    p_dir_name=$P_DIR/$(basename $dir_name).csv
    if [[ ! -e $dir_name  ]]; then mkdir $dir_name; fi
    bsub -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p 1000 -d $N_DRIVERS -u "$mutation_rate"e-9 -v 26e-6 -G $G -s $fitness -R $R; Rscript Scripts/process_simulation.R $dir_name 1 50000 1 3 $p_dir_name"
  done
done

for mutation_rate in $mrh
  do
  mutation_rate=$(echo "scale=2; $mutation_rate*13" | bc)
  mutation_rate=0$(echo "scale=3; $mutation_rate/100" | bc)
  for fitness in $fh
  do
    fitness=0$(echo "scale=4; $fitness/1000*13" | bc)
    dir_name=$DIR/hsc_"$fitness"_"$mutation_rate"
    p_dir_name=$P_DIR/$(basename $dir_name).csv
    if [[ ! -e $dir_name  ]]; then mkdir $dir_name; fi
    bsub -M 2000 -n 1 -o /dev/null -e /dev/null "$CLONEX_PATH -n $N -N $N -g $g -f $dir_name -t 0 -p 1000 -d $N_DRIVERS -u "$mutation_rate"e-9 -v 26e-6 -G $G -s $fitness -R $R; Rscript Scripts/process_simulation.R $dir_name 1 50000 1 3 $p_dir_name"
  done
done
