## 200K 13

output_dir=hsc_output_200k_13_processed
input_dir=hsc_output_200k_13
mkdir -p $output_dir

for folder in $input_dir/*/
do
    output_file=$output_dir/$(basename $folder).csv
    bsub -n 2 -M 2000 -o /dev/null -e /dev/null Rscript Scripts/process_simulation.R $folder 13 200000 5 40 $output_file
done


## 100K 5

output_dir=hsc_output_100k_5_processed
input_dir=hsc_output_100k_5
mkdir -p $output_dir

for folder in $input_dir/*/
do
    output_file=$output_dir/$(basename $folder).csv
    bsub -n 2 -M 2000 -o /dev/null -e /dev/null Rscript Scripts/process_simulation.R $folder 5 100000 5 15 $output_file
done

## 50K 1

output_dir=hsc_output_50k_1_processed
input_dir=hsc_output_50k_1
mkdir -p $output_dir

for folder in $input_dir/*/
do
    output_file=$output_dir/$(basename $folder).csv
    bsub -n 2 -M 2000 -o /dev/null -e /dev/null Rscript Scripts/process_simulation.R $folder 1 50000 1 3 $output_file
done
