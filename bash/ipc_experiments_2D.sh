#!/bin/bash
path2_datasets="/home/slamemix/Data/RobustSLAM/2D/"
dataset_list="CSAIL FR079 FRH INTEL"
outliers="10 20 30 40 50 60 70 80 90 100"
monte_runs="00 01 02 03 04 05 06 07 08 09"

# G2O related solutiona
g2o_opt="G2O_IPC"
exp_date="260202"

date 
for dataset in ${dataset_list}
do
    cfg_file=${path2_datasets}${dataset}"/params.yaml"
    for out in ${outliers}
    do
        for run in ${monte_runs}
        do
            input_file=${path2_datasets}${dataset}"/SPOILED_DATA/"${out}"/"${run}".g2o"
            output_traj=${path2_datasets}${dataset}"/EXP/"${exp_date}"/"${g2o_opt}"/"${out}"/"${run}".TRJ"
            tmp_yaml="./"${g2o_opt}"_"${out}"_"${run}".yaml"
            cp ${cfg_file} ${tmp_yaml}
            yq -i ".dataset=\"$input_file\"" ${tmp_yaml} 
            yq -i ".output=\"$output_traj\"" ${tmp_yaml}
            yq -i ".s_factor=10.0" ${tmp_yaml}
            yq -i ".use_recovery=true" ${tmp_yaml}
            yq -i ".fast_reject_th=0.58" ${tmp_yaml}
            yq -i ".fast_reject_iter_base=25" ${tmp_yaml}
            yq -i ".slow_reject_th=6.251" ${tmp_yaml}
            yq -i ".slow_reject_iter_base=50" ${tmp_yaml}
            ../build/thirdParty/IPC/ipc_tester_2D -c ${tmp_yaml} &
        done
        jobs
        wait
        rm ./*.yaml
    done
    echo "Finished "${dataset}
done
date
