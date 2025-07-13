#!/bin/bash
path2_datasets="/home/slam-emix/Datasets/BACK_END/3D/"
dataset_list="TUM_FR1_DESK KITTI_05 KITTI_00"
outliers="10 20 30 40 50 60 70 80 90 100"
monte_runs="00 01 02 03 04 05 06 07 08 09"

# G2O related solutiona
# exp_date="301124"
# g2o_opt="G2O_IPC_REC"

exp_date="110125"
g2o_opt="G2O_IPC_REC_1"

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
            yq -i ".s_factor=50.0" ${tmp_yaml}
            yq -i ".use_best_k_buddies=false" ${tmp_yaml}
            yq -i ".use_recovery=true" ${tmp_yaml}
            ../build/ipc_tester_3D -c ${tmp_yaml} &
        done
        jobs
        wait
        rm ./*.yaml
    done
    echo "Finished "${dataset}
done
date
