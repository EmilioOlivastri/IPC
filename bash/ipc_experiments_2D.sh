#!/bin/bash
path2_datasets="/home/slam-emix/Datasets/BACK_END/2D/"
#dataset_list="CSAIL FR079 FRH MIT INTEL M3500"
dataset_list="CSAIL FR079 FRH MIT INTEL"
outliers="10 20 30 40 50 60 70 80 90 100"
monte_runs="00 01 02 03 04 05 06 07 08 09"

# G2O related solutiona
#g2o_opt="G2O_IPC_K2"
#g2o_opt="G2O_IPC_REC_20"
g2o_opt="G2O_IPC_QUAD_GATED_M_EST"
#exp_date="210823"
#exp_date="311024"
exp_date="301025"

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
            yq -i ".m_estimator_delta=0.584" ${tmp_yaml}
            yq -i ".k_buddies=2" ${tmp_yaml}
            yq -i ".use_best_k_buddies=false" ${tmp_yaml}
            yq -i ".use_recovery=false" ${tmp_yaml}
            yq -i ".fast_reject_iter_base=20" ${tmp_yaml}
            yq -i ".slow_reject_iter_base=50" ${tmp_yaml}
            yq -i ".fast_reject_th=11.345" ${tmp_yaml}
            yq -i ".slow_reject_th=12.838" ${tmp_yaml}
            #yq -i ".fast_reject_th=10.64" ${tmp_yaml}
            #yq -i ".slow_reject_th=10.64" ${tmp_yaml}
            ../build/ipc_tester_2D -c ${tmp_yaml} &
        done
        jobs
        wait
        rm ./*.yaml
    done
    echo "Finished "${dataset}
done
date
