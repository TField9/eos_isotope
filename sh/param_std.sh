#!/bin/bash

# 简化版Be数据处理任务生成脚本
# 参数直接内置在脚本中

# 设置输出文件
output_file="/eos/ams/user/s/selu/mdst/tianye/sh/param_std.txt"
# 内置参数设置
beam_types="7 9 10"
train_start_num=1
train_end_num=5000
# test_start_num=6
# test_end_num=10
increment=10

# 创建输出目录
inputdir="/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_train_og"
outdir="/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_stdcut"
# 清空或创建输出文件
> "$output_file"

# 生成任务命令并写入文件
for beam in $beam_types; do
    current_start=$train_start_num
    while [[ $current_start -le $train_end_num ]]; do
        current_end=$(( current_start + increment - 1 ))
        
        # 确保不超过结束序号
        if [[ $current_end -gt $train_end_num ]]; then
            current_end=$train_end_num
        fi
        
        # 构造完整命令
        # cmd="be$beam /cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc/Be.B1403/be$beam.pl1.l1.4400.6_05.txt $current_start $current_end /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_train/sdst_be${beam}_${current_start}_${current_end}.root"
        cmd="$inputdir/sdst_be${beam}_${current_start}_${current_end}.root $outdir/sdst_be${beam}_${current_start}_${current_end}_std.root"
        
        # 写入文件并显示
        echo "$cmd" >> "$output_file"
        echo "已生成: $cmd"
        
        # 更新起始序号
        current_start=$(( current_end + 1 ))
    done
done

# for beam in $beam_types; do
#     current_start=$test_start_num
#     while [[ $current_start -le $test_end_num ]]; do
#         current_end=$(( current_start + increment - 1 ))
        
#         # 确保不超过结束序号
#         if [[ $current_end -gt $test_end_num ]]; then
#             current_end=$test_end_num
#         fi
        
#         # 构造完整命令
#         # cmd="be$beam /cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc/Be.B1403/be$beam.pl1.l1.4400.6_05.txt $current_start $current_end /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_test/sdst_be${beam}_${current_start}_${current_end}.root"

#         cmd="be$beam /cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc/Be.B1403/be$beam.pl1.l1.4400.6_05.txt $current_start $current_end sdst_be${beam}_${current_start}_${current_end}.root"

        
#         # 写入文件并显示
#         echo "$cmd" >> "$output_file"
#         echo "已生成: $cmd"
        
#         # 更新起始序号
#         current_start=$(( current_end + 1 ))
#     done
# done

echo "任务命令已保存到: $output_file"