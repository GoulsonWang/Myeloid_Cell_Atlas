#!/bin/bash
#SBATCH --job-name=Prediction  #设置作业名称
#SBATCH --partition=q04    #指定作业运行分区
#SBATCH --get-user-env     #继承环境变量
#SBATCH --chdir=/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Prediction/code
#SBATCH --output=0.Exact_ESCC_TestData.R.out
#SBATCH --error=0.Exact_ESCC_TestData.R.error
#SBATCH --mem=160G
#SBATCH --cpus-per-task=10      #每个任务分配十个cpu

Rscript --no-save  0.Exact_ESCC_TestData.R

# 创建完成标志文件
# ecscript 的退出状态
if [ $? -eq 0 ]; then
    # 正常结束，创建成功标志文件
    echo "作业完成" > /home/user/wangguosheng/singlecell/job_completed.txt
else
    # 出现错误，创建错误标志文件
    echo "作业失败" > /home/user/wangguosheng/singlecell/job_failed.txt
fi
