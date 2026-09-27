#!/bin/bash
#SBATCH --job-name=2Merge_Integration  #设置作业名称
#SBATCH --partition=q06    #指定作业运行分区
#SBATCH --get-user-env     #继承环境变量
#SBATCH --chdir=/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/code/TotalData
#SBATCH --output=Merge_Integration2.R.out
#SBATCH --error=Merge_Integration2.R.error
#SBATCH --mem=170G
#SBATCH --cpus-per-task=10      #每个任务分配十个cpu

Rscript --no-save  Merge_Data_2.R &&/
Rscript --no-save  Integration_Total_2.R

# 创建完成标志文件
# ecscript 的退出状态
if [ $? -eq 0 ]; then
    # 正常结束，创建成功标志文件
    echo "作业完成" > /home/user/wangguosheng/singlecell/job_completed.txt
else
    # 出现错误，创建错误标志文件
    echo "作业失败" > /home/user/wangguosheng/singlecell/job_failed.txt
fi
