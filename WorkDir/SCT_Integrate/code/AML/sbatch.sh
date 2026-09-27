#!/bin/bash
#SBATCH --job-name=processing  #设置作业名称
#SBATCH --partition=q07    #指定作业运行分区
#SBATCH --get-user-env     #继承环境变量
#SBATCH --chdir=/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/AML
#SBATCH --output=modified_processing.R.out
#SBATCH --error=modified_processing.R.out
#SBATCH --mem=40G
#SBATCH --mincpus=10

Rscript --no-save  modified_processing.R

# 创建完成标志文件
# ecscript 的退出状态
if [ $? -eq 0 ]; then
    # 正常结束，创建成功标志文件
    echo "作业完成" > /home/user/wangguosheng/singlecell/job_completed.txt
else
    # 出现错误，创建错误标志文件
    echo "作业失败" > /home/user/wangguosheng/singlecell/job_failed.txt
fi
