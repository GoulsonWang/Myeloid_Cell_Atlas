#!/bin/bash
#SBATCH --job-name=QC  #设置作业名称
#SBATCH --partition=q07    #指定作业运行分区
#SBATCH --get-user-env     #继承环境变量
#SBATCH --chdir=/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code
#SBATCH --output=QC2.r.out
#SBATCH --error=QC2.r.out
#SBATCH --mem=40G
#SBATCH --mincpus=1

Rscript --no-save  QC2.R

# 创建完成标志文件
# ecscript 的退出状态
if [ $? -eq 0 ]; then
    # 正常结束，创建成功标志文件
    echo "作业完成" > /home/user/wangguosheng/singlecell/job_completed.txt
else
    # 出现错误，创建错误标志文件
    echo "作业失败" > /home/user/wangguosheng/singlecell/job_failed.txt
fi
