#!/bin/bash
#SBATCH --job-name=build_shiny       # 作业名称
#SBATCH --partition=q03              # 分区（q08不存在，改用q02，190GB内存）
#SBATCH --get-user-env               # 继承环境变量
#SBATCH --mem=160G                   # 内存（和你的其他作业一样）
#SBATCH --cpus-per-task=10            # CPU 核心数
#SBATCH --output=build_shiny.out     # 标准输出日志
#SBATCH --error=build_shiny.err      # 错误输出日志
#SBATCH --chdir=/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Shiny

# 加载 conda 环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate shiny

# 运行构建脚本
Rscript --no-save build_shiny_app.R

# 检查退出状态并创建标志文件
if [ $? -eq 0 ]; then
    echo "============================================"
    echo "ShinyCell2 App 构建成功！"
    echo "产出目录: Shiny/shinyApp/"
    echo ""
    echo "在登录节点启动 App:"
    echo "  cd /home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Shiny/shinyApp"
    echo "  R -e \"shiny::runApp(host='0.0.0.0', port=3838)\""
    echo "============================================"
    # 创建成功标志文件
    echo "ShinyCell2 App 构建完成" > /home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Shiny/build_shiny_success.txt
else
    echo "============================================"
    echo "构建失败！请查看 build_shiny.err 日志。"
    echo "============================================"
    # 创建失败标志文件
    echo "构建失败，请查看 build_shiny.err" > /home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Shiny/build_shiny_failed.txt
    exit 1
fi
