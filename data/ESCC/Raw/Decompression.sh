#!/bin/bash

# --- 配置 ---
# TAR 文件名的公共前缀和后缀
PREFIX="OMIX005710-"
SUFFIX=".tar"

# TAR 文件编号范围
START_NUM=3
END_NUM=49

# 需要跳过的编号
SKIP_NUM=33

# 目标目录前缀
DATA_PREFIX="data"

# --- 脚本主体 ---
# 初始化目录计数器
counter=1

# 循环处理指定范围内的每个编号
for ((num=START_NUM; num<=END_NUM; num++))
do
  # 如果当前编号是需要跳过的编号，则跳过
  if [[ $num -eq $SKIP_NUM ]]; then
    echo "Skipping number $num as configured."
    continue
  fi

  # 格式化编号为两位数 (例如 3 -> 03)
  formatted_num=$(printf "%02d" $num)
  
  # 构造 tar 文件名
  tar_file="${PREFIX}${formatted_num}${SUFFIX}"
  
  # 构造目标目录名
  data_dir="${DATA_PREFIX}${counter}"
  
  # 检查 tar 文件是否存在
  if [[ -f "$tar_file" ]]; then
    echo "Processing $tar_file..."
    
    # 创建目标目录
    mkdir -p "./RawData(decompression)/$data_dir"
    
    # 解压 tar 文件到目标目录
    # 使用 --strip-components=1 可以去除 tar 内部的一层目录结构（如果需要）
    # tar -xf "$tar_file" -C "$data_dir" --strip-components=1
    tar -xf "$tar_file" -C "./RawData(decompression)/$data_dir"
    
    if [[ $? -eq 0 ]]; then
      echo "  Successfully extracted $tar_file to ./RawData(decompression)/$data_dir"
    else
      echo "  Error extracting $tar_file" >&2
    fi
  else
    echo "Warning: File $tar_file not found. Skipping." >&2
  fi
  
  # 增加目录计数器 (只有在处理了文件后才增加)
  if [[ -f "$tar_file" ]]; then
      ((counter++))
  fi
done

echo "Processing complete. Created $((counter-1)) data directories."