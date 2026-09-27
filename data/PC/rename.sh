#!/bin/bash

# 定义目标文件名
target_filename="barcodes.tsv.gz"

# 遍历当前目录及其子目录中的所有文件
find . -type f -name "*$target_filename" | while read -r file; do
    # 获取文件所在的目录
    dir=$(dirname "$file")
    
    # 移动并重命名文件
    mv "$file" "$dir/$target_filename"
    
    echo "Renamed: $file to $dir/$target_filename"
done

target_filename="features.tsv.gz"

# 遍历当前目录及其子目录中的所有文件
find . -type f -name "*$target_filename" | while read -r file; do
    # 获取文件所在的目录
    dir=$(dirname "$file")
    
    # 移动并重命名文件
    mv "$file" "$dir/$target_filename"
    
    echo "Renamed: $file to $dir/$target_filename"
done

target_filename="matrix.mtx.gz"

# 遍历当前目录及其子目录中的所有文件
find . -type f -name "*$target_filename" | while read -r file; do
    # 获取文件所在的目录
    dir=$(dirname "$file")
    
    # 移动并重命名文件
    mv "$file" "$dir/$target_filename"
    
    echo "Renamed: $file to $dir/$target_filename"
done

echo "Rename operation complete."