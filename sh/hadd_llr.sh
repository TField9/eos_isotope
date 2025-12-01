#!/bin/bash

# 设置基础路径
BASE_DIR="/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_llr_reweight"

# 定义同位素配置
# declare -A isotopes=(
#     ["he3"]="He." 
#     ["he4"]="He."
#     ["li6"]="Li." 
#     ["li7"]="Li."
#     ["be7"]="Be." 
#     ["be9"]="Be." 
#     ["be10"]="Be."
#     ["b10"]="B." 
#     ["b11"]="B."
# )

declare -A isotopes=(
    ["li6"]="Li." 
    ["li7"]="Li."
    ["be10"]="Be."
)

# declare -A isotopes=(
#     ["be7"]="Be." 
#     ["be9"]="Be." 
#     ["be10"]="Be."
# )

# 生成索引数组
temp_idx=()
target_idx=()

# 生成完整的数字序列 (1, 11, 21, ..., 1991)
all_numbers=()
for i in {1..390}; do
    all_numbers+=($(( (i - 1) * 10 + 1 )))
done

# 随机打乱整个序列
# 设置固定种子（可重现）
# 设置固定种子（可重现）
FIXED_SEED=42
shuffled_numbers=($(printf '%s\n' "${all_numbers[@]}" | awk -v seed="$FIXED_SEED" 'BEGIN{srand(seed);} {print rand() "\t" $0}' | sort -k1,1n | cut -f2-))

# 方法2：如果方法1还有问题，使用系统随机源
# shuffled_numbers=($(shuf -e "${all_numbers[@]}"))

echo "打乱后的数组长度: ${#shuffled_numbers[@]}"
# shuffled_numbers=($(shuf -e "${all_numbers[@]}"))

# 从打乱的序列中随机选择两组不重叠的数字
temp_idx=("${shuffled_numbers[@]:0:195}")      # 前195个数字
target_idx=("${shuffled_numbers[@]:195:195}")  # 后195个数字

# 清理旧文件函数 - 重新添加
cleanup_old_files() {
    local isotope=$1
    local element_dir=$2
    
    # 删除临时文件和目标文件
    local target_file="${BASE_DIR}/${element_dir}/${isotope}_llr_target.root"
    local temp_file="${BASE_DIR}/${element_dir}/${isotope}_llr_temp.root"
    
    echo "清理 ${isotope} 的旧文件..."
    
    if [ -f "$target_file" ]; then
        echo "删除目标文件: $target_file"
        rm -f "$target_file"
    fi
    
    if [ -f "$temp_file" ]; then
        echo "删除临时文件: $temp_file"
        rm -f "$temp_file"
    fi
}

# 合并文件函数（带调试信息）
merge_files() {
    local isotope=$1
    local element_dir=$2
    local output_suffix="${@: -1}"
    
    local num_args=$#
    local idx_array=("${@:3:$((num_args-3))}")
    
    echo "=== 处理 ${isotope} (${output_suffix}) ==="
    echo "预期文件数: ${#idx_array[@]}"
    
    # 使用find命令查找文件
    local found_files=()
    local missing_indices=()
    
    for idx in "${idx_array[@]}"; do
        # 尝试多种文件名模式
        local file=$(find "${BASE_DIR}/${element_dir}" -maxdepth 1 -name "${isotope}_${idx}_*llr.root" -o -name "${isotope}*${idx}*llr.root" 2>/dev/null | head -1)
        
        if [ -n "$file" ] && [ -f "$file" ]; then
            found_files+=("$file")
        else
            missing_indices+=("$idx")
        fi
    done
    
    # 输出调试信息
    echo "找到文件数: ${#found_files[@]}"
    if [ ${#missing_indices[@]} -gt 0 ]; then
        echo "缺失索引: ${missing_indices[*]:0:10}"  # 只显示前10个缺失索引
    fi
    
    if [ ${#found_files[@]} -eq 0 ]; then
        echo "错误: 没有找到任何文件"
        return 1
    fi
    
    # 执行合并
    local output_file="${BASE_DIR}/${element_dir}/${isotope}_llr_${output_suffix}.root"
    echo "开始合并 ${#found_files[@]} 个文件到 ${output_file}"
    
    hadd -f "$output_file" "${found_files[@]}"
    local result=$?
    
    if [ $result -eq 0 ]; then
        echo "合并成功!"
    else
        echo "合并失败!"
    fi
    
    return $result
}

# 主执行逻辑
main() {
    echo "开始清理旧文件..."
    for isotope in "${!isotopes[@]}"; do
        element_dir="${isotopes[$isotope]}"
        cleanup_old_files "$isotope" "$element_dir"
    done
    
    echo "开始合并文件..."
    # 合并临时文件
    for isotope in "${!isotopes[@]}"; do
        element_dir="${isotopes[$isotope]}"
        merge_files "$isotope" "$element_dir" "${temp_idx[@]}" "temp"
        echo "---"
    done
    
    # 合并目标文件
    for isotope in "${!isotopes[@]}"; do
        element_dir="${isotopes[$isotope]}"
        merge_files "$isotope" "$element_dir" "${target_idx[@]}" "target"
        echo "---"
    done
    
    echo "所有文件处理完成!"
}

# 错误处理
set -e  # 遇到错误时退出

# 执行主函数
main "$@"



