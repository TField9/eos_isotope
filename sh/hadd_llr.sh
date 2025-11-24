#!/bin/bash

# 设置基础路径
BASE_DIR="/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_llr"

# 定义同位素配置
declare -A isotopes=(
    ["he3"]="He." 
    ["he4"]="He."
    ["li6"]="Li." 
    ["li7"]="Li."
    ["be7"]="Be." 
    ["be9"]="Be." 
    ["be10"]="Be."
    ["b10"]="B." 
    ["b11"]="B."
)

# 生成索引数组
temp_idx=()
target_idx=()

# 生成完整的数字序列 (1, 11, 21, ..., 1991)
all_numbers=()
for i in {1..200}; do
    all_numbers+=($(( (i - 1) * 10 + 1 )))
done

# 随机打乱整个序列
shuffled_numbers=($(shuf -e "${all_numbers[@]}"))

# 从打乱的序列中随机选择两组不重叠的数字
temp_idx=("${shuffled_numbers[@]:0:100}")      # 前100个数字
target_idx=("${shuffled_numbers[@]:100:100}")  # 后100个数字

# 清理旧文件函数
cleanup_old_files() {
    local isotope=$1
    local element_dir=$2
    
    # 删除临时文件和目标文件
    rm -f "${BASE_DIR}/${element_dir}/${isotope}_llr_target.root"
    rm -f "${BASE_DIR}/${element_dir}/${isotope}_llr_temp.root"
}

# 合并文件函数
merge_files() {
    local isotope=$1
    local element_dir=$2
    local output_suffix="${@: -1}"  # 获取最后一个参数
    
    # 获取除最后一个参数外的所有参数作为索引数组
    local num_args=$#
    local idx_array=("${@:3:$((num_args-3))}")
    
    # 构建文件模式数组
    local file_patterns=()
    for idx in "${idx_array[@]}"; do
        file_patterns+=("${BASE_DIR}/${element_dir}/${isotope}_${idx}_*_llr.root")
    done
    
    # 检查是否有匹配的文件
    local found_files=()
    for pattern in "${file_patterns[@]}"; do
        if compgen -G "$pattern" > /dev/null; then
            found_files+=($pattern)
        fi
    done
    
    if [ ${#found_files[@]} -eq 0 ]; then
        echo "警告: 没有找到 ${isotope} 的匹配文件，跳过合并"
        return 1
    fi
    
    # 执行合并
    local output_file="${BASE_DIR}/${element_dir}/${isotope}_llr_${output_suffix}.root"
    echo "合并 ${isotope} 文件到 ${output_file} (${#found_files[@]} 个文件)"
    
    hadd -f "$output_file" "${found_files[@]}"
    return $?
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
        echo "处理 ${isotope} (temp)..."
        merge_files "$isotope" "$element_dir" "${temp_idx[@]}" "temp"
    done
    
    # 合并目标文件
    for isotope in "${!isotopes[@]}"; do
        element_dir="${isotopes[$isotope]}"
        echo "处理 ${isotope} (target)..."
        merge_files "$isotope" "$element_dir" "${target_idx[@]}" "target"
    done
    
    echo "所有文件处理完成!"
}

# 错误处理
set -e  # 遇到错误时退出

# 执行主函数
main "$@"