#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_path> <output_path>"
    exit 1
fi

# 获取输入和输出路径
input_path=$1
output_path=$2

# 创建输出路径（如果不存在）
mkdir -p "$output_path"

# 清空输出目录（如果需要）
# rm -f "${output_dir}"/*.bam "${output_dir}"/*.txt

# 遍历所有以 .sorted.bam 结尾的文件
for bam_file in "${input_dir}"/*.sorted.bam; do
    # 获取文件名（不带路径和扩展名）
    base_name=$(basename "$bam_file" .sorted.bam)
    
    # 1. 过滤 BAM 文件并生成新的 BAM 文件
    samtools view -@ 16 -hf 0x100 "$bam_file" -bS > "${output_dir}/${base_name}_sec.bam"
    
    # 2. 使用 bedtools 转换 BAM 为 BED 格式
    bedtools bamtobed -cigar -i "${output_dir}/${base_name}_sec.bam" | sed -e s/\\//\ /g | awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' > "${output_dir}/${base_name}_sec.txt"
    
    # 3. 生成 AB010_sec_uid.txt 文件（假设需要手动创建或已有）
    # 这里是一个示例，您需要根据实际情况提供 UID 列表
    # echo -e "UID1\nUID2\nUID3" > "${output_dir}/AB010_sec_uid.txt"
    awk '{print $4}' "${output_dir}/${base_name}_sec.txt" | sort | uniq > "${output_dir}/${base_name}_sec_uid.txt"
    # 4. 使用 UID 文件生成新的 BAM 文件
    samtools view -h -N "${output_dir}/${base_name}_sec_uid.txt" "$bam_file" > "${output_dir}/${base_name}_spc_paired.bam"
done

# 提示完成
echo "Batch processing completed."
