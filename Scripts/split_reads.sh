#!/bin/bash

# # input path and output path
#input_dir="/public/home/eccDNA/script_test/sample"
#output_dir="/public/home/spc_eccDNA/test"
input_dir="your path"
output_dir="your path"

# Clear the output directory (if necessary)
# rm -f "${output_dir}"/*.bam "${output_dir}"/*.txt

# Iterate over all files ending with .sorted.bam
for bam_file in "${input_dir}"/*.sorted.bam; do
    
    base_name=$(basename "$bam_file" .sorted.bam)
    
  
    samtools view -@ 16 -hf 0x100 "$bam_file" -bS > "${output_dir}/${base_name}_sec.bam"
    
    
    bedtools bamtobed -cigar -i "${output_dir}/${base_name}_sec.bam" | sed -e s/\\//\ /g | awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' > "${output_dir}/${base_name}_sec.txt"
    
 
   
    # echo -e "UID1\nUID2\nUID3" > "${output_dir}/AB010_sec_uid.txt"
    awk '{print $4}' "${output_dir}/${base_name}_sec.txt" | sort | uniq > "${output_dir}/${base_name}_sec_uid.txt"
    
    samtools view -h -N "${output_dir}/${base_name}_sec_uid.txt" "$bam_file" > "${output_dir}/${base_name}_spc_paired.bam"
done

# 提示完成
echo "Batch processing completed."
