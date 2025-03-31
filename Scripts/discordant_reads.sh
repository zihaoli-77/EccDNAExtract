#!/bin/bash


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_path> <output_path>"
    exit 1
fi

# input path and output path
input_path=$1
output_path=$2

# Create the output path if it does not exist
mkdir -p "$output_path"

# Iterate over all files ending with .sorted.bam
for bam_file in "$input_path"/*sorted.bam; do
    
    base_name=$(basename "$bam_file" .sorted.bam)

    echo "Processing $bam_file..."

    
    samtools view -h "$bam_file" | grep "SA:Z:" | awk '{print $1}' | sort -u > "$output_path/remove_${base_name}_split_qname.txt"

    
    samtools view -h "$bam_file" | grep -v -Ff "$output_path/remove_${base_name}_split_qname.txt" | samtools view -h -b -o "$output_path/no_split_${base_name}.bam"

    
    samtools view -h -F 0xE "$output_path/no_split_${base_name}.bam" -b -o "$output_path/no_flag248_${base_name}.bam"

    
    samtools sort -n "$output_path/no_flag248_${base_name}.bam" -o "$output_path/sorted_no_flag248_${base_name}.bam"

    
    samtools view "$output_path/sorted_no_flag248_${base_name}.bam" | awk '{print $1}' | sort | uniq -c | awk '$1 == 2 {print $2}' > "$output_path/twice_no_flag248_${base_name}_qnames.txt"

    
    samtools view -h -N "$output_path/twice_no_flag248_${base_name}_qnames.txt" "$output_path/sorted_no_flag248_${base_name}.bam" -b -o "$output_path/twice_sorted_no_flag248_${base_name}.bam"

    
    samtools view -h "$output_path/twice_sorted_no_flag248_${base_name}.bam" | awk '
    BEGIN { OFS="\t" }
    {
    # 
    if ($0 ~ /^@/) {
        print;
        next;
    }
    
    # 
    if ((and($2 , 16)) && !(and($2 , 32))) {  
        r1 = $0;                     
        r1_name = $1;                
        r1_chr = $3;                 
        r1_pos = $4;                 
        
         
        if (getline > 0) {
           
            if ($1 == r1_name) {    
                if (!(and($2 , 16)) && (r1_chr == $3) && (r1_pos < $4)) {
                    print r1;      
                    print $0;     
                }
                else if (r1_chr != $3) {
                    print r1;      
                    print $0;      
                }
            }
        }
      }
    }' | samtools view -bS - | samtools sort -o "$output_path/discordant_${base_name}.bam"

    
    samtools sort -n "$output_path/discordant_${base_name}.bam" -o "$output_path/${base_name}_sorted_discordant.bam"

    
    rm "$output_path/remove_${base_name}_split_qname.txt" \
       "$output_path/no_split_${base_name}.bam" \
       "$output_path/no_flag248_${base_name}.bam" \
       "$output_path/sorted_no_flag248_${base_name}.bam" \
       "$output_path/twice_no_flag248_${base_name}_qnames.txt" \
       "$output_path/twice_sorted_no_flag248_${base_name}.bam" \
       "$output_path/discordant_${base_name}.bam"

    echo "$bam_file processed and intermediate files removed."
done

