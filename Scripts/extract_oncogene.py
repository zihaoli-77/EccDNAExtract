import pandas as pd
import subprocess
import os
import glob

# 定义输入路径和输出路径
input_dirs = [
    "/public/home/hpc216085/eccDNA/useSample/combineBam_sort/bamfiles/temp/",

]

output_dirs = [
    "/public/home/hpc216085/eccDNA/useSample/combineBam_sort/immunomodulatoryGene_bamfiles/temp/",

]

# Read gene position file
gene_positions_file = "genes.xlsx"  # Contains information files
gene_data = pd.read_excel(gene_positions_file,)

# Iterate through each input directory
for input_dir, output_dir in zip(input_dirs, output_dirs):

    os.makedirs(output_dir, exist_ok=True)

    bam_files = glob.glob(os.path.join(input_dir, "*.bam"))

    for bam_file in bam_files:

        base_name = os.path.splitext(os.path.basename(bam_file))[0]

        base_name = base_name.replace("combine.sorted", "").strip(".")

        for index, row in gene_data.iterrows():
            gene_name = row['geneName']
            chromosome = row['chr']
            start = row['start']
            end = row['end']
            extracted_bam = f"{output_dir}/{gene_name}.{base_name}.bam"
            sorted_bam = f"{output_dir}/{gene_name}.{base_name}.sorted.bam"
            cmd_extract = ["samtools", "view", "-b", bam_file, f"{chromosome}:{start}-{end}"]
            with open(extracted_bam, 'wb') as out_file:
                subprocess.run(cmd_extract, stdout=out_file)

            print(f"Extracted segments for {gene_name} from {bam_file} into {extracted_bam}")

            cmd_sort = ["samtools", "sort", extracted_bam, "-o", sorted_bam]
            subprocess.run(cmd_sort)

            print(f"Sorted segments for {gene_name} into {sorted_bam}")

            cmd_index = ["samtools", "index", sorted_bam]
            subprocess.run(cmd_index)

            print(f"Indexed sorted file for {gene_name}")

            os.remove(extracted_bam)
            print(f"Deleted extracted BAM file: {extracted_bam}")

print("All operations completed.")
