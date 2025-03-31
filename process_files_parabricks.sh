#!/bin/bash

# Activate PyCharm virtual python environment to manage dependencies. This will look different in a different runtime environment
source ./venv/bin/activate

# Clear files and logs from prior runs
rm -rf ./F_*.fq ../F_*.bam ../F_*.bai gatk_*.log prep_files_*.log

# Clear intermediate processing files
rm -rf ../processing/1-fastqc/F* ../processing/2-trim/F* ../processing/3-fastqc/F* ../processing/4-bwa/F*

## declare an array variable
declare -a arr=("ND6W8N_1_Alu_SB" "ND6W8N_2_Alu_Cout" "ND6W8N_3_Fat_SB" "ND6W8N_4_FAT_Cout")

# Start processing each sample
for i in "${arr[@]}"
do
  echo "$i"
  prep_job_id=$(sbatch --parsable --job-name="prep_files" --partition=LocalQ --time=24:00:00 --output="prep_files_%j.log")
  # Submit alignment job as a dependent job
  sbatch --dependency=afterok:$prep_job_id GFP_SEQ_TrmAlgn_human.hg38.parabricks.bash "$i"
done