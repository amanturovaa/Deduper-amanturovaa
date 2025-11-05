#!/bin/bash
#dedup

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=4                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=10:00:00 
#SBATCH --job-name=dedup_final_3.12py_                 #optional: job name
#SBATCH --output=dedup_final_3.12py_%j.out              #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=dedup_final_3.12py_%j.err               #optional: file to store stderr from job, %j adds the assigned jobID

#created environment with python 3.12 : conda activate deduper_env

out_dir="/projects/bgmp/amant/bioinfo/Deduper-amanturovaa"

input="/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam"
sorted_bam="${out_dir}/C1_SE_uniqAlign_sorted.bam"
sorted_sam="${out_dir}/C1_SE_uniqAlign_sorted.sam"
out_dedup="${out_dir}/dedup_final_3.12py_.sam"
umi="${out_dir}/STL96.txt"

/usr/bin/time -v samtools sort -@ 4 -o "${sorted_bam}" "${input}" 
samtools view -h "${sorted_bam}" > "${sorted_sam}"

/usr/bin/time -v python manturova_deduper.py \
    -f "${sorted_sam}" \
    -o "${out_dedup}" \
    -u "${umi}"