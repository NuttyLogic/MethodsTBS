#!/bin/bash
#$ -pe shared 2
#$ -l h_rt=10:00:00
#$ -l h_data=1G


while getopts f:d:o:s:r:i:a: option
do
 case "${option}"
 in
 f) files_to_process=${OPTARG};;
 d) file_directory=${OPTARG};;
 o) output_folder=${OPTARG};;
 s) scratch=$OPTARG;;
 esac
done

cd ${output_folder}

. /u/local/Modules/default/init/modules.sh
module load python/3.7.2

sample_name=$(cat $files_to_process | head -${SGE_TASK_ID} | tail -1 )

echo ~/.local/bin/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --pair-filter=any -m 40 -o ${output_folder}/${sample_name}_1.fastq -p ${output_folder}/${sample_name}_2.fastq ${file_directory}/${sample_name}_L001_R1_001.fastq.gz ${file_directory}/${sample_name}_L001_R2_001.fastq.gz

~/.local/bin/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --pair-filter=any -m 40 -o ${output_folder}/${sample_name}_1.fastq -p ${output_folder}/${sample_name}_2.fastq ${file_directory}/${sample_name}_L001_R1_001.fastq.gz ${file_directory}/${sample_name}_L001_R2_001.fastq.gz > ${output_folder}/${sample_name}.cutadapt.log
