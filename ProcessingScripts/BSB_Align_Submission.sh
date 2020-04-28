#!/bin/bash
#$ -pe shared 8
#$ -l h_rt=2:00:00
#$ -l h_data=2G


while getopts f:d:o:i: option
do
 case "${option}"
 in
 f) files_to_process=${OPTARG};;
 d) file_directory=${OPTARG};;
 o) output_folder=${OPTARG};;
 i) index=${OPTARG};;

 esac
done

mkdir -p $output_folder

cd ${output_folder}

. /u/local/Modules/default/init/modules.sh
module load python/3.7.2

sample=$(cat $files_to_process | head -${SGE_TASK_ID} | tail -1 )

python3 -m BSBolt Align -t 8 -DB ${index} -F1 ${file_directory}${sample}_1.fastq.gz -F2 ${file_directory}${sample}_2.fastq.gz -O ${output_folder}${sample} -T 110  > ${output_folder}${sample}.log
