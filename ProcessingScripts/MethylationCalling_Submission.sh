#!/bin/bash
#$ -pe shared 8
#$ -l h_rt=4:00:00
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
module load samtools

sample=$(cat $files_to_process | head -${SGE_TASK_ID} | tail -1 )

samtools fixmate -p -m ${file_directory}${sample}_masked.bam ${file_directory}${sample}_masked.fixmates.bam
samtools sort -@ 8 -o ${file_directory}${sample}_masked.sorted.bam ${file_directory}${sample}_masked.fixmates.bam
rm ${file_directory}${sample}_masked.fixmates.bam
samtools markdup ${file_directory}${sample}_masked.sorted.bam ${file_directory}${sample}_masked.dup.bam
rm ${file_directory}${sample}_masked.sorted.bam
samtools index ${file_directory}${sample}_masked.dup.bam
python3 -m BSBolt CallMethylation -t 8 -BQ 5 -MQ 0 -IO -min 1 -CG -O ${file_directory}${sample}_masked -I ${file_directory}${sample}_masked.dup.bam -DB ${index} -IO > ${file_directory}${sample}_masked.meth.log
