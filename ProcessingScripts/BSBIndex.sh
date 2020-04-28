#!/bin/bash
#$ -pe shared 8
#$ -l h_rt=20:00:00
#$ -l h_data=4G

. /u/local/Modules/default/init/modules.sh
module load python/3.7.2

wd=working_directory

python3 -m BSBolt Index -G hg38_lambda.fa.gz -DB ${wd}hg38_index