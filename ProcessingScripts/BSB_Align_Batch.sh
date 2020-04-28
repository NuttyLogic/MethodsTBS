files_to_process=samples.txt

# get the number of lines in txt file
number_files=$(cat $files_to_process | wc -l)

wd=working_directory

# set variables to pass to alignemnt script
file_directory=${wd}/trimmed_fastqs/
output_folder=${wd}/alignments/
index=${wd}/hg38_masked


# of jobs to process simultaneously
JOBS=40


# submit jobs

qsub -M colinpatfarrell@g.ucla.edu -m a -t 1-$number_files -tc $JOBS BSB_Align_Submission.sh -f $files_to_process -d $file_directory -o $output_folder  -i $index
