# aggregate matrix 

wd=/alignment_directory/

python3 -m BSBolt AggregateMatrix -F ${wd}alignments/S10.CGmap.gz,${wd}alignments/S11.CGmap.gz,${wd}alignments/S12.CGmap.gz,${wd}alignments/S13.CGmap.gz,${wd}alignments/S14.CGmap.gz,${wd}alignments/S15.CGmap.gz,${wd}alignments/S16.CGmap.gz,${wd}alignments/S17.CGmap.gz,${wd}alignments/S18.CGmap.gz,${wd}alignments/S19.CGmap.gz,${wd}alignments/S1.CGmap.gz,${wd}alignments/S20.CGmap.gz,${wd}alignments/S21.CGmap.gz,${wd}alignments/S22.CGmap.gz,${wd}alignments/S23.CGmap.gz,${wd}alignments/S24.CGmap.gz,${wd}alignments/S25.CGmap.gz,${wd}alignments/S26.CGmap.gz,${wd}alignments/S27.CGmap.gz,${wd}alignments/S28.CGmap.gz,${wd}alignments/S29.CGmap.gz,${wd}alignments/S2.CGmap.gz,${wd}alignments/S30.CGmap.gz,${wd}alignments/S31.CGmap.gz,${wd}alignments/S32.CGmap.gz,${wd}alignments/S33.CGmap.gz,${wd}alignments/S34.CGmap.gz,${wd}alignments/S35.CGmap.gz,${wd}alignments/S36.CGmap.gz,${wd}alignments/S37.CGmap.gz,${wd}alignments/S38.CGmap.gz,${wd}alignments/S39.CGmap.gz,${wd}alignments/S3.CGmap.gz,${wd}alignments/S40.CGmap.gz,${wd}alignments/S41.CGmap.gz,${wd}alignments/S42.CGmap.gz,${wd}alignments/S43.CGmap.gz,${wd}alignments/S44.CGmap.gz,${wd}alignments/S45.CGmap.gz,${wd}alignments/S46.CGmap.gz,${wd}alignments/S47.CGmap.gz,${wd}alignments/S48.CGmap.gz,${wd}alignments/S4.CGmap.gz,${wd}alignments/S5.CGmap.gz,${wd}alignments/S6.CGmap.gz,${wd}alignments/S7.CGmap.gz,${wd}alignments/S8.CGmap.gz,${wd}alignments/S9.CGmap.gz -S S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S1,S20,S21,S22,S23,S24,S25,S26,S27,S28,S29,S2,S30,S31,S32,S33,S34,S35,S36,S37,S38,S39,S3,S40,S41,S42,S43,S44,S45,S46,S47,S48,S4,S5,S6,S7,S8,S9 -CG -O ${wd}methylation_matrix.txt -min-coverage 5 -min-sample 0.7

# impute missing samples

python3 -m BSBolt Impute -M methylation_matrix.txt -W 3000 -k 3 -t 8 -O TBS_methylation_matrix.tsv

gzip TBS_methylation_matrix.tsv

