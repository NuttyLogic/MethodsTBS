# import libraries
import os
import subprocess

import joblib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
from tqdm.notebook import tqdm


# use latex formatting for figures, latex must be on system path for this to work 
rc('text', usetex=False)

# set environment plotting params
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)


# set working directory
wd = '~/working_directory/'


samples = []
with open('samples.txt', 'r') as sample_list:
    for sample in sample_list:
        samples.append(sample.strip())


def parse_log(log_file):
    alignment_log = {}
    if not os.path.exists(log_file):
        return log_file.split('/')[-1], alignment_log
    # open log file as plain text file 
    with open(log_file, 'r') as log:
        for line in log:
            if 'Alignment Complete: Time ' in line:
                alignment_time = line.replace('Alignment Complete: Time ', '').strip().split(':')
                alignment_seconds = int(alignment_time[0]) * 60 * 60 + int(alignment_time[1]) * 60 + int(alignment_time[2])
                alignment_log['AlignmentTimeSec'] = alignment_seconds
            elif 'Total Reads: ' in line:
                alignment_log['TotalReads'] = int(line.replace('Total Reads: ', '').strip())
            elif 'Mappability' in line:
                alignment_log['Mappability'] = float(line.replace('Mappability: ', '').replace('get_ipython().run_line_magic("',", " '').strip())")
            elif 'Watson_C2T' in line:
                alignment_log['WatsonC2T'] = int(line.split(':')[1].strip())
            elif 'Watson_G2A' in line:
                alignment_log['WatsonG2A'] = int(line.split(':')[1].strip())
            elif 'Crick_C2T' in line:
                alignment_log['CrickC2T'] = int(line.split(':')[1].strip())
            elif 'Crick_G2A' in line:
                alignment_log['CrickG2A'] = int(line.split(':')[1].strip())
            elif 'Unmapped' in line:
                alignment_log['Unmapped'] = int(line.split(':')[1].strip())
            elif 'Ambiguous' in line:
                alignment_log['BSAmbiguous'] = int(line.split(':')[1].strip())
    return log_file.split('/')[-1], alignment_log


# parse all logs and save data in dict

sample_stats = {}

for sample in tqdm(samples):
    _, log = parse_log(f'{wd}alignments/{sample}.log')
    sample_stats[sample] = log


# transform dict into pandas dataframe and plot
stats_df = pd.DataFrame(sample_stats).T


fig, ax = plt.subplots(figsize=(12,12))

plot_cats = ['TotalReads','WatsonC2T','CrickC2T','WatsonG2A','CrickG2A','Unmapped','BSAmbiguous']
sns.boxplot(data=stats_df[plot_cats], ax=ax, palette='Set2', linewidth=2.5)
ax.set_title('BSBolt Log Stats')
#plt.savefig('bsb_log_stats.png', dpi=200, bbox_inches='tight')
plt.show()


def get_coverage(sample_path, target_bed, sample_name, include_dups=False):
    sample_coverage = []
    bed_args = ['bedtools', 'multicov', '-bams', sample_path, '-bed', target_bed]
    if include_dups:
        bed_args = ['bedtools', 'multicov', '-D', '-bams', sample_path, '-bed', target_bed]
    with subprocess.Popen(args=bed_args, stdout=subprocess.PIPE, universal_newlines=True) as bed_process:
        while True:
            bed_line = bed_process.stdout.readline()
            if not bed_line:
                break
            line_split = bed_line.replace('\n', '').split('\t')
            sample_coverage.append(dict(site=f'{line_split[0]}:{line_split[1]}-{line_split[2]}', coverage=int(line_split[4]), 
                                        probe_region=line_split[3], sample=sample_name, dup='De-Duplicated' if not include_dups else 'Duplicates'))
    return sample_coverage


target_bed = os.getcwd() + '/combined_target.bed.gz'


formatted_coverage = []

for sample in samples:
    formatted_coverage.append([f'{alignment_dir}{sample}.dup.bam', target_bed, sample])
    formatted_coverage.append([f'{alignment_dir}{sample}.dup.bam', target_bed, sample, True])


coverages = joblib.Parallel(n_jobs=16, verbose=10)(joblib.delayed(get_coverage)(*run) for run in formatted_coverage)


coverages = list(coverages)


dup_cov = {}
dedup_cov = {}

for coverage in coverages:
    sample = coverage[0]['sample']
    sample_coverage = {}
    for site in coverage:
        sample_coverage[site['site']] = site['coverage']
    if coverage[0]['dup'] == 'De-Duplicated':
        dedup_cov[sample] = sample_coverage
    else:
        dup_cov[sample] = sample_coverage


dup_cov_df = pd.DataFrame(dup_cov)


dedup_cov_df = pd.DataFrame(dedup_cov)


coverage_df = pd.DataFrame([dup_cov_df.mean(axis=1), dedup_cov_df.mean(axis=1)]).T


coverage_df.columns = ['Dup Mean', 'DeDup Mean']


coverage_df['Log2_cov Dup'] = np.log2(coverage_df['Dup Mean'] + 1)
coverage_df['Log2_cov DeDup'] = np.log2(coverage_df['DeDup Mean'] + 1)


fig, ax = plt.subplots(figsize=(12,12))

sns.distplot(coverage_df['Log2_cov Dup'].values, label='Dup. Reads Included', ax=ax)
sns.distplot(coverage_df['Log2_cov DeDup'].values, label='Dup. Reads Removed', ax=ax)

plt.legend(fontsize=18)
ax.set_xlabel(f'$Log_2(coverage + 1)$', fontsize=24)
ax.set_ylabel(f'$Density$', fontsize=24)

plt.show()


# run samtools flagstat on each sample and parse output 
# command will fail if samtools is not on path
def get_flagstats(sample_path, sample_name):
    sam_stats = {}
    sam_args = ['samtools', 'flagstat', sample_path]
    with subprocess.Popen(args=sam_args, stdout=subprocess.PIPE, universal_newlines=True) as sam_process:
        for count, line in enumerate(iter(sam_process.stdout)):
            if count == 0:
                sam_stats['TotalReads'] = int(line.split(' ')[0])
            elif count == 1:
                sam_stats['SecondaryAlignment'] = int(line.split(' ')[0])
            elif count == 2:
                sam_stats['Supplementary']= int(line.split(' ')[0])
            elif count == 3:
                sam_stats['Duplicate']= int(line.split(' ')[0])
            elif count == 4:
                sam_stats['Mapped']= int(line.split(' ')[0])
            elif count == 5:
                sam_stats['Paired']= int(line.split(' ')[0])
            elif count == 6:
                sam_stats['Read1']= int(line.split(' ')[0])
            elif count == 7:
                sam_stats['Read2']= int(line.split(' ')[0])
            elif count == 8:
                sam_stats['ProperPair']= int(line.split(' ')[0])
            elif count == 9:
                sam_stats['MateMapped']= int(line.split(' ')[0])
            elif count == 10:
                sam_stats['Singletons']= int(line.split(' ')[0])
            elif count == 11:
                sam_stats['MateMappedToDiffChrom']= int(line.split(' ')[0])
            elif count == 12:
                sam_stats['MateMappedToDiffChromMapq>=5']= int(line.split(' ')[0])
    return sample_name, sam_stats


# use multiple processing threads to parallelize
flag_stats = joblib.Parallel(n_jobs=14, verbose=10)(joblib.delayed(get_flagstats)(*[f'{wd}alignments/{sample}.dup.bam', sample]) for sample in samples)


# format flagstat results

processed_flag_stats = {}

for sample, f_stats in flag_stats:
    processed_flag_stats[sample] = f_stats


# make pandas data frame 
flag_df = pd.DataFrame(processed_flag_stats).T


# add duplication rate statistic
flag_df['DuplicationRate'] = flag_df['Duplicate'] / flag_df['Mapped'] 


## Get number of reads mapping to coverage region
def get_on_target_coverage(sample_name):
    dup_removed_count =  dedup_cov_df[sample_name].sum()
    dup_count = dup_cov_df[sample_name].sum()
    return dup_count, dup_removed_count


on_target_count_dedup = []
on_target_count_dup = []

for sample in flag_df.index:
    dup_count, dup_removed_count = get_on_target_coverage(sample)
    on_target_count_dup.append(dup_count)
    on_target_count_dedup.append(dup_removed_count)


flag_df['OnTarget-All'] = on_target_count_dup
flag_df['OnTarget-DeDuplicated'] = on_target_count_dedup
flag_df['Mapped-Deduplicated'] = flag_df['Mapped'] - flag_df['Duplicate']


fig, ax = plt.subplots(figsize=(12,12))

plot_cats = ['TotalReads', 'Duplicate', 'Mapped', 'Mapped-Deduplicated', 'Paired', 'ProperPair', 'OnTarget-All', 'OnTarget-DeDuplicated']
sns.boxplot(data=flag_df[plot_cats], ax=ax, palette='Set2', linewidth=2.5)
ax.set_title('Samtools Flagstat', fontsize=20)
ax.set_xticklabels(ax.get_xticklabels(),rotation=30)
#plt.savefig('flagstats.png', dpi=200, bbox_inches='tight')
plt.show()


fig, ax = plt.subplots(figsize=(12,12))

plot_cats = ['DuplicationRate']
sns.swarmplot(data=flag_df[plot_cats], ax=ax, palette='Set2', linewidth=2.5)
ax.set_title('Observed Duplication')
#plt.savefig('duplication.png', dpi=200, bbox_inches='tight')
plt.show()
