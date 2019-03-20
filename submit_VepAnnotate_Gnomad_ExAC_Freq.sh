#!/bin/sh
module load python/3.5.2
module load R/3.2.0

#
work_dir=`pwd`
snakefile="VepAnnotate_Gnomad_ExAC_Freq.snakefile"
## Make project result directories
#mkdir FASTQ
#mkdir BAM
mkdir LOGFILES
#mkdir Raw_gVCF
#mkdir Recal_Data

chmod -R 770 $work_dir

## Execute SnakeMake , GATK pipeline

snakemake \
 --nolock \
 --jobname 'SM.{rulename}.{jobid}' \
 --directory $work_dir \
 --snakefile $snakefile \
 -k -p -w 10 \
 --rerun-incomplete \
 --stats $work_dir/GATK_genotypeGVCF_MultiSample.stats \
 -j 50 \
 --timestamp \
 --cluster "qsub -l vmem=60g -l walltime=23:00:00 \
 -e LOGFILES -o LOGFILES" >& VepAnnotate_Gnomad_ExAC_Freq.log
