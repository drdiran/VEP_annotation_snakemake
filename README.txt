########################################################
##### README file for VEP_annotation        ############
##### Developer: Oyediran Akinrinade, PhD ##############
########################################################

The pipeline runs vep92 (variant effects predictor) on WGS vcf file.
It first splits the genome-wide vcf file by chromosome and submits
chromosome level jobs to the cluster for annotation.
Each chromosome level annotated file is written to the VepAnnotation folder
together with merged annotated vcf file for downstream analyses.

To use, run:

bash submit_VepAnnotate_Gnomad_ExAC_Freq.sh &

NB: check file 'VepAnnotate_Gnomad_ExAC_Freq.log' for progress report on the
submitted jobs.
