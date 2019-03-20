import glob, os


# Declare the resources being used
dbsnp_150 = "~/references/homo_sapiens/hg19/dbSNP_150/GATK/All_20170403.vcf.gz"
plugins_dir = "~/references/homo_sapiens/VEP/Plugins"
FATHMM_MKL_db = "~/references/homo_sapiens/VEP/Plugins/fathmm-MKL_Current.tab.gz"
cache_dir = "~/references/vep/92/cache"
Human_Gnome_ref = "~/references/ucsc.hg19.fasta"
ref = "~/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta"
fasta_refseq = "~/references/homo_sapiens/VEP/homo_sapiens/92_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
dbNSFP2_9_2 = "~/references/homo_sapiens/hg19/dbNSFP_2.9.2/dbNSFP2.9.2.txt.gz"
Human_Ref = "~/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta"
Human_Ref_bwa = "~/references/homo_sapiens/v37_decoy/bwa/human_g1k_v37_decoy_mod.fasta"
INDELS = "~/references/homo_sapiens/v37_decoy/gatk/1000G_phase1.indels.b37_mod.vcf"
MILLS = "~/references/homo_sapiens/v37_decoy/gatk/Mills_and_1000G_gold_standard.indels.b37_mod.vcf"
DBSNP = "~/references/homo_sapiens/v37_decoy/gatk/dbsnp_138.b37_mod.vcf"
HAPMAP = "~/references/homo_sapiens/v37_decoy/gatk/hapmap_3.3.b37_mod.vcf"
OMNI = "~/references/homo_sapiens/v37_decoy/gatk/1000G_omni2.5.b37_mod.vcf"
SNPS = "~/references/homo_sapiens/v37_decoy/gatk/1000G_phase1.snps.high_confidence.b37_mod.vcf"
dbsnp_147 = "~/references/homo_sapiens/hg19/dbSNP_147/All_20160601.vcf.gz"
GNOMAD_Genome = "~/references/homo_sapiens/hg19/gnomAD/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz"
GNOMAD_Exome = "~/references/homo_sapiens/hg19/gnomAD/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz"
GNOMAD_Genome_decomposed = "~/references/homo_sapiens/hg19/gnomAD/gnomad.genomes.r2.0.1.sites.noVEP_decomposed.vcf.gz"
GNOMAD_Exome_decomposed = "~/references/homo_sapiens/hg19/gnomAD/gnomad.exomes.r2.0.1.sites.noVEP_decomposed.vcf.gz"

# Tools directory
picard_dir = "/hpf/tools/centos6/picard-tools/2.5.0"
gatk_dir = "/hpf/tools/centos6/gatk/3.7.0"
snpEff="~/TOOLS/snpEff/4.3"


SAMP = [os.path.basename(sample) for sample in glob.glob('*.vcf.gz')]
SAMPLES = [sample.split('.vcf.gz')[0] for sample in SAMP]

CHROMOSOMES = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'.split()

#rule all:
#     input:"VEP_FILTERED/VepAnnotated_merged_GenomeWide_Rare1pcnt.vcf",
#           "VEP_FILTERED/VepAnnotated_merged_GenomeWide.vcf"

rule all:
     input:expand("VepAnnotation/{sample}_Vep92Ann_merged_genomewide.vcf",sample=SAMPLES)

rule splitByChr:
     input: "{sample}.vcf.gz"
     output: "SplitByChrom/{sample}_{chr}.vcf.gz"
     params: "{chr}"
     log:"LOGFILES/splitVCFbyChrom.log"
     shell:"""
           module load snpEff/4.3
           module load bcftools/1.3.1
           bcftools filter -r {params} -O z -o {output} {input} 2> {log}
           """

rule vep92_Ann:
     input: "SplitByChrom/{sample}_{chr}.vcf.gz"
     params: cache=cache_dir,
             Gno_Gen=GNOMAD_Genome,
             Gno_Exo=GNOMAD_Exome,
             Gno_Gen_decomp=GNOMAD_Genome_decomposed,
             Gno_Exo_decomp=GNOMAD_Exome_decomposed,
             Ref=ref
     log: "LOGFILES/{sample}_vep92_Ann_{chr}.log"
     output: "VepAnnotation/{sample}_{chr}_vep92Ann.vcf"

     shell: """
            module load vep/92
            module load java/1.8.0_91
            module load snpEff/4.3
            module load perl/5.20.1
            java -Xmx10G -jar $snpEff/SnpSift.jar annotate -v -noId -info AF,AC,Hom,AN,AS_FilterStatus -name gnomADgenome_ {params.Gno_Gen_decomp} {input}\
            | java -Xmx10G -jar $snpEff/SnpSift.jar annotate -v -noId -info AF,AC,Hom,AN,AS_FilterStatus -name gnomADexome_ {params.Gno_Exo_decomp} -\
            | vep --assembly GRCh37 --fasta {params.Ref} --no_progress --flag_pick --cache --dir {params.cache} --refseq --offline \
            -o {output} --everything --vcf --no_stats --format vcf --af_gnomad -custom {params.Gno_Gen},gnomADg,vcf,exact,0,AF,AC,AN,Hom 2> {log}
            """

rule joinVCF:
     input: expand("VepAnnotation/{sample}_{chr}_vep92Ann.vcf", chr=CHROMOSOMES)
     output: merged="VepAnnotation/{sample}_Vep92Ann_merged_genomewide.vcf"
     log: "LOGFILES/joinVCF.log"
     shell:"""
           module load snpEff/4.3
           module load java/1.8.0_91
           module load vep/90
           module load bcftools/1.4
           module load perl/5.20.1
           java -jar -Xmx20g /hpf/tools/centos6/snpEff/4.3/SnpSift.jar split -j {input} > {output.merged} 2> {log}
           perl -pi -e "s/^##INFO=<ID=gnomADg,Number=.,Type=String/##INFO=<ID=gnomADg,Number=A,Type=Integer/g" {output.merged}
           perl -pi -e "s/^##INFO=<ID=gnomADg_AF,Number=.,Type=String/##INFO=<ID=gnomADg,Number=A,Type=Integer/g" {output.merged}
           perl -pi -e "s/^##INFO=<ID=gnomADg_AC,Number=.,Type=String/##INFO=<ID=gnomADg,Number=A,Type=Integer/g" {output.merged}
           perl -pi -e "s/^##INFO=<ID=gnomADg_AN,Number=.,Type=String/##INFO=<ID=gnomADg,Number=A,Type=Integer/g" {output.merged}
            """
