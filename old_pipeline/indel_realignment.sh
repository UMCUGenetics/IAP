#!/bin/bash



/hpc/cog_bioinf/common_scripts/GATK_v2/bundle/b37
/hpc/cog_bioinf/GENOMES/Mus_musculus_GRCm38_GATK_illumina_bwa075/Mus_musculus_GRCm38.fasta
#IndelTargetCreator

java -Xmx2g -jar /hpc/cog_bioinf/common_scripts/GATK_v2/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
        -R ref.fasta \
        -I input.bam \
        -o forIndelRealigner.intervals \


	
    #IndelRealignment	
    java -Xmx4g -jar  /hpc/cog_bioinf/common_scripts/GATK_v2//GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-I dedup\_$i \
	-R /data/GENOMES/human_GATK_illumina_GRCh37/human_GATK_illumina_GRCh37.fasta \
	-targetIntervals /data/common_scripts/GATK_v2/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf.gz.intervals \
	-o realigned_dedup_$i \
	-known /data/common_scripts/GATK_v2/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
	--consensusDeterminationModel KNOWNS_ONLY \
	-LOD 0.4 \
	-L /data/ENRICH/SS_exome_50mb_V4_nochr.bed \
	2>&1 | tee -a $i.log
	