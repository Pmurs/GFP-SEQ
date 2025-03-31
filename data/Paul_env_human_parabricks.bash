#!/bin/bash

# Set environment variables
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
export PATH=$PATH:/usr/local/bin

# Reference genome paths (update these paths as necessary)
export refgenome=/human/Homo_sapiens_assembly38.fasta
export refvar=/human
export dbsnp=/human/Homo_sapiens_assembly38.dbsnp138.vcf
export dbsnp_gz=/human/Homo_sapiens_assembly38.dbsnp138.vcf.gz
export chr_file=/human/chr_hg38.txt

export docker_refvar=/human:/human

# Make Picard available
export picard=/usr/local/bin/picard.jar

# Make GATK docker image available
export gatk="broadinstitute/gatk:latest"