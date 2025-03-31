#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -J nanopore_mapping
#SBATCH -o nanopore_mapping_%A.log
#SBATCH -n 8
#SBATCH --mem=32G

echo "Nanopore DNA mapping pipeline started on:"; date
echo "========================================="

# === LOAD ENVIRONMENT ===
source Paul_env_human_parabricks.bash  # <-- loads refgenome, refvar, etc.

# === USER PARAMETERS ===
samplename=$1         # Sample name (expects ./data/${samplename}.fq)
np=8                  # Number of threads

echo "Sample: $samplename"
echo "Reference genome: $refgenome"
echo "Ref var path: $refvar"

# === DIRECTORY SETUP ===
mkdir -p processing/1-fastqc
mkdir -p processing/2-filtered
mkdir -p processing/3-alignment
mkdir -p processing/4-variant

# === 1. FASTQC (optional quality check) ===
echo "Running FastQC..."; date
fastqc -o processing/1-fastqc ./data/${samplename}.fq

# === 2. FILTERING (optional) ===
echo "Filtering reads with Filtlong..."; date
filtlong --min_length 400 ./data/${samplename}.fq > processing/2-filtered/${samplename}.filtered.fq

# === 3. ALIGNMENT with Minimap2 ===
echo "Aligning with Minimap2..."; date
minimap2 -t $np -ax map-ont $refgenome processing/2-filtered/${samplename}.filtered.fq \
  | samtools view -b - \
  | samtools sort -@ $np -o processing/3-alignment/${samplename}.sorted.bam

samtools index processing/3-alignment/${samplename}.sorted.bam

# === 4. VARIANT CALLING with Medaka ===
echo "Running Medaka variant calling..."; date
medaka_consensus \
  -i processing/2-filtered/${samplename}.filtered.fq \
  -d $refgenome \
  -o processing/4-variant/medaka_${samplename} \
  -t $np \
  -m r941_min_high_g360  # Change model if needed

# === DONE ===
echo "Pipeline complete for $samplename"
echo "Finished on:"; date
