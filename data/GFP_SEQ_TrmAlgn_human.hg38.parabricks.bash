#!/bin/bash
#SBATCH -t 2-00:00:00
##SBATCH -p normal
#SBATCH -J gatk
#SBATCH -o gatk_%A.log
#SBATCH -n 8
##SBATCH --mem=32GB
##SBATCH --exclude=cpu-752


echo "DNA mapping for human start on"; date; echo "====="

# global parameters 
sampleroot=..
samplename=$1
np=8

source Paul_env_human_parabricks.bash

cd $sampleroot
echo "sampleroot is:"
pwd
echo "samplename=$samplename;np=$np"
echo "The reference will be use $refgenome"

# set the work path
mkdir -p ./processing/1-fastqc
mkdir -p ./processing/2-trim
mkdir -p ./processing/3-fastqc
mkdir -p ./processing/4-bwa


# 1. fastqc
echo "1th fastqc start on"; date; echo "====="

fastqc -o ./processing/1-fastqc \
./data/$samplename.fq \

echo "====="; echo "1th fastqc end on"; date


# 2. trim, manual correct rrbs enzyme parameters
echo "2th trim galore start on"; date; echo "====="

trim_galore --illumina -o ./processing/2-trim \
./data/$samplename.fq \

mv ./processing/2-trim/"${samplename}"_trimmed.fq ./processing/2-trim/"${samplename}".trimmed.fq

echo "====="; echo "2th trim galore end on"; date


# 3. fastqc
echo "3th fastqc start on"; date; echo "====="

fastqc -o ./processing/3-fastqc \
./processing/2-trim/$samplename.trimmed.fq

echo "====="; echo "3th fastqc end on"; date

# 4. bwa
echo "4th bwa start on"; date; echo "====="

bwa mem -t $np -M -S $refgenome \
./processing/2-trim/$samplename.trimmed.fq \
> ./processing/4-bwa/$samplename-pe.sam

echo "====="; echo "4th bwa end on"; date

# 5. picard
cd processing/4-bwa

mkdir tmp
chmod 777 tmp

# 5.1 convert sam to bam; then sort and remove duplication read
samtools view -Sb -@ $np $samplename-pe.sam > $samplename-pe.pre.bam # Parabricks outputs aligned file as a bam already
samtools sort -@ $np $samplename-pe.pre.bam > $samplename-pe.sorted.bam
samtools rmdup $samplename-pe.sorted.bam $samplename-pe.sorted.rmdup.bam
#ln -s $samplename-pe.sorted.bam $samplename-pe.sorted.rmdup.bam

# 5.2 Picard
java -Xmx20g -jar $picard CleanSam \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=`pwd`/tmp \
I=$samplename-pe.sorted.rmdup.bam \
O=$samplename-pe.bam

java -Xmx20g -jar $picard AddOrReplaceReadGroups \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=`pwd`/tmp \
I=$samplename-pe.bam \
O=$samplename-pe.rd.bam \
RGID=$samplename RGLB=$samplename RGPL=illumina RGPU=$samplename RGSM=$samplename

java -Xmx20g -jar $picard FixMateInformation \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=`pwd`/tmp \
I=$samplename-pe.rd.bam \
O=$samplename-pe.rd.fixed.bam \
SO=coordinate

java -Xmx20g -jar $picard SortSam \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=`pwd`/tmp \
I=$samplename-pe.rd.fixed.bam \
O=$samplename-pe.rd.fixed.sorted.bam \
SO=coordinate \
CREATE_INDEX=true

java -Xmx20g -jar $picard ReorderSam \
TMP_DIR=`pwd`/tmp \
I=$samplename-pe.rd.fixed.sorted.bam \
O=$samplename-pe.rd.fixed.sorted.reorder.bam \
CREATE_INDEX=true \
REFERENCE_SEQUENCE=$refgenome \
SEQUENCE_DICTIONARY="${refgenome%.fasta}.dict"
#
## 6.1 Indel Realigner
#sudo docker run --rm -v /media/sf_SMM_SEQ/human:/data -w /data $gatk \
#gatk \
#--java-options "-Xmx20g" \
#RealignerTargetCreator \
#-R $refgenome \
#-O $samplename-pe.rd.fixed.sorted.reorder.bam.intervals \
#-I $samplename-pe.rd.fixed.sorted.reorder.bam \
#--known $refvar/Homo_sapiens_assembly38.known_indels.vcf.gz \
#--known $refvar/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#
#sudo docker run --rm -v /media/sf_SMM_SEQ/human:/data -w /data $gatk \
#gatk \
#--java-options "-Xmx20g" \
#IndelRealigner \
#-R $refgenome \
#-I $samplename-pe.rd.fixed.sorted.reorder.bam \
#--targetIntervals $samplename-pe.rd.fixed.sorted.reorder.bam.intervals \
#--known $refvar/Homo_sapiens_assembly38.known_indels.vcf.gz \
#--known $refvar/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
#-O $samplename-pe.rd.fixed.sorted.reorder.realign.bam

# 6.2 Base Recalibration

# BaseRecalibrator
sudo docker run --rm -v /human:/human -v $(pwd)/:/data -w /data $gatk \
  gatk BaseRecalibrator \
  --java-options "-Xmx20g" \
  -R $refgenome \
  -I $samplename-pe.rd.fixed.sorted.reorder.bam \
  --known-sites $refvar/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites $refvar/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites $refvar/Homo_sapiens_assembly38.dbsnp138.vcf \
  -O $samplename-pe.rd.fixed.sorted.reorder.realign.bam.recal.grp >> gatk_%A.log
# currently using an older version of the dbSNP because that's what's available for free on S3

# ApplyBQSR
sudo docker run --rm -v /human:/human -v $(pwd)/:/data -w /data $gatk \
  gatk ApplyBQSR \
  --java-options "-Xmx20g" \
  --max-cycle 1000 \
  -R $refgenome \
  -I $samplename-pe.rd.fixed.sorted.reorder.bam \
  --bqsr-recal-file $samplename-pe.rd.fixed.sorted.reorder.realign.bam.recal.grp \
  -O $samplename-pe.rd.fixed.sorted.reorder.realign.recal.bam >> gatk_%A.log

ln -s $samplename-pe.rd.fixed.sorted.reorder.realign.recal.bam $samplename.gatk.bam
samtools index $samplename.gatk.bam

echo "DNA mapping for human end on"; date; echo "====="

# Useful for local testing

ln -s ./processing/4-bwa/$samplename.gatk.bam ../../$samplename.gatk.bam
ln -s ./processing/4-bwa/$samplename.gatk.bam.bai ../../$samplename.gatk.bam.bai

# 7. mpileup prepared for varscan2
# samtools mpileup -f $refgenome2 $samplename.gatk.bam > $samplename.gatk.mpu
# echo "mpileup for human end on"; date; echo "====="


