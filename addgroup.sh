#!/bin/bash
#SBATCH --partition=cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=40G
#SBATCH --cpus-per-task=5
#SBATCH --time=5-00:00:00
#SBATCH --output=serial_test_%j.log

gatk_bin="/gpfs/data/igorlab/software/GenomeAnalysisTK/gatk-4.1.9.0/gatk"
list="/gpfs/data/davolilab/projects/dbGAP/phs001572/bam/list"
bam_path="/gpfs/data/davolilab/projects/dbGAP/phs001572/bam/"
output="/gpfs/data/davolilab/projects/dbGAP/phs001572/bam/addgroup"
while read -r sample;
do
echo $sample
$gatk_bin AddOrReplaceReadGroups \
      I=$bam_path"/"$sample.bam \
      O=$output"/"$sample.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=$sample
done < $list

