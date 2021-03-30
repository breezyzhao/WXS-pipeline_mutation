# WGS-pipeline (mutation)
### 1. Alignment (BWA) ---bash :http://bio-bwa.sourceforge.net/bwa.shtml
module add bwa/0.7.17 <br />
module add samtools/1.9-new <br />

im1_files=(/gpfs/data/davolilab/DNAseq_Jan2020/WGS/rawdata/~*~/~*~.fastq) <br />
im3_files=(/gpfs/data/davolilab/DNAseq_Jan2020/WGS/rawdata/~*~/) <br />
im4_files=(/gpfs/data/davolilab/DNAseq_Jan2020/WGS/rawdata/~*~/aln-se.sam) <br />

for ((i=0;i<=${#im1_files[@]};i++)); <br />
do <br />
bwa mem /gpfs/scratch/zhaox12/referencce/ucsc-hg38/hg38.fa "${im1_files[i]}" > "${im3_files[i]}"aln-se.sam <br />
samtools view -S "${im3_files[i]}"$i.aln.sam -b | samtools sort -o "${im3_files[i]}"$i.sorted.bam <br />
samtools rmdup "${im3_files[i]}"$i.sorted.bam "${im3_files[i]}"$i.sorted.rmdup.bam ; samtools index "${im3_files[i]}"$i.sorted.rmdup.bam <br />
done <br />

### 2.sh mutect2-funco.sh https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-





