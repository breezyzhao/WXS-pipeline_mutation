#!/bin/bash
#SBATCH --partition=cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=8
#SBATCH --time=5-00:00:00
#SBATCH --output=serial_test_%j.log


bed="/gpfs/data/davolilab/backup/zhaox12/referencce/b37/genome.bed"
ref_fasta="/gpfs/data/davolilab/backup/zhaox12/referencce/b37/references_b37_Homo_sapiens_assembly19.fasta"
#hg38 
#ref_fasta="/gpfs/data/davolilab/backup/zhaox12/referencce/ucsc-hg38/hg38.order/genome.fa"
list="/gpfs/data/davolilab/projects/dbGAP/test/list"
bam_path="/gpfs/data/davolilab/projects/dbGAP/test/bam/"
proj_dir="/gpfs/data/davolilab/projects/dbGAP/test/GATK"

module add samtools/1.9-new

while read -r sample;
do
echo $sample
# call variants with Mutect2 (2.1 part of GATK 4)
# check for correct number of arguments
# arguments
name_t="${sample}.T"
bam_t="${bam_path}${sample}.T.bam"
sample_t="$(samtools view -H $bam_t | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)"
bam_n="${bam_path}${sample}.N.bam"
name_n="${sample}.N"
sample_n="$(samtools view -H $bam_n | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" |uniq)"
threads=10


#########################


# settings and files

sample="${sample_t}:${sample_n}"

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
contamination_csv="${summary_dir}/${sample_t}-${sample_n}.VCF-Mutect2-contamination.csv"

vcf_dir="${proj_dir}/VCF-Mutect2"
mkdir -p "$vcf_dir"

mutect_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$mutect_logs_dir"

vcf_unfiltered="${mutect_logs_dir}/${sample_t}-${sample_n}.unfiltered.vcf"
idx_unfiltered="${vcf_unfiltered}.idx"

vcf_filtered="${vcf_dir}/${sample_t}-${sample_n}.original.vcf.gz"
idx_filtered="${vcf_filtered}.tbi"

vcf_fixed="${vcf_dir}/${sample_t}-${sample_n}.vcf"

pileup_table_t="${mutect_logs_dir}/${sample_t}-${sample_n}.t.pileup.txt"
pileup_table_n="${mutect_logs_dir}/${sample_t}-${sample_n}.n.pileup.txt"
contamination_table="${mutect_logs_dir}/${sample_t}-${sample_n}.contamination.txt"

anno_dir="${proj_dir}/annotation"
mkdir -p "$anno_dir"
annotation="${anno_dir}/${name_t}.funcotator.maf"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# check for output

# skip to annotation if final VCF exists already
if [ -s "$vcf_fixed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n CMD: $annot_cmd \n"
	($annot_cmd)
	exit 0
fi

# delete unfiltered VCF (likely incomplete since the final fixed VCF was not generated)
if [ -s "$vcf_unfiltered" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT VCF $vcf_unfiltered EXISTS \n" >&2
	rm -fv "$vcf_unfiltered"
	rm -fv "$idx_unfiltered"
fi

# delete filtered VCF (likely incomplete since the final fixed VCF was not generated)
if [ -s "$vcf_filtered" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT VCF $vcf_filtered EXISTS \n" >&2
	rm -fv "$vcf_filtered"
	rm -fv "$idx_filtered"
fi


#########################


# GATK settings

# known variants
germline_resource_arg="--germline-resource /gpfs/data/davolilab/backup/zhaox12/referencce/b37/af-only-gnomad.raw.sites.vcf"
pileup_variants="/gpfs/data/davolilab/backup/zhaox12/referencce/b37/somatic-b37-small_exac_common_3.vcf"
#/gpfs/data/igorlab/ref/hg38/gatk-bundle/af-only-gnomad.hg38.vcf.gz
#/gpfs/data/igorlab/ref/hg38/gatk-bundle/small_exac_common_3.hg38.vcf.gz


#########################


# GATK Mutect2

module add python/cpu/3.6.5
module add samtools/1.9-new
# command
gatk_bin="/gpfs/data/igorlab/software/GenomeAnalysisTK/gatk-4.1.9.0/gatk"

echo
echo " * GATK: $(readlink -f $gatk_bin) "
echo " * GATK version: $($gatk_bin Mutect2 --version 2>&1 | grep "Version") "
echo " * sample T : $sample_t "
echo " * BAM T : $bam_t "
echo " * sample N : $sample_n "
echo " * BAM N : $bam_n "
echo " * intervals BED: $bed "
echo " * VCF unfiltered: $vcf_unfiltered "
echo

# --native-pair-hmm-threads: how many threads should a native pairHMM implementation use
# --max-reads-per-alignment-start: maximum number of reads to retain per alignment start position (50)
# --dont-use-soft-clipped-bases: do not analyze soft clipped bases in the reads
# --germline-resource: population vcf of germline sequencing containing allele fractions

mutect_cmd="
$gatk_bin --java-options \"-Xms8G -Xmx8G\" Mutect2 \
--seconds-between-progress-updates 600 \
--native-pair-hmm-threads $threads \
--reference $ref_fasta \
$germline_resource_arg \
--dont-use-soft-clipped-bases \
--max-reads-per-alignment-start 100 \
--intervals $bed \
--interval-padding 10 \
--input $bam_t \
--input $bam_n \
--tumor-sample $sample_t \
--normal-sample $sample_n \
--output $vcf_unfiltered \
"
echo -e "\n CMD: $mutect_cmd \n"
eval "$mutect_cmd"


#########################


# check that output generated

# check if VCF file is present
if [ ! -s "$vcf_unfiltered" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_unfiltered NOT GENERATED \n" >&2
	exit 1
fi

# check if VCF index is present (should be present if VCF is complete)
if [ ! -s "$idx_unfiltered" ] ; then
	echo -e "\n $script_name ERROR: VCF IDX $idx_unfiltered NOT GENERATED \n" >&2
	# delete VCF since something went wrong and it might be corrupted
	rm -fv "$vcf_unfiltered"
	exit 1
fi


#########################


# estimate cross-sample contamination for genomes with available resources

# GetPileupSummaries: tabulates read counts that support ref and alt alleles for the sites in the resource
# CalculateContamination: calculate the fraction of reads coming from cross-sample contamination

if [ -n "$pileup_variants" ] ; then

	# tabulates pileup metrics for T
	pileup_t_cmd="
	$gatk_bin --java-options \"-Xms8G -Xms8G\" GetPileupSummaries \
	--verbosity WARNING \
	--variant $pileup_variants \
	--intervals $bed \
	--interval-padding 100 \
	--input $bam_t \
	--output $pileup_table_t \
	"
	echo -e "\n CMD: $pileup_t_cmd \n"
	eval "$pileup_t_cmd"

	# tabulates pileup metrics for N
	pileup_n_cmd="
	$gatk_bin --java-options \"-Xms8G -Xms8G\" GetPileupSummaries \
	--verbosity WARNING \
	--variant $pileup_variants \
	--intervals $bed \
	--interval-padding 100 \
	--input $bam_n \
	--output $pileup_table_n \
	"
	echo -e "\n CMD: $pileup_n_cmd \n"
	eval "$pileup_n_cmd"

	# the resulting table provides the fraction contamination, one line per sample
	contamination_cmd="
	$gatk_bin --java-options \"-Xms8G -Xms8G\" CalculateContamination \
	--verbosity WARNING \
	--input $pileup_table_t \
	--matched-normal $pileup_table_n \
	--output $contamination_table \
	"
	echo -e "\n CMD: $contamination_cmd \n"
	eval "$contamination_cmd"

	# check if contamination table is present
	if [ ! -s "$contamination_table" ] ; then
		echo -e "\n $script_name ERROR: CONTAMINATION TABLE $contamination_table NOT GENERATED \n" >&2
		# delete pileups since something went wrong and they might be corrupted
		rm -fv "$pileup_table_t"
		rm -fv "$pileup_table_n"
		exit 1
	fi

	# set the parameter for next step
	contamination_arg="--contamination-table $contamination_table"

	# contamination summary
	contamination_score=$(cat "$contamination_table" | grep -v 'contamination' | head -1 | cut -f 2)
	echo "contamination original: $contamination_score"
	contamination_score=$(printf "%.5f" "$contamination_score")
	echo "contamination rounded: $contamination_score"
	contamination_header="#SAMPLE,contamination"
	echo "${contamination_header}" > "$contamination_csv"
	echo "${sample},${contamination_score}" >> "$contamination_csv"
	sleep 5
	cat ${summary_dir}/*.VCF-Mutect2-contamination.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.VCF-Mutect2-contamination.csv"

fi


#########################


# GATK FilterMutectCalls

echo
echo " * VCF unfiltered: $vcf_unfiltered "
echo " * VCF filtered: $vcf_filtered "
echo

# --unique-alt-read-count: filter a variant if fewer than this many unique reads supporting the alternate allele (0)

mutect_filter_cmd="
$gatk_bin --java-options \"-Xms8G -Xms8G\" FilterMutectCalls \
--verbosity WARNING \
--reference $ref_fasta \
--unique-alt-read-count 5 \
--variant $vcf_unfiltered \
$contamination_arg \
--output $vcf_filtered \
"
echo -e "\n CMD: $mutect_filter_cmd \n"
eval "$mutect_filter_cmd"


#########################


# check that output generated

# check if VCF file is present
if [ ! -s "$vcf_filtered" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_filtered NOT GENERATED \n" >&2
	exit 1
fi

# check if VCF index is present (should be present if VCF is complete)
if [ ! -s "$idx_filtered" ] ; then
	echo -e "\n $script_name ERROR: VCF IDX $idx_filtered NOT GENERATED \n" >&2
	# delete VCF since something went wrong and it might be corrupted
	rm -fv "$vcf_filtered"
	exit 1
fi


#########################


# adjust VCF for ANNOVAR compatibility (http://annovar.openbioinformatics.org/en/latest/articles/VCF/)

module add samtools/1.9

echo
echo " * samtools: $(readlink -f $(which samtools)) "
echo " * samtools version: $(samtools --version | head -1) "
echo " * VCF filtered: $vcf_filtered "
echo " * VCF fixed: $vcf_fixed "
echo

# 1) keep header and only passing variants
# 2) split multi-allelic variants calls into separate lines (uses VCF 4.2 specification)
# 3) perform indel left-normalization (start position shifted to the left until it is no longer possible to do so)

fix_vcf_cmd="
zcat $vcf_filtered \
| grep -E '^#|PASS' \
| bcftools norm --multiallelics -both --output-type v - \
| bcftools norm --fasta-ref $ref_fasta --output-type v - \
> $vcf_fixed
"
echo -e "\n CMD: $fix_vcf_cmd \n"
eval "$fix_vcf_cmd"


#########################


# check that output generated

if [ ! -s "$vcf_fixed" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_fixed NOT GENERATED \n" >&2
	exit 1
fi


#########################


# clean up

rm -fv "$vcf_unfiltered"
rm -fv "$idx_unfiltered"
rm -fv "$pileup_table_t"
rm -fv "$pileup_table_n"


#########################
# funcotator
funcotator_cmd="
$gatk_bin --java-options \"-Xms8G -Xms8G\" Funcotator \
     --variant $vcf_fixed \
     --reference $ref_fasta \
     --ref-version hg19 \
     --data-sources-path /gpfs/data/davolilab/backup/zhaox12/referencce/GATK/funcotator_dataSources.v1.6.20190124s \
     --output $annotation \
     --output-file-format MAF \
     --add-output-vcf-command-line false \
     --force-b37-to-hg19-reference-contig-conversion true \
     --annotation-default normal_barcode:$name_n \
     --annotation-default tumor_barcode:$name_t \
     --annotation-default Center:broad.mit.edu \
     --add-output-sam-program-record false
"
echo -e "\n CMD: $funcotator_cmd \n"
eval "$funcotator_cmd"


#########################

done < $list

# end
