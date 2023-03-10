# Bovine Transcriptome Assembly and Annotation

## Bioinformatics Methods and Analysis Steps

### RNA-seq Data Analysis

* **Trim Adapters and Low-quality Bases**
```
module load trimgalore

rnaseq_dir='/usr/RNA-seq'
mkdir trimmed
cd ${rnaseq_dir}/trimmed
for i in $(ls ${rnaseq_dir}/*.fastq)
	do
		trim_galore ${rnaseq_dir}/${i} --quality 20 --length 20;
done;
cd ${rnaseq_dir}

```

* **Map High Quality Reads to Genome (ARS-UCD1.2)**
```
module load STAR

for i in $(ls ${rnaseq_dir}/trimmed/*.fastq)
	do
		tissue=$(echo ${i} | sed 's/\.fastq/g')
		STAR --runThreadN 40 --genomeDir /usr/ARS-UCD1.2.RepeadMasked.STAR/ \
		--readFilesIn ${i} --outFileNamePrefix ${tissue} \
		--outFilterMatchNminOverLread 0.9 --outFilterMismatchNoverLmax 0.05
		--outSAMtype BAM SortedByCoordinate
done;
```

* **Assemble Transcripts**
``` 
# normalise reads to maximum coverage of 50
module load trinity
mkdir transcriptomes
cd ${transcriptomes}
for i in $(ls ${rnaseq_dir}/trimmed/*.fastq)
	do
		tissue=$(echo ${i} | sed 's/\.fastq/g')
		insilico_read_normalization.pl --seqType fq --single ${i} --JM 350G \
		--max_cov 50  --output ${tissue}_normalised
		Trinity --seqType fq --SS_lib_type F --single ${tissue}_normalised/normalised_reads.fq \
		--CPU 36  --max_memory 180G --output ${id}_trinity --bflyGCThreads 31 \
		--bflyCPU 31 --no_distributed_trinity_exec
done;

# download HpcGridRunner
wget https://github.com/HpcGridRunner/HpcGridRunner/archive/refs/tags/v1.0.1.tar.gz
tar -xf v1.0.1.tar.gz

# denovo transcriptome assembly
module purge
module load miniconda
source activate trinityenv
module load parallel
for i in $(ls ${transcriptomes} | grep trinity$)
	do
		tissue=$(echo ${i} | sed 's/\_trinity//g')
		perl ${transcriptomes}/HpcGridRunner-master/hpc_cmds_GridRunner.pl --grid_conf grid.conf \
		-c ${transcriptomes}/${tissue}_trinity/recursive_trinity.cmds
		find ${transcriptomes}/${tissue}_trinity/read_partitions/ -name '*inity.fasta'  \
		| parallel --pipe -N1000 partitioned_trinity_aggregator.pl \
		TRINITY_DN > ${transcriptomes}/${tissue}_trinity//Trinity.tmp
		mv ${transcriptomes}/${tissue}_trinity/Trinity.tmp \
		${transcriptomes}/${tissue}_trinity/Trinity.fasta
done;

# mapp assembled transcript reads to bovine genome
module purge
module load gmap-gsnap-legacy
for i in $(ls ${transcriptomes} | grep trinity$)
	do
		tissue=$(echo ${i} | sed 's/\_trinity//g')
		gmap -D /usr/ARS-UCD1.2.RepeadMasked.GMAP -d GMAP_bostaurus \
		-f samse --min-trimmed-coverage 0.90 --min-identity 0.95 \
		${transcriptomes}/${tissue}_trinity/Trinity.fasta \
		>${transcriptomes}/${tissue}_trinity/Trinity.sam 2 \
		>${transcriptomes}/${tissue}_trinity/Trinity.sam.log
		sort -k 3,3 -k 4,4n Trinity.sam >Trinity_sorted.sam
done;

# collapse & group reads into putative gene models
module purge
module load smrtlink
module load cufflinks
module load bedtools2
genome_dir="/usr/ARS-UCD1.2.RepeadMasked.GMAP"
for i in $(ls ${transcriptomes} | grep trinity$)
	do
		tissue=$(echo ${i} | sed 's/\_trinity//g')
		collapse_isoforms_by_sam.py --max_fuzzy_junction 15 \
		--min-coverage 0.90 --min-identity 0.95 --input \
		${transcriptomes}/${tissue}_trinity/Trinity.fasta \
		-s ${transcriptomes}/${tissue}_trinity/Trinity_sorted.sam \
		-o ${tissue}
		sed -i 's/chr//g' ${tissue}.collapsed.gff
		gffread -w ${tissue}_transcripts.fa -g ${genome_dir}/bt_ref_ARS-UCD1.2.RepeadMasked.fa \
		${tissue}.collapsed.gff
done;

# check splice-junction validity
# check transcript coverage
# filter retained intron transcripts and non-canonical splice-junctions
# filter genomic DNA reads
module purge
module load tophat	
for i in $(ls ${transcriptomes} | grep trinity$)
	do
		tissue=$(echo ${i} | sed 's/\_trinity//g')
		tophat --library-type fr-firststrand -p 36 -o ${tissue}_tophat \
		-G ${tissue}.collapsed.gff --transcriptome-only -j junctions.juncs \
		/usr/ARS-UCD1.2.RepeadMasked.bowtie2/ARS-UCD1.2 ${tissue}.fastq		
done;
for i in $(ls ${transcriptomes} | grep trinity$)
	do
		#check splice-junction validity
		module purge
		module load  py-pandas/0.21.1-py2-326uzkn
		module load py-biopython/1.65-py2-fgq2vqv
		tissue=$(echo ${i} | sed 's/\_trinity//g')
		python get_introns_bed.py ${tissue}.collapsed.gff
		sed 's/\,/\t/g' ${tissue}_tophat/junctions.bed | \
		awk 'BEGIN { OFS = "\t" } { $17 = $2 + $13 +1} 1' | \
		awk 'BEGIN { OFS = "\t" } { $18 = $3 - $14} 1' | \
		awk 'NR>1{print $1, $17, $18, $5}' | tr ' ' '\t' | \
		awk '$4>=1' >${tissue}_valid_splice_tmp
		python merger.py ${tissue}_introns.bed ${tissue}_valid_splice_tmp
  		python merger.py ${tissue}_introns.bed ${tissue}_valid_splice_tmp
  		awk 'NF==8'  ${tissue}_output_tmp >${tissue}_output_tmp_2
  		awk 'NR>1 {print $4}' ${tissue}_output_tmp_2 | sed 's/intron.*$//g' | \ sort \
  			| uniq -c >${tissue}_info_tmp
  		awk '{print $4}' ${tissue}_introns.bed | sed 's/intron.*//g' | sort | uniq -c \
  			>${tissue}_a_tmp
  		join -j 2 -o 1.2 1.1 2.1 ${tissue}_info_tmp ${tissue}_a_tmp >${tissue}_c_tmp
  		awk '$2==$3 {print $1}' ${tissue}_c_tmp >${tissue}_valid_splice_transcripts  		
		#check transcript coverage
  		module purge
  		module load parallel
  		module load singularity
  		module load bedtools2
  		module load samtools/1.9-k6deoga
  		module load mosdepth		
  		awk '$3="transcript" {print $12}' ${tissue}.collapsed.gff | sed 's/\"\;//g; s/\"//g' \
  			| sort | uniq>${tissue}_tmp_transcripts
  		awk '$3=="exon" {print $12}' ${tissue}.collapsed.gff | sed 's/\"\;//g; s/\"//g' | \
  			sort | uniq -c>${tissue}_tmp_transcript.exons
  		awk '$3=="exon"' ${tissue}.collapsed.gff | awk '{print $1, $4,$5, $12, 99,$7, $12}' \
  			| sed 's/ /\t/g; s/\"\;//g; s/\"//g' >${tissue}_tmp_exons.bed
  		awk '$6=="+"' ${tissue}_tmp_exons.bed >${tissue}_tmp_forward.exons.bed
  		awk '$6=="-"' ${tissue}_tmp_exons.bed >${tissue}_tmp_reverse.exons.bed
  		samtools index ${tissue}_Aligned.out.bam
		mosdepth --by ${tissue}_tmp_forward.exons.bed ${tissue}_tmp_forward.coverage-output \
			${tissue}_Aligned.out.bam
		mosdepth --by ${tissue}_tmp_reverse.exons.bed ${tissue}_tmp_reverse.coverage-output \
			${tissue}_Aligned.out.bam	
		gunzip ${tissue}_tmp_forward.coverage-output.per-base.bed.gz \
			${tissue}_tmp_reverse.coverage-output.per-base.bed.gz
  		bedtools intersect -a ${tissue}_tmp_forward.exons.bed -b \
  			${tissue}_tmp_forward.coverage-output.per-base.bed -wao>\
  			${tissue}_tmp_forward.coverage.results
  		bedtools intersect -a ${tissue}_tmp_reverse.exons.bed -b \
  			${tissue}_tmp_reverse.coverage-output.per-base.bed -wao>\
  			${tissue}_tmp_reverse.coverage.results
		cat ${tissue}_tmp_forward.coverage.results ${tissue}_tmp_reverse.coverage.results\
			>${tissue}_tmp_coverage.results
		awk '$11>=1 {print $4}' ${tissue}_tmp_coverage.results >${tissue}_tmp_a1
		sort ${tissue}_tmp_a1 | uniq >${tissue}_tmp_more_than1
		awk '$11<1 {print $4}' ${tissue}_tmp_coverage.results >${tissue}_tmp_b1
		sort ${tissue}_tmp_b1 | uniq >${tissue}_tmp_less_than1
		cat ${tissue}_tmp_more_than1 | parallel --pipe -N 1000 grep -v -Fw -f \
			${tissue}_tmp_less_than1 >${tissue}_valid_coverage_transcripts_1
		awk '$11>=3 {print $4}' ${tissue}_tmp_coverage.results >${tissue}_tmp_a3
		sort ${tissue}_tmp_a3 | uniq >${tissue}_tmp_more_than3
		awk '$11<3 {print $4}' ${tissue}_tmp_coverage.results >${tissue}_tmp_b3
		sort ${tissue}_tmp_b3 | uniq >${tissue}_tmp_less_than3
		cat ${tissue}_tmp_more_than3 | parallel --pipe -N 1000 grep -v -Fw -f \
			${tissue}_tmp_less_than3 >${tissue}_valid_coverage_transcripts_3
		#filter genomic DNA reads
		module purge
		module load parallel
		module load bedtools2
		module load emboss
		parallel -a ${tissue}_tmp_transcripts --pipepart ./get_valid3end.sh \
			${tissue}.collapsed.gff
		sort ${tissue}_valid3end-transcripts | uniq >${tissue}_valid_3end_transcripts 
		# filter retained intron transcripts and non-canonical splice-junctions		   		
  		module purge
  		module load py-pandas/0.21.1-py2-326uzkn
  		cat ${tissue}.collapsed.gff | grep -Fw -f ${tissue}_valid_splice_transcripts \
  			>${tissue}_tmp_input.gff
  		time python get_premature_mRNAs.py ${tissue}_tmp_input.gff
  		#filter transcripts
		module load parallel
  		awk '$3=="exon" {print $12}' ${tissue}.collapsed.gff | sed 's/\"//g; s/\;//g' \
  			|  sort | uniq -c | awk '$1>1 {print $2}' >${tissue}_tmp_spliced_transcripts
  		awk '$3=="exon" {print $12}' ${tissue}.collapsed.gff | sed 's/\"//g; s/\;//g' \
  			|  sort | uniq -c | awk '$1<2 {print $2}' >${tissue}_tmp_un_spliced_transcripts
  		cat ${tissue}_valid_coverage_transcripts_3 | parallel --pipe -N1000 grep -Fw -f \
  			${tissue}_valid_splice_transcripts >${tissue}_tmp_valid-coverage_splice_transcripts
  		cat ${tissue}_tmp_valid-coverage_splice_transcripts | grep -v -Fw -f \
  			${tissue}_valid_splice_premature_mRNAs>${tissue}_tmp_valid_spliced_transcripts
  		cat ${tissue}_valid_coverage_transcripts_3 | parallel --pipe -N1000 grep -Fw -f \
  			${tissue}_tmp_un_spliced_transcripts >\
  		${tissue}_tmp_valid-coverage_unsplice_transcripts
  		cat ${tissue}_valid_3end_transcripts  | parallel --pipe -N1000 grep -Fw -f \
  			${tissue}_tmp_valid-coverage_unsplice_transcripts \
  			>${tissue}_tmp_valid_unspliced_transcripts
  		cat ${tissue}_tmp_valid_spliced_transcripts \
  			${tissue}_tmp_valid_unspliced_transcripts >${tissue}_tmp_valid_transcripts
  		cat ${tissue}.collapsed.group.txt | grep -Fwf ${tissue}_tmp_valid_transcripts \
  			| awk '{print $2}' | tr '\,' '\n' >${tissue}_tmp_valid_reads
		#re-clustering
 		module purge
  		module load samtools
  		module load picard
  		module load tofu2
  		module load seqtk
		samtools view -Sb -@ 36 ${transcriptomes}/${tissue}_trinity/Trinity.sam \
			>${transcriptomes}/${tissue}_trinity/Trinity.bam 
  		picard FilterSamReads I=${transcriptomes}/${tissue}_trinity/Trinity.bam \
  			O=${transcriptomes}/${tissue}_trinity/Trinity_valid_reads.bam \
  			READ_LIST_FILE=${tissue}_tmp_valid_reads FILTER=includeReadList 
  		samtools view -h -o ${transcriptomes}/${tissue}_trinity/Trinity_valid_reads.sam \
  			${transcriptomes}/${tissue}_trinity/Trinity_valid_reads.bam 
  		sort -k 3,3 -k 4,4n ${transcriptomes}/${tissue}_trinity/Trinity_valid_reads.sam \
  			>${transcriptomes}/${tissue}_trinity/Trinity_valid_reads_sorted.sam
  		seqtk subseq ${transcriptomes}/${tissue}_trinity/${tissue}-formatted-Trinity.fasta \
  			${tissue}_tmp_valid_reads > ${tissue}_tmp_valid_reads.fa
  		collapse_isoforms_by_sam.py --max_fuzzy_junction 15 --min-coverage 0.90 \
  			--min-identity 0.95 --input ${tissue}_tmp_valid_reads.fa -s \
  			${transcriptomes}/${tissue}_trinity/Trinity_valid_reads_sorted.sam \
  			-o ${tissue}_final
		ls | grep tmp | grep ${tissue} | xargs -I "{}" rm -r {}  
done;

```

* **WTTS-seq Data Analysis**
``` 
# trim adapters and T-rich stretches located at the 5’ end of each read
# trim_T-rich_regions_from_reads_5end.pl were obtained from PMC4896187 study
module purge
module load perl
wtts_data="/usr/WTTS-seq"
cd ${wtts_data}
flag=1
for i in $(ls | grep fastq$);
	do
		id=$(echo ${i} | sed 's/\.fastq//g')
		perl trim_T-rich_regions_from_reads_5end.pl ${i} ${id}.5end_trimmed.fq 12
		echo $flag "file is done"
		flag=$(expr $flag + 1)
done;

# error-correction, quality filter, and mapping
for i in $(ls | grep 5end_trimmed.fq$);
	do
		id=$(ech ${i} | sed 's/\.5end*$//g')
		#error-correction
		coral -fq ${i} -o ${id}.5end_trimmed_errorCorrected.fq \
			-p 36 -mr 2 -mm 2 -g 3 -e 0.07 
		#length filtering
		bioawk -cfastx 'length($seq) <=300 {print "@"$name"\n"$seq"\n+\n"$qual}' \
			${id}.5end_trimmed_errorCorrected.fq >${id}_5end_trimmed_length_filtered.fq
		#quality filtering
		fastq_quality_filter -v -q 20 -p 50 -i ${id}_5end_trimmed_length_filtered.fq \
			-o ${id}_5end_trimmed_length_QC_filtered.fq
		#mapping
		STAR --runThreadN 36 --genomeDir /usr/ARS-UCD1.2.RepeadMasked.STAR/ --readFilesIn \
			${id}_5end_trimmed_length_QC_filtered.fq  --outFileNamePrefix \
			${id}_5end_trimmed_length_QC_filtered_ --outFilterMatchNminOverLread 0.9 \
			--outFilterMismatchNoverLmax 0.05 --outSAMtype BAM SortedByCoordinate		
done;
```
