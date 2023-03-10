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

# combine tissue replicates
module purge
module load picard
flag=1
for tissue in $(cat tissues);
	do		
		grep -Fw ${tissue} file_annotation | awk '{print $2}' | \
			sed 's/fastq/5end\_trimmed\_errorCorrected\_filtered\.bam/g' >tmp
		cat tmp | tr '\n' ' ' | sed 's/Ion/\I=Ion/g' | xargs -I "{}" \
			sed "s/INPUTFILES/{}/g; s/TISSUE/${tissue}/g" picard_command.sh >tmp.sh\
		bash tmp.sh
		flag=$(expr $flag + 1 )
done;

# get bed files
module load bedtools2
for i in  $(ls | grep clean_WTTS_reads.bam$);
	do
		id=$(echo $i | sed 's/\.bam//g')
		bamToBed -i ${i} >${id}.bed
done;

```

* **RAMPAGE Data Analysis**
```
# RAMPAGE bed files were obtained from PMC8015843 study
rampage_dir='/usr/RAMPAGE'

```
* **Combine Tissue Transcripts**
```
module load miniconda3/4.3.30-qdauveb
cd ${rnaseq_dir}
source activate /usr/.conda/envs/py-libs
ls | grep gff$ | grep final> input_files
python combine-gtf_files.py input_files --combine ---wtts_dir '/usr/WTTS-seq'\
 	---atac_dir '/usr/ATAC-seq'

# get transcript sequences
module load cufflinks
gffread -w rna_seq_transcriptome.fa -g \
	${genome_dir}/bt_ref_ARS-UCD1.2.RepeadMasked.fa rna_seq_transcriptome.gtf
```

* **Quantify Tissue Transcripts/Genes**
```
module load trinity
cd ${rnaseq_dir}
for i in $(ls ${rnaseq_dir}/trimmed/*.fastq)
	do
		tissue=$(echo ${i} | sed 's/\.fastq/g')
		echo tissue>>tissues
done;

for tissue in $(cat tissues)
	do
		align_and_estimate_abundance.pl --transcripts rna_seq_transcriptome.fa \
			--seqType fq --single ${rnaseq_dir}/trimmed/${tissue}.fastq \
			--aln_method bowtie --est_method RSEM --SS_lib_type RF --thread_count 36 \
			--gene_trans_map gene_trans_map --thread_count 36 \
			--output_dir ${tissue}-RSEM_quantification
done;

```

* **Transcript Biotypes**
```
# get transcript's ORFs
wget https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz
gunzip ORFfinder.gz
./ORFfinder -in rna_seq_transcriptome.fa -ml 44 -s 1 -strand plus \
	-out rna_seq_transcriptome.orfs -outfmt 0

# select 3 longest ORFs
samtools faidx rna_seq_transcriptome.orfs
cut -f1-2  rna_seq_transcriptome.orfs.fai >orf.length
awk 'BEGIN { FS = "_" } ; { print $3}' orf.length >b
awk 'BEGIN { FS = "_" } ; { print $1, $2}' orf.length | sed 's/ /\_/g' >a
paste a b>orf.length2
# R script
>data=read.table(file="orf.length2",header=T)
>library(dplyr)
> res=data %>%
  group_by(isoform) %>%
  top_n(n = 3, wt = orf_length)
> write.table(res, file = "rna_seq_transcriptome.selected_orfs", col.names = T,\
  row.names = F,quote=FALSE,sep="\t")

# blastp against vertebrates peptide sequences from Uniprot database
module load ncbi-toolkit
export BLASTDB=/usr/unire90_vertebrates_db
blastp -query rna_seq_transcriptome.selected_orfs -db unire90_vertebrates_db \
	-out rna_seq_transcriptome.selected_orfs.Blastp.txt -evalue 1e-6 -outfmt \
	"7 qseqid qlen qstart qend sseqid slen sstart send identity pident length qcovs mismatch evalue" \
	-max_target_seqs 1 -num_threads 40

# filter for identity & coverage
awk '$9>=95 && $12>=60' rna_seq_transcriptome.selected_orfs.Blastp.txt \
	>rna_seq_transcriptome.selected_orfs.Blastp_filtered.txt

# gene-level homology
pyhton get_longest_transcript_per_gene.py rna_seq_transcriptome.gtf rna_seq_transcriptome.fa
#get_longest_transcript_per_gene.py is available https://github.com/hamidbeiki/Cattle-Genome
export BLASTDB=/usr/unire90_vertebrates_db
blastx -query genes_longest_transcripts.fa -db unire90_vertebrates_db \
	-out gene-level_homology.Blastx.txt -evalue 1e-6 -outfmt \
	"7 qseqid qlen qstart qend sseqid slen sstart send identity pident length qcovs mismatch evalue" \
	-max_target_seqs 1 -num_threads 40
awk '$9>=95 && $12>=95' gene-level_homology.Blastx.txt \
	>gene-level_homology.Blastx_filtered.txt
# get putative ncRNAs sequences
module load cufflinks
grep '>' rna_seq_transcriptome.selected_orfs | sed 's/>//g'>ids
grep -Fwfv ids rna_seq_transcriptome.gtf | awk '$3=="transcript"' >putative_ncrnas.gtf
gffread -w putative_ncrnas.fa -g \
	${genome_dir}/bt_ref_ARS-UCD1.2.RepeadMasked.fa putative_ncrnas.gtf
# check the coding-potential of putative ncRNAs
module load cpc2
CPC2.py -i putative_ncrnas.fa -o putative_ncrnas_cpc2_results

# get transcripts/genes biotypes
module load pyhton
python get_transcript-gene_biotypes.py rna_seq_transcriptome.gtf \
	   rna_seq_transcriptome.selected_orfs.Blastp.txt putative_ncrnas_cpc2_results
# check homology of ncRNAs with known-ncRNAs
module purge
module load cufflinks
module load ncbi-toolkit
grep ncRNA rna_seq_transcriptome_transcript_biotypes | awk '{print $1}' >ids
grep -Fwf ids rna_seq_transcriptome.gtf >ncRNAs.gtf
gffread -w ncRNAs.fa -g \
	${genome_dir}/bt_ref_ARS-UCD1.2.RepeadMasked.fa ncRNAs.gtf
export BLASTDB=/usr/ncbi_ensemble_release_102_ncrnas_db
blastn -query ncRNAs.fa -db ncbi_ensemble_release_102_ncrnas_db \
	-out ncRNAs.Blastn.txt -evalue 1e-6 -outfmt \
	"7 qseqid qlen qstart qend sseqid slen sstart send identity pident length qcovs mismatch evalue" \
	-max_target_seqs 100 -num_threads 40

```

* **Alternative-splicing Analysis**
```
AS_dir='/usr/RNA-seq/AS'
cd ${AS_dir}
module load suppa
suppa.py generateEvents -i rna_seq_transcriptome.gtf -o SUPPA_output -e FL -f ioe
suppa.py generateEvents -i rna_seq_transcriptome.gtf -o SUPPA_output -e SE -f ioe
suppa.py generateEvents -i rna_seq_transcriptome.gtf -o SUPPA_output -e SS -f ioe
suppa.py generateEvents -i rna_seq_transcriptome.gtf -o SUPPA_output -e MX -f ioe
suppa.py generateEvents -i rna_seq_transcriptome.gtf -o SUPPA_output -e RI -f ioe

# get unique-splice exons
module load python
python get_unique_splice_exon_AS_events.py

```

* **ATAC-seq Data Analysis**
```
atacseq_data="/usr/ATAC-seq"  
cd ${atacseq_data}
# trim reads
module load trimgalore
for i in $(ls *_1.fastq)
	do
		id=$(echo ${i} | sed 's/\_1\.fastq//g')
		trim_galore --paired ${id}_1.fastq ${id}_2.fastq;
done;

# alignment
module load bowtie
module load samtools
for i in $(ls 1_val_1)
	do
		id=$(echo ${i} | sed 's/\_1\_val\_1//g')
		treat=$(echo ${i} | sed 's/\_\_.*$//g')
		tissue=$(echo ${i} | sed 's/^.*\_\_$//g' | sed 's/\_1\_val\_1//g')
		bowtie -x /usr/ARS-UCD1.2.RepeadMasked.bowtie/ARS --fr -1 \
		${treat}__${tissue}_1_val_1.fastq -2 ${treat}__${tissue}_2_val_2.fastq \
		-S ${tissue}_${treat}.sam
		samtools view -bSq 30 ${tissue}_${treat}.sam >${tissue}_${treat}_unique.sam 
		samtools sort ${tissue}_${treat}_unique.sam>${tissue}_${treat}_unique_sorted.sam	
done;

# filter duplicated reads
module purge
module load picard
module load samtools
for i in $(ls unique_sorted.sam)
	do
		id=$(echo ${i} | sed 's/\_unique\_sorted.sam//g')
		tissue=$(echo ${id} | sed 's/\_.*$//g')
		treat=$(echo ${id} | sed 's/^.*\_//g')
		picard MarkDuplicates REMOVE_DUPLICATES=true \
		METRICS_FILE=${tissue}_Duplicates.txt INPUT=${i} \
		OUTPUT=${tissue}_${treat}_DuplicatesRemoved.sam
		samtools view -bS ${tissue}_${treat}_DuplicatesRemoved.sam \
		>${tissue}_${treat}_DuplicatesRemoved.bam	
done;

# peak calling
module load macs2
for i in $(ls DuplicatesRemoved.bam$ | grep -v input)
	do
		id=$(echo ${i} | sed 's/\_DuplicatesRemoved\.bam//g')
		tissue=$(echo ${id} | sed 's/\_.*$//g')
		treat=$(echo ${id} | sed 's/^.*\_//g')		
		macs2 callpeak -t ${i} -c ${tissue}_input_DuplicatesRemoved.bam  --broad  \
		--gsize 2.80e9 --broad-cutoff 0.1 --outdir ${tissue}_macs2.OUTPUT -n ${tissue}  \
		-B --nomodel -q 0.0
done;

```



