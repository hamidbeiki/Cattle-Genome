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

### WTTS-seq Data Analysis

* **Trim Adapters and T-rich Stretches Located at The 5’ End of Each Read**
``` 
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

```

* **Error-correction, Quality Filter, and Mapping**
```
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
* **Combine Tissue Replicates**
```
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

```

* **Get Bed Files**
```
module load bedtools2
for i in  $(ls | grep clean_WTTS_reads.bam$);
	do
		id=$(echo $i | sed 's/\.bam//g')
		bamToBed -i ${i} >${id}.bed
done;

```

### RAMPAGE Data Analysis

```
#RAMPAGE bed files were obtained from PMC8015843 study
rampage_dir='/usr/RAMPAGE'

```

### Combine Asembled Tissues RNA-seq-based Transcripts

```
module load miniconda3/4.3.30-qdauveb
cd ${rnaseq_dir}
source activate /usr/.conda/envs/py-libs
ls | grep gff$ | grep final> input_files
python combine-gtf_files.py input_file ---wtts_dir '/usr/WTTS-seq'\
 	---atac_dir '/usr/ATAC-seq'

# get transcript sequences
module load cufflinks
gffread -w rna_seq_transcriptome.fa -g \
	${genome_dir}/bt_ref_ARS-UCD1.2.RepeadMasked.fa rna_seq_transcriptome.gtf
```

### Quantify Asembled Tissues RNA-seq-based Transcripts/Genes

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

### Transcript Biotypes

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

```

### Check ncRNAs Homology

```
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

### Alternative-splicing Analysis

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

### ATAC-seq Data Analysis

* **Trim Adapters and Low-quality Bases**
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

```

* **Alignment**
```
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

```

* **Filter Duplicated Reads**
```
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
```

* **Peak Calling**
```
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

### ChIP-seq Data Analysis

```
#### ChIP-seq peaks (bed files) were received from PMC7988148 study
chip_dir='/usr/ChIP-seq'

```

### PacBio Iso-seq Data Analysis

* **Convert Bax to Bam**
```
### NOTE: the matched RNA-seq data were processed as described in "RNA-seq pre-processing",
### "RNA-seq mapp to genome" and "Assemble RNA-seq-based Transcripts" sections

module load smrtlink/7.0.0  		
module load samtools
module load bamtools
pacbio_data='/usr/pacbio'
cd ${pacbio_data}

for i in $(cat libraries);
	do 
		bax2bam -o $lib $lib.1.bax.h5 $lib.2.bax.h5 $lib.3.bax.h5 
done;  		
```

* **Get CCS Reads**
```
for i in $(cat libraries);
	do
		ccs --noPolish --minLength=300 --minPasses=1 --minZScore=-999 --maxDropFraction=0.8\
	 		--minPredictedAccuracy=0.8 --minSnr=4 --reportFile ${lib}_report.txt \
	 		${lib}.subreads.bam ${lib}.ccs.bam
done;

for i in $(ls | grep ccs.bam$)
	do
   		id=$(echo $i | rev | cut -c 5- | rev)
   		dataset create --type ConsensusReadSet ${id}.XML ${i}
done;

```

* **Get Full length (FL) Reads**
```
for i in $(ls | grep ccs.bam$)
	do
		id=$(echo $i | rev | cut -c 9- | rev)
		lima ${id}.ccs.bam primers.fasta ${id}.ccs.primerCleaned.bam --isoseq --no-pbi
		isoseq3 refine ${id}.ccs.primerCleaned.primer_5p--primer_3p.bam primers.fasta \
			${id}.refinedFl.bam
		isoseq3 cluster ${id}.refinedFl.bam ${id}.unpolishedFl.bam --verbose
done;

```

* **Get Tissue's FL Reads**
```
for i in $(awk '{print $3}' tissue-lib-info |  sort | uniq);
	do 
		awk -v tissue="${i}" '$3==tissue {print $1}' tissue-lib-info>tmp
     		ls | grep refinedFl.bam$ | grep -f tmp>tmp_list
     		bamtools merge -list tmp_list -out ${i}.refinedFl.bam
     		isoseq3 cluster ${i}.refinedFl.bam ${i}.unpolishedFl.bam --verbose
done;

```

* **Get Full Length Non-Chimeric (FLNC) reads**
```
for i in $(ls | grep ccs.bam$)
	do
  		id=$(echo $i | rev | cut -c 9- | rev)
  		samtools bam2fq ${i} | seqtk seq -A > ${id}.ccs.fasta
  		pbtranscript classify ${id}.ccs.fasta ${id}_draft.fasta --primer \
  			primers_formatted.fasta --cpus 70 --flnc=${id}.flnc.fasta \
  			--nfl=${id}.nfl.fasta --outDir ${id}_tmp
  		rm -r ${id}_tmp
done;

#get tissue's FLNC reads
for i in $(awk '{print $3}' tissue-lib-info |  sort | uniq);
 do
	echo "working on $i"
   	awk -v tissue="${i}" '$3==tissue {print $1}' tissue-lib-info | \
   		awk '{print $1"_s1_p0.flnc.fasta"}'>tmp
   	{ xargs cat < tmp; } >${i}.flnc.fasta
   	rm -r tmp
   	echo "$i is done"
done;

```

* **Genome-Guided Error-Correction with FMLRC**
```
module purge
module load fmlrc/1.0.0
for tissue in $(cat pacbio_tissue_list):
	do 
		fmlrc -p 70 ${pacbio_data}/${rna_seq}/${tissue}-norm/comp_msbwt.npy \
			${tissue}.flnc.fasta ${tissue}_fmlrc_corrected_flnc.fa
done;

```

* **denovo Error-Correction with Proovread**
```
module load miniconda
source activate proovreadenv
for tissue in $(cat pacbio_tissue_list):
	do
		SeqChunker -s 20M -o ${tissue}-%03d.fa ${tissue}.flnc.fasta  
		proovread  --long-reads=${tissue}-001.fa \
			-s ${pacbio_data}/${rna_seq}/${tissue}-norm/${tissue}.normalized_reads.fq \
			-p ${tissue}_proovread_corrected_flnc   -t 70 --no-sampling  --coverage=20 --ignore-sr-length
done;		

```

* **Combine Proovread, and FMLRC Error-corrected reads***
```
cat ${tissue}_proovread_corrected_flnc.fa ${tissue}_fmlrc_corrected_flnc.fa\
	> ${tissue}_fmlrc-proovread_corrected_flnc.fa

```

* **Mapp Error-corrected PacBio Reads to Bovine Genome**
```
module purge
module load gmap-gsnap-legacy
for tissue in $(cat pacbio_tissue_list)
	do
		tissue=$(echo ${i} | sed 's/\_trinity//g')
		gmap -D /usr/ARS-UCD1.2.RepeadMasked.GMAP -d GMAP_bostaurus \
			-f samse --min-trimmed-coverage 0.90 --min-identity 0.95 \
			${tissue}_fmlrc-proovread_corrected_flnc.fa \
			>${tissue}_fmlrc-proovread_corrected_flnc.sam 2 \
			>${tissue}_fmlrc-proovread_corrected_flnc.sam.log
		sort -k 3,3 -k 4,4n ${tissue}_fmlrc-proovread_corrected_flnc.sam \
			>${tissue}_fmlrc-proovread_corrected_flnc_sorted.sam
done;

```

* **Collapse & Group Reads Into Putative Gene Models**
```
module purge
module load smrtlink
module load cufflinks
module load bedtools2
genome_dir="/usr/ARS-UCD1.2.RepeadMasked.GMAP"
for tissue in $(cat pacbio_tissue_list)
	do
		collapse_isoforms_by_sam.py --max_fuzzy_junction 0 \
			--min-coverage 0.90 --min-identity 0.95 --input \
			${tissue}_fmlrc-proovread_corrected_flnc.fa \
			-s ${tissue}_fmlrc-proovread_corrected_flnc_sorted.sam \
			-o ${tissue}
		sed -i 's/chr//g' ${tissue}.collapsed.gff
		gffread -w ${tissue}_transcripts.fa -g ${genome_dir}/bt_ref_ARS-UCD1.2.RepeadMasked.fa \
			${tissue}.collapsed.gff
done;

### NOTE: The following steps were performed similar to what described in\
### "Assemble RNA-seq-based Transcripts" section to "get ${tissue}_final.gff" files:\
### "check splice-junction validity"; "check transcript coverage";\
### "filter retained intron transcripts and non-canonical splice-junctions";
### "filter genomic DNA reads"

```

* **Combine Tissue Transcripts**
```
module load miniconda3/4.3.30-qdauveb
source activate /home/beiki/.conda/envs/py-libs
ls | grep gff$ | grep final> input_files
cd ${pacbio_data}
python combine-gtf_files.py input_files

```

### Oxford Nanopore Data Analysis

```
### transcript coordinate file (gtf files) were received from PMC8173071 study
nanopore_dir='/usr/nanopore'

```

### Transcript Structure Validation

```
module load miniconda3/4.3.30-qdauveb
source activate /home/beiki/.conda/envs/py-libs
cd ${rnaseq_dir}
cat input_files
	/usr/RNA-seq/rna_seq_transcriptome.gtf	
	/usr/pacbio/pacbio_transcriptome.gtf
	/usr/nanopore_dir/nanopore_transcriptome.gtf	
	/usr/Ensembl_Bos_taurus.ARS-UCD1.2.gtf
	/usr/NCBI_Bos_taurus.ARS-UCD1.2.gtf		
python compare-gtf_files.py input_files ---wtts_dir '/usr/WTTS-seq'\
 	---atac_dir '/usr/ATAC-seq'

# get known gene's border extension
module load python
python get_gene_extensions.py

# get the effect of gene extension on expression
python relate_gene_extension_to_expression.py

```

### Transcript Support With Epigenetic Data

```
module load bedtools2
cd ${rnaseq_dir}

# get transcript bed file   		
awk '$3=="transcript"{print $1,$4,$5,$12,".",$7,$12}' rna_seq_transcriptome.gtf \
	sed 's/\"//g;s/\;//g' | tr ' ' '\t' > transcripts.bed

# get promoters bed file
awk '$6=="+"' transcripts.bed >plus
awk '$6=="-"' transcripts.bed >revs
awk 'BEGIN { OFS = "\t" } { $8 = $2 - 500, $9 = $2 + 300} 1' plus \
	| awk '{print $1,$8,$9,$7,".",$6,$7}' | tr ' ' '\t' > a
awk 'BEGIN { OFS = "\t" } { $8 = $3 + 500, $9 = $3 - 300} 1' revs \
	| awk '{print $1,$9,$8,$7,".",$6,$7}' | tr ' ' '\t' > b
cat a b >transcripts_promoters.bed  		

# get transcript supported by ATAC-seq
for i in $(ls ${atacseq_data} | grep macs2.OUTPUT$)
	do
		tissue=$(echo ${i} | se 's/\macs2.*$//g')
		bedtools intersect -s -a transcripts_promoters.bed \
			-b ${atacseq_data}/${tissue}_macs2.OUTPUT/peaks.bed | awk {print $4} | sort | uniq\
			>>atac_supported_transcripts
done;

```

### Transcript Border Validation

```
module load bedtools2

# get transcript's 5' offset bed file
awk '$6=="+"' transcripts.bed >plus
awk '$6=="-"' transcripts.bed >revs
awk 'BEGIN { OFS = "\t" } { $8 = $2 - 30, $9 = $2 + 10} 1' plus \
	| awk '{print $1,$8,$9,$7,".",$6,$7}' | tr ' ' '\t' > a
awk 'BEGIN { OFS = "\t" } { $8 = $3 + 30, $9 = $3 - 10} 1' revs \
	| awk '{print $1,$9,$8,$7,".",$6,$7}' | tr ' ' '\t' > b
cat a b >TSS_offset.bed   
	
# get transcript supported by RAMPAGE		
for i in $(ls ${rampage_dir} | grep bed$)
	do
		bedtools intersect -s -a TSS_offset.bed \
			-b ${rampage_dir}/${i} | awk {print $4} | sort | uniq \
			>>rampage_supported_transcripts
done;

# get transcript's 3' offset bed file
awk '$6=="+"' transcripts.bed >plus
awk '$6=="-"' transcripts.bed >revs
awk 'BEGIN { OFS = "\t" } { $8 = $3 - 10, $9 = $3 + 165} 1' plus \
	| awk '{print $1,$8,$9,$7,".",$6,$7}' | tr ' ' '\t' > a
awk 'BEGIN { OFS = "\t" } { $8 = $2 - 165, $9 = $2 + 10} 1' revs \
	| awk '{print $1,$9,$8,$7,".",$6,$7}' | tr ' ' '\t' > b
cat a b >TTS_offset.bed   
	
### gets transcript supported by RAMPAGE		
for i in $(ls ${wtts_dir} | grep bed$)
	do
		bedtools intersect -s -a TTS_offset.bed \
			-b ${wtts_dir}/${i} | awk {print $4} | sort | uniq \
			>>wtts_supported_transcripts
done;

```

### miRNA Data Analysis

* **Trim Adapters and Low-quality Bases**
```
mirna_dir="/usr/miRNA-seq"
cd ${mirna_dir}
module load cutadapt
module load fastqc
for i in $(ls | grep fq.gz$)
	do
    		trim_galore --gzip --cores 70 --quality 20 --fastqc --length 16 \
    			--max_length 30 -a AACTGTAGGCACCATCAAT ${i}
done;

```

* **Mapp Trimmed Readss**
```
module load mirdeep2

for i in $(cat trimmed_file_list)
  do
    	id="$(echo $i | rev | cut -c 15- | rev)"
    	pigz -p 70 -c -d ${i}>${id}_trimmed.fq
    	mapper.pl  ${id}_trimmed.fq -e -h -q -j -l 16 -o 70 -r 1 -m \
    		-p /usr/ARS-UCD1.2.RepeadMasked.bowtie2/ARS-UCD1.2 -s reads_collapsed.${id}.fa \
    		-t reads_vs_genome.${id}.arf -v -n |&tee>${id}-info
    grep "#desc" -A2 ${id}-info >${id}-reads_mapping_info
    mv bowtie.log ${id}-UNIQUE_reads_mapping_info
    rm -r ${id}_trimmed.fq ${id}-info
done;

```

* **Quantify miRNAs**
```
* get known miRNA sequences fro, mirBase
wget https://www.mirbase.org/ftp/CURRENT/mature.fa.gz
wget https://www.mirbase.org/ftp/CURRENT/hairpin.fa.gz
gunzip -c mature.fa.gz | grep "Bos taurus" -A1 >bos_taurus-mature.fa
gunzip -c hairpin.fa.gz | grep "Bos taurus" -A1 >bos_taurus-hairpin.fa

#quantification
module load mirdeep2
for i in $(ls | grep reads_collapsed | grep .fa$)
    do
     	id=$(echo $i | sed 's/reads\_collapsed\.//g; s/\.fa//g')
	 miRDeep2.pl reads_collapsed.${id}.fa \
	 	${genome}/bt_ref_ARS-UCD1.2.RepeadMasked_miRNA_quantification.fa \
	 	reads_vs_genome.${id}.arf ${supplementary}/bos_taurus-mature.fa \
	 	bos_taurus-mature.fa bos_taurus-hairpin.fa -t bta -c -v 2>report_filter.${id}.log
done;

# summarize miRNA expression results
module load python
python get_miRNAs_summary.py

```

### Traits similarity Network

```
module load python
python trait_similarity_network.py

```

### Additional Analysis
```
module load python
python compare_transcript_biotype_expression.py
python compare_gene_biotype_expression.py
python get_number_of_UTRs_per_gene.py
python get_relation_betwee_protein_coding_genes_and_aberrant_transcripts.py
python get_gene_biotype_changes_across_tissues.py
python get_tissue_specificity_info.py

```

### The constructed bovine trait similarity network is publicly available through the Animal Genome database:
https://www.animalgenome.org/host/reecylab/a

### The constructed transcriptome and related sequences are publicly available in the Open Science Framework database:
https://osf.io/jze72/?view_only=d2dd1badf37e4bafae1e12731a0cc40d

### 


#### Citation:


#### contact:
beiki.h.m@gmail.com




	

