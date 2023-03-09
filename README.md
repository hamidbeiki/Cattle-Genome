# Bovine Transcriptome Assembly and Annotation

## Bioinformatics methods and Analysis Steps

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


```
