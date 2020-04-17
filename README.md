# Apodemus_research
Analysing the carryover sequences of Apodemus mice sequences obtained from RAD-seq data                                         
Programmes needed: Seqtk, blast+, bwa-mem, parallel (can all be installed using conda)

# Obtaining the unmapped reads 

Burrow Wheeler Alignmer (BWA) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to a few megabases. BWA-MEM and BWA-SW share similar features such as the support of long reads and chimeric alignment, but BWA-MEM, which is the latest, is generally recommended as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

The latest version (0.7.17) can be installed through conda: https://anaconda.org/bioconda/bwa
```
conda install -c bioconda bwa 
```
**Documentation**

Detailed documentation can be found at http://bio-bwa.sourceforge.net/bwa.shtml

# BWA-mem analysis 

We mapped the sequences using BWA-mem against the draft Apodemus sylvaticus from ncbi (university of Liverpool, 498341) 
https://www.ncbi.nlm.nih.gov/assembly/GCA_001305905.1/

in order to obtain all the unmapped reads                                                                                               
1a: index the reference sequence, obtained from ncbi 
```
bwa index GCA_001305905.1_ASM130590v1_genomic.fna.gz
```

1b: Map all sequences against the reference genome and convert the files from sam to bam
```
cd ~
cd Research 
cd 80_sequences_fq
for a in *.fq.gz
do
bwa mem GCA_001305905.1_ASM130590v1_genomic.fna.gz $a | samtools sort -o /Volumes/HaeyamData/Mapping/output/$a.bam -
done
```

1c: filter out all the unmapped reads
```
cd /Volumes/HaeyamData/Mapping/output
for a in *.fq.gz.bam
do 
samtools view -b -f 4 $a > unmapped_bam/$a.unmapped.bam
done
```
1d: convert all files to fasta 
```
cd /Volumes/HaeyamData/Mapping/output/unmapped_bam
for a in *.fq.gz.bam.unmapped.bam;
do
samtools fasta $a > unmapped_fa/$a.fa
done
```

**Blast the unmapped reads to the whole ncbi nt_v5 database**                                                                    
2a: blast all the unmapped reads on the cluster, need to copy all the input files and the nt_v5 database on the cluster 
```
qsub orion_job.txt
"ls *.txt |parallel -j 8 'blastn -task megablast -db ../blastdb_v5/nt_v5 -query {} -dust no -max_target_seqs 1 -perc_identity 75 qcov_hsp_perc 50 -outfmt "6 qseqid sseqid evalue pident stitle" -out output/{.}.txt'"
```
2b: copy all outputs onto desktop

