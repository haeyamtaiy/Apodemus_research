# Apodemus_research
Analysing the carryover sequences of Apodemus mice sequences obtained from RAD-seq data                                         
Programmes needed: bwa-mem, samtools, blast+, parallel (can all be installed using conda)

**Overview**
Step |Procedure | outputs
------------ | ------------ | -------------
1 | ddRAD-seq | A.sylvaticus and A.flavicollis FASTQ files
2 | BWA-mem | Unmapped reads 
3 | BLAST+ 2.2.31+  | Blasted results


# Obtaining the unmapped reads 

Burrow Wheeler Alignmer (BWA) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to a few megabases. BWA-MEM and BWA-SW share similar features such as the support of long reads and chimeric alignment, but BWA-MEM, which is the latest, is generally recommended as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

The latest version (0.7.17) can be installed through [conda] https://anaconda.org/bioconda/bwa
```
conda install -c bioconda bwa 
```
**Documentation**

Detailed documentation can be found at http://bio-bwa.sourceforge.net/bwa.shtml

# BWA-mem analysis 

We mapped the sequences using BWA-mem against the draft *Apodemus sylvaticus* from ncbi (university of Liverpool, 498341) 
https://www.ncbi.nlm.nih.gov/assembly/GCA_001305905.1/
                                                                                            
**1a:** Before alignment, indexing is necessary therefore index the reference sequence ie the *Apodemus Sylvaticus* sequence downloaded from ncbi
```
bwa index GCA_001305905.1_ASM130590v1_genomic.fna.gz
```
This step will take a while and by the end there will be 5 output files which will be saved where the original reference sequence file is stored: 
1. <reference_file_name>.fa.amb -> text file
2. <reference_file_name>.fa.ann -> text file
3. <reference_file_name>.fa.bwt -> binary file
4. <reference_file_name>.fa.pac -> binary file
5. <reference_file_name>.fa.sa -> binary file

**1b:** Map all sequences against the reference genome and convert the files from sam to bam
In order to convert the mapped output sam files (The Sequence Alignment/Map) to bam files (binary file version of the sam files) need to instal 'samtools'

**Samtools** 

Samtools is a set of utilities that manipulate alignments in the BAM format. It imports from and exports to the SAM (Sequence Alignment/Map) format, does sorting, merging and indexing, and allows to retrieve reads in any regions swiftly.

The latest version (1.10) can be installed through [conda] https://anaconda.org/bioconda/samtools
```
conda install -c bioconda samtools
```

**Documentation**

Detailed documentation can be found at http://www.htslib.org/doc/samtools.html

**Loop**

The loop created combines both the alignment step and the conversion of the output sam files to bam files.   
This is carried out in the folder which contains the input files to be alignmed as well as the index reference files and the output is directed to an output folder created
```
for a in *.fq.gz
do
bwa mem GCA_001305905.1_ASM130590v1_genomic.fna.gz $a | samtools sort -o /Volumes/HaeyamData/Mapping/output/$a.bam -
done
```

**1c:** In our case we are only interested in the unmapped reads therefore we need to filter out all the mapped reads so we are left with only the unmapped. 
The output files have the prefix "unmapped.bam" and are directed into another folder called <unmapped_bam> 
```
cd /Volumes/HaeyamData/Mapping/output
for a in *.fq.gz.bam
do 
samtools view -b -f 4 $a > unmapped_bam/$a.unmapped.bam
done
```
**1d:** These files are still in the bam format. This is a binary file (non-human readable) therefore, need to convert all the files to FASTA format. (text-based format for representing either nucleotide sequences or amino acid (protein) sequences). This step will also require samtools. 
```
cd /Volumes/HaeyamData/Mapping/output/unmapped_bam
for a in *.fq.gz.bam.unmapped.bam;
do
samtools fasta $a > unmapped_fa/$a.fa
done
```
# BLAST+

The NCBI provides a suite of command-line tools to run BLAST called BLAST+. This allows users to perform BLAST searches on their own server without size, volume and database restrictions. BLAST+ can be used with a command line so it can be integrated directly into your workflow.

[Documentation and installation] https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

**Blast the unmapped reads to the whole ncbi nt_v5 database**   

nt_v5 database: This is a nucleotide database which can be obtained through the [ncbi website] https://www.ncbi.nlm.nih.gov/books/NBK537770/

**2a:** blast all the unmapped reads on a high performance computer (hpc) as it requires high memory storage and RAM usage. 

First need to copy all the input files and the nt_v5 database onto the hpc cluster
```
qsub orion_job.txt 
"ls *.txt |parallel -j 8 'blastn -task megablast -db ../blastdb_v5/nt_v5 -query {} -dust no -max_target_seqs 1 -perc_identity 75 qcov_hsp_perc 50 -outfmt "6 qseqid sseqid evalue pident stitle" -out output/{.}.txt'"
```
This script contains specific paramaters: 

**2b:** copy all outputs onto desktop using scp -r

