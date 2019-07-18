# Running jobs on Spider

1. [Job preparation](#prepare-jobs)

### <a name="spider-jobs"></a> 1. Job preparation

We will download a set of trimmed FASTQ files to work with. These are small subsets of our real trimmed data we prepared earlier, and will enable us to run our variant calling workflow quite quickly.

cd /project/spidercourse/Data/ecoli-analysis/data
mkdir trimmed_fastq_small
cd trimmed_fastq_small
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz

Our variant calling workflow has the following steps:

1. Index the reference genome for use by bwa and samtools  
2. Align reads to reference genome  
3. Convert the format of the alignment to sorted BAM, with some intermediate steps  
4. Calculate the read coverage of positions in the genome  
5. Detect the single nucleotide polymorphisms (SNPs)  
6. Filter and report the SNP variants in VCF (variant calling format) 


 ```sh
 ls /project 
 ```
> **_Food for brain:_**
>
> * Do you know what project you belong to? What all access does it provide to you?
> * What are each of ther project directories for? What is public/private? Do you have read write permissions on all spaces?




