# Running jobs on Spider

We will download a set of trimmed FASTQ files to work with. These are small subsets of our real trimmed data we prepared earlier, and will enable us to run our variant calling workflow quite quickly. Later on if you have more time, you can try using the full data.

```sh
#As data manager
cd /project/spidercourse/Data/ecoli-analysis/data

# As a regular user
cd $HOME/ecoli-analysis/data

mkdir trimmed_fastq_small
cd trimmed_fastq_small
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
mv sub/* .
```

Our variant calling workflow has the following steps:

1. Index the reference genome for use by bwa and samtools  
2. Align reads to reference genome  
3. Convert the format of the alignment to sorted BAM, with some intermediate steps  
4. Calculate the read coverage of positions in the genome  
5. Detect the single nucleotide polymorphisms (SNPs)  
6. Filter and report the SNP variants in VCF (variant calling format) 

```sh
#As data manager
cd /project/spidercourse/Data/ecoli-analysis/

# As a regular user
cd $HOME/ecoli-analysis/

mkdir results
wget https://raw.githubusercontent.com/sara-nl/2019-08-28-htdp-elixir/gh-pages/_episodes/scripts/job-submit-variant-calling.sh
```

Let us inspect the contents of the script that will run the job of variant calling

```sh
cat job-submit-variant-calling.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake
bash /project/spidercourse/Data/ecoli-analysis/run-variant-calling.sh 
```

The job script in turn calls another script that will run the variant calling. Let us dwonload that script first

```sh
wget https://raw.githubusercontent.com/sara-nl/2019-08-28-htdp-elixir/gh-pages/_episodes/scripts/run-variant-calling.sh
```

Let us submit the job first and then inspect the steps while the job runs

```sh
sbatch --job-name=var-call -J 'var-call' --output=%x-%j.out job-submit-variant-calling.sh
squeue -u $USER

cat run-variant-calling.sh

#!/bin/bash
set -e
ecolipath=/project/spidercourse/Data/ecoli-analysis

cd $ecolipath/results

genome=$ecolipath/data/ref_genome/ecoli_rel606.fasta

bwa index $genome

mkdir -p sam bam bcf vcf

for fq1 in $ecolipath/data/trimmed_fastq_small/*_1.trim.sub.fastq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.sub.fastq)
    echo "base name is $base"

    fq1=$ecolipath/data/trimmed_fastq_small/${base}_1.trim.sub.fastq
    fq2=$ecolipath/data/trimmed_fastq_small/${base}_2.trim.sub.fastq
    sam=$ecolipath/results/sam/${base}.aligned.sam
    bam=$ecolipath/results/bam/${base}.aligned.bam
    sorted_bam=$ecolipath/results/bam/${base}.aligned.sorted.bam
    raw_bcf=$ecolipath/results/bcf/${base}_raw.bcf
    variants=$ecolipath/results/bcf/${base}_variants.vcf
    final_variants=$ecolipath/results/vcf/${base}_final_variants.vcf 

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done
```

Let us see if the job is running and what it is doing. You can inspect the output log file even if the job is not completed. 

```sh
squeue -u $USER
cat var-call-jobid.out #replace the jobid with your jobid 
```

You probably received an error that says

```sh
[bwa_index] Pack FASTA... [bns_fasta2bntseq] fail to open file '/project/spidercourse/Data/ecoli-analysis/data/ref_genome/ecoli_rel606.fasta.pac' : Permission denied
```

> **_Food for brain:_**
>
> * This error indicates that it failed to open a file, do you know why? Hint: check if such a file exists in this path
> * The project Data folder path is provided in the script, check the path $HOME/ecoli-analysis/data/ref_genome/ and you can see that no such file exists. So what is going on? Why is it trying to open this file?

The bwa tool is trying to create the file ecoli_rel606.fasta.pac in the Data project space where as you know you do not have write permissions. How can you fix this? Try the following:

```sh
In the job-submit-variant-calling.sh script replace the following line 

bash /project/spidercourse/Data/ecoli-analysis/run-variant-calling.sh 

to

bash $HOME/ecoli-analysis/run-variant-calling.sh 

-----------------
In the run-variant-calling.sh script replace the path

ecolipath=/project/spidercourse/Data/ecoli-analysis

to

ecolipath=$HOME/ecoli-analysis
```

And run the job again

```sh
sbatch --job-name=var-call -J 'var-call' --output=%x-%j.out job-submit-variant-calling.sh
squeue -u $USER
```
So did the job run properly this time? Check the log file

```sh
squeue -u $USER
cat var-call-jobid.out #replace the jobid with your jobid 
```

To do - The above will fail as the output will still be written to Data folders. Introduce the Shared space or make them do it in home. Introduce the 'overwrite in share' space errors and then indicate the correct paths

extras
--check if tmpdir performs with a significant difference
--they can run the full example if they are ready

