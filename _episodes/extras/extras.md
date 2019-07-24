# Running jobs with additional features

1. [Using local scratch on worker nodes](#job-tmpdir)

### <a name="sjob-tmpdir"></a> 1. Using local scratch on worker nodes

In all the previous jobs, the jobs had input as well as output in your project space (on CephFS; Ceph File System). The 
Spider worker nodes have a large scratch area on local SSD, particularly efficient for large I/O. Here we will run a job where
you can copy input/output to/from the local scratch.

```sh
cd $HOME/ecoli-analysis
wget https://raw.githubusercontent.com/sara-nl/2019-08-28-htdp-elixir/gh-pages/_episodes/scripts/job-submit-variant-calling-tmpdir.sh
wget https://raw.githubusercontent.com/sara-nl/2019-08-28-htdp-elixir/gh-pages/_episodes/scripts/run-variant-calling-tmpdir.sh
```
The workflow is the same as we ran in the previous example, except first we copy the files and scripts to the local scratch space on the worker node where your job lands. Let us inspect these scripts and see the difference with the previous example.

```sh
cat job-submit-variant-calling-tmpdir.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis/run-variant-calling-tmpdir.sh .

time bash run-variant-calling-tmpdir.sh 
```
Here we first created a directory with the help of a globally defined variable $TMPDIR. This directory will be created at the start of the job on the local scratch space and removed when the job is done. We copy the variant calling script to this directory and run it. To compare the performance with jobs that ran on the project spaces, we will 'time' the job - this will tell us how long it took for the full job to finish.

```sh
cat run-variant-calling-tmpdir.sh

#!/bin/bash
set -e
ecolipath=$PWD

mkdir -p data/ref_genome
cp $HOME/ecoli-analysis/data/ref_genome/ecoli_rel606.fasta data/ref_genome/
ls data/ref_genome

mkdir data/trimmed_fastq_small
cp $HOME/ecoli-analysis/data/trimmed_fastq_small/*fastq data/trimmed_fastq_small/
ls data/trimmed_fastq_small

mkdir results
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

cp -r $PWD/results $HOME/ecoli-analysis/
```
Here we copy the input data to the $TMPDIR. The parent paths are redefined and hence the rest of the workflow remains the same. In the end we copy the output to our $HOME directory as the $TMPDIR is removed after thew job finishes amd we will lose our results. You can run this example and compare if the performance was better/worse/equivalent to the performance with the jobs when the data is in project spaces.

```sh
#Make sure the path to store the results in the variant calling script does not already have the results

sbatch --job-name=var-call-tmpdir -J 'var-call-tmpdir' --output=%x-%j.out job-submit-variant-calling-tmpdir.sh
```

> **_Food for brain:_**
>
> * What does the time command do? How do you interpret the output?
> * You need to rerun the previous example with data in the project space by adding the 'time' command.
> * Does the $TMPDIR example have better performance? When is it advantageous to use it?
