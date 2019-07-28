# Running jobs on Spider


1. [Set up job enviromment](#job-setup)
2. [Run variant calling](#run-var-call)
3. [Handling errors](#error-handling)
4. [Sharing results in a project](#share-data)

### <a name="spider-spaces"></a> 1. Set up job environment

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
### <a name="run-var-call"></a> 2. Run variant calling jobs

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

#index the reference genome
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
    
    #Align reads to reference genome
    bwa mem $genome $fq1 $fq2 > $sam
    
    #Convert from sam to bam format
    samtools view -S -b $sam > $bam

    #Sort the bam files    
    samtools sort -o $sorted_bam $bam 
    
    #Calculate the read coverage of positions in the genome
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    
    #Detect the single nucleotide polymorphisms (SNPs)
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    
    #Filter the SNPs for the final output in VCF format
    vcfutils.pl varFilter $variants > $final_variants
   
    done
```

Let us see if the job is running and what it is doing. You can inspect the output log file even if the job is not completed. 

```sh
squeue -u $USER
cat var-call-jobid.out #replace the jobid with your jobid 
```
### <a name="error-handling"></a> 3. Handling errors

You probably received an error that says

```sh
[bwa_index] Pack FASTA... [bns_fasta2bntseq] fail to open file '/project/spidercourse/Data/ecoli-analysis/data/ref_genome/ecoli_rel606.fasta.pac' : Permission denied
```

> **_Food for brain:_**
>
> * This error indicates that it failed to open a file, do you know why? Hint: check if such a file exists in this path
> * The project Data folder path is provided in the script, check the path $HOME/ecoli-analysis/data/ref_genome/ and you can see that no such file exists. So what is going on? Why is it trying to open this file?

The bwa tool is trying to create the file ecoli_rel606.fasta.pac in the Data project space where as you know you do not have write permissions. How can you fix this? 

1. Fix paths in job-submit-variant-calling.sh

```sh
#In the job-submit-variant-calling.sh script replace the following line 

bash /project/spidercourse/Data/ecoli-analysis/run-variant-calling.sh 

#to

bash $HOME/ecoli-analysis/run-variant-calling.sh 
```

2. Fix paths in run-variant-calling.sh 

```sh

#In the run-variant-calling.sh script replace the path

ecolipath=/project/spidercourse/Data/ecoli-analysis

#to

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
You can see that now the job runs properly which is great.

> **_Food for brain:_**
>
> * Often in a project you want to share results with your colleagues. Your $HOME is not accessible to other members in the project so how would you share the results? Hint: check the folders in your /project/spidercourse directory
> * Do you have write access to such a folder? Do all project members have read and write access to a common folder?

### <a name="share-data"></a> 4. Sharing results in a project

You have already been introduced to the Data and Software project spaces and the associated roles. A shared Data folder reduces unnecessary replicas of the same data, and having a manager to handle the data reduces the risk of accidental removal. However, you also need a space where the results can be shared and the benefits are twofold - not everyone needs to run the same analysis (and co-ordinate on how this can be done uniformly) and the results can be used for post processing without the hassle of transfering the data to each other. This functionality is provided by the Shared project space.

```sh
cd /project/spidercourse/Share
mkdir my-results
ls -l
```

> **_Food for brain:_**
>
> * So you created this folder but looks like someone else became the owner. Any idea why? Hint: This is a "Shared" directory.
> * Can you delete the directories/files not owned by you?

So the etiquette to keep in mind when using the Share folder is that you can acidentally delete/overwrite someone else's
results and vice versa, or the entire results of your project. So be careful!
 
```sh
cd /project/spidercourse/Share
mkdir $USER-results
ls -l
```

You may copy your results to the folder you created above to share the results with your colleagues. Please check first if your job is finished, else you will copy partial results only.

```sh
cd /project/spidercourse/Share/$USER-results
cp -r $HOME/ecoli-analysis/results .
```

You can also avoid having to copy results back and forth by running the analysis in the Share folder. 

> **_Food for brain:_**
>
> * What paths need to be changed for the analysis to write the output in Share folder?
> * How can you make sure you don't accidentally overwrite someone else's/your own results? 

Once you have the results, lets see how many variants are in the vcf files

```sh
cd $HOME/ecoli-analysis
grep -v "#" results/vcf/SRR2589044_final_variants.vcf | wc -l
```

> **_Food for brain:_**
>
> * How many variants are detected in each file? 
> * Solution: For SRR2589044 from generation 5000 there were 10 mutations, for SRR2584863 from generation 15000 there were 25 mutations, and SRR2584866 from generation 766 mutations. In the last generation, a hypermutable phenotype had evolved, causing this strain to have more mutations.

Now that you finished your analysis for a small dataset, shall we see how much resources did you use on the cluster? You can use the following command to find this information:

```sh
sacct --format=JobID,JobName,AveCPU,Elapsed  # sacc displays accounting data for all jobs
```

> **_Food for brain:_**
>
> * What do these numbers indicate?
> * How can you find more information about your jobs e.g., the maximum memory? Hint: add the flag "MaxRSS" in your query
> * How much time do you estimate your typical workflow will take to run on Spider - at the moment it has 240 cores?
> * How can your work benefit from running on Spider e.g., with project spaces, public views? In the enxt session we will demonstrate yet more features on Spider

There are some extra examples [here](https://github.com/sara-nl/2019-08-28-htdp-elixir/blob/gh-pages/_episodes/extras/extras.md) if you have already finished the previous examples and would like to explore more.
