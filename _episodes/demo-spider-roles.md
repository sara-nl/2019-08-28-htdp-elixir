# Spider project spaces and roles

1. [Project environment](#spider-spaces)
2. [Software management](#spider-sm)
3. [Data management](#spider-dm)

### <a name="spider-spaces"></a> 1. Spider project environment

#### 1.1 Project environment

Familiarize yourself with your environment :

 ```sh
 ls /project 
 ```
> **_Food for brain:_**
>
> * Do you know what project you belong to? What all access does it provide to you?
> * What are each of ther project directories for? What is public/private? Do you have read write permissions ion all spaces?

### <a name="job-submit"></a> 2. Software management

#### 2.1 Software manager Role

 ```sh
 id $USER
 ```
 
> **_Food for brain:_**
>
> * Are you a software manager? If not, do you know who is the software manager?

 ```sh
 getent group spidercourse-sw
 getent group spidercourse-user
 ```
 
#### 2.2 Miniconda installation

To install software for the project users, you should be a software manager. We will use Miniconda which is a package manager that simplifies the installation process. Please first install miniconda3 and then proceed to the installation of individual tools.

 ```sh
 wget https://repo.continuum.io/miniconda/Miniconda2-4.6.14-Linux-x86_64.sh
 bash Miniconda2-4.6.14-Linux-x86_64.sh
 ```

It will ask you for an installation path. If you are not a software manager you cannot use the project Software space. But a regular user can also install the software in their $HOME. 

 ```sh
 #Please provide the following path for installation as a software manager
 /project/surfadvisors/Software/ecoli-analysis-software/miniconda3 

 or 

 #Please provide the following path for installation as a regular user
 $HOME/ecoli-analysis-software

 exit 
 ```

Login again to Spider and inspect what environment variables have been set up

 ```sh
 cat $HOME/.bashrc
 ```

#### 2.3 Variant calling tools installation

Follow the further instructions for the installation of individual tools

 ```sh
 conda install -c bioconda fastqc=0.11.7=5

 conda install -c bioconda trimmomatic=0.38=0

 conda install -c bioconda bwa=0.7.17=ha92aebf_3

 conda install -c bioconda samtools=1.9=h8ee4bcc_1

 conda install -c bioconda bcftools=1.8=h4da6232_3 
 ```
 
> **_Food for brain:_**
>
> * Do you know where the above pacakages are installed? How do you test if the installation was successful?


### <a name="spider-dm"></a> 3. Data management

#### 3.1 Data manager Role

 ```sh
 id $USER
 ```
 
> **_Food for brain:_**
>
> * Are you a data manager? If not, do you know who is the data manager?

 ```sh
 getent group spidercourse-data
 getent group spidercourse-user
 ```
 
#### 3.2 Data download

Let us download some data that we will use later to run jobs on the cluster. The data we are going to use is part of a long-term evolution experiment led by [Richard Lenski](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment) to assess adaptation in E. coli. A population was propagated for more than 50,000 generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations. 

Let us download the paired-end data from [European Nucleotide Archive](https://www.ebi.ac.uk/ena).

 ```sh
 cd /project/spidercourse/Data/ecoli-analysis/data
 mkdir untrimmed_fastq/
 cd untrimmed_fastq

 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz 
 ```

Let us also download the reference genome for E. coli REL606.

 ```sh
 cd /project/spidercourse/Data/ecoli-analysis/data
 mkdir ref_genome
 cd ref_genome
 curl -L -o ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
 gunzip ecoli_rel606.fasta.gz
 ```

You may also download the above data in your $HOME directory instead of project Data space. However, please note that for running the examples today the paths defined in the workflows expect the data to be in the project Data space. You will need to modify all the paths in the further scripts if you download the data in your $HOME.

#### 3.3 Data cleanup

Now that you have the raw data, we will assess the quality of the sequence reads contained in our fastq files and run filtering.

```sh
cd /project/surfadvisors/Data/ecoli-analysis
cat job-submit--datacleaning.sh

#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --constraint=skylake
time bash /project/spidercourse/Data/ecoli-analysis/read_qc.sh 
```

This script in turns call another script which we can inspect further

```sh
cat read_qc.sh 

set -e
ecolipath=/project/surfadvisors/Data/ecoli-analysis

mkdir -p $ecolipath/data/fastqc_untrimmed_reads

cd $ecolipath/data/fastqc_untrimmed_reads/

echo "Running FastQC ..."
fastqc $ecolipath/data/untrimmed_fastq/*.fastq* -o ./ 

cd $ecolipath/data/untrimmed_fastq
cp /project/surfadvisors/Software/dc-genomics/miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/NexteraPE-PE.fa .
echo "Running trimmomatic"
for infile in *_1.fastq.gz
do
   base=$(basename ${infile} _1.fastq.gz)
   trimmomatic PE ${infile} ${base}_2.fastq.gz \
               ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
               ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
               SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
done

mkdir $ecolipath/data/trimmed_fastq	
mv *.trim* $ecolipath/data/trimmed_fastq

echo "done"
```

#### Acknowledgements 
This example was adopted from https://datacarpentry.org/wrangling-genomics/ 

 
