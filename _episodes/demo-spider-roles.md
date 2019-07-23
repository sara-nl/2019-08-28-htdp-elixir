# Spider project spaces and roles

1. [Project environment](#spider-spaces)
2. [Data management](#spider-dm)
3. [Software management](#spider-sm)
4. [Data cleanup](#data-cleanup)

### <a name="spider-spaces"></a> 1. Spider project environment

#### 1.1 Project environment

Familiarize yourself with your environment :

 ```sh
 ls /project 
 ```
> **_Food for brain:_**
>
> * Do you know what project you belong to? Hint: id $USER
> * What all access does it provide to you?
> * What are each of the project directories for? What is public/private? Do you have read write permissions on all spaces?

### <a name="spider-dm"></a> 2. Data management

#### 2.1 Data manager Role

 ```sh
 id $USER
 ```
 
> **_Food for brain:_**
>
> * Are you a data manager? If not, do you know who is the data manager?

 ```sh
 getent group spidercourse-data
 getent group spidercourse-user
 mkdir /project/spidercourse/Data/mydata
 ```
 
What happened? As you are not a data manager, you do not have write permissions in the project's Data directory. However, all users still have read permissions. This means data can be shared within the project without worrying about someone accidentally deleting/overwriting the data for the project.


#### 2.2 Data download

Let us download some data that we will use later to run jobs on the cluster. The data we are going to use is part of a long-term evolution experiment led by [Richard Lenski](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment) to assess adaptation in E. coli. A population was propagated for more than 50,000 generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations. 

Let us download the paired-end data from [European Nucleotide Archive](https://www.ebi.ac.uk/ena).

 ```sh
 
 #As data manager
 mkdir -p /project/spidercourse/Data/ecoli-analysis/data/untrimmed_fastq/
 cd /project/spidercourse/Data/ecoli-analysis/data/untrimmed_fastq/
 
 #As a regular user
 mkdir -p $HOME/ecoli-analysis/data/untrimmed_fastq/
 cd $HOME/ecoli-analysis/data/untrimmed_fastq/

 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz 
 ```

Let us also download the reference genome for E. coli REL606.

 ```sh
 #As data manager
 cd /project/spidercourse/Data/ecoli-analysis/data
 mkdir ref_genome
 cd ref_genome
 
 #As a regular user
 cd $HOME/ecoli-analysis/data
 mkdir ref_genome
 cd ref_genome
 
 curl -L -o ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
 gunzip ecoli_rel606.fasta.gz
 ```

### <a name="job-submit"></a> 3. Software management

#### 3.1 Software manager Role

 ```sh
 id $USER
 ```
 
> **_Food for brain:_**
>
> * Are you a software manager? If not, do you know who is the software manager?

 ```sh
 getent group spidercourse-sw
 getent group spidercourse-user
 mkdir /project/spidercourse/Software/mysoftware
 ```
 
What happened? As you are not a Software manager, you do not have write permissions in the project's Software directory. However, all users still have read permissions. This means that you or one of your colleagues can install complicated software and the dependencies, maintain it, and all project users can use that uniformly. This makes life easier and is crucial for reproducability of your results. Also apart from project wide software, individual users can use their own software (or different versions). 

#### 3.2 Miniconda installation

We will use Miniconda; it is a package manager. We shall first install miniconda2 and then proceed to the installation of individual tools. Below are two set of instructions - one for Sfotware manager and one for regular user. 

 ```sh
 #As Software manager
 cd /project/spidercourse/Software/ 
 
 #As a regular user
 cd $HOME
 
 wget https://repo.continuum.io/miniconda/Miniconda2-4.6.14-Linux-x86_64.sh
 bash Miniconda2-4.6.14-Linux-x86_64.sh
 ```

It will ask you for an installation path. (please answer yes)

 ```sh
 yes 
 ```

It will ask you for an installation path. 

 ```sh
 #As a software manager
 /project/spidercourse/Software/ecoli-analysis-software/miniconda2 

 #As a regular user
 $HOME/ecoli-analysis-software/miniconda2 

 exit 
 ```

We will see later how Software installed in the Software project space can be used by all users. If you installed it in your $HOME directoy, other users in the project cannot access it (unless you open up permissions which is NOT recommended). Login again to Spider and inspect what environment variables have been set up

 ```sh
 cat $HOME/.bashrc
 ```

#### 3.3 Variant calling tools installation

Follow the further instructions for the installation of individual tools

 ```sh
  #As Software manager
 cd /project/spidercourse/Software/ecoli-analysis-software 
 
 #As a regular user
 cd $HOME
 
 conda install -c bioconda fastqc=0.11.7=5

 conda install -c bioconda trimmomatic=0.38=0

 conda install -c bioconda bwa=0.7.17=ha92aebf_3

 conda install -c bioconda samtools=1.9=h8ee4bcc_1

 conda install -c bioconda bcftools=1.8=h4da6232_3 
 ```
 
> **_Food for brain:_**
>
> * Do you know where the above pacakages are installed? How do you test if the installation was successful? #Hint - try running one of the above installed software

```sh
fastqc -h       (runs quality control checks on raw sequence data)

trimmomatic     (data trimming tool)

bwa             (maps DNA sequences against reference genome)

samtools.       (utility for manipulating data in the SAM format)

bcftools        (variant calling tool)
```

For the last two tools you can see that a library is missing! This will very often be the case that not every installation will be completely successful on every system. If the Software manager resolves these issues, users can freely use the software and avoid hassles of resolving software dependencies. The Softare manager has already resolved this for you and instead of troubleshooting you can already start using this installation by performing the following steps:

```sh
nano $HOME/.bashrc

#In the conda initialize set up replace all the paths 

# >>> conda initialize >>>

current path in your file 

$HOME/ecoli-analysis-software/miniconda2/bin/conda 

replace this with the following path instead

/project/spidercourse/Software/ecoli-analysis-software/miniconda2/bin/conda

save the changes

exit 
```
Now you are using the software environment that was set up by the Software manager.

```sh
samtools
bcftools

echo $PATH
```
As you can see the error is resolved and you can proceed to running the analysis. Just sop you know in this particular case it was a simple solution, the missing library was downloaded to Spider and copied as follows: 

```sh
cp /scratch/libcrypto.so.1.0.0 ecoli-analysis-software/miniconda2/lib/
```

### <a name="data-cleanup"></a> 4. Data cleanup

Now that you have the raw data and the software installed, we will assess the quality of the sequence reads contained in our fastq files and run filtering.

```sh
#As data manager
cd /project/spidercourse/Data/ecoli-analysis/

# As a regular user
cd $HOME/ecoli-analysis

wget https://raw.githubusercontent.com/sara-nl/2019-08-28-htdp-elixir/gh-pages/_episodes/scripts/job-submit-datatrimming.sh
 
cat job-submit-datatrimming.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

#As data manager uncomment the following line
#bash /project/spidercourse/Data/ecoli-analysis/data_qc.sh 

#As a regular user uncomment the following line
bash $HOME/ecoli-analysis/data_qc.sh 
```

Let us inspect what steps we follow in the data trimming

```sh
wget https://raw.githubusercontent.com/sara-nl/2019-08-28-htdp-elixir/gh-pages/_episodes/scripts/data_qc.sh

cat data_qc.sh 

#!/bin/bash
set -e

#if you are a data manager, uncomment the following line (remove the #)
#ecolipath=/project/spidercourse/Data/ecoli-analysis

#if you are a regular user, uncomment the following line (remove the #)
#ecolipath=$HOME/ecoli-analysis

mkdir -p $ecolipath/data/fastqc_untrimmed_reads

cd $ecolipath/data/fastqc_untrimmed_reads/

echo "Running FastQC ..."
fastqc $ecolipath/data/untrimmed_fastq/*.fastq* -o ./ 

cd $ecolipath/data/untrimmed_fastq
cp /project/spidercourse/Software/ecoli-analysis-software/miniconda2/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/NexteraPE-PE.fa .
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
Now let's submit the job

```sh
sbatch --job-name=data-trim -J 'data-trim' --output=%x-%j.out job-submit-datatrimming.sh
squeue -u $USER
```
This job will create trimmed data file and files in .html format that can be used to display the summary report.

---Depending on whether the public view issue is resolved or not, make them copy the html files to their laptop and view it, or use the Public view - this is the short step to introduce public data folder as a nice functionality.
Also if this is getting longer thasn sn hour skip the fastqc steps - or do them and use ir in the later session

#### Acknowledgements 
This example was adopted from https://datacarpentry.org/wrangling-genomics/ 

 
