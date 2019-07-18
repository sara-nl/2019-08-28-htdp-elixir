# Various roles in Spider project spaces

1. [Spider_roles](#spider-roles)
2. [Data management](#spider-dm)
3. [Software manager role](#spider-sm)

### <a name="spider-roles"></a> 1. Spider roles

#### 1.1 Project environment

Familiarize yourself with your environment :

 ```sh
 ls /project 
 ```
> **_Food for brain:_**
>
> * Do you know what project you belong to? What all access does it provide to you?
> * What are each of ther project directories for? What is public/private? Do you have read write permissions ion all spaces?
  
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
 getent group spidercourse-sw
 getent group spidercourse-user
 ```
 
#### 2.2 Data download

Let us download some data that we will use later to run jobs on the cluster. The data we are going to use is part of a long-term evolution experiment led by [Richard Lenski](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment) to assess adaptation in E. coli. A population was propagated for more than 50,000 generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations. 

Let us download the paired-end data from [European Nucleotide Archive](https://www.ebi.ac.uk/ena).

```sh
cd /project/spidercourse/Data
mkdir -p ecoli-analysis/data/untrimmed_fastq/
cd ecoli-analysis/data/untrimmed_fastq

curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz 
```

Let us also download the reference genome for E. coli REL606.

```sh
mkdir -p data/ref_genome
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip data/ref_genome/ecoli_rel606.fasta.gz
```

### <a name="job-submit"></a> 3. Software manager role

  

#### Acknowledgements 
This example was adopted from https://datacarpentry.org/wrangling-genomics/ 

 
