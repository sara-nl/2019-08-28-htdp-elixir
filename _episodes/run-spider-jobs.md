# Running jobs on Spider

Running a life science example on spider with the dataq in data folder and preferably software in software folder 

We will also download a set of trimmed FASTQ files to work with. These are small subsets of our real trimmed data, and will enable us to run our variant calling workflow quite quickly.

cd /project/surfadvisors/Data/ecoli-analysis/data
mkdir trimmed_fastq_small
cd trimmed_fastq_small
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
