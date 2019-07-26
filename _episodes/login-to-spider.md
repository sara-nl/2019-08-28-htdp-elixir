# Logging onto Spider

Here are your first steps:

1. [Login to Spider](#spider-login)
2. [Get familiar with the login node](#spider-env)
3. [Submit jobs](#job-submit)

### <a name="spider-login"></a> 1. Login to Spider

The login credentials (e.g., spidercourse-g01) will be  provided to you at the start of the session. Open a terminal in your laptop 
    
 ```sh
 ssh username@spider.surfsara.nl #replace `username` with the username assigned to you
 ```
  
### <a name="cartesius-env"></a> 2. Get familiar with the login node

Familiarize yourself with your environment :

 ```sh
 whoami
 id $USER
 pwd
 sinfo
```

> **_Food for brain:_**
>
> * What does the output of the command "id $USER" tell you?
> * How many worker nodes are available on the cluster?


### <a name="job-submit"></a> 3. Submit your first job

#### 3.1 Submit a simple job 

Let us download our first job script and inspect it:
  
 ```sh
 wget https://raw.githubusercontent.com/sara-nl/2019-08-28-htdp-elixir/gh-pages/_episodes/scripts/my-first-job.sh

 cat my-first-job.sh
 
 #!/bin/bash
 #SBATCH -t 10:00
 #SBATCH -c 1
 #SBATCH --constraint=skylake
 sinfo
 squeue
 echo "You just ran your first job on" $HOSTNAME " with a job ID " $SLURM_JOBID
 ```
 -t: max total run time of the job, here it is 10 minutes  
 -c: 1 core requested 
 --constraint: here we request nodes with skylake processors and ssd local scratch disk space
 
Now that you have inspected the script that will submit your job, let's submit it running the following command:
  
 ```sh
 sbatch my-first-job.sh  #This command will submit a job and give you a job ID in return
 squeue -u $USER  #Check the status of your job
 ls  #Check if the output,slurm-jobID.out, is present
 cat slurm-yourjobid.out
 ```
 
