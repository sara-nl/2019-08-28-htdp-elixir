# Running jobs with additional features


1. [Using local scratch on worker nodes](#job-tmpdir)
2. [Accessing external storage](#dcache)


### <a name="sjob-tmpdir"></a> 1. Using local scratch on worker nodes

In all the previous jobs, the jobs had input as well as output in your project space (on CephFS; Ceph File System). The 
Spider worker nodes have a large scratch area on local SSD, particularly efficient for large I/O. Here we will run a job where
you can copy input/output to/from the local scratch.

```sh
cd $HOME/ecoli-analysis
wget https://raw.githubusercontent.com/sara-nl/2019-08-28-htdp-elixir/gh-pages/_episodes/scripts/job-submit-variant-calling-tmpdir.sh
wget https://raw.githubusercontent.com/sara-nl/2019-08-28-htdp-elixir/gh-pages/_episodes/scripts/run-variant-calling-tmpdir.sh
```
