## Getting started


Begin by logging into VM:

```
ssh -X ubuntu@137.205.69.49
```

May need this:
```
source .profile 
```

Discussion point environment variables and configuration files.

Some useful repositories should already be on the system (do not run the following):


~~mkdir ~/repos~~

~~cd repos~~

~~git clone https://github.com/chrisquince/WorkshopSept2017.git~~


and another two we will need:


~~git clone https://github.com/chrisquince/MAGAnalysis.git~~
~~git clone https://github.com/chrisquince/StrainMetaSim.git~~

## Data sets

All the tutorial data sets are available in the Data directory go and have a look:

```
cd Data
ls
```

Have a look in each sub directory. Then go into the AD folder and look at the meta data file:

```
cd AD
more Meta.csv 
```

And the reads

```
cd ReadsSub
```

### Fastq file format

The reads are stored as pairs of fastq files.
```
head -n 10 S102_Sub_R1.fastq
```
Sometimes these will be zipped.

Lets have a look at the wikipedia page on [fastq format](https://en.wikipedia.org/wiki/FASTQ_format).

Count up number of reads in a fastq file:

```
cat S102_Sub_R1.fastq | echo $((`wc -l`/4))
```

Try to understand the anatomy of this command. What does 'cat' do and the '|' what about 'echo' and 'wc'? 


What will the number be in S102_Sub_R2.fastq?

Discussion point: What are paired end reads?

In fact these reads have been subsampled to make the tutorial (barely) tractable. This 
was done with the following commands (do not run):

~~cd ~/Projects/AD~~
~~mkdir ReadsSub~~
~~for file in Reads/*R1*fastq~~
~~do
    ~~base=${file##*/}~~
    ~~stub=${base%_R1.fastq}~~
    ~~echo $stub~~
    ~~seqtk sample -s100 $file 1000000 > ReadsSub/${stub}_Sub_R1.fastq&~~
    ~~seqtk sample -s100 Reads/${stub}_R2.fastq 1000000 > ReadsSub/${stub}_Sub_R2.fastq&~~
~~done~~


### Now lets look at the gut samples

```
cd ~/Data/Gut
```

These are in different sequence format, what is it?

```
head -n 10 Reads/C13_R1.fasta
```

This is a convenient way to count reads in a fasta file:
```
grep -c ">" Reads/C13_R1.fasta
```


### And we have one ancient DNA sample

```
cd ~/Data/Ragna
```

What do you notice about this file format. Can you figure out a command to count the number of reads in this file the key is the command 'zcat'

### And some synthetic samples

```
cd ~/Data/Synthetic
ls
```

What is the 'Genomes.tar.gz'? Here is how we can extract it:

```
tar -xvzf Genomes.tar.gz
```

## Sequence quality control

The first thing you will want to do with raw sequence is check its quality. 

Go to the home directory:

```
cd ~
```

and make a directory called 'Projects' and go into it:

```
cd Projects
```

Lets make a sub-directory for the AD analysis:

```
mkdir AD
cd AD
```

And a directory for the AD analysis. Now we want to conveniently access the Data but we 
want that directory kept clean. We can use a symbolic link for this:

```
ln -s ~/Data/AD/ReadsSub ReadsSub
```

Now we run fastqc on one of the samples:
```
fastqc ./ReadsSub/S102_Sub_R1.fastq
fastqc ./ReadsSub/S102_Sub_R2.fastq
```

Look at the output files (probably best to download to local computer):
```
scp ubuntu@137.205.69.49:~/Projects/AD/ReadsSub/*fastqc*html .
```

You should see something like this:

[S102_Sub_R1_trimmed_fastqc.html](https://github.com/chrisquince/CopenhagenWorkshop/blob/master/Results/FastQC/S102_Sub_R1_trimmed_fastqc.html)


[S102_Sub_R2_trimmed_fastqc.html](https://github.com/chrisquince/CopenhagenWorkshop/blob/master/Results/FastQC/S102_Sub_R2_trimmed_fastqc.html)

