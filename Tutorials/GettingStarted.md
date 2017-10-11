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

In fact these reads have been subsampled to make the tutorial (barely) tractable.

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

What is the 'Genomes.tar.gz' here is how we can extract it:

```
tar -xvzf Genomes.tar.gz
```


###Subsampling reads


Lets counts up reads in all files:
```
for file in *_R1.fastq; do cat $file | echo $((`wc -l`/4)); done > Counts.txt
```

And plot histogram, median read number is 4,917,354:
```
R
>library(ggplot2)
>Counts <- read.csv('Counts.txt',header=FALSE)
>summary(Counts$V1)
>pdf("Counts.pdf")
>qplot(Counts$V1, geom="histogram") 
>dev.off()
>q()
```

Then we make ourselves a Projects directory:

```
mkdir ~/Projects
mkdir ~/Projects/AD
cd ~/Projects/AD
mkdir Reads
cd Rea


Now we run fastqc on one of the samples:
```
fastqc S102_R1.fastq
```

Look at the output files:
```
ls
firefox S102_R1_fastqc.html 