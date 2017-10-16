## Getting started


Begin by logging into VM:

```
ssh -X ubuntu@137.205.69.49
```

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

Explore results any problems?

In actual fact some of the reads have transposon adaptors left on. We can remove these with trim_galore:

```
trim_galore ReadsSub/S102_Sub_R1.fastq 
```

Rerun fastqc - is it better. Lets tidy up and trim all samples...


```
mkdir ReadsTrim
cd ReadsTrim
for file in ../*ReadsSub/*fastq; do trim_galore $file ; done
```

Try to understand bash for loop structure it is very useful!

## SeqQC for ancient DNA sample

I want you to run fastqc and trim_galore on the Ragna sample by yourselves. Put the results in a directory ~/Projects/Ragna and ~/Projects/Ragna/Trim.

## Removing human contamination 

The easiest way to remove any contaminating DNA is by simply mapping onto that genome. Obviously this assumes that you have the genome. The Ragna sample will contain human DNA
that may interfere with downstream analysis. We can remove the human reads as follows:

```
cd ~/Projects/Ragna
```

We have saved a copy of the human genome on the server:

```
more ~/Databases/Human/human_g1k_v37.fasta
```

This has been indexed previously:

~~bwa index human_g1k_v37.fasta~~

Enabling us to map reads against it with BWA:

```
mkdir MapHuman
cd MapHuman
bwa mem -t 8 ~/Databases/Human/human_g1k_v37.fasta ../Trim/sk152_dentine_L1.1_trimmed.fq.gz  > sk152_dentine_trimmed_maphuman.sam
```

What does the '-t 8' flag do?

This generates a 'sam file' alignment of the trimmed reads against the human genome. Have a look at it? Make sense?


We convert this sam to bam...

```
samtools view -bS sk152_dentine_trimmed_maphuman.sam > sk152_dentine_trimmed_maphuman.bam
```

and then filter out the *unmapped* reads only:

```
cd ..
mkdir NotHuman
samtools view -b -f 4 MapHuman/sk152_dentine_trimmed_maphuman.bam > NotHuman/${stub2}.bam
```

Now we use a third party tool to turn the bam alignment back into fastq:

```
bamToFastq -i NotHuman/sk152_dentine_nothuman.bam -fq NotHuman/sk152_dentine_nothuman.fq
```

How many reads did not map to human?

Also separate of those that did map...

```
samtools view -b -F 4 MapHuman/sk152_dentine_trimmed_maphuman.bam > MapHuman/sk152_dentine_human.bam

bamToFastq -i MapHuman/sk152_dentine_human.bam -fq MapHuman/sk152_dentine_human.fq

```

Just for fun lets calculate the coverage of human chromosomes:

```
samtools sort MapHuman/sk152_dentine_human.bam > MapHuman/sk152_dentine_human.sorted.bam

bedtools genomecov -ibam MapHuman/sk152_dentine_human.sorted.bam -g ~/Databases/Human/human_g1k_v37.len > MapHuman/sk152_dentine_human_cov.txt
```

Calculates coverage histograms across positions use some awk to calculate mean chromosome coverages:
```
awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' MapHuman/sk152_dentine_human_cov.txt > MapHuman/sk152_dentine_human_cov.csv
```

What sex was Ragna?
```
more MapHuman/sk152_dentine_human_cov.csv
```