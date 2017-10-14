# Strain Resolution

## Getting started

This tutorial describes a complete walkthrough for CONCOCT2 and DESMAN on the STAMPS servers using 4 threads. Some of of the pre-processing steps are time consuming and not that informative. These have been pre-run and you will just copy across the output.
We will begin by constructing a new working directory in our home directories:

```bash
cd ~/Projects
mkdir Synthetic
cd Synthetic
```

The tutorial data consists of 32 synthetic samples of 500,000 2X150 bp paired end reads taken from 8 species and 15 genomes. Copy across the simulation config file just to have a look at it with the less command:

```bash
cp ~/Data/Synthetic/coverage.tsv genome_coverage.tsv
less genome_coverage.tsv
```

This has format "Sample\tSpecies\Strain\Coverage\Relative frequency". We will now use CONCOCT/DESMAN to resolve these *de novo*.

The simulated reads and genomes are in the Synthetic Data folder

```
ln -s ~/Data/Synthetic/Reads Reads
ln -s ~/Data/Synthetic/Genomes Genomes
```

## Running CONCOCT2 and DESMAN on the 32 sample mock

We now describe how to perform a complete analysis binning and resolving strains on 
this synthetic community. Some of these steps are time consuming so we recommend 
that you copy across the output instead. This assumes that DESMAN and CONCOCT are installed and their paths set to the variables DESMAN and CONCOCT respectively e.g. (changing paths to your system set-up).

The first step in the analysis is to assemble the reads. 

### Coassembly

I previously co-assembled the reads using MEGAHIT 1.1.1 and default parameters:

```
ls Reads/*r1*gz | tr "\n" "," | sed 's/,$//' > r1.csv
ls Reads/*r2*gz | tr "\n" "," | sed 's/,$//' > r2.csv
megahit -1 $(<r1.csv) -2 $(<r2.csv) -t 4 -o Assembly 
```

This should take approximately 20 minutes. 

Evaluate the assembly quality with the script provided:

```bash
contig-stats.pl < Assembly/final.contigs.fa 
```

Results in output:
```bash
sequence #: 8602	total length: 31455569	max length: 1071063	N50: 15733	N90: 1419
```

Do you think this is a good co-assembly?

### Calculating contig coverages

The next step in CONCOCT is to cut up contigs and map the reads back onto them. Again 
__do not run__ any of the following:


```
cd Assembly
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m final.contigs.fa > final_contigs_c10K.fa
bwa index final_contigs_c10K.fa
cd ..
```

Then we perform mapping onto contigs:

```
mkdir Map

for file in Reads/*.r1.fq.gz
do

    stub=${file%.r1.fq.gz}

    base=${stub##*/}

    echo $base

    bwa mem -t 8 Assembly/final_contigs_c10K.fa $file ${stub}.r2.fq.gz > Map/$base.sam

done
```

And calculate coverages:

```
python $DESMAN/scripts/Lengths.py -i Assembly/final_contigs_c10K.fa | tr " " "\t" > Assembly/Lengths.txt

for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub  
    (samtools view -h -b -S $file > ${stub}.bam; samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam; samtools sort -m 1000000000 ${stub}.mapped.bam -o ${stub}.mapped.sorted.bam; bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g Assembly/Lengths.txt > ${stub}_cov.txt)
done
```



Collate coverages together:

```
for i in Map/*_cov.txt 
do 
   echo $i
   stub=${i%_cov.txt}
   stub=${stub#Map\/}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > Map/${stub}_cov.csv
done

$DESMAN/scripts/Collate.pl Map > Coverage.csv
sed -i 's/\.\.\/CDTutorial\/Map\///g' Coverage.csv 

```
This just contains a header and a row for each contig givings its coverage (**how is this defined?**) in each sample. If you want to use CONCOCT with mappers other than bwa this is the input that is required.

### Contig binning

We are going to use a more efficient version of CONCOCT, CONCOCT2. I will not show you how to do the refinement step discussed in the presentation. If that interests you we can work through it in an open lab.
Now we can run CONCOCT:

```bash
mkdir Concoct

mv Coverage.csv Concoct

cd Concoct

tr "," "\t" < Coverage.csv > Coverage.tsv

concoct --coverage_file Coverage.tsv --composition_file ../Assembly/final_contigs_c10K.fa -t 4

```

As the program runs it outputs iteration number and fit and change in fit i.e.:

```
137,-171867.297247,0.000349
138,-171867.297516,0.000269
139,-171867.297723,0.000207
140,-171867.297881,0.000158
141,-171867.298002,0.000121
142,-171867.298094,0.000092
CONCOCT Finished, the log shows how it went.
```
as the clustering converges the improvement in fit should got to zero at which point the program terminates. 

We can see the arguments that CONCOCT takes as follows:
```bash
concoct 
```

The only two that generally need adjusting are the the number of initial clusters '-c' which defaults at 400 but should be at least twice the number of genomes in your assembly. See how to determine that below and the number of threads to use '-t' which ensures maximum performance on a given machine. The key output file is the cluster assignments.

```bash
head -n 10 clustering_gt1000.csv 
```

This gives for each contig the bin it is placed in as an integer index. To get the number of contigs in each of the 21 clusters that should result:

```bash
sed '1d' clustering_gt1000.csv | cut -d"," -f 2 | sort | uniq -c 
```

Other files generated by the program are:

*  original_data_gt1000.csv: the original data after normalisation    

* PCA_components_data_gt1000.csv: the PCA components

* PCA_transformed_data_gt1000.csv: the data after PCA transformation

We can also visualise the clusters in the PCA space, run:
```bash
Rscript ~/bin/ClusterPlot2.R -c clustering_gt1000.csv -p PCA_transformed_data_gt1000.csv -o clustering_gt1000_PCA.pdf
```
and download the resulting pdf it should look like this:

![PCA clusters](../Figures/clustering_gt1000_PCA.pdf)

### Contig annotation

To determine how good our binning is we need to annotate the contigs with clusters of orthologous genes (COGs). Find genes using prodigal:

```
    cd ..
    
    mkdir Annotate

    cd Annotate/

    python $DESMAN/scripts/LengthFilter.py ../Assembly/final_contigs_c10K.fa -m 1000 >     final_contigs_gt1000_c10K.fa

    prodigal -i final_contigs_gt1000_c10K.fa -a final_contigs_gt1000_c10K.faa -d     final_contigs_gt1000_c10K.fna  -f gff -p meta -o final_contigs_gt1000_c10K.gff 
```

Assign COGs change the -c flag which sets number of parallel processes appropriately:
```
    export COGSDB_DIR=/home/ubuntu/Databases/CDTutorial/rpsblast_cog_db/
    $CONCOCT/scripts/RPSBLAST.sh -f final_contigs_gt1000_c10K.faa -p -c 4 -r 1
```


We will also write out lists of cogs and genes in each contig. These will be useful later:
```
python $DESMAN/scripts/ExtractCogs.py -b final_contigs_gt1000_c10K.out --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv -g final_contigs_gt1000_c10K.gff > final_contigs_gt1000_c10K.cogs
python $DESMAN/scripts/ExtractGenes.py -g final_contigs_gt1000_c10K.gff > final_contigs_gt1000_c10K.genes
```


We will also need contig lengths later on:
```
python $DESMAN/scripts/Lengths.py -i final_contigs_gt1000_c10K.fa > final_contigs_gt1000_c10K.len
```

The annotation files are just simple comma separated lists, have a look at them:

```bash
head final_contigs_gt1000_c10K.cogs 
head final_contigs_gt1000_c10K.genes
```

Format is:

* gene_id,cog_id,startpos,endpos,cogstart,cogend,strand
* gene_id,contig_id,startpos,endpos,strand


### Determine number of complete genomes

We are going to determine the number of complete and pure genomes in our binning using single-core gene frequencies. Before doing that it is useful to estimate how many genomes 
we think are in our assembly. We can get this from the assembly by just calculating the 
median number of single-copy core genes:

```bash
?cd $METASIMWD
?cp /class/stamps-shared/CDTutorial/scgs.txt .
$CDSCRIPTS/CountSCGs.pl scgs.txt < Annotate/final_contigs_gt1000_c10K.cogs | cut -d"," -f2 | awk -f $CDSCRIPTS/median.awk
```

This should give **8** which since that is our species number is reassuring.
Now we calculate scg frequencies on the CONCOCT clusters:

```
cd Concoct
python $CONCOCT/scripts/COG_table.py -b ../Annotate/final_contigs_gt1000_c10K.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c clustering_gt1000.csv  --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_gt1000_scg.tsv
```

This produces a tab separated file of SCG frequencies, which we can then visualise in R:
```
Rscript $CONCOCT/scripts/COGPlot.R -s clustering_gt1000_scg.tsv -o clustering_gt1000_scg.pdf
```
You may need to download this pdf off the class servers to visualise.
We should see 8 nearly complete bins and some fragmentary ones:

![SCG Results](../Figures/clustering_gt1000_scg.pdf)

We want to store a list of 75% complete and pure for DESMAN analysis:
```
python $CDSCRIPTS/CompleteClusters.py clustering_gt1000_scg.tsv > Cluster75.txt
```

### How to calculate cluster coverages


In real studies it is important to know the coverage of each bin in each sample. This can then be transformed into relative copy numbers for correlation with meta-data:

```
python $CDSCRIPTS/ClusterMeanCov.py Coverage.csv clustering_gt1000.csv ../Assembly/final_contigs_c10K.fa > Cluster_Cov.csv
```

We can quickly use R to sum up the cluster coverages:

```
R
>Cluster_Cov <- read.csv("Cluster_Cov.csv",header=TRUE,row.names=1)
>sort(rowSums(Cluster_Cov))
>q()
```

We have two clusters, Cluster7 and Cluster16 that have a total coverage 
of > 100 and are 75% pure and complete, this is typically the minimum coverage necessary for strain resolution.

### Comparing to known reference genome assignments

These reads were generated synthetically _in silico_ from known genomes. It is therefore possible to compare the CONCOCT binning results to the assignments of contigs to genomes. This will not be possible for a genuine experimental metagenome. The assignments of contigs to genomes have been precomputed for you:

```bash
cd $METASIMWD
cp /class/stamps-shared/CDTutorial/AssignContigs.tar.gz .
tar -xvzf AssignContigs.tar.gz
```

Then move into this directory and run the following CONCOCT script which calculates accuracy statistics for a comparison of bins to reference genomes:

```bash
cd AssignContigs
$CONCOCT/scripts/Validate.pl --cfile=../Concoct/clustering_gt1000.csv --sfile=Contig_Species.csv --ffile=../Assembly/final_contigs_c10K.fa  
```

The results should looks as follows:
```
N	M	TL	S	K	Rec.	Prec.	NMI	Rand	AdjRand
5648	5642	2.7103e+07	8	21	0.765935	0.999235	0.869097	0.928564	0.695141
```
This confirms SCG plots we have some good bins high precision, but they are fragmented with lower recall.

We can also visualise the resulting confusion matrix.

```bash
Rscript $CONCOCT/scripts/ConfPlot.R -c Conf.csv -o Conf.pdf 
```

Which should give something like...

![Confusion Matrix](Figures/Conf.pdf)

## Run the DESMAN pipeline to resolve strains in each high quality bin

### Getting core variant frequencies

First we separate out the contigs, cogs, and genes into the separate bins:
```bash
cd ..
mkdir Split
cd Split
$DESMAN/scripts/SplitClusters.pl ../Annotate/final_contigs_gt1000_c10K.fa ../Concoct/clustering_gt1000.csv
$CDSCRIPTS/SplitCOGs.pl ../Annotate/final_contigs_gt1000_c10K.cogs ../Concoct/clustering_gt1000.csv
$CDSCRIPTS/SplitGenes.pl ../Annotate/final_contigs_gt1000_c10K.genes ../Concoct/clustering_gt1000.csv
cd ..
```

Then we select the SCGS for each cluster:
```bash
cd $METASIMWD
while read -r cluster 
do
    echo $cluster
    $DESMAN/scripts/SelectContigsPos.pl scgs.txt < Split/${cluster}/${cluster}.cog > Split/${cluster}/${cluster}_core.cogs
done < Concoct/Cluster75.txt
```

The first step in pre-processing for DESMAN would be to split up the bam files by each cluster in turn. **Do not** run this yourselves:

```
mkdir SplitBam

while read -r cluster 
do
    grep ">" Split/${cluster}/${cluster}.fa | sed 's/>//g' > Split/${cluster}/${cluster}_contigs.txt
    $CDSCRIPTS/AddLengths.pl Annotate/final_contigs_gt1000_c10K.len < Split/${cluster}/${cluster}_contigs.txt > Split/${cluster}/${cluster}_contigs.tsv
    mkdir SplitBam/${cluster}

    for bamfile in Map/*.mapped.sorted.bam
    do
        stub=${bamfile#Map\/}
        stub=${stub%.mapped.sorted.bam}
        
        samtools view -bhL Split/${cluster}/${cluster}_contigs.tsv $bamfile > SplitBam/${cluster}/${stub}_Filter.bam&

    done 
    wait    
done < Concoct/Cluster75.txt 
```
**do not run the above**

and use a third party program bam-readcount to get base frequencies at each position on each contig:

```
module load bam-readcount

while read line
do
    mag=$line

    echo $mag

    cd SplitBam
    cd ${mag}

    cp ../../Split/${mag}/${mag}_contigs.tsv ${mag}_contigs.tsv
    samtools faidx ../../Split/${mag}/${mag}.fa

    echo "${mag}_contigs.tsv"
    mkdir ReadcountFilter
    for bamfile in *_Filter.bam
    do
        stub=${bamfile%_Filter.bam}
            
        echo $stub

        (samtools index $bamfile; bam-readcount -w 1 -q 20 -l "${mag}_contigs.tsv" -f ../../Split/${mag}/${mag}.fa $bamfile 2> ReadcountFilter/${stub}.err > ReadcountFilter/${stub}.cnt)&
    
     done
     cd ..
     cd ..
     wait
    
done < Concoct/Cluster75.txt

``` 
**do not run the above**


Then we can get the base counts on the core cogs:

```

mkdir Variants
while read -r cluster 
do
    echo $cluster
    cd ./SplitBam/${cluster}/ReadcountFilter
    gzip *cnt
    cd ../../..
    python $DESMAN/scripts/ExtractCountFreqGenes.py Split/${cluster}/${cluster}_core.cogs ./SplitBam/${cluster}/ReadcountFilter --output_file Variants/${cluster}_scg.freq > Variants/${cluster}log.txt&

done < Concoct/Cluster75.txt 
``` 
**do not run the above**

Instead just copy the following tar archive into your working directory:
```bash
cp /class/stamps-shared/CDTutorial/Variants.tar.gz $METASIMWD
tar -xvzf Variants.tar.gz
```

The directory contains 8 .freq files one for each cluster. If we look at one:
```bash
head -n 10 Variants/Cluster2_scg.freq 
```
We see it comprises a header, plus one line for each core gene position, giving base frequencies in the order A,C,G,T. This is the input required by DESMAN.

### Detecting variants on core genes

This may be necessary to get DESMAN accessory scripts to work:
```bash
export PATH=/class/stamps-software/desman/lib/python2.7/site-packages/desman-0.1.dev0-py2.7-linux-x86_64.egg/desman:$PATH
```

First we detect variants on both clusters that were identified as 75% pure and complete and had a coverage of greater than 100:

```bash

cd $METASIMWD

mkdir SCG_Analysis

for file in ./Variants/Cluster7_scg.freq ./Variants/Cluster16_scg.freq
do
    echo $file

    stub=${file%.freq}
    stub=${stub#./Variants\/}

    echo $stub

    mkdir SCG_Analysis/$stub
    
    cp $file SCG_Analysis/$stub
    cd SCG_Analysis/$stub    

    Variant_Filter.py ${stub}.freq -o $stub -m 1.0 -f 25.0 -c -sf 0.80 -t 2.5 
    
    cd ../..
done

```

The Variant_Filter.py script produces an output file ${stub}sel_var.csv which lists those positions that are identified as variants by the log-ratio test with FDR < 1.0e-3. We can compare variant frequencies in the two clusters:

```bash
cd SCG_Analysis
wc */*sel_var.csv
```

We can also go into the Cluster 16 directory and look at the output files:
```bash
cd Cluster16_scg
more Cluster16_scgsel_var.csv 
```

The other important file is:
```bash
more Cluster16_scgtran_df.csv
```

This is an estimate of base error rates which is used as a starting point for the haplotype inference.


### Inferring haplotypes

So accounting for the header line we observe 17 and 0 variants in Clusters 16 and 7 respectively. For only Cluster 16 then can we attempt to actually resolve haplotypes. Using the desman executable:

```
cd $METASIMWD/SCG_Analysis/Cluster16_scg

varFile='Cluster16_scgsel_var.csv'

eFile='Cluster16_scgtran_df.csv'
    
for g in 1 2 3 4  
do
    echo $g
    for r in 0 1 2 3 4
    do
	    echo $r
        desman $varFile -e $eFile -o Cluster16_${g}_${r} -g $g -s $r -m 1.0 -i 100 
    done
done
```

```
cat */fit.txt | cut -d"," -f2- > Dev.csv
sed -i '1iH,G,LP,Dev' Dev.csv 
```

Which we can then visualise:
```
$DESMAN/scripts/PlotDev.R -l Dev.csv -o Dev.pdf
```

![Posterior deviance](Figures/Dev.pdf)

There are clearly two haplotypes. We can also run the heuristic to determine haplotype number:

```bash
python $DESMAN/scripts/resolvenhap.py Cluster16
```

This should output:
```
2,2,3,0.0,Cluster16_2_3/Filtered_Tau_star.csv
```

This has the format:
```
No of haplotypes in best fit, No. of good haplotypes in best fit, Index of best fit, Average error rate, File with base predictions
```

Have a look at the prediction file:
```
more Cluster16_2_3/Filtered_Tau_star.csv
```

The position encoding is ACGT so what are the base predictions at each variant position? 
We can turn these into actual sequences with the following commands:

```bash
    cut -d"," -f 1 < Cluster16_scgcogf.csv | sort | uniq | sed '1d' > coregenes.txt

    mkdir SCG_Fasta_2_3
    
    python $DESMAN/scripts/GetVariantsCore.py ../../Annotate/final_contigs_gt1000_c10K.fa ../..//Split/Cluster16/Cluster16_core.cogs Cluster16_2_3/Filtered_Tau_star.csv coregenes.txt -o SCG_Fasta_2_3/
```

This generates one fasta sequence file for each gene with the two strains in:

```bash
ls SCG_Fasta_2_3
```

In this case we know Cluster16 maps onto species 2095, which reassuringly has two strains present. I have pre-calculated the true variants between these two strains.

```bash 
cp  /class/stamps-shared/CDTutorial/DesmanValidation/Cluster16_core_tauRF.csv .

python $DESMAN/scripts/validateSNP2.py Cluster16_2_3/Filtered_Tau_star.csv Cluster16_core_tauRF.csv
```

This gives distance matrices between the true variants and the predictions in terms of SNV and fractions:
```bash
Intersection: 16
[[ 0 16]
 [16  0]]
[[ 0.  1.]
 [ 1.  0.]]
```
Each predicted haplotype should match onto a reference strain with no errors.

### Accessory Gene assignment

The next step would be to calculate the accessory gene presence and absences for these strains. That too will be left to an open lab if there is interest.




```
cd $METASIMWD
cp /class/stamps-shared/CDTutorial/Cluster16.tar.gz .
tar -xvzf Cluster16.tar.gz
cd Cluster16
python $DESMAN/desman/GeneAssign.py ${cluster}_coremean_sd_df.csv $SEL_RUN/Gamma_star.csv ../${cluster}/${cluster}_gene_cov.csv $SEL_RUN/Eta_star.csv -m 20 -v ../${cluster}/${cluster}M0sel_var.csv -
o ${cluster} --assign_tau
```

```