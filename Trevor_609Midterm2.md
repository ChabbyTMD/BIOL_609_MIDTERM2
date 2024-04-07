# BIOL 609 MIDTERM 2

Trevor Mugoya.

## Sequence Import

Import sequence reads. Using gunzipped sequence reads from our current directory, the command creates sample names using the string prefix before the first "_", which in our case is "zr13074"


```
make.file(inputdir=., type=fastq, prefix=midterm609)
```

The next command is performed to address PCR and sequencing errors. It takes the midterm609.files as input and creates contigs to be used downstream from the forward and reverse reads from the zr13074 sample.
```
make.contigs(file=midterm609.files)
```

The next command generates summary statistics on the contig generation process in the previous step. This gives us some information we can use in following steps. For one we can see that we have 284398 sequences with variation between 320 and 468 bps. Given the length of the base pair at the 97.5 percentile, we can approximate a cutoff that we can use to screen for reads that did not assemble well.

```
summary.seqs(fasta=midterm609.trim.contigs.fasta, count=midterm609.contigs.count_table)
```
Here, we are removing ambiguous base calls and any reads that are longer than 490. We obatined this threshold using the 97.5 percentile summary statistic from the summary table in the previous command.

```
screen.seqs(fasta=midterm609.trim.contigs.fasta, count=midterm609.contigs.count_table, maxambig=0, maxlength=490, maxhomop=8)
```

This gives us a summary of the filtering reads longer than 490bp. We have 273560 unique sequences that survived the trimming process.
```
summary.seqs(fasta=current, count=current)
```
The next command will merge duplicates to lessen computational load in the next steps.
```
unique.seqs(fasta=midterm609.trim.contigs.good.fasta, count=midterm609.contigs.good.count_table)
```

Get summary of unique sequences. 219783 unique

```
summary.seqs(count=midterm609.trim.contigs.good.count_table)
```

Create an alignment based on our region of interest(V3V4). Here we are using the silva bacteria database.

> **SIDENOTE**: Figuring out the coordinates for the V3V4 region.(https://mothur.org/blog/2016/Customization-for-your-region/). Downloaded the silav.seed align file from the mothur site.
```
align.seqs(fasta=midterm609.trim.contigs.good.unique.fasta, reference=silva.bacteria/silva.seed_v138_1.align)

summary.seqs(fasta=midterm609.trim.contigs.good.unique.align)
```
From the output of the summary command, I opted to use the median Start(6387) and end(25318) coordinates for input parameters for pcr.seq.

Sanity checking with full reference. Everything seems the same i.e. coordinates.
```
align.seqs(fasta=midterm609.trim.contigs.good.unique.fasta, reference=silva.bacteria/silva.bacteria.fasta)

summary.seqs(fasta=midterm609.trim.contigs.good.unique.align)
```

Subsetting the reference database using coordinates of our region of interest
```
pcr.seqs(fasta=silva.bacteria/silva.bacteria.fasta, start=6387, end=25318, keepdots=F)
```

> **NOTE:** The files silva.bacteria.8mer and silva.bacteria.pcr.fasta were coppied to the current directory prior to the execution of the next steps.

Renaming the reference database to use later on to construct alignment
```
rename.file(input=silva.bacteria.pcr.fasta, new=silva.v3v4.fasta)
```

Checkpoint.
```
summary.seqs(fasta=silva.v3v4.fasta)
```

Constructing alignment using modified reference database.
```
align.seqs(fasta=midterm609.trim.contigs.good.unique.fasta, reference=silva.v3v4.fasta)

summary.seqs(fasta=midterm609.trim.contigs.good.unique.align, count=midterm609.trim.contigs.good.count_table)
```

Extract sequences that align best to our region of interest.
```
screen.seqs(fasta=midterm609.trim.contigs.good.unique.align, count=midterm609.trim.contigs.good.count_table, start=6387, end=25318)
```
Filtering any sequence overhangs of our region of interest.
```
filter.seqs(fasta=midterm609.trim.contigs.good.unique.good.align, vertical=T, trump=.)
```

Identify and merge duplicate sequences in our alignment.
```
unique.seqs(fasta=midterm609.trim.contigs.good.unique.good.filter.fasta, count=midterm609.trim.contigs.good.good.count_table)
```
Split sequences by group and subsequently sort by abundance.
```
pre.cluster(fasta=midterm609.trim.contigs.good.unique.good.filter.unique.fasta, count=midterm609.trim.contigs.good.unique.good.filter.count_table, diffs=2)
```
Identify and remove all chimeric sequences from our data.
```
chimera.vsearch(fasta=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
```

Using Bayesian classifier to classify sequences.
First download and specify paths to reference files from https://mothur.org/wiki/rdp_reference_files/

```
classify.seqs(fasta=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, reference=trainset/trainset19_072023.pds/trainset19_072023.pds.fasta, taxonomy=trainset/trainset19_072023.pds/trainset19_072023.pds.tax)
```

Removal of undesirable eukaryotic, mitochondrial and cholorplast sequences.

```
remove.lineage(fasta=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, taxonomy=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
```

## Error rate assement

This code lets us assess the sequencing error rate using a mock bacterial community as a benchmark. Given that our data does not have a labelled mock community we did not find nor subsequently filter out mock samples.
```
get.groups(count=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, fasta=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta, groups=Mock)
```

This command failed due to the file `midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.fasta` being blank.

```
seq.error(fasta=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.fasta, count=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, reference=HMP_MOCK.v35.fasta, aligned=F)
```

Renaming the fasta, count and taxonomy files.
```
rename.file(fasta=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta, count=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=midterm609.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.pick.taxonomy, prefix=final)
```

## OTU Clustering

Here where calculate the distance and cluster our sequences with the following commands.

Distance Calculation.
```
dist.seqs(fasta=final.fasta, cutoff=0.03)
```

Clustering

```
cluster(column=final.dist, count=final.count_table)
```

```
cluster.split(fasta=final.fasta, count=final.count_table, taxonomy=final.taxonomy, taxlevel=4, cutoff=0.03)
```

```
make.shared(list=final.opti_mcc.list, count=final.count_table, label=0.03)

```

```
make.shared(count=final.count_table)
```

```
classify.otu(list=final.asv.list, count=final.count_table, taxonomy=final.taxonomy, label=ASV)

```

## Phylogenetic methods
```
dist.seqs(fasta=final.fasta, output=lt)

clearcut(phylip=final.phylip.dist)
```

Encountered error: 
> Clearcut: 
Incorrect number of values in the distance matrix. Expected 284924630, and found 1297268539.
Clearcut: Syntax error in distance matrix at offset 61272225246. 

Mothur subsequently exited after this error.


We then calculate how many sequences are present in each sample. Given our data, we only have one sample. Therefore all reads will be assigned to one sample.
```
count.groups(shared=final.opti_mcc.shared)

Result:
Size of smallest group: 223975.

Total seqs: 223975.
```

Sub sampling our data prior to rarefaction. The SOP called for subsampling based on the sample that had the smallest number of sequences. In this case since we only had one sample, we use the value from the `count.groups` command above.
```
sub.sample(shared=final.opti_mcc.shared, size=223975)
```
### Alpha Diversity Analysis

Using the command below we generate rarefaction curves describing the number of OTU's. This output was used to generate rarefaction plot in R.
```
rarefaction.single(shared=final.opti_mcc.shared, calc=sobs, freq=100)
```
Plot generation code adapted from [this repo.](https://github.com/kahoughton/Mothur-commands-MiSeq-16S/blob/master/Creating%20rarefaction%20curve%20from%20mothur%20output)
```R

# Rarefaction output from mothur.
library(tidyverse)

fl <- "final.opti_mcc.groups.rarefaction"


# read data, change name of samples
dat <- read.table(fl, header = T) %%
  gather('var', 'val', -numsampled) %%
  mutate(
    samp = gsub('.*(NDS_GN[0-9]*)$', '\\1', var),
    var = ifelse(grepl('hci', var), 'hi', ifelse(grepl('lci', var), 'lo', 'val'))
    ) %%
  spread(var, val) %%
  arrange(samp, numsampled)

# Rarefaction plot
ggplot(dat, aes(x = numsampled, y = val, group = samp)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = 'lightgrey', fill = 'lightgrey') +
  geom_line(aes(colour = samp)) + 
  theme_bw() + 
  scale_x_continuous('Number of reads sampled') +  
  scale_y_continuous('Number of OTUs observed')
```
![
](https://)

Using our subsampling parameter of 223975 above the command below will subsequently attempt to randomly sample all sequences from our singular sample.
```
summary.single(shared=final.opti_mcc.shared, calc=nseqs-coverage-sobs-invsimpson)
```
The following commands failed to execute.
```
dist.shared(shared=final.opti_mcc.shared, calc=braycurtis-jclass)

Error:
Using 16 processors.
0.03 
[ERROR]: You have not provided enough valid groups.  I cannot run the command.
```

```
pcoa(phylip=final.opti_mcc.braycurtis.0.03.lt.ave.dist)

Error:
Unable to open final.opti_mcc.braycurtis.0.03.lt.ave.dist. Trying mothur's executable directory final.opti_mcc.braycurtis.0.03.lt.ave.dist.
Unable to open final.opti_mcc.braycurtis.0.03.lt.ave.dist.
Unable to open final.opti_mcc.braycurtis.0.03.lt.ave.dist
[ERROR]: did not complete pcoa.
```