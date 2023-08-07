---
layout: post
title:  "RNA-seq using Salmon (under development)"
date:   2023-07-26 15:24:42 -0400
---
# Preamble
<img src="{{site.baseurl}}/assets/img/RNAseq_flowcell_pixelart.png">
This tutorial will cover a boilerplate RNAseq analysis using Salmon psuedo-alignment of illumina generated paired-end reads to the *Arabidopsis thaliana* transcriptome.  

Reads are from a publically available dataset for testing the [transcriptomic response of *Arabidopsis thaliana* to dopamine](https://doi.org/10.3390/stresses3010026). 

We will cover:  
- Data cleaning (shell and other tools)
- Alignment of reads (Salmon)
- Statistical analysis of alignment data (R and edgeR)
- Annotation, heatmap production and cluster analysis
- Gene Set Enrichment Analysis (R)

The data cleaning step of this tutorial is broadly applicable to any workflow and errs on the conservative side with multiple quality checks. Once the data is cleaned, there are many ways to proceed with alignment. here we use Salmon, an efficient psuedo-alignment program.  

## Setting up your environment
This tutorial uses open-source bioinformatics tools and R.

### Installing the tools
We will install our tools using conda. I have included a conda environment yml, to get things up and running quickly. Skip this step if these programs are already available to you or if you wish to install them nativly.

packages include: 

- [FastQC](https://github.com/s-andrews/FastQC)
- [Bowtie2](https://github.com/BenLangmead/bowtie2.git)
- [Salmon](https://github.com/COMBINE-lab/salmon.git)
- [rCorrector](https://github.com/mourisl/Rcorrector.git)
- [TrimGalore](https://github.com/FelixKrueger/TrimGalore.git)
- [SRA-Tools](https://github.com/ncbi/sra-tools.git)

### R
If you don't have R and R-studio installed, get it installed. If you are in a hurry, I have included a conda environment configureation which will automatically deploy R-studio with all the R packages needed for this tutorial preconfigured.

### one last thing to download...
Finally, there is a python script from the Harvard Bioinformatics team which is useful here. Go to the link below, download the python script and place it in the working directory you have chosen to use for this tutorial. 

[FilterUncorrectabledPEfastq.py](https://github.com//TranscriptomeAssemblyTools/blob/dfe258636088c11eb60d3ce69da2fd5cd00ef5b3/FilterUncorrectabledPEfastq.py)

 <img src="{{site.baseurl}}/assets/img/HBTdwnl.jpeg">

### Variables {cName}, {tName} and {threads} 
If you wish to make all of this code run without editing it, you can simply set the $cName, $tName and $threads variables in your bash session. These variables are for your control read names, treatment read names and the number of threads you want to use for the analysis. Here, we will simply cName our data cDA2 and tName tDA2. The thread number is set to 20 here, for use on 24 core machines to give you some headroom for using the computer while the analysis is running.  

{% highlight ruby %}
cName='cDA2'
tName='tDA2'
threads='20'
{% endhighlight %}

# Quantification of reads using Salmon
<img src="{{site.baseurl}}/assets/img/salmon.png">
## Get data (or use your own)
First we need to download our data. For this we will use Sequence Reach Archive Tools.

{% highlight ruby %}

{% endhighlight %}

## Clean Data
Assess the quality of the raw data using FASTQC ([How to read and interperate FastQC reports](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html))

mkdir qcRaw

{% highlight ruby %}
fastqc --threads $threads --outdir ./qcRaw ./${cName}_1.fq.gz
fastqc --threads $threads --outdir ./qcRaw ./${cName}_2.fq.gz
fastqc --threads $threads --outdir ./qcRaw ./${tName}_1.fq.gz
fastqc --threads $threads --outdir ./qcRaw ./${tName}_2.fq.gz
{% endhighlight %}

### Run rCorrector
rCorrector will repair read pairs which are low quality.
{% highlight ruby %}
mkdir ./rCorr
{% endhighlight %}

{% highlight ruby %}
rcorrector \
    -ek 20000000000 \
    -t $threads \
    -od ./rCorr \
    -1 ../${cName}_1.fq \
    -2 ../${cName}_2.fq

rcorrector \
    -ek 20000000000 \
    -t $threads \
    -od ./rCorr \
    -1 ../${tName}_1.fq \
    -2 ../${tName}_2.fq
{% endhighlight %}

{% highlight ruby %}
cd rCorr
python ../FilterUncorrectabledPEfastq.py -1 ${cName}_1.cor.fq -2 ${cName}_2.cor.fq -s ${cName}
python ../FilterUncorrectabledPEfastq.py -1 ${tName}_1.cor.fq -2 ${tName}_2.cor.fq -s ${tName}
{% endhighlight %}


### Remove unfixable reads
{% highlight ruby %}
trim_galore \
    -j 8 \
    --paired \
    --retain_unpaired \
    --phred33 \
    --output_dir ./trimmed_reads \
    --length 36 \
    -q 5 \
    --stringency 1 \
    -e 0.1
    unfixrm_${cName}_1.cor.fq unfixrm_${cName}_2.cor.fq

trim_galore \
    -j 8 \
    --paired \
    --retain_unpaired \
    --phred33 \
    --output_dir ./trimmed_reads \
    --length 36 \
    -q 5 \
    --stringency 1 \
     -e 0.1 \
     unfixrm_${tName}_1.cor.fq unfixrm_${tName}_2.cor.fq
{% endhighlight %}


### Remove reads originating from ribsomal RNA
Most modern RNAseq library prep methods are very good at minimizing rRNA contamination, but there is always *some*.  

To correct for this, we are going to map our reads to the SILVA rRNA blacklist, and output those reads which don't align. This will remove rRNA contamination. 

First, lets create a directory to keep things organized
{% highlight ruby %}
mkdir riboMap   
cd riboMap  
{% endhighlight %}

Download the SILVA rRNA blacklists for the Large and Small subunits, which we will subsequently use to clean our reads of rRNA contamination.
{% highlight ruby %}
wget "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz"
wget "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
{% endhighlight %}

Unzip the files, concatonate them and replace "U" in the fasta sequences with "T"

{% highlight ruby %}
gunzip ./*
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta > SILVArRNAdb.fa.temp
cat SILVA_138.1_SSURef_NR99_tax_silva.fasta > SILVArRNAdb.fa.temp
sed '/^[^>]/s/U/T/g' SILVArRNAdb.fa.temp > SILVAcDNAdb.fa
rm SILVArRNAdb.fa.temp
{% endhighlight %}

Align the reads using Bowie2 to the rRNA blacklist. Those reads which *do not* align to the blacklist will be output as clean_{cName}_1.fq.gz and clean_{cName}_2.fq.gz

{% highlight ruby %}
bowtie2 --quiet --very-sensitive-local --phred33  \
    -x SILVAcDNAdb \
    -1 ../rCorr/trimmed_reads/unfixrm_${cName}_1.cor_val_1.fq \
    -2 ../rCorr/trimmed_reads/unfixrm_${cName}_2.cor_val_2.fq \
    --threads $threads \
    --met-file ${cName}_bowtie2_metrics.txt \
    --al-conc-gz blacklist_paired_aligned_${cName}.fq.gz \
    --un-conc-gz clean_${cName}.fq.gz  \
    --al-gz blacklist_unpaired_aligned_${cName}.fq.gz \
    --un-gz blacklist_unpaired_unaligned_${cName}.fq.gz

bowtie2 --quiet --very-sensitive-local --phred33  \
    -x SILVAcDNAdb \
    -1 ../rCorr/trimmed_reads/unfixrm_${tName}_1.cor_val_1.fq \
    -2 ../rCorr/trimmed_reads/unfixrm_${tName}_2.cor_val_2.fq \
    --threads $threads \
    --met-file ${tName}_bowtie2_metrics.txt \
    --al-conc-gz blacklist_paired_aligned_${tName}.fq.gz \
    --un-conc-gz clean_${tName}.fq.gz  \
    --al-gz blacklist_unpaired_aligned_${tName}.fq.gz \
    --un-gz blacklist_unpaired_unaligned_${tName}.fq.gz
{% endhighlight %}

Navigate to your primary working directory, and create a directory for your FastQC outputs of your cleaned data
{% highlight ruby %}
cd ..
mkdir ./qcClean
cd ./qcClean
{% endhighlight %}

Reasses the quality of your data using FastQC. 
{% highlight ruby %}
fastqc --threads $threads --outdir ./ ../riboMap/clean_${cName}_1.fq
fastqc --threads $threads --outdir ./ ../riboMap/clean_${cName}_2.fq
fastqc --threads $threads --outdir ./ ../riboMap/clean_${tName}_1.fq
fastqc --threads $threads --outdir ./ ../riboMap/clean_${tName}_2.fq
{% endhighlight %}

At this point, your data should be in good shape. 

## Generate count data using Salmon

### Make *Arabidopsis thaliana* index
First, we need to download the *Arabidopsis thaliana* reference transcriptome and index it for Salmon. There are a number of options for indexing in Salmon ([docs](https://salmon.readthedocs.io/en/latest/index.html)). You may notice the documentation recommends running in *decoy-aware* mode. This can be safely ignored in Arabidopsis, as the reference transcripome is very good. For organisms with less robust reference transcriptomes decoy awareness can avoid spurious mapping of your reads ([to, for example transcribed psuedogenes](https://www.biostars.org/p/456231/")).   

{% highlight ruby %}
wget -O athal.fa.gz "ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz"
salmon index -t athal.fa.gz -i athal_index
{% endhighlight %}

### Run Salmon
Soft-link your cleaned FQ files. Not really necissary, but makes the run command tidier. 
{% highlight ruby %}
mkdir salmon_out
cd salmon_out
ln -s ../riboMap/clean_${cName}_1.fq ./clean_${cName}_1.fq
ln -s ../riboMap/clean_${cName}_2.fq ./clean_${cName}_2.fq
ln -s ../riboMap/clean_${cName}_1.fq ./clean_${tName}_1.fq
ln -s ../riboMap/clean_${cName}_2.fq ./clean_${tName}_2.fq
{% endhighlight %}

{% highlight ruby %}
salmon quant \
    -i ../athal_index \
    -l A -1 clean_${cName}_1.fq \
    -2 clean_${cName}_2.fq \
    --gcBias \
    -p 20 \
    --validateMappings \
    -o ${cName}_transQ

salmon quant \
    -i ../athal_index \
    -l A \
    -1 clean_${tName}_1.fq \
    -2 clean_${tName}_2.fq \
    --gcBias \
    -p 20 \
    --validateMappings \
    -o ${tName}_transQ
{% endhighlight %}


# Statistical Analysis of count data using R
<img src="{{site.baseurl}}/assets/img/Rstudio.png">

First we need to download the gene features file for the Arabidopsis transcriptome. This will the production of gene-level abundence estemates using [tximport](https://f1000research.com/articles/4-1521). This step dramatically simplifies analysis and the interpritation of results. 

{% highlight ruby %}
wget "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"
{% endhighlight %}

Now it's time to leave the bash environment and fire up R studio. 

## install needed libraries and load them up
{% highlight ruby %}
setwd("./")
library('edgeR')
library('RColorBrewer')
library("GenomicFeatures")
library("tximport")
library("readr")
library("tximport")
library("edgeR")
{% endhighlight %}
## make TxDB for TAIR

{% highlight ruby %}
gtffile <- "./TAIR10_GFF3_genes.gff"
file.exists(gtffile)
txdb <- makeTxDbFromGFF(gtffile, format = "gff3", circ_seqs = character())
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
write.csv(tx2gene, file = "./TAIR10tx2gene.gencode.v27.csv")
head(tx2gene)
{% endhighlight %}

## add variables to avoid having to change the file names
{% highlight ruby %}
TP = "2H"
{% endhighlight %}

## Import data using tximport
{% highlight ruby %}
#Read in tx2gene database
tx2gene <- read_csv("./TAIR10tx2gene.gencode.v27.csv")
{% endhighlight %}

## Read in sample list (file names and factors)
In order to set up the analysis, we need to create a file containing the file names and factors. 
{% highlight ruby %}
samples <- read.table(file = paste0(TP, "/samples.txt"), header = T)
head(samples)
{% endhighlight %}

## Retrieve file paths for salmon .sf files
{% highlight ruby %}
files <- file.path(TP, paste0(samples$samples,"_quant.sf") )
names(files) <- paste0("sample", 1:6)
all(file.exists(files))
head(files)
{% endhighlight %}

## Import salmon .sf files, and summarize to gene using tx2gene database
{% highlight ruby %}
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
{% endhighlight %}

## Prep Data for EdgeR 
{% highlight ruby %}
cts <- txi$counts
head(cts)
seqDataGroups <- c(paste0(TP,"_0", TP, "_0", TP, "_0", TP,"_50", TP, "_50", TP, "_50"))
seqDataGroups
{% endhighlight %}

# use edgeR function DGEList to make read count list
{% highlight ruby %}
d <- DGEList(counts=cts,group=factor(seqDataGroups))
head(d)
{% endhighlight %}
## Data Filtering
{% highlight ruby %}
dim(d)
d.full <- d # keep in case things go awry
head(d$counts)
head(cpm(d))

#total counts/sample
apply(d$counts,2,sum)
{% endhighlight %}

## Trim lowly-abundant genes. cpm number can be changed based on library size, and detected tags.
{% highlight ruby %}
keep <- rowSums(cpm(d)>100) >= 2
keep
d <- d[keep,]
dim (d)

d$samples$lib.size <- colSums(d$counts)
d$samples

d <- calcNormFactors(d, logratioTrim = 0)
d

plotMDS(d, method="bcv", col=as.numeric(d$samples$group))

d1 <- estimateCommonDisp(d, verbose=T)
names(d1)

d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")

# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
topTags(et12,n=100)


#de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
#for GSEA, uncomment above for standard analysis

de1 <- decideTestsDGE(et12, adjust.method="BH")
summary(de1)

# differentially expressed tags from the naive method in d1
de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
design.mat
fit <- glmFit(d2, design.mat)
summary(de1)

DE <- topTags(et12,n=50000)
as.data.frame(DE)
write.csv(DE, paste0("./edgeRout_", TP, "_LRT0.csv"))
{% endhighlight %}