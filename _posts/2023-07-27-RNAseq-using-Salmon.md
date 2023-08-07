---
layout: post
title:  "RNA-seq using Salmon (under development)"
date:   2023-07-26 15:24:42 -0400
---



This tutorial will cover a boilerplate RNAseq analysis using Salmon psuedo-alignment of illumina generated paired-end reads to the *Arabidopsis thaliana* transcriptome.  

Reads are from a publically available dataset for testing the [transcriptomic response of *Arabidopsis thaliana* to dopamine](https://doi.org/10.3390/stresses3010026). 

We will cover:  
- Data cleaning (shell and other tools)
- Alignment of reads (Salmon)
- Statistical analysis of alignment data (R and DeSeq2)
- Annotation, heatmap production and cluster analysis
- Gene Set Enrichment Analysis (R)


## Preamble
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

### one last thing...
Finally, there is a python script from the Harvard Bioinformatics team which is useful here. Go to the link below, download the python script and place it in the working directory you have chosen to use for this tutorial. 

[FilterUncorrectabledPEfastq.py](https://github.com//TranscriptomeAssemblyTools/blob/dfe258636088c11eb60d3ce69da2fd5cd00ef5b3/FilterUncorrectabledPEfastq.py)

 <img src="{{site.baseurl}}/assets/img/HBTdwnl.jpeg">

## Data cleaning

### Running rCorrector
rCorrector will repair read pairs which are low quality. The unfixable read pairs will also be removed during this process. 
{% highlight ruby %}
rcorrector -ek 20000000000 -t 24 -od ./ -1 ${2}_1.fq -2 ${2}_2.fq
{% endhighlight %}

{% highlight ruby %}
python ./FilterUncorrectabledPEfastq.py -1 ${2}_1.cor.fq -2 ${2}_2.cor.fq -s $2
{% endhighlight %}


{% highlight ruby %}
mkdir ./trimmed_reads
{% endhighlight %}

{% highlight ruby %}
trim_galore -j 8 --paired --retain_unpaired --phred33 --output_dir ./trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 unfixrm_${2}_1.cor.fq unfixrm_${2}_2.cor.fq
{% endhighlight %}

{% highlight ruby %}
mkdir riboMap   
cd riboMap  
{% endhighlight %}

### Remove reads originating from rRNA
Most modern RNAseq library prep methods are very good at minimizing rRNA contamination, but there is always *some*.  

To correct for this, we are going to map our reads to the SILVA rRNA blacklist, and output those reads which don't align. This will remove rRNA contamination. 

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


{% highlight ruby %}
bowtie2 --quiet --very-sensitive-local --phred33  -x SILVAcDNAdb -1 ../rCorr/trimmed_reads/unfixrm_${2}_1.cor_val_1.fq -2 ../rCorr/trimmed_reads/unfixrm_${2}_2.cor_val_2.fq --threads 22 --met-file ${2}_bowtie2_metrics.txt --al-conc-gz blacklist_paired_aligned_${2}.fq.gz --un-conc-gz clean_blacklist_paired_unaligned_${2}.fq.gz  --al-gz blacklist_unpaired_aligned_${2}.fq.gz --un-gz blacklist_unpaired_unaligned_${2}.fq.gz
{% endhighlight %}

### Run fastQC

Create the directory structure for your outputs
{% highlight ruby %}
mkdir ./qcClean
cd ./qcClean
{% endhighlight %}

soft-link your data, which is now cleaned from rRNA contamination
{% highlight ruby %}
ln -s ../riboMap/clean_${2}_1.fq ./clean_${2}_1.fq
ln -s ../riboMap/clean_${2}_2.fq ./clean_${2}_2.fq
{% endhighlight %}

Reasses the quality of your data. 
{% highlight ruby %}
fastqc --threads 24 --outdir ./ ./clean_$2_1.fq
fastqc --threads 24 --outdir ./ ./clean_$2_2.fq
{% endhighlight %}

{% highlight ruby %}
rm clean_${2}_1.fq
rm clean_${2}_2.fq
{% endhighlight %}

{% highlight ruby %}
mkdir ./rcorrected
cd ./rcorrected
{% endhighlight %}

{% highlight ruby %}
#ln -s ../../rCorr/${2}_1.cor.fq ./${2}_1.cor.fq
#ln -s ../../rCorr/${2}_2.cor.fq ./${2}_2.cor.fq
{% endhighlight %}

{% highlight ruby %}
#fastqc --threads 24 --outdir ./ ./$2_1.cor.fq
#fastqc --threads 24 --outdir ./ ./$2_2.cor.fq
{% endhighlight %}

### Make *Arabidopsis thaliana* index
First, we need to download the *Arabidopsis thaliana* reference transcriptome and index it for Salmon. There are a number of options for indexing in Salmon ([docs](https://salmon.readthedocs.io/en/latest/index.html)). You may notice the documentation recommends running in *decoy-aware* mode. This can be safely ignored in Arabidopsis, as the reference transcripome is very good. For organisms with less robust reference transcriptomes decoy awareness can avoid spurious mapping of your reads ([to, for example transcribed psuedogenes](https://www.biostars.org/p/456231/").   

{% highlight ruby %}
wget -O athal.fa.gz "ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz"
salmon index -t athal.fa.gz -i athal_index
{% endhighlight %}

### Run Salmon
Soft-link your cleaned FQ files. Not really necissary, but makes the run command tidier. 
{% highlight ruby %}
mkdir salmon_out
cd salmon_out
ln -s ../riboMap/clean_${2}_1.fq ./clean_${2}_1.fq
ln -s ../riboMap/clean_${2}_2.fq ./clean_${2}_2.fq
{% endhighlight %}

{% highlight ruby %}
salmon quant -i ../athal_index -l A -1 clean_${2}_1.fq -2 clean_${2}_2.fq --gcBias -p 20 --validateMappings -o ${2}_transQ
{% endhighlight %}

 <img src="{{site.baseurl}}/assets/img/2H.png">

{% highlight ruby %}
{% endhighlight %}