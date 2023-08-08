#!/usr/bin env bash

#declare variables. 
cName='cDA2'
tName='tDA2'
threads='20'

#Set this variable as an array of your replicate measurements. For nearly every command in the following section, we will be using a for loop to run each command on all of the replicates. In our example, we will have 3.
declare -a reps=("R1" "R2" "R3")


# Data Acquisition & Cleaning

## Get data (or use your own)
wget "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"

fasterq-dump --fasta --split-files SRX18321587
mv SRX18321587_1* cDA2_R1_1.fq
mv SRX18321587_2* cDA2_R1_2.fq
fasterq-dump --fasta --split-files SRX18321588
mv SRX18321588_1* cDA2_R2_1.fq
mv SRX18321588_2* cDA2_R2_2.fq
fasterq-dump --fasta --split-files SRX18321591
mv SRX18321591_1* cDA2_R3_1.fq
mv SRX18321591_2* cDA2_R3_2.fq
fasterq-dump --fasta --split-files SRX18321595
mv SRX18321595_1* tDA2_R1_1.fq
mv SRX18321595_2* tDA2_R1_2.fq
fasterq-dump --fasta --split-files SRX18321596
mv SRX18321596_1* tDA2_R2_1.fq
mv SRX18321596_2* tDA2_R2_2.fq
fasterq-dump --fasta --split-files SRX18321597
mv SRX18321597_1* tDA2_R3_1.fq
mv SRX18321597_2* tDA2_R3_2.fq
gzip *.fq

## Clean Data

mkdir qcRaw

for i in "${reps[@]}"
do
    fastqc --threads $threads --outdir ./qcRaw ./${cName}_{i}_1.fq.gz
    fastqc --threads $threads --outdir ./qcRaw ./${cName}_{i}_2.fq.gz
    fastqc --threads $threads --outdir ./qcRaw ./${tName}_{i}_1.fq.gz
    fastqc --threads $threads --outdir ./qcRaw ./${tName}_{i}_2.fq.gz
done

## Run rCorrector

mkdir ./rCorr

for i in "${reps[@]}"
do
    rcorrector \
        -ek 20000000000 \
        -t $threads \
        -od ./rCorr \
        -1 ../${cName}_{i}_1.fq \
        -2 ../${cName}_{i}_2.fq
    rcorrector \
        -ek 20000000000 \
        -t $threads \
        -od ./rCorr \
        -1 ../${tName}_{i}_1.fq \
        -2 ../${tName}_{i}_2.fq
done

cd rCorr

for i in "${reps[@]}"
do
    python ../FilterUncorrectabledPEfastq.py -1 ${cName}_{i}_1.cor.fq -2 ${cName}_{i}_2.cor.fq -s ${cName}_{i}
    python ../FilterUncorrectabledPEfastq.py -1 ${tName}_{i}_1.cor.fq -2 ${tName}_{i}_2.cor.fq -s ${tName}_{i}
done

## Remove unfixable reads

for i in "${reps[@]}"
do
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
        unfixrm_${cName}_{i}_1.cor.fq unfixrm_${cName}_{i}_2.cor.fq
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
         unfixrm_${tName}_{i}_1.cor.fq unfixrm_${tName}_{i}_2.cor.fq
done

## Remove reads originating from ribsomal RNA
mkdir riboMap   
cd riboMap  

wget "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz"
wget "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"

gunzip ./*
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta > SILVArRNAdb.fa.temp
cat SILVA_138.1_SSURef_NR99_tax_silva.fasta > SILVArRNAdb.fa.temp
sed '/^[^>]/s/U/T/g' SILVArRNAdb.fa.temp > SILVAcDNAdb.fa
rm SILVArRNAdb.fa.temp

for i in "${reps[@]}"
do
    bowtie2 --quiet --very-sensitive-local --phred33  \
        -x SILVAcDNAdb \
        -1 ../rCorr/trimmed_reads/unfixrm_${cName}_{i}_1.cor_val_1.fq \
        -2 ../rCorr/trimmed_reads/unfixrm_${cName}_{i}_2.cor_val_2.fq \
        --threads $threads \
        --met-file ${cName}_{i}_bowtie2_metrics.txt \
        --al-conc-gz blacklist_paired_aligned_${cName}_{i}.fq.gz \
        --un-conc-gz clean_${cName}_{i}.fq.gz  \
        --al-gz blacklist_unpaired_aligned_${cName}_{i}.fq.gz \
        --un-gz blacklist_unpaired_unaligned_${cName}_{i}.fq.gz
    bowtie2 --quiet --very-sensitive-local --phred33  \
        -x SILVAcDNAdb \
        -1 ../rCorr/trimmed_reads/unfixrm_${tName}_{i}_1.cor_val_1.fq \
        -2 ../rCorr/trimmed_reads/unfixrm_${tName}_{i}_2.cor_val_2.fq \
        --threads $threads \
        --met-file ${tName}_{i}_bowtie2_metrics.txt \
        --al-conc-gz blacklist_paired_aligned_${tName}_{i}.fq.gz \
        --un-conc-gz clean_${tName}_{i}.fq.gz  \
        --al-gz blacklist_unpaired_aligned_${tName}_{i}.fq.gz \
        --un-gz blacklist_unpaired_unaligned_${tName}_{i}.fq.gz
done

cd ..
mkdir ./qcClean
cd ./qcClean

## Reasses the quality of your data using FastQC. 

for i in "${reps[@]}"
do
    fastqc --threads $threads --outdir ./ ../riboMap/clean_${cName}_{i}_1.fq
    fastqc --threads $threads --outdir ./ ../riboMap/clean_${cName}_{i}_2.fq
    fastqc --threads $threads --outdir ./ ../riboMap/clean_${tName}_{i}_1.fq
    fastqc --threads $threads --outdir ./ ../riboMap/clean_${tName}_{i}_2.fq
done

# Generate Count Data Using Salmon

wget -O athal.fa.gz "ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz"
salmon index -t athal.fa.gz -i athal_index

mkdir salmon_out
cd salmon_out
for i in "${reps[@]}"
do
    ln -s ../riboMap/clean_${cName}_{i}_1.fq ./clean_${cName}_{i}_1.fq
    ln -s ../riboMap/clean_${cName}_{i}_2.fq ./clean_${cName}_{i}_2.fq
    ln -s ../riboMap/clean_${cName}_{i}_1.fq ./clean_${tName}_{i}_1.fq
    ln -s ../riboMap/clean_${cName}_{i}_2.fq ./clean_${tName}_{i}_2.fq
    salmon quant \
        -i ../athal_index \
        -l A -1 clean_${cName}_{i}_1.fq \
        -2 clean_${cName}_{i}_2.fq \
        --gcBias \
        -p 20 \
        --validateMappings \
        -o ${cName}_{i}_transQ
    salmon quant \
        -i ../athal_index \
        -l A \
        -1 clean_${tName}_{i}_1.fq \
        -2 clean_${tName}_{i}_2.fq \
        --gcBias \
        -p 20 \
        --validateMappings \
        -o ${tName}_{i}_transQ
done

