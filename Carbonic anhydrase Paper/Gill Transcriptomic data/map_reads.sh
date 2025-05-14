#!/bin/bash

fastqdir=resources/samples
mapdir=mapped
mkdir $mapdir


hisat2 -p 8 --dta -x index/ON_tran -1 $fastqdir/brain_1.fastq.gz -2 $fastqdir/brain_2.fastq.gz -S $mapdir/brain.sam

hisat2 -p 8 --dta -x index/ON_tran -1 $fastqdir/heart_1.fastq.gz -2 $fastqdir/heart_2.fastq.gz -S $mapdir/heart.sam

hisat2 -p 8 --dta -x index/ON_tran -1 $fastqdir/kidney_1.fastq.gz -2 $fastqdir/kidney_2.fastq.gz -S $mapdir/kidney.sam

hisat2 -p 8 --dta -x index/ON_tran -1 $fastqdir/liver_1.fastq.gz -2 $fastqdir/liver_2.fastq.gz -S $mapdir/liver.sam


