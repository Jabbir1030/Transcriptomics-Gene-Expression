#!/bin/bash

mapdir=mapped

samtools sort -@ 8 -o $mapdir/SRR15563640.bam $mapdir/SRR15563640.sam
samtools sort -@ 8 -o $mapdir/SRR15563641.bam $mapdir/SRR15563641.sam
samtools sort -@ 8 -o $mapdir/SRR15563642.bam $mapdir/SRR15563642.sam
samtools sort -@ 8 -o $mapdir/SRR15563648.bam $mapdir/SRR15563648.sam
 

