#!/bin/bash

mapdir=mapped

samtools sort -@ 8 -o $mapdir/brain.bam $mapdir/brain.sam
samtools sort -@ 8 -o $mapdir/heart.bam $mapdir/heart.sam
samtools sort -@ 8 -o $mapdir/kidney.bam $mapdir/kidney.sam
samtools sort -@ 8 -o $mapdir/liver.bam $mapdir/liver.sam
 

