#!/bin/bash

gtf=resources/ON.gtf
assembly=assembly
mapdir=mapped
mkdir $assembly

echo "assembling brain"
stringtie $mapdir/brain.bam -l brain -p 8 -G $gtf -o $assembly/brain.gtf

echo "assembling heart"
stringtie $mapdir/heart.bam -l heart -p 8 -G $gtf -o $assembly/heart.gtf

echo "assembling kidney"
stringtie $mapdir/kidney.bam -l kidney -p 8 -G $gtf -o $assembly/kidney.gtf

echo "assembling liver"
stringtie $mapdir/liver.bam -l liver -p 8 -G $gtf -o $assembly/liver.gtf

