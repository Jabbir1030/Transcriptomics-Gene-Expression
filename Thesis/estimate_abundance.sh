#!/bin/bash

abundancedir=abundance
mapdir=mapped

echo "Estimating abundance"

stringtie -e -B -p 8 -G merged.gtf -o $abundancedir/SRR15563640/SRR15563640.gtf $mapdir/SRR15563640.bam

stringtie -e -B -p 8 -G merged.gtf -o $abundancedir/SRR15563641/SRR15563641.gtf $mapdir/SRR15563641.bam

stringtie -e -B -p 8 -G merged.gtf -o $abundancedir/SRR15563642/SRR15563642.gtf $mapdir/SRR15563642.bam

stringtie -e -B -p 8 -G merged.gtf -o $abundancedir/SRR15563648/SRR15563648.gtf $mapdir/SRR15563648.bam

echo "Job completed"
