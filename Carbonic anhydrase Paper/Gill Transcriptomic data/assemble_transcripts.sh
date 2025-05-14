#!/bin/bash

gtf=resources/ON.gtf
assembly=assembly
mapdir=mapped
mkdir $assembly

echo "assembling SRR15563640"
stringtie $mapdir/SRR15563640.bam -l SRR15563640 -p 8 -G $gtf -o $assembly/SRR15563640.gtf

echo "assembling SRR15563641"
stringtie $mapdir/SRR15563641.bam -l SRR15563641 -p 8 -G $gtf -o $assembly/SRR15563641.gtf

echo "assembling SRR15563642"
stringtie $mapdir/SRR15563642.bam -l SRR15563642 -p 8 -G $gtf -o $assembly/SRR15563642.gtf

echo "assembling SRR15563648"
stringtie $mapdir/SRR15563648.bam -l SRR15563648 -p 8 -G $gtf -o $assembly/SRR15563648.gtf

