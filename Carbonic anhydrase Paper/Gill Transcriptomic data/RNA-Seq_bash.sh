#!/bin/bash

# Define paths
REF_GENOME="ON.fna"
REF_GTF="ON.gtf"
OUTDIR="rnaseq_results"
mkdir -p $OUTDIR

# HISAT2 Index (if not already built)
hisat2-build $REF_GENOME $OUTDIR/genome_index

# Process each sample
for FASTQ in SRR1556364{0..3}.fastq.gz SRR1556364{8,9}.fastq.gz; do
    SAMPLE=$(basename $FASTQ .fastq.gz)
    echo "Processing $SAMPLE..."
    
    # HISAT2 Alignment
    hisat2 \
        -x $OUTDIR/genome_index \
        -U $FASTQ \
        -S $OUTDIR/${SAMPLE}.sam \
        --dta-cufflinks 2> $OUTDIR/${SAMPLE}_align_stats.txt
    
    # SAM to BAM (Samtools)
    samtools sort -@ 4 -o $OUTDIR/${SAMPLE}.bam $OUTDIR/${SAMPLE}.sam
    samtools index $OUTDIR/${SAMPLE}.bam
    
    # StringTie Quantification
    stringtie \
        $OUTDIR/${SAMPLE}.bam \
        -G $REF_GTF \
        -o $OUTDIR/${SAMPLE}.gtf \
        -p 4 \
        -e  # Enable expression estimation
    
    # Cleanup
    rm $OUTDIR/${SAMPLE}.sam
done

# Merge transcripts (optional)
stringtie --merge \
    -G $REF_GTF \
    -o $OUTDIR/merged_transcripts.gtf \
    $OUTDIR/*.gtf

echo "Pipeline complete! Results in $OUTDIR/"