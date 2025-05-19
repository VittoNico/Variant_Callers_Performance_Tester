#!/bin/bash

# Read SRR samples from sample.txt
mapfile -t SAMPLES < sample.txt

# Reference file
REF="reference.fasta"
THREADS=12

for SAMPLE in "${SAMPLES[@]}"; do
    echo "🔽 Processing $SAMPLE..."
    mkdir -p "$SAMPLE"
    cd "$SAMPLE" || exit

    # ========== DATA DOWNLOAD ==========
    prefetch "$SAMPLE" --progress
    fasterq-dump --threads "$THREADS" --outdir ./ "$SAMPLE"
    
    # Combine FASTQs
    find . -type f -name "${SAMPLE}*.fastq" -exec cat {} + > "${SAMPLE}_combined.fastq"
    RAW_READS="${SAMPLE}_combined.fastq"
    rm -rf "$SAMPLE" *.sra
    find . -type f -name "${SAMPLE}*.fastq" ! -name "${RAW_READS}" -delete

    # ========== QC ==========
    mkdir -p nanoplot_${SAMPLE}
    NanoPlot --fastq "${RAW_READS}" -o nanoplot_${SAMPLE} --threads $THREADS

    # ========== PROCESSING ==========
    porechop -i "${RAW_READS}" -o "${SAMPLE}_trimmed.fastq" --threads $THREADS
    cat "${SAMPLE}_trimmed.fastq" | NanoFilt -q 7 > "${SAMPLE}_filtered.fastq"
    FILTERED_READS="${SAMPLE}_filtered.fastq"

    # ========== MAPPING ==========
    minimap2 -x map-ont -A 2 -B 8 -O 4,24 -E 2,1 -m 5000 -z 200,100 --secondary=no -s 30 -a ../"$REF" "$FILTERED_READS" -t $THREADS > "${SAMPLE}_mapped.sam"
    
    # Convert to BAM and sort
    samtools view -bS "${SAMPLE}_mapped.sam" > "${SAMPLE}_mapped.bam"
    samtools sort "${SAMPLE}_mapped.bam" -o "${SAMPLE}_mapped.sort.bam"
    samtools index "${SAMPLE}_mapped.sort.bam"
    samtools depth -a "${SAMPLE}_mapped.sort.bam" > "${SAMPLE}_coverage.txt"

    # ========== MD FILE FOR SNIFFLES ==========
    echo "🛠 Creating temporary MD file for Sniffles..."
    samtools calmd -b "${SAMPLE}_mapped.sort.bam" ../"$REF" > "${SAMPLE}_mapped.md.bam"
    samtools index "${SAMPLE}_mapped.md.bam"

    # ========== SV CALLING ==========
    # Sniffles (with MD file)
    echo "⏱️ Starting Sniffles for $SAMPLE"
    /usr/bin/time -v -o "${SAMPLE}_sniffles.time" sniffles -m "${SAMPLE}_mapped.md.bam" -v "${SAMPLE}_sniffles_tmp.vcf" \
        --allow-overwrite --min_support 5 --min_seq_size 50
    bcftools view -i 'SVTYPE="DUP" & QUAL>=20' "${SAMPLE}_sniffles_tmp.vcf" > "${SAMPLE}_sniffles.vcf"
    
    # Remove temporary files
    rm -f "${SAMPLE}_mapped.md.bam" "${SAMPLE}_mapped.md.bam.bai" "${SAMPLE}_sniffles_tmp.vcf"

    # [Resto della pipeline identico...]
    # cuteSV, DeBreak, SVIM come nello script precedente
    
    # ========== COMPARATIVE ANALYSIS ==========
    mkdir -p comparative_analysis_${SAMPLE}
    
    # Create consensus
    truvari collapse -i "${SAMPLE}_sniffles.vcf" "${SAMPLE}_cutesv.vcf" svim_${SAMPLE}/tandem_dups.vcf debreak_${SAMPLE}/variants.vcf \
        -o comparative_analysis_${SAMPLE}/consensus.vcf \
        -f ../"$REF" --passonly --chain -r 1000 -p 0.0 -P 0.7 -O 0.2 -k first

    # Benchmark against consensus
    for TOOL in sniffles cutesv svim debreak; do
        truvari bench -b comparative_analysis_${SAMPLE}/consensus.vcf \
            -c "${SAMPLE}_${TOOL}.vcf" \
            -f ../"$REF" \
            -o comparative_analysis_${SAMPLE}/${TOOL}_vs_consensus \
            --passonly --gtcomp -r 1000 -p 0.0 -P 0.7 -O 0.2
    done

    # ========== CLEANUP ==========
    rm -f "${RAW_READS}" "${SAMPLE}_trimmed.fastq" "${SAMPLE}_mapped."* "${SAMPLE}_"*.time
    cd ..
    echo "✅ Finished $SAMPLE"
done

# Generate reports...
