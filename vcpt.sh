#!/bin/bash

# Legge i sample SRR da un file chiamato sample.txt
mapfile -t SAMPLES < sample.txt

# Reference file
REF="reference.fasta"
THREADS=24
TRUTH="../truth_set.vcf.gz"  # ← path al tuo truth set

# Inizializza Conda
source "$(conda info --base)/etc/profile.d/conda.sh"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "🔽 Processing $SAMPLE..."

    mkdir -p "$SAMPLE"
    cd "$SAMPLE" || exit

    # Scarica la run SRA (Nanopore)
    prefetch "$SAMPLE" --progress
    fasterq-dump --threads "$THREADS" --outdir ./ "$SAMPLE"

    # Unisci tutti i fastq trovati nella cartella
    find . -type f -name "${SAMPLE}*.fastq" -exec cat {} + > "${SAMPLE}_combined.fastq"
    RAW_READS="${SAMPLE}_combined.fastq"

    # Rimuovi tutto tranne il fastq combinato
    rm -rf "$SAMPLE"
    rm -rf *.sra
    find . -type f -name "${SAMPLE}*.fastq" ! -name "${RAW_READS}" -delete

    # QC iniziale con NanoPlot
    mkdir -p nanoplot_${SAMPLE}
    NanoPlot --fastq "${RAW_READS}" -o nanoplot_${SAMPLE} --threads $THREADS

    # Trimming con Porechop
    conda activate sv-env
    porechop -i "${RAW_READS}" -o "${SAMPLE}_trimmed.fastq" --threads $THREADS

    # Filtro di qualità con NanoFilt (es. Q≥10)
    cat "${SAMPLE}_trimmed.fastq" | NanoFilt -q 10 > "${SAMPLE}_filtered.fastq"
    FILTERED_READS="${SAMPLE}_filtered.fastq"

    # Mapping
    minimap2 -x map-ont -A 2 -B 8 -O 4,24 -E 2,1 -m 5000 -z 200,100 --secondary=no -s 30 -a ../"$REF" "$FILTERED_READS" -t $THREADS > "${SAMPLE}_mapped.sam"

    # BAM processing
    samtools view -bS "${SAMPLE}_mapped.sam" > "${SAMPLE}_mapped.bam"
    samtools sort "${SAMPLE}_mapped.bam" -o "${SAMPLE}_mapped.sort.bam"
    samtools index "${SAMPLE}_mapped.sort.bam"
    samtools depth -a "${SAMPLE}_mapped.sort.bam" > "${SAMPLE}_coverage.txt"

    # SV Calling: Sniffles (filter only DUP) with timing
    echo "⏱️ Starting Sniffles for $SAMPLE"
    /usr/bin/time -v -o "${SAMPLE}_sniffles.time" sniffles --input "${SAMPLE}_mapped.sort.bam" --vcf "${SAMPLE}_sniffles_tmp.vcf" --allow-overwrite
    bcftools view -i 'SVTYPE="DUP"' "${SAMPLE}_sniffles_tmp.vcf" > "${SAMPLE}_sniffles.vcf"
    rm -f "${SAMPLE}_sniffles_tmp.vcf"

    # SV Calling: cuteSV (filter only DUP) with timing
    echo "⏱️ Starting cuteSV for $SAMPLE"
    /usr/bin/time -v -o "${SAMPLE}_cutesv.time" cuteSV "${SAMPLE}_mapped.sort.bam" ../"$REF" "${SAMPLE}_cutesv_tmp.vcf" cute --max_size 1100000 --min_support 5
    bcftools view -i 'SVTYPE="DUP"' "${SAMPLE}_cutesv_tmp.vcf" > "${SAMPLE}_cutesv.vcf"
    rm -f "${SAMPLE}_cutesv_tmp.vcf"

    # DeBreak + filter DUP with timing
    echo "⏱️ Starting DeBreak for $SAMPLE"
    conda activate debreak
    mkdir -p debreak_${SAMPLE}
    /usr/bin/time -v -o "${SAMPLE}_debreak.time" debreak --bam "${SAMPLE}_mapped.sort.bam" --outpath debreak_${SAMPLE}
    bcftools view -i 'SVTYPE="DUP"' debreak_${SAMPLE}/variants.vcf > debreak_${SAMPLE}/dup.vcf
    mv debreak_${SAMPLE}/dup.vcf debreak_${SAMPLE}/variants.vcf
    conda deactivate

    # SVIM - call only tandem duplications with timing
    echo "⏱️ Starting SVIM for $SAMPLE"
    conda activate svim
    mkdir -p svim_${SAMPLE}
    /usr/bin/time -v -o "${SAMPLE}_svim.time" svim alignment svim_${SAMPLE} "${SAMPLE}_mapped.sort.bam" --reference ../"$REF" \
        --tandem-duplications --skip-insertions --skip-deletions --skip-inversions --skip-translocations
    conda deactivate

    # Truvari collapse (now including DeBreak VCF)
    conda activate sv-env
    mkdir -p truvari_${SAMPLE}
    truvari collapse \
        -i "${SAMPLE}_sniffles.vcf" "${SAMPLE}_cutesv.vcf" svim_${SAMPLE}/variants.vcf debreak_${SAMPLE}/variants.vcf \
        -o truvari_${SAMPLE}/${SAMPLE}_merged.vcf \
        -f ../"$REF" \
        --passonly --chain -r 1000 -p 0.0 -P 0.7 -O 0.2 -k first

    # Truvari benchmarking
    mkdir -p truvari_bench_${SAMPLE}
    truvari bench \
        -b "$TRUTH" \
        -c truvari_${SAMPLE}/${SAMPLE}_merged.vcf \
        -f ../"$REF" \
        -o truvari_bench_${SAMPLE} \
        --passonly --gtcomp -r 1000 -p 0.0 -P 0.7 -O 0.2

    # Generate performance report
    echo "📊 Generating performance report for $SAMPLE"
    {
        echo "Caller,Elapsed Time (s),CPU Percentage,Max Memory (KB),User Time (s),System Time (s)"
        grep -E 'Elapsed|CPU|Maximum' "${SAMPLE}_sniffles.time" | awk -v caller="Sniffles" 'BEGIN{OFS=","} \
            /Elapsed \(wall clock\)/{split($8,a,":"); time=a[1]*60 + a[2]} \
            /CPU/{cpu=$NF} \
            /Maximum resident set size/{mem=$NF} \
            /User time/{utime=$4} \
            /System time/{stime=$4} \
            END{print caller,time,cpu,mem,utime,stime}'
        grep -E 'Elapsed|CPU|Maximum' "${SAMPLE}_cutesv.time" | awk -v caller="cuteSV" 'BEGIN{OFS=","} \
            /Elapsed \(wall clock\)/{split($8,a,":"); time=a[1]*60 + a[2]} \
            /CPU/{cpu=$NF} \
            /Maximum resident set size/{mem=$NF} \
            /User time/{utime=$4} \
            /System time/{stime=$4} \
            END{print caller,time,cpu,mem,utime,stime}'
        grep -E 'Elapsed|CPU|Maximum' "${SAMPLE}_debreak.time" | awk -v caller="DeBreak" 'BEGIN{OFS=","} \
            /Elapsed \(wall clock\)/{split($8,a,":"); time=a[1]*60 + a[2]} \
            /CPU/{cpu=$NF} \
            /Maximum resident set size/{mem=$NF} \
            /User time/{utime=$4} \
            /System time/{stime=$4} \
            END{print caller,time,cpu,mem,utime,stime}'
        grep -E 'Elapsed|CPU|Maximum' "${SAMPLE}_svim.time" | awk -v caller="SVIM" 'BEGIN{OFS=","} \
            /Elapsed \(wall clock\)/{split($8,a,":"); time=a[1]*60 + a[2]} \
            /CPU/{cpu=$NF} \
            /Maximum resident set size/{mem=$NF} \
            /User time/{utime=$4} \
            /System time/{stime=$4} \
            END{print caller,time,cpu,mem,utime,stime}'
    } > "${SAMPLE}_performance_report.csv"

    # Cleanup
    rm -f "${RAW_READS}" "${SAMPLE}_trimmed.fastq" "${SAMPLE}_mapped.sam" "${SAMPLE}_mapped.bam" \
          "${SAMPLE}_mapped.sort.bam" "${SAMPLE}_mapped.sort.bam.bai" \
          "${SAMPLE}_sniffles.time" "${SAMPLE}_cutesv.time" "${SAMPLE}_debreak.time" "${SAMPLE}_svim.time"

    cd ..
    echo "✅ Finished $SAMPLE"
done

# Generate combined performance report for all samples
echo "📊 Generating combined performance report"
{
    echo "Sample,Caller,Elapsed Time (s),CPU Percentage,Max Memory (KB),User Time (s),System Time (s)"
    for SAMPLE in "${SAMPLES[@]}"; do
        if [ -f "${SAMPLE}/${SAMPLE}_performance_report.csv" ]; then
            # Skip header and add sample name to each line
            tail -n +2 "${SAMPLE}/${SAMPLE}_performance_report.csv" | sed "s/^/${SAMPLE},/"
        fi
    done
} > "combined_performance_report.csv"

echo "🏁 All samples processed! Performance report saved to combined_performance_report.csv"
