 """
 #!/bin/bash

# ✅ Legge i sample SRR da un file chiamato sample.txt
mapfile -t SAMPLES < sample.txt

#Reference file
REF="reference.fasta"
THREADS=24

# Inizializza Conda
source "$(conda info --base)/etc/profile.d/conda.sh"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "🔽 Processing $SAMPLE..."

    mkdir -p "$SAMPLE"
    cd "$SAMPLE" || exit

    #Scarica la run SRA (Nanopore)
    prefetch "$SAMPLE" --progress
    fasterq-dump --threads "$THREADS" --outdir ./ "$SAMPLE"

    #Unisci tutti i fastq trovati nella cartella
    find . -type f -name "${SAMPLE}*.fastq" -exec cat {} + > "${SAMPLE}_combined.fastq"
    RAW_READS="${SAMPLE}_combined.fastq"

    # Rimuovi tutto tranne il fastq combinato
    rm -rf "$SAMPLE"
    rm -rf *.sra
    find . -type f -name "${SAMPLE}*.fastq" ! -name "${RAW_READS}" -delete

    #QC iniziale con NanoPlot
    mkdir -p nanoplot_${SAMPLE}
    NanoPlot --fastq "${RAW_READS}" -o nanoplot_${SAMPLE} --threads $THREADS

    # Trimming con Porechop
    conda activate sv-env
    porechop -i "${RAW_READS}" -o "${SAMPLE}_trimmed.fastq" --threads $THREADS

    # 🧹 Filtro di qualità con NanoFilt (es. Q≥10)
    cat "${SAMPLE}_trimmed.fastq" | NanoFilt -q 10 > "${SAMPLE}_filtered.fastq"
    FILTERED_READS="${SAMPLE}_filtered.fastq"

    # 🧬 Mapping
    minimap2 -x map-ont -A 2 -B 8 -O 4,24 -E 2,1 -m 5000 -z 200,100 --secondary=no -s 30 -a ../"$REF" "$FILTERED_READS" -t $THREADS > "${SAMPLE}_mapped.sam"

    # 🧪 BAM processing
    samtools view -bS "${SAMPLE}_mapped.sam" > "${SAMPLE}_mapped.bam"
    samtools sort "${SAMPLE}_mapped.bam" -o "${SAMPLE}_mapped.sort.bam"
    samtools index "${SAMPLE}_mapped.sort.bam"
    samtools depth -a "${SAMPLE}_mapped.sort.bam" > "${SAMPLE}_coverage.txt"

    #SV Calling
    sniffles --input "${SAMPLE}_mapped.sort.bam" --vcf "${SAMPLE}_sniffles.vcf" --allow-overwrite
    cuteSV "${SAMPLE}_mapped.sort.bam" ../"$REF" "${SAMPLE}_cutesv.vcf" cute --max_size 1100000 --min_support 5

    #DeBreak
    conda activate debreak
    mkdir -p debreak_${SAMPLE}
    debreak --bam "${SAMPLE}_mapped.sort.bam" --outpath debreak_${SAMPLE}
    conda deactivate

    #SVIM
    conda activate svim
    mkdir -p svim_${SAMPLE}
    svim alignment svim_${SAMPLE} "${SAMPLE}_mapped.sort.bam" --reference ../"$REF"
    conda deactivate

    #Cleanup (mantieni solo output utile e log)
    rm -f "${RAW_READS}" "${SAMPLE}_trimmed.fastq" "${SAMPLE}_mapped.sam" "${SAMPLE}_mapped.bam" "${SAMPLE}_mapped.sort.bam" "${SAMPLE}_mapped.sort.bam.bai"

    cd ..
    echo "✅ Finished $SAMPLE"
done
"""
