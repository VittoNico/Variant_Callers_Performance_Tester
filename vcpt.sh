#!/bin/bash

# Read SRR samples from sample.txt
mapfile -t SAMPLES < sample.txt

# Reference file
REF="Reference.fasta"
THREADS=12

log_time() {
	local LABEL=$1
	local START_TIME=$2
	local END_TIME=$(date +%s)
	local ELAPSED=$((END_TIME - START_TIME))
	echo "🕒 $LABEL took ${ELAPSED}s"
}

for SAMPLE in "${SAMPLES[@]}"; do
	echo "🔽 Processing $SAMPLE..."
	mkdir -p "$SAMPLE"
	cd "$SAMPLE" || exit

	# ========== DATA DOWNLOAD ==========
	prefetch "$SAMPLE" --progress --max-size 100G
	fasterq-dump --threads "$THREADS" --outdir ./ "$SAMPLE"

	# Combine FASTQs
	find . -type f -name "${SAMPLE}*.fastq" -exec cat {} + > "${SAMPLE}_combined.fastq"
	RAW_READS="${SAMPLE}_combined.fastq"
	rm -rf "$SAMPLE" *.sra
	find . -type f -name "${SAMPLE}*.fastq" ! -name "${RAW_READS}" -delet



	# ========== QC ==========
	mkdir -p nanoplot_${SAMPLE}
	NanoPlot --fastq "${RAW_READS}" -o nanoplot_${SAMPLE} --threads $THREADS

	# ========== PROCESSING ==========
	porechop -i "${RAW_READS}" -o "${SAMPLE}_trimmed.fastq" --threads $THREADS
	cat "${SAMPLE}_trimmed.fastq" | NanoFilt -q 7 > "${SAMPLE}_filtered.fastq"
	FILTERED_READS="${SAMPLE}_filtered.fastq"



	# ========== MAPPING TO REFERENCE ==========
	minimap2 -x map-ont --MD -Y -L -a "$REF" "$FILTERED_READS" -t $THREADS > "${SAMPLE}_mapped.sam"
	samtools view -bS "${SAMPLE}_mapped.sam" > "${SAMPLE}_mapped.bam"
	samtools sort "${SAMPLE}_mapped.bam" -o "${SAMPLE}_mapped.sort.bam"
	samtools index "${SAMPLE}_mapped.sort.bam"
	samtools depth -a "${SAMPLE}_mapped.sort.bam" > "${SAMPLE}_coverage.txt"


	# ========== SV CALLING ==========
	source "$(conda info --base)/etc/profile.d/conda.sh"
	# Existing tools
	# Sniffles v1
	echo "Running Sniffles v1..."
	START=$(date +%s)
	conda activate sniffles1
	/usr/bin/time -v -o "${SAMPLE}_sniffles1.time" \
	sniffles -m "${SAMPLE}_mapped.sort.bam" -v "${SAMPLE}_sniffles1.vcf" -t $THREADS --cluster --genotype
	conda deactivate
	log_time "Sniffles v1 for $SAMPLE" "$START"

	# Sniffles2
	echo "Running Sniffles2..."
	START=$(date +%s)
	/usr/bin/time -v -o "${SAMPLE}_sniffles.time" \
	sniffles --threads $THREADS --input "${SAMPLE}_mapped.sort.bam" --vcf "${SAMPLE}_sniffles.vcf" --allow-overwrite
	log_time "Sniffles2 for $SAMPLE" "$START"

	# cuteSV
	echo "Running cuteSV..."
	START=$(date +%s)
	mkdir -p cute
	/usr/bin/time -v -o "${SAMPLE}_cutesv.time" \
	cuteSV "${SAMPLE}_mapped.sort.bam" "$REF" "${SAMPLE}_cutesv.vcf" cute --max_cluster_bias_INS 100 -t $THREADS --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --max_size 1000000000 --diff_ratio_merging_DEL 0.3
	log_time "cuteSV for $SAMPLE" "$START"


	# DeBreak
	echo "Running DeBreak..."
	START=$(date +%s)
	conda activate debreak
	mkdir -p debreak_${SAMPLE}
	/usr/bin/time -v -o "${SAMPLE}_debreak.time" \
	debreak --bam "${SAMPLE}_mapped.sort.bam" --outpath debreak_${SAMPLE} --rescue_large_ins --rescue_dup --poa -r "$REF" -t $THREADS
	conda deactivate
	log_time "DeBreak for $SAMPLE" "$START"

	# dysgu
	echo "Running dysgu..."
	START=$(date +%s)
	conda activate dysgu
	mkdir -p dysgu_${SAMPLE}
	/usr/bin/time -v -o "${SAMPLE}_dysgu.time" \
	dysgu run --overwrite --mode nanopore --min-size 30 "$REF" dysgu_${SAMPLE} "${SAMPLE}_mapped.sort.bam" > "${SAMPLE}_dysgu.vcf"
	conda deactivate
	log_time "dysgu for $SAMPLE" "$START"

	# SVIM
	echo "Running SVIM..."
	START=$(date +%s)
	conda activate svim
	mkdir -p svim_${SAMPLE}
	/usr/bin/time -v -o "${SAMPLE}_svim.time" \
	svim alignment svim_${SAMPLE} "${SAMPLE}_mapped.sort.bam" "$REF" --min_sv_size 30 --max_sv_size 1000000000
	conda deactivate
	log_time "SVIM for $SAMPLE" "$START"
    
	# pbsv
	echo "Running pbsv..."
	START=$(date +%s)
	conda activate pbsv
	samtools sort -n "${SAMPLE}_mapped.bam" -o "${SAMPLE}_mapped.sort_name.bam"
	mkdir -p pbsv_${SAMPLE}
	samtools addreplacerg -r "ID:${SAMPLE}\tSM:${SAMPLE}" -o "${SAMPLE}_mapped.sort_name.RG.bam" "${SAMPLE}_mapped.sort_name.bam"
	pbsv discover "${SAMPLE}_mapped.sort_name.RG.bam" "pbsv_${SAMPLE}/${SAMPLE}.svsig.gz"
	/usr/bin/time -v -o "${SAMPLE}_pbsv.time" \
	pbsv call "$REF" "pbsv_${SAMPLE}/${SAMPLE}.svsig.gz" "pbsv_${SAMPLE}/${SAMPLE}_pbsv.vcf" -j $THREADS --max-dup-length 1000000000
	conda deactivate
	log_time "pbsv for $SAMPLE" "$START"
    

	# SVDSS
	echo "🔍 Running SVDSS..."
	START=$(date +%s)
	conda activate svdss
	mkdir -p svdss_${SAMPLE}
	SVDSS index --threads $THREADS -d "$REF" > "${SAMPLE}.fmd"
	SVDSS smooth --threads $THREADS --reference "$REF" --bam "${SAMPLE}_mapped.sort.bam"  > "${SAMPLE}_mapped.sort.smooth.bam"
	samtools index "${SAMPLE}_mapped.sort.smooth.bam"
	SVDSS search --threads $THREADS --index "${SAMPLE}.fmd" --bam "${SAMPLE}_mapped.sort.smooth.bam" > "svdss_${SAMPLE}/specifics.txt"
	/usr/bin/time -v -o "${SAMPLE}_svdss.time" \
	SVDSS call --threads  $THREADS --reference "$REF" --bam "${SAMPLE}_mapped.sort.smooth.bam" --sfs "svdss_${SAMPLE}/specifics.txt" > "svdss_${SAMPLE}/${SAMPLE}_svdss.vcf"
	conda deactivate
	log_time "SVDSS for $SAMPLE" "$START"


	# ========== PERFORMANCE REPORT ==========
	echo "Tool,Elapsed(s),CPU%,MaxMemory(KB),UserTime(s),SysTime(s)" > "${SAMPLE}_performance_report.csv"

	for TOOL in sniffles1 sniffles cutesv debreak svim dysgu svdss pbsv; do
    	if [[ -f "${SAMPLE}_${TOOL}.time" ]]; then
        	ELAPSED=$(grep "Elapsed (wall clock)" "${SAMPLE}_${TOOL}.time" | awk '{split($8,t,":"); print (t[1]*60 + t[2])}')
        	CPU=$(grep "Percent of CPU" "${SAMPLE}_${TOOL}.time" | awk '{print $NF}')
        	MEM=$(grep "Maximum resident set size" "${SAMPLE}_${TOOL}.time" | awk '{print $NF}')
        	UTIME=$(grep "User time" "${SAMPLE}_${TOOL}.time" | awk '{print $4}')
        	STIME=$(grep "System time" "${SAMPLE}_${TOOL}.time" | awk '{print $4}')
        	echo "$TOOL,$ELAPSED,$CPU,$MEM,$UTIME,$STIME" >> "${SAMPLE}_performance_report.csv"
    	fi
	done

	# ========== CLEANUP ==========
	rm -f "${RAW_READS}" "${SAMPLE}_trimmed.fastq" "${SAMPLE}.fastq" "${SAMPLE}_mapped.sam" "${SAMPLE}_mapped.bam" "${SAMPLE}_mapped.sort.bam" "${SAMPLE}_filtered.fastq" *.bam "dysgu_${SAMPLE}/"
	cd ..
	echo "Finished $SAMPLE"
done


