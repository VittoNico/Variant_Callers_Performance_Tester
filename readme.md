# Nanopore SV Duplication Detection Pipeline

This pipeline is designed for structural variant (SV) discovery focused specifically on **duplications (SVTYPE=DUP)** using Oxford Nanopore sequencing data. It incorporates QC, trimming, mapping, SV calling (Sniffles, cuteSV, DeBreak, SVIM), and benchmarking using Truvari.

---

## 🔧 Requirements

* Linux/Unix system
* Miniconda (with conda-forge and bioconda channels configured)
* Tools (installed in appropriate conda environments):

  * `sniffles`
  * `cuteSV`
  * `debreak`
  * `svim`
  * `truvari`
  * `NanoPlot`
  * `NanoFilt`
  * `porechop`
  * `fasterq-dump`, `prefetch` (from `sra-tools`)
  * `minimap2`, `samtools`

---

## 📁 Input

* `sample.txt`: Plaintext file with one SRR ID per line.
* `reference.fasta`: Reference genome.

---

## 🧪 Pipeline Steps

### 1. Read Samples

```bash
mapfile -t SAMPLES < sample.txt
```

### 2. For Each Sample

```bash
for SAMPLE in "${SAMPLES[@]}"; do
```

### 3. Download Raw Nanopore Data

```bash
prefetch "$SAMPLE" --progress
fasterq-dump --threads 24 --outdir ./ "$SAMPLE"
```

### 4. Merge All FASTQ Files

```bash
find . -type f -name "${SAMPLE}*.fastq" -exec cat {} + > "${SAMPLE}_combined.fastq"
```

### 5. Clean Intermediate Files

```bash
rm -rf "$SAMPLE" *.sra
find . -type f -name "${SAMPLE}*.fastq" ! -name "${SAMPLE}_combined.fastq" -delete
```

### 6. QC - NanoPlot

```bash
NanoPlot --fastq "${SAMPLE}_combined.fastq" -o nanoplot_${SAMPLE} --threads 24
```

### 7. Adapter Trimming - Porechop

```bash
porechop -i "${SAMPLE}_combined.fastq" -o "${SAMPLE}_trimmed.fastq" --threads 24
```

### 8. Read Filtering - NanoFilt

```bash
cat "${SAMPLE}_trimmed.fastq" | NanoFilt -q 10 > "${SAMPLE}_filtered.fastq"
```

### 9. Mapping with Minimap2

```bash
minimap2 -x map-ont -A 2 -B 8 -O 4,24 -E 2,1 -m 5000 -z 200,100 --secondary=no -s 30 -a reference.fasta "${SAMPLE}_filtered.fastq" -t 24 > "${SAMPLE}_mapped.sam"
```

### 10. Convert SAM to BAM + Sorting + Indexing

```bash
samtools view -bS "${SAMPLE}_mapped.sam" > "${SAMPLE}_mapped.bam"
samtools sort "${SAMPLE}_mapped.bam" -o "${SAMPLE}_mapped.sort.bam"
samtools index "${SAMPLE}_mapped.sort.bam"
samtools depth -a "${SAMPLE}_mapped.sort.bam" > "${SAMPLE}_coverage.txt"
```

### 11. SV Calling: Sniffles

```bash
sniffles --input "${SAMPLE}_mapped.sort.bam" --vcf "${SAMPLE}_sniffles.vcf" --allow-overwrite --svtypes DUP
```

### 12. SV Calling: cuteSV

```bash
cuteSV "${SAMPLE}_mapped.sort.bam" reference.fasta "${SAMPLE}_cutesv.vcf" cute --max_size 1100000 --min_support 5 --genotype --sv_type DUP
```

### 13. SV Calling: DeBreak

```bash
conda activate debreak
debreak --bam "${SAMPLE}_mapped.sort.bam" --outpath debreak_${SAMPLE} --sv_types DUP
conda deactivate
```

### 14. SV Calling: SVIM

```bash
conda activate svim
mkdir -p svim_${SAMPLE}
svim alignment svim_${SAMPLE} "${SAMPLE}_mapped.sort.bam" --reference reference.fasta --skip_unspecified --sv_types DUP
conda deactivate
```

### 15. Truvari Benchmarking (Optional)

Ensure you have a merged or collapsed truth VCF (`truth.vcf`) available for comparison.

```bash
truvari bench -b ../truth.vcf -c ${SAMPLE}_cutesv.vcf -o truvari_${SAMPLE}_cutesv --includebed target_regions.bed -p 0 -P 0.5 -r 500 -O 0 --passonly --gtcomp --pctseq 0 --multimatch
```

### 16. Cleanup (Keep QC, Coverage, BAM index, and VCFs)

```bash
rm -f "${SAMPLE}_combined.fastq" "${SAMPLE}_trimmed.fastq" "${SAMPLE}_mapped.sam" "${SAMPLE}_mapped.bam" "${SAMPLE}_mapped.sort.bam" "${SAMPLE}_mapped.sort.bam.bai"
```

---

## 📤 Output (per sample)

* Quality plots (`nanoplot_*`)
* Filtered FASTQ (`*_filtered.fastq`)
* VCF files: `*_sniffles.vcf`, `*_cutesv.vcf`, `debreak_*/`, `svim_*/`
* Coverage: `*_coverage.txt`
* Benchmark (optional): `truvari_*`

---

## ✅ Focus: Duplications Only

Each caller is configured with parameters that restrict SV detection to **DUP** events to streamline downstream analysis and reduce irrelevant output.

---

## 🧠 Tips

* Be sure all tools are in their respective `conda` environments.
* Always check basecalling quality before investing compute time.
* You can parallelize this pipeline using GNU Parallel or a job scheduler (e.g. SLURM).

---

## 📚 References

* Sniffles2: [https://github.com/schwarzlab/sniffles](https://github.com/schwarzlab/sniffles)
* cuteSV: [https://github.com/tjiangHIT/cuteSV](https://github.com/tjiangHIT/cuteSV)
* SVIM: [https://github.com/eldariont/svim](https://github.com/eldariont/svim)
* DeBreak: [https://github.com/chhylp123/DeBreak](https://github.com/chhylp123/DeBreak)
* Truvari: [https://github.com/ACEnglish/truvari](https://github.com/ACEnglish/truvari)
* NanoPlot: [https://github.com/wdecoster/NanoPlot](https://github.com/wdecoster/NanoPlot)
* Porechop: [https://github.com/rrwick/Porechop](https://github.com/rrwick/Porechop)
* NanoFilt: [https://github.com/wdecoster/nanofilt](https://github.com/wdecoster/nanofilt)
