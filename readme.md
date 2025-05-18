# Structural Variant Detection with Nanopore Long Reads (Duplication-Focused)

This pipeline performs a full analysis of structural variants (SVs), specifically targeting **duplications (SVTYPE=DUP)**, using long-read sequencing data from Oxford Nanopore Technologies. It automates read preprocessing, alignment, SV calling, and benchmarking using `Sniffles`, `cuteSV`, `SVIM`, `DeBreak`, and `Truvari`. The pipeline supports multiple samples listed in a `sample.txt` file and runs each analysis in an isolated environment.

---

## 🔧 Stage 0: Requirements & Installation

### Software Requirements

* `conda`
* `minimap2`
* `samtools`
* `sniffles`
* `cuteSV`
* `debreak`
* `svim`
* `truvari`
* `bcftools`
* `NanoPlot`, `NanoFilt`, `Porechop`

Ensure all tools are available in appropriate `conda` environments:

```bash
conda create -n sv-env minimap2 samtools sniffles cutadapt bcftools porechop nanoplot nanofilt -c bioconda -c conda-forge
conda create -n debreak debreak -c bioconda
conda create -n svim svim -c bioconda
conda create -n truvari truvari -c bioconda
```

---

## 🧬 Stage 1: Sample Preparation

### Input File

A file named `sample.txt` should list the SRA Run Accessions (e.g., SRRXXXXXXX) one per line.

### Reference Genome

Place the reference genome in the current working directory and name it `reference.fasta`. Index it using minimap2 and samtools as needed.

---

## 🚀 Stage 2: Pipeline Execution

For each sample listed in `sample.txt`, the following steps are performed:

### 1. Download & Preprocessing

* **Download**: Run is downloaded using `prefetch` and converted to FASTQ with `fasterq-dump`.
* **Merge**: All FASTQ files are merged into one combined file.
* **QC**: Initial QC is done with `NanoPlot`.
* **Trimming**: Adapters are removed with `Porechop`.
* **Filtering**: Quality-filtered with `NanoFilt` (Q >= 10).

### 2. Alignment

* **minimap2**: Long reads are aligned to the reference with specified scoring parameters.
* **samtools**: BAM file is sorted and indexed. Coverage is extracted with `samtools depth`.

### 3. SV Calling (Only DUPs)

* **Sniffles**: SV calling → filtered with `bcftools` for `SVTYPE=DUP`
* **cuteSV**: Same as above.
* **SVIM**: Called with `--tandem-duplications` only.
* **DeBreak**: Full call → filtered for duplications.

### 4. Benchmarking

* **Truvari bench**: Evaluates precision, recall, F1, genotyping accuracy against a `truth_set.vcf.gz` file filtered to duplications using `truvari collapse`.

---

## 📂 Output Organization

Each sample has its own folder:

* `nanoplot_<sample>`: Initial QC reports
* `<sample>_coverage.txt`: Coverage profile
* `<tool>_<sample>.vcf`: Filtered VCF for duplications
* `<tool>_<sample>` folder: Caller-specific results
* `truvari_<sample>`: Benchmarking result folder

---

## 📋 Example Summary Script

Use this to summarize recall, precision, F1:

```bash
grep -A 5 "Overall" truvari_<sample>/summary.txt
```

Or automate for all:

```bash
for d in truvari_*/summary.txt; do echo "$d"; grep -A 5 Overall "$d"; done > summary_report.txt
```

---

## 📌 Notes

* This pipeline focuses **only on duplications**. Other SV types are filtered out to reduce noise.
* Make sure the reference and VCF files are indexed (`samtools faidx`, `bgzip`, `tabix`).
* Use `conda deactivate` between environment switches.

---

## ✅ Summary

This pipeline provides an automated, scalable solution for high-confidence detection of duplications using long reads. It supports preprocessing, trimming, alignment, SV calling with 4 different tools, and benchmarking. The output is well-structured and suitable for publication-ready results.

---

## 🔗 Reference Tools

* [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
* [cuteSV](https://github.com/tjiangHIT/cuteSV)
* [SVIM](https://github.com/eldariont/svim)
* [DeBreak](https://github.com/PacificBiosciences/debreak)
* [Truvari](https://github.com/ACEnglish/truvari)
* [NanoPlot](https://github.com/wdecoster/NanoPlot)
* [Porechop](https://github.com/rrwick/Porechop)
* [NanoFilt](https://github.com/wdecoster/nanofilt)

---

