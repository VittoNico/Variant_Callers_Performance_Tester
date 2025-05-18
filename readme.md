# Nanopore SV Duplication Detection Pipeline

This pipeline is designed for structural variant (SV) discovery focused specifically on **duplications (SVTYPE=DUP)** using Oxford Nanopore sequencing data. It incorporates QC, trimming, mapping, SV calling (Sniffles, cuteSV, DeBreak, SVIM, along with a report of their performance in terms of CPU usage and time), and benchmarking using Truvari. 

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

# Structural Variant Detection Pipeline

This pipeline processes Nanopore sequencing data to detect tandem duplications (DUP SVs) using multiple callers, then merges and evaluates the results.

## Pipeline Overview

### 1. Data Download and Processing
- **Input**: Takes SRA accession numbers from `sample.txt`
- **Download**: Retrieves Nanopore reads using `prefetch` and `fasterq-dump`
- **QC**: Performs initial quality assessment with `NanoPlot`
- **Processing**:
  - Adapter trimming with `Porechop`
  - Quality filtering (Q≥10) using `NanoFilt`

### 2. Mapping
- Aligns filtered reads to reference genome using `minimap2` with optimized parameters for Nanopore data
- Processes alignments with `samtools`:
  - Converts SAM to BAM
  - Sorts and indexes BAM files
  - Generates coverage statistics

### 3. Variant Calling (Tandem Duplications Only)
Runs four specialized SV callers in parallel:
- **Sniffles**: Sensitive breakpoint detection
- **cuteSV**: Precise SV calling with size filtering (≤1.1Mb)
- **DeBreak**: Specialized for breakpoint refinement
- **SVIM**: Tuned specifically for tandem duplications  
All callers are configured to report only duplication-type SVs (DUP)

### 4. Truvari Analysis
- **Collapse**: Merges variant calls from all four methods using `truvari collapse` with:
  - 1000bp reciprocal overlap threshold
  - 70% sequence similarity requirement
  - Priority given to callers in specified order
- **Benchmarking**: Evaluates merged calls against truth set with `truvari bench`
  - Generates precision/recall metrics
  - Produces comprehensive output reports


### 📤 Output (per sample)

* Quality plots (`nanoplot_*`)
* VCF files: `*_sniffles.vcf`, `*_cutesv.vcf`, `debreak_*/`, `svim_*/`
* Benchmark (optional): `truvari_*`

---

## ✅ Focus: Duplications Only

Each caller is configured with parameters that restrict SV detection to **DUP** events to streamline downstream analysis and reduce irrelevant output.

---

## 🧠 Warnings

* Be sure all tools are in their respective `conda` environments.
* Always check basecalling quality before investing compute time.
* You can parallelize this pipeline using GNU Parallel or a job scheduler (e.g. SLURM).
* The .sh file of the tools has been tested and it is ready to use. The python version needs testing

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
