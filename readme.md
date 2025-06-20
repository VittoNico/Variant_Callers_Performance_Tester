# 🧬 Nanopore SV Duplication Detection Pipeline

A complete pipeline for **structural variant discovery**, focusing on **duplications (SVTYPE=DUP)**, using **Oxford Nanopore** data. The workflow includes **read preprocessing, mapping, SV calling (8 tools), and performance reporting**. Designed for **large-scale SVs**, with optimized parameters for duplications.

---

## 🚀 Overview

This pipeline:

1. Downloads Nanopore reads using `SRA` tools.
2. Trims and filters reads for quality.
3. Aligns reads to a reference genome with `minimap2`.
4. Calls structural variants using **eight SV callers**:
   - Sniffles1
   - Sniffles2
   - cuteSV
   - DeBreak
   - Dysgu
   - SVIM
   - pbsv
   - SVDSS
5. Measures performance (time, memory, CPU) for each tool.
6. Produces VCF files and indexed BAMs.
7. Supports downstream benchmarking with `Truvari`.

---

## 🛠 Requirements

- Unix/Linux system
- [`conda`](https://docs.conda.io/en/latest/) (with `conda-forge` + `bioconda`)
- Tools (installed in dedicated conda environments):

| Task         | Tool(s) Used                                  |
|--------------|-----------------------------------------------|
| Download     | `prefetch`, `fasterq-dump` (from `sra-tools`) |
| QC           | `NanoPlot`, `NanoFilt`, `Porechop`            |
| Mapping      | `minimap2`, `samtools`                        |
| SV Calling   | `sniffles`, `cuteSV`, `debreak`, `svim`, `pbsv`, `dysgu`, `SVDSS` |
| Evaluation   | `Truvari`, `bcftools`                         |

---

## 📁 Inputs

- `sample.txt`: One **SRA accession** per line.
- Reference genome: FASTA format (`*.fasta` or `*.fna`).
  - Make sure it's indexed if needed (`samtools faidx`, `minimap2 -d`).
- `*.svsig.gz` files are reused for **pbsv** to avoid repeating `discover`.

---

## 🧬 Pipeline Structure

### Step 1: Data Download and Preprocessing

- Download reads with `prefetch`, convert with `fasterq-dump`.
- Combine and QC reads (`NanoPlot`, `Porechop`, `NanoFilt`).
- Trimming: adapters + low Q-score filtering.

### Step 2: Mapping

- Align reads using `minimap2 -x map-ont`.
- Convert, sort, index BAM files.
- Generate per-base coverage.

### Step 3: Structural Variant Calling

The pipeline uses 8 state-of-the-art SV callers, each with unique algorithms and strengths:

| Tool        | Description |
|-------------|-------------|
| **Sniffles1** | One of the first SV callers for long reads; identifies SVs using split-read and coverage signals. |
| **Sniffles2** | Improved version of Sniffles; supports better genotyping and is optimized for high-throughput data. |
| **cuteSV**    | Fast and memory-efficient SV caller that clusters signals and refines SV boundaries. |
| **DeBreak**   | Accurate SV caller using partial order alignment and duplication recovery for noisy long reads. |
| **Dysgu**     | Versatile SV tool supporting many sequencing types; efficient with nanopore data. |
| **SVIM**      | Detects a wide range of SVs using a signal-based approach from long-read alignments. |
| **pbsv**      | PacBio’s official SV caller, works well with aligned reads including ONT; supports split-read analysis. |
| **SVDSS**     | High-resolution caller using smoothing and suffix-filtering strategies to detect SVs with precision. |

All tools support **multi-threading** and are configured for maximum SV size where applicable.

### Step 4: Performance Logging

Each SV caller logs:
- Elapsed wall time
- CPU usage
- Peak memory
- User/system time

Results are saved in `<SAMPLE>_performance_report.csv`.

---

## 📤 Output Per Sample

| File/Folder | Description |
|-------------|-------------|
| `*_mapped.sort.bam` + `.bai` | Aligned BAM and index |
| `*_coverage.txt`             | Per-base coverage |
| `*_caller.vcf`               | VCF output from each SV tool |
| `<tool>_<SAMPLE>/`           | Tool-specific outputs (e.g., `svim/`, `pbsv/`) |
| `*_performance_report.csv`   | Time, memory, CPU stats per caller |
| `nanoplot_<SAMPLE>/`         | Quality metrics and plots |

---

## 🔍 Post-processing (optional)

- Run [`truvari`](https://github.com/ACEnglish/truvari) to compare output to golden standard.
- Filter only `SVTYPE=DUP` from VCFs using `bcftools view -i 'INFO/SVTYPE=="DUP"'`.

---

## 📚 References

- [Sniffles](https://github.com/schwarzlab/sniffles)
- [cuteSV](https://github.com/tjiangHIT/cuteSV)
- [DeBreak](https://github.com/chhylp123/DeBreak)
- [SVIM](https://github.com/eldariont/svim)
- [pbsv](https://github.com/PacificBiosciences/pbsv)
- [Dysgu](https://github.com/simon-ritchie/dysgu)
- [SVDSS](https://github.com/mucole/svdss)
- [Truvari](https://github.com/ACEnglish/truvari)
- [NanoPlot](https://github.com/wdecoster/NanoPlot)
- [Porechop](https://github.com/rrwick/Porechop)
- [NanoFilt](https://github.com/wdecoster/nanofilt)

---

## 💡 Tips

- Run on an HPC cluster or use `GNU parallel` to scale samples.
- Ensure your reference genome and read quality are appropriate.
- Tune parameters (e.g. `--max_size`) depending on expected SV length.

---

## 🧪 Status

✅ Bash version: fully tested  
⚠️ Python version: still under testing
