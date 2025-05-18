#!/usr/bin/env python3
import os
import subprocess
import time
import psutil
from pathlib import Path
import json

def run_command(cmd, check=True, stage_name=""):
    """Helper function to run shell commands with timing"""
    print(f"Running {stage_name}: {' '.join(cmd)}")
    start_time = time.time()
    cpu_start = psutil.cpu_percent(interval=None)
    mem_start = psutil.virtual_memory().used
    
    process = subprocess.run(cmd, check=check)
    
    elapsed = time.time() - start_time
    cpu_usage = psutil.cpu_percent(interval=None) - cpu_start
    mem_usage = (psutil.virtual_memory().used - mem_start) / (1024 ** 2)  # MB
    
    return {
        "command": ' '.join(cmd),
        "time_sec": round(elapsed, 2),
        "cpu_percent": round(cpu_usage, 1),
        "memory_mb": round(mem_usage, 1)
    }

def process_sample(sample, ref, truth, threads):
    """Process a single sample through the entire pipeline"""
    print(f"🔽 Processing {sample}...")
    metrics = {}
    sample_dir = Path(sample)
    sample_dir.mkdir(exist_ok=True)
    os.chdir(sample_dir)
    
    try:
        # 1. Data Download and Processing
        metrics["prefetch"] = run_command(
            ["prefetch", sample, "--progress"], stage_name="prefetch")
        metrics["fasterq-dump"] = run_command(
            ["fasterq-dump", "--threads", str(threads), "--outdir", ".", sample],
            stage_name="fasterq-dump")
        
        # Combine FASTQs
        combined_fastq = f"{sample}_combined.fastq"
        with open(combined_fastq, "wb") as outfile:
            for f in Path(".").glob(f"{sample}*.fastq"):
                with open(f, "rb") as infile:
                    outfile.write(infile.read())
        
        # Cleanup
        for f in Path(".").glob(f"{sample}*.fastq"):
            if f.name != combined_fastq:
                f.unlink()
        Path(sample).rmdir()
        for f in Path(".").glob("*.sra"):
            f.unlink()
        
        # QC with NanoPlot
        nanoplot_dir = f"nanoplot_{sample}"
        Path(nanoplot_dir).mkdir(exist_ok=True)
        metrics["nanoplot"] = run_command(
            ["NanoPlot", "--fastq", combined_fastq, "-o", nanoplot_dir, "--threads", str(threads)],
            stage_name="NanoPlot")
        
        # 2. Read Processing
        trimmed_fastq = f"{sample}_trimmed.fastq"
        metrics["porechop"] = run_command(
            ["porechop", "-i", combined_fastq, "-o", trimmed_fastq, "--threads", str(threads)],
            stage_name="porechop")
        
        filtered_fastq = f"{sample}_filtered.fastq"
        with open(filtered_fastq, "w") as outfile:
            metrics["nanofilt"] = run_command(
                ["NanoFilt", "-q", "10"], 
                stdin=open(trimmed_fastq), 
                stdout=outfile,
                stage_name="NanoFilt")
        
        # 3. Mapping
        sam_file = f"{sample}_mapped.sam"
        with open(sam_file, "w") as sam_out:
            metrics["minimap2"] = run_command([
                "minimap2", "-x", "map-ont", "-A", "2", "-B", "8", "-O", "4,24", "-E", "2,1",
                "-m", "5000", "-z", "200,100", "--secondary=no", "-s", "30", "-a", 
                f"../{ref}", filtered_fastq, "-t", str(threads)
            ], stdout=sam_out, stage_name="minimap2")
        
        # 4. BAM Processing
        bam_file = f"{sample}_mapped.bam"
        sorted_bam = f"{sample}_mapped.sort.bam"
        metrics["samtools_view"] = run_command(
            ["samtools", "view", "-bS", sam_file, "-o", bam_file],
            stage_name="samtools view")
        metrics["samtools_sort"] = run_command(
            ["samtools", "sort", bam_file, "-o", sorted_bam],
            stage_name="samtools sort")
        metrics["samtools_index"] = run_command(
            ["samtools", "index", sorted_bam],
            stage_name="samtools index")
        metrics["samtools_depth"] = run_command(
            ["samtools", "depth", "-a", sorted_bam, "-o", f"{sample}_coverage.txt"],
            stage_name="samtools depth")
        
        # 5. Variant Calling with Metrics
        # Sniffles
        sniffles_tmp = f"{sample}_sniffles_tmp.vcf"
        sniffles_vcf = f"{sample}_sniffles.vcf"
        metrics["sniffles"] = run_command(
            ["sniffles", "--input", sorted_bam, "--vcf", sniffles_tmp, "--allow-overwrite"],
            stage_name="sniffles")
        with open(sniffles_vcf, "w") as outfile:
            run_command(["bcftools", "view", "-i", 'SVTYPE="DUP"', sniffles_tmp], stdout=outfile)
        Path(sniffles_tmp).unlink()
        
        # cuteSV
        cutesv_tmp = f"{sample}_cutesv_tmp.vcf"
        cutesv_vcf = f"{sample}_cutesv.vcf"
        metrics["cuteSV"] = run_command([
            "cuteSV", sorted_bam, f"../{ref}", cutesv_tmp, "cute",
            "--max_size", "1100000", "--min_support", "5", "--threads", str(threads)
        ], stage_name="cuteSV")
        with open(cutesv_vcf, "w") as outfile:
            run_command(["bcftools", "view", "-i", 'SVTYPE="DUP"', cutesv_tmp], stdout=outfile)
        Path(cutesv_tmp).unlink()
        
        # DeBreak
        debreak_dir = Path(f"debreak_{sample}")
        debreak_dir.mkdir(exist_ok=True)
        metrics["debreak"] = run_command(
            ["debreak", "--bam", sorted_bam, "--outpath", str(debreak_dir), "--threads", str(threads)],
            stage_name="debreak")
        with open(debreak_dir/"dup.vcf", "w") as outfile:
            run_command(["bcftools", "view", "-i", 'SVTYPE="DUP"', debreak_dir/"variants.vcf"], stdout=outfile)
        (debreak_dir/"variants.vcf").unlink()
        (debreak_dir/"dup.vcf").rename(debreak_dir/"variants.vcf")
        
        # SVIM
        svim_dir = Path(f"svim_{sample}")
        svim_dir.mkdir(exist_ok=True)
        metrics["svim"] = run_command([
            "svim", "alignment", str(svim_dir), sorted_bam, "--reference", f"../{ref}",
            "--tandem-duplications", "--skip-insertions", "--skip-deletions",
            "--skip-inversions", "--skip-translocations", "--threads", str(threads)
        ], stage_name="svim")
        
        # 6. Truvari Analysis
        truvari_dir = Path(f"truvari_{sample}")
        truvari_dir.mkdir(exist_ok=True)
        metrics["truvari_collapse"] = run_command([
            "truvari", "collapse",
            "-i", sniffles_vcf, cutesv_vcf, str(svim_dir/"variants.vcf"), str(debreak_dir/"variants.vcf"),
            "-o", str(truvari_dir/f"{sample}_merged.vcf"),
            "-f", f"../{ref}",
            "--passonly", "--chain", "-r", "1000", "-p", "0.0", "-P", "0.7", "-O", "0.2", "-k", "first"
        ], stage_name="truvari collapse")
        
        truvari_bench_dir = Path(f"truvari_bench_{sample}")
        truvari_bench_dir.mkdir(exist_ok=True)
        metrics["truvari_bench"] = run_command([
            "truvari", "bench",
            "-b", truth,
            "-c", str(truvari_dir/f"{sample}_merged.vcf"),
            "-f", f"../{ref}",
            "-o", str(truvari_bench_dir),
            "--passonly", "--gtcomp", "-r", "1000", "-p", "0.0", "-P", "0.7", "-O", "0.2"
        ], stage_name="truvari bench")
        
        # 7. Generate Performance Report
        generate_performance_report(metrics, sample)
        
        print(f"✅ Finished {sample}")
        
    finally:
        os.chdir("..")

def generate_performance_report(metrics, sample):
    """Generate JSON report with timing and resource usage"""
    report = {
        "sample": sample,
        "total_time": round(sum(m["time_sec"] for m in metrics.values()), 2),
        "caller_metrics": {
            "sniffles": metrics["sniffles"],
            "cuteSV": metrics["cuteSV"],
            "debreak": metrics["debreak"],
            "svim": metrics["svim"]
        },
        "resource_summary": {
            "max_cpu": max(m["cpu_percent"] for m in metrics.values()),
            "max_memory": max(m["memory_mb"] for m in metrics.values()),
            "total_cpu_hours": round(sum(
                m["time_sec"] * m["cpu_percent"] / 100 / 3600 
                for m in metrics.values()), 3)
        }
    }
    
    with open(f"performance_{sample}.json", "w") as f:
        json.dump(report, f, indent=2)
    
    # Print summary to console
    print(f"\n⏱️ Performance Summary for {sample}:")
    print(f"Total time: {report['total_time']} seconds")
    print("Caller Times:")
    for caller, data in report["caller_metrics"].items():
        print(f"  {caller:<8}: {data['time_sec']}s (CPU: {data['cpu_percent']}%, Mem: {data['memory_mb']}MB)")

def main():
    with open("sample.txt") as f:
        samples = [line.strip() for line in f if line.strip()]
    
    ref = "reference.fasta"
    truth = "../truth_set.vcf.gz"
    threads = 24
    
    for sample in samples:
        process_sample(sample, ref, truth, threads)

if __name__ == "__main__":
    main()
