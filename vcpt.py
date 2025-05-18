vcpy.py
#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path

def run_command(cmd, check=True):
    """Helper function to run shell commands"""
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=check)

def process_sample(sample, ref, truth, threads):
    """Process a single sample through the entire pipeline"""
    print(f" Processing {sample}...")
    
    # Create sample directory
    sample_dir = Path(sample)
    sample_dir.mkdir(exist_ok=True)
    os.chdir(sample_dir)
    
    try:
        # Download SRA run
        run_command(["prefetch", sample, "--progress"])
        run_command(["fasterq-dump", "--threads", str(threads), "--outdir", ".", sample])
        
        # Combine FASTQ files
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
        run_command(["NanoPlot", "--fastq", combined_fastq, "-o", nanoplot_dir, "--threads", str(threads)])
        
        # Trimming with Porechop
        trimmed_fastq = f"{sample}_trimmed.fastq"
        run_command(["porechop", "-i", combined_fastq, "-o", trimmed_fastq, "--threads", str(threads)])
        
        # Quality filtering with NanoFilt
        filtered_fastq = f"{sample}_filtered.fastq"
        with open(filtered_fastq, "w") as outfile:
            run_command(["NanoFilt", "-q", "10"], stdin=open(trimmed_fastq), stdout=outfile)
        
        # Mapping with minimap2
        sam_file = f"{sample}_mapped.sam"
        with open(sam_file, "w") as sam_out:
            run_command([
                "minimap2", "-x", "map-ont", "-A", "2", "-B", "8", "-O", "4,24", "-E", "2,1",
                "-m", "5000", "-z", "200,100", "--secondary=no", "-s", "30", "-a", 
                f"../{ref}", filtered_fastq, "-t", str(threads)
            ], stdout=sam_out)
        
        # BAM processing
        bam_file = f"{sample}_mapped.bam"
        sorted_bam = f"{sample}_mapped.sort.bam"
        run_command(["samtools", "view", "-bS", sam_file, "-o", bam_file])
        run_command(["samtools", "sort", bam_file, "-o", sorted_bam])
        run_command(["samtools", "index", sorted_bam])
        run_command(["samtools", "depth", "-a", sorted_bam, "-o", f"{sample}_coverage.txt"])
        
        # SV Calling: Sniffles (DUP only)
        sniffles_tmp = f"{sample}_sniffles_tmp.vcf"
        sniffles_vcf = f"{sample}_sniffles.vcf"
        run_command(["sniffles", "--input", sorted_bam, "--vcf", sniffles_tmp, "--allow-overwrite"])
        with open(sniffles_vcf, "w") as outfile:
            run_command(["bcftools", "view", "-i", 'SVTYPE="DUP"', sniffles_tmp], stdout=outfile)
        Path(sniffles_tmp).unlink()
        
        # SV Calling: cuteSV (DUP only)
        cutesv_tmp = f"{sample}_cutesv_tmp.vcf"
        cutesv_vcf = f"{sample}_cutesv.vcf"
        run_command([
            "cuteSV", sorted_bam, f"../{ref}", cutesv_tmp, "cute",
            "--max_size", "1100000", "--min_support", "5"
        ])
        with open(cutesv_vcf, "w") as outfile:
            run_command(["bcftools", "view", "-i", 'SVTYPE="DUP"', cutesv_tmp], stdout=outfile)
        Path(cutesv_tmp).unlink()
        
        # DeBreak (DUP only)
        debreak_dir = Path(f"debreak_{sample}")
        debreak_dir.mkdir(exist_ok=True)
        run_command(["debreak", "--bam", sorted_bam, "--outpath", str(debreak_dir)])
        with open(debreak_dir/"dup.vcf", "w") as outfile:
            run_command(["bcftools", "view", "-i", 'SVTYPE="DUP"', debreak_dir/"variants.vcf"], stdout=outfile)
        (debreak_dir/"variants.vcf").unlink()
        (debreak_dir/"dup.vcf").rename(debreak_dir/"variants.vcf")
        
        # SVIM (tandem duplications only)
        svim_dir = Path(f"svim_{sample}")
        svim_dir.mkdir(exist_ok=True)
        run_command([
            "svim", "alignment", str(svim_dir), sorted_bam, "--reference", f"../{ref}",
            "--tandem-duplications", "--skip-insertions", "--skip-deletions",
            "--skip-inversions", "--skip-translocations"
        ])
        
        # Truvari collapse
        truvari_dir = Path(f"truvari_{sample}")
        truvari_dir.mkdir(exist_ok=True)
        run_command([
            "truvari", "collapse",
            "-i", sniffles_vcf, cutesv_vcf, str(svim_dir/"variants.vcf"),
            "-o", str(truvari_dir/f"{sample}_merged.vcf"),
            "-f", f"../{ref}",
            "--passonly", "--chain", "-r", "1000", "-p", "0.0", "-P", "0.7", "-O", "0.2", "-k", "first"
        ])
        
        # Truvari benchmarking
        truvari_bench_dir = Path(f"truvari_bench_{sample}")
        truvari_bench_dir.mkdir(exist_ok=True)
        run_command([
            "truvari", "bench",
            "-b", truth,
            "-c", str(truvari_dir/f"{sample}_merged.vcf"),
            "-f", f"../{ref}",
            "-o", str(truvari_bench_dir),
            "--passonly", "--gtcomp", "-r", "1000", "-p", "0.0", "-P", "0.7", "-O", "0.2"
        ])
        
        # Cleanup
        for f in [combined_fastq, f"{sample}_trimmed.fastq", sam_file, bam_file, sorted_bam, f"{sorted_bam}.bai"]:
            if Path(f).exists():
                Path(f).unlink()
        
        print(f"✅ Finished {sample}")
        
    finally:
        os.chdir("..")

def main():
    # Read samples from file
    with open("sample.txt") as f:
        samples = [line.strip() for line in f if line.strip()]
    
    # Configuration
    ref = "reference.fasta"
    truth = "../truth_set.vcf.gz"
    threads = 24
    
    # Initialize conda
    conda_sh = os.path.join(os.path.dirname(os.path.realpath(subprocess.check_output(["which", "conda"]).decode().strip())), 
                           "etc/profile.d/conda.sh")
    run_command(["source", conda_sh])
    
    # Process each sample
    for sample in samples:
        process_sample(sample, ref, truth, threads)

if __name__ == "__main__":
    main()
