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
    print(f"🔽 Processing {sample}...")
    
    # Create sample directory
    sample_dir = Path(sample)
    sample_dir.mkdir(exist_ok=True)
    os.chdir(sample_dir)
    
    try:
        # [Previous steps remain exactly the same until the Truvari collapse step]
        
        # Truvari collapse (now including DeBreak VCF)
        truvari_dir = Path(f"truvari_{sample}")
        truvari_dir.mkdir(exist_ok=True)
        run_command([
            "truvari", "collapse",
            "-i", sniffles_vcf, cutesv_vcf, str(svim_dir/"variants.vcf"), str(debreak_dir/"variants.vcf"),
            "-o", str(truvari_dir/f"{sample}_merged.vcf"),
            "-f", f"../{ref}",
            "--passonly", "--chain", "-r", "1000", "-p", "0.0", "-P", "0.7", "-O", "0.2", "-k", "first"
        ])
        
        # [Remaining steps remain the same]
        
    finally:
        os.chdir("..")

def main():
    # [Main function remains the same]
    
if __name__ == "__main__":
    main()
