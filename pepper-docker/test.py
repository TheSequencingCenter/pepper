import os
import subprocess

# Set base directory
print("Set base directory...")
base = "/home/ubuntu/pepper-docker/ont-case-study"

# Set up input data
print("Set up input data...")
input_dir = os.path.join(base, "input/data")
ref = "GRCh38_no_alt.chr20.fa"
bam = "HG002_pass_2_GRCh38.R10.4_q20.chr20.bam"

# Set the number of CPUs to use
print("Set number of cpu's...")
threads = "64"

# Set up output directory
print("Set up output directory...")
output_dir = os.path.join(base, "output")
output_prefix = "HG002_ONT_R10_Q20_2_GRCh38_PEPPER_Margin_DeepVariant.chr20"
output_vcf = output_prefix + ".vcf.gz"

# Create local directory structure
print("Create local directory structure...")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(input_dir, exist_ok=True)

# Download the data to input directory
print("Download data to input directory...")
subprocess.run(["wget", "-P", input_dir, "https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_pass_2_GRCh38.R10.4_q20.chr20.bam"])
subprocess.run(["wget", "-P", input_dir, "https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/HG002_pass_2_GRCh38.R10.4_q20.chr20.bam.bai"])
subprocess.run(["wget", "-P", input_dir, "https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr20.fa"])
subprocess.run(["wget", "-P", input_dir, "https://storage.googleapis.com/pepper-deepvariant-public/usecase_data/GRCh38_no_alt.chr20.fa.fai"])

# Pull the docker image
print("Pull docker image...")
subprocess.run(["sudo", "docker", "pull", "kishwars/pepper_deepvariant:r0.8"])

# Run PEPPER-Margin-DeepVariant
print("Run pepper-margin-deepvariant...")
subprocess.run([
    "sudo", "docker", "run",
    "-v", f"{input_dir}:{input_dir}",
    "-v", f"{output_dir}:{output_dir}",
    "kishwars/pepper_deepvariant:r0.8",
    "run_pepper_margin_deepvariant", "call_variant",
    "-b", f"{input_dir}/{bam}",
    "-f", f"{input_dir}/{ref}",
    "-o", output_dir,
    "-p", output_prefix,
    "-t", threads,
    "--ont_r10_q20"
])

print("Done.")
