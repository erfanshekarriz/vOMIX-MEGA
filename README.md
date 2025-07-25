# vOMIX-MEGA 

vOMIX-MEGA is a reproducible, scalable, and fast viral metagenomic pipeline with rigorously benchmarked backing on its results. It is built on a snakemake backend, can be containerized, and is ready for cloud deployment.


# Quick Start 

**1.1 Install the vOMIX-MEGA base environment:**

```bash
# Set channel priority to strict before running vOMIX-MEGA to ensure reproducibility [IMPORTANT]
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Update conda for snakemake compatibility
conda update -n base -c defaults conda

# Install base environment
conda create -n vomix -c conda-forge snakemake=8.25.5 biopython=1.84 pip=25.1.1 -y # does not include cluster execution plugs. See more at https://snakemake.github.io/snakemake-plugin-catalog/index.html
conda activate vomix

# Verify the snakemake backend is running
snakemake -v
```

**1.2 Download and install vOMIX-MEGA from the GitHub repository:**

```bash
# Clone from GitHub
git clone https://github.com/holab-hku/vOMIX-MEGA
cd vOMIX-MEGA

# Ensure conda env is activated
conda info --envs | grep '^vomix' # you do not want to pip install into your main environment

# Install vOMIX-MEGA using pip 
pip install .

# Verify installation 
vomix -h 
```

**1.3 Test Viral Contig Identification using Sample Data**
```bash
vomix viral-identify --outdir test_res --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --splits 8 -j 4 --latency-wait 20
```


# Wiki

For the full documentation on inputs, outputs, configurations, and modules of vOMIX-MEGA, please visit our Wiki page on github at https://github.com/holab-hku/vOMIX-MEGA/wiki ! 




