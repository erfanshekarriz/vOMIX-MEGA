conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
# conda update -n base -c defaults conda (taking a long time to finish, so skipping for now)
conda create -n vomix -c conda-forge snakemake=8.25.5 biopython=1.84 -y 
source activate base
# conda activate vomix