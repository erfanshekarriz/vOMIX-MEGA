from vomix.module import Module

class PreProcessingModule(Module):
    # snakemake --config module="preprocess" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "preprocess"
    def __init__(self, decontam_host=False, hasOptions=False, dwnld_params=None, pigz_paramss=None, fastp_params=None, hostile_params=None, hostile_aligner=None, aligner_params=None, index_path=None):
        self.decontam_host = decontam_host
        self.hasOptions = hasOptions
        self.dwnld_params = dwnld_params
        self.pigz_params = pigz_paramss
        self.fastp_params = fastp_params
        self.hostile_params = hostile_params 
        self.hostile_aligner = hostile_aligner
        self.aligner_params = aligner_params
        self.index_path = index_path
        # dwnld-only

class AssemblyCoAssemblyModule(Module):
    # snakemake --config module="assembly" assembler="megahit" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "assembly"
    def __init__(self, assembler=False, hasOptions=False, megahit_minlen=300, megahit_params=None, spades_params=None, spades_memory=250):
        self.assembler = assembler
        self.hasOptions = hasOptions
        self.megahit_minlen = megahit_minlen
        self.megahit_params = megahit_params
        self.spades_params = spades_params
        self.spades_memory = spades_memory
    
class ViralIdentifyModule(Module):
    # snakemake --config module="viral-identify" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "viral-identify"
    def __init__(self, hasOptions=False, contig_minlen=0, genomad_db=None, genomad_minlen=1500, genomad_params=None, genomad_cutoff=0.7, checkv_original=False, checkv_params=None, checkv_database=None, clustering_fast=True, cdhit_params=None, vOTU_ani=95, vOTU_targetcov=85, vOTU_querycov=0):
        self.hasOptions = hasOptions
        self.contig_minlen = contig_minlen
        self.genomad_db = genomad_db
        self.genomad_minlen = genomad_minlen
        self.genomad_params = genomad_params
        self.genomad_cutoff = genomad_cutoff
        self.checkv_original = checkv_original
        self.checkv_params = checkv_params
        self.checkv_database = checkv_database
        self.clustering_fast = clustering_fast
        self.cdhit_params = cdhit_params
        self.vOTU_ani = vOTU_ani
        self.vOTU_targetcov = vOTU_targetcov
        self.vOTU_querycov = vOTU_querycov
    
class ViralTaxonomyModule(Module):
    # snakemake --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "viral-taxonomy"
    def __init__(self, hasOptions=False, viphogs_hmmeval=0.01, viphogs_prop=0.06, PhaBox2_db=None, phagcn_minlen=1500, phagcn_params=None, diamond_params=None, genomad_db=None, genomad_params=None):
        self.hasOptions = hasOptions
        self.viphogs_hmmeval = viphogs_hmmeval
        self.viphogs_prop = viphogs_prop
        self.PhaBox2_db = PhaBox2_db
        self.phagcn_minlen = phagcn_minlen
        self.phagcn_params = phagcn_params
        self.diamond_params = diamond_params
        self.genomad_db = genomad_db
        self.genomad_params = genomad_params

class ViralHostModule(Module):
    # snakemake --config module="viral-host" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "viral-host"
    def __init__(self, hasOptions=False, CHERRY_params=None, PhaTYP_params=None, iphop_cutoff=90, iphop_params=None):
        self.hasOptions = hasOptions
        self.CHERRY_params = CHERRY_params
        self.PhaTYP_params = PhaTYP_params
        self.iphop_cutoff = iphop_cutoff
        self.iphop_params = iphop_params

class ViralCommunityModule(Module):
    # snakemake --config module="viral-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
    name = "viral-community"
    def __init__(self, hasOptions=False, mpa_indexv=None, mpa_params=None):
        self.hasOptions = hasOptions
        self.mpa_indexv = mpa_indexv
        self.mpa_params = mpa_params    


class ViralAnnotateModule(Module):
    # snakemake --config module="viral-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "viral-annotate"
    def __init__(self, hasOptions=False, eggNOG_params=None, PhaVIP_params=None):
        self.hasOptions = hasOptions
        self.eggNOG_params = eggNOG_params
        self.PhaVIP_params = PhaVIP_params

class ProkaryoticCommunityModule(Module):
    # snakemake --config module="prok-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
    name = "prok-community"
    def __init__(self, hasOptions=False, mpa_params=None, mpa_indexv=None):
        self.hasOptions = hasOptions
        self.mpa_params = mpa_params
        self.mpa_indexv = mpa_indexv

# TBD module
class ProkaryoticBinningModule(Module):
    name = "prok-binning"

class ProkaryoticAnnotateModule(Module):
    # snakemake --config module="prok-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "prok-annotate"
    def __init__(self, hasOptions=False):
        self.hasOptions = hasOptions

class EndToEndModule(Module):
    # snakemake --config module="end-to-end" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4
    name = "end-to-end"
    def __init__(self, hasOptions=False):
        self.hasOptions = hasOptions

class ClusterFastModule(Module):
    # snakemake --config module="cluster-fast" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "cluster-fast"
    def __init__(self, hasOptions=False, clustering_fast=True, cdhit_params=None, vOTU_ani=95, vOTU_targetcov=85, vOTU_querycov=0):
        self.hasOptions = hasOptions
        self.clustering_fast = clustering_fast
        self.cdhit_params = cdhit_params
        self.vOTU_ani = vOTU_ani
        self.vOTU_targetcov = vOTU_targetcov
        self.vOTU_querycov = vOTU_querycov

class CheckVPyHMMERModule(Module):
    # snakemake --config module="checkv-pyhmmer" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "checkv-pyhmmer"
    def __init__(self, hasOptions=False, checkv_original=False, checkv_params=None, checkv_database=None):
        self.hasOptions = hasOptions
        self.checkv_original = checkv_original
        self.checkv_params = checkv_params
        self.checkv_database = checkv_database

class SetupDatabaseModule(Module):
    # snakemake --config module="setup-database" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "setup-database"
    def __init__(self, hasOptions=False, PhaBox2_db=None, genomad_db=None, checkv_db=None, eggNOG_db=None, eggNOG_db_params=None, virsorter2_db=None, iphop_db=None, humann_db=None):
        self.hasOptions = hasOptions
        self.PhaBox2_db = PhaBox2_db
        self.genomad_db = genomad_db
        self.checkv_db = checkv_db
        self.eggNOG_db = eggNOG_db
        self.eggNOG_db_params = eggNOG_db_params
        self.virsorter2_db = virsorter2_db
        self.iphop_db = iphop_db
        self.humann_db = humann_db
