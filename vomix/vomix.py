import click
import sys
import logging
import os
from vomix.snakemakeFlags import SnakemakeFlags
from vomix.vomix_actions import vomix_actions
from vomix.modules import PreProcessingModule, AssemblyCoAssemblyModule, ViralIdentifyModule, ViralTaxonomyModule, ViralHostModule, ViralCommunityModule, ViralAnnotateModule, ProkaryoticCommunityModule, ProkaryoticBinningModule, ProkaryoticAnnotateModule, EndToEndModule, ClusterFastModule, CheckVPyHMMERModule, SetupDatabaseModule

logging.basicConfig(level=logging.INFO)

modules_list = ['assembly','checkv-pyhmmer','checkv'
,'clustering-fast','clustering-sensitive','host-cherry','host','preprocess-decontam','preprocess','prok-annotate','prok-binning','prok-community','refilter-genomad'
,'setup-database','symlink','viral-annotate','viral-benchmark','viral-binning','viral-community','viral-host','viral-identify','viral-refilter','viral-taxonomy']

def useLastOptionsCheck(ctx, param, value):
    i = 2
    totalParams = len(ctx.command.params)
    if value == True:
        while i < totalParams:
            ctx.command.params[i].prompt_required = False
            ctx.command.params[i].required = False
            i += 1

# common options decorator
def common_options(function):
    function = click.option('--workdir', default=None, required=False, help = 'Set the working directory for Snakefile (We recommend not changing this)')(function)
    function = click.option('--outdir', default=None, required=False, help = 'Select the output directory for hierarchal results formatting || default: "./results"')(function)
    function = click.option('--datadir', default=None, required=False, help = '')(function)
    function = click.option('--samplelist', default=None, required=False, help = '')(function)
    function = click.option('--fasta', default=None, required=False, help = '')(function)
    function = click.option('--fastadir', default=None, required=False, help = '')(function)
    function = click.option('--sample-name', default=None, required=False, help = '')(function)
    function = click.option('--assembly-ids', default=None, required=False, help = '')(function)
    function = click.option('--latest-run', default=None, required=False, help = '')(function)
    function = click.option('--splits', default=0, required=False, help = 'Splits data into N chunks to reduce memory usage wherever possible || default: 0')(function)
    function = click.option('--viral-binning', is_flag=True, default=False, required=False, help = '')(function)
    function = click.option('--intermediate', is_flag=True, default=False, required=False, help = 'Flag to keep LARGE intermediate files generated during analysis || default: False')(function)
    function = click.option('--setup-database', is_flag=True, default=False, required=False, help = '')(function)
    function = click.option('--max-cores', default=4, required=False, help = '')(function)
    function = click.option('--email', default=None, required=False, help = '')(function)
    function = click.option('--NCBI-API-key', default=None, required=False, help = '')(function)
    function = click.option('--custom-config', default=None, required=False, help = 'Path to your custom config.yml')(function)
    return function

# common snakemake options decorator
def snakemake_options(function):
    function = click.option('--dry-run', '--dryrun', '-n', required=False, default=False, flag_value=True, help = '[Snakemake]')(function)
    function = click.option('--forceall', '-F', required=False, default=False, flag_value=True, help = '[Snakemake]')(function)
    function = click.option('--configfile', default=None, required=False, help = '[Snakemake]')(function)
    function = click.option('--unlock', required=False, default=False, flag_value=True, help = '[Snakemake]')(function)
    function = click.option('--cores', '-c', default=1, required=False, help = '[Snakemake]')(function)
    function = click.option('--jobs', '-j', default=4, required=False, help = '[Snakemake]')(function)
    function = click.option('--latency-wait', default=20, required=False, help = '[Snakemake]')(function)
    function = click.option('--rerun-incomplete', '-ri', required=False, default=False, flag_value=True, help = '[Snakemake]')(function)
    function = click.option('--rerun-triggers', required=False, default=None, help = '[Snakemake]')(function)
    function = click.option('--sdm', required=False, default=None, help = '[Snakemake]')(function)
    function = click.option('--executor', '-e', required=False, default=None, help = '[Snakemake]')(function)
    function = click.option('--quiet', required=False, default=False, flag_value=True, help = '[Snakemake]')(function)

    function = click.option('--add-args', required=False, default=None, help = '[Snakemake]')(function)

    return function


def setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config):
    module_obj.workdir = workdir
    module_obj.outdir = outdir 
    module_obj.datadir = datadir
    module_obj.samplelist = samplelist
    module_obj.fasta = fasta
    module_obj.fastadir = fastadir
    module_obj.sample_name = sample_name
    module_obj.assembly_ids = assembly_ids
    module_obj.latest_run = latest_run
    module_obj.splits = splits
    module_obj.viral_binning = viral_binning
    module_obj.intermediate = intermediate
    module_obj.setup_database = setup_database
    module_obj.max_cores = max_cores
    module_obj.email = email
    module_obj.NCBI_API_key = ncbi_api_key
    module_obj.custom_config = custom_config

    return module_obj

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    """
    vOMIX-MEGA - a reproducible, scalable, and fast viral metagenomic pipeline
    """

# @cli.command(
#         'activate',
#         context_settings=dict(ignore_unknown_options=True),
#         short_help='Set up vomix conda environment'
#     )
# def activate():
#     output = vomix_actions.env_setup_script()
#     logging.info(f"Output: {output}")
#     message = "Completed vomix environment creation"
#     logging.info(message)


@cli.command(
    'preprocess',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Pre-processing module'
)
@common_options
@click.option('--decontam-host', default=True, required=False, help='')
@click.option('--dwnld-params', required=False, default=None, help='Parameters for fasterq-dump for downloading from sra tools https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump || default: ""')
@click.option('--pigz-params', required=False, default=None, help='Parameters of pigz for compressing downloaded fastq files https://github.com/madler/pigz || default: ""')
@click.option('--fastp-params', required=False, default=None, help='Parameters to pass on fastp software https://github.com/OpenGene/fastp || default: ""')
@click.option('--hostile-params', required=False, default=None, help='Parameters for hostile decontamination https://github.com/bede/hostile|| default: ""')
@click.option('--hostile-aligner', required=False, default=None, help='Which mapper to use for host decontamination- bowtie2 or minimap2 (recommended) || default: "minimap2"')
@click.option('--aligner-params', required=False, default=None, help='PLEASE DO NOT change the -x sr for minimap2 to make sure it can accurately map short reads || default: "-x sr"')
@click.option('--index-path', required=False, default=None, help='Path to host contamination || default: "./workflow/database/hostile/human-t2t-hla.fa.gz"')
@snakemake_options
def run_preprocess(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, decontam_host, dwnld_params, pigz_params, fastp_params, hostile_params, hostile_aligner, aligner_params, index_path, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: preprocess")
        logging.info(f"decontamHost: {decontam_host}, outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        # logging.info(f"Run with snakemake flags: " + "dry_run:" + str(dry_run) + ", forceall:" + str(forceall) + ", configfile:" + str(configfile))

        module_obj = PreProcessingModule()
        module_obj.name = "preprocess"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        module_obj.decontam_host = decontam_host

        # optional params
        if dwnld_params:
            module_obj.dwnld_params = dwnld_params
            module_obj.hasOptions = True
        if pigz_params:
            module_obj.pigz_params = pigz_params  
            module_obj.hasOptions = True
        if fastp_params:
            module_obj.fastp_params = fastp_params 
            module_obj.hasOptions = True
        if hostile_params:
            module_obj.hostile_params = hostile_params
            module_obj.hasOptions = True
        if hostile_aligner:
            module_obj.hostile_aligner = hostile_aligner
            module_obj.hasOptions = True
        if aligner_params:
            module_obj.aligner_params = aligner_params
            module_obj.hasOptions = True
        if index_path:
            module_obj.index_path = index_path
            module_obj.hasOptions = True

        # snakemake options 
        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("preprocess", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'assembly',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Assembly & Co-assembly module'
)
@common_options
@click.option('--assembler', default=None, required=False, help = '')
@click.option('--megahit-minlen', required=False, default=None, help = 'Minimum length for MEGAHIT to use for contig building || default: 300')
@click.option('--megahit-params', required=False, default=None, help = 'Extra parameters to hand off to MEGAHIT software https://github.com/voutcn/megahit || default: "--prune-level 3"')
@click.option('--spades-params', required=False, default=None, help = 'Parameters to pass on fastp software https://github.com/OpenGene/fastp || default: "--meta"')
@click.option('--spades-memory', required=False, default=None, help = 'Parameters for hostile decontamination https://github.com/bede/hostile|| default: 250 NUM')
@snakemake_options
def run_assembly(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, assembler, megahit_minlen, megahit_params, spades_params, spades_memory, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: assembly")
        logging.info(f"assembler: {assembler}, outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = AssemblyCoAssemblyModule()
        module_obj.name = "assembly"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        module_obj.assembler = assembler

        if megahit_minlen:
            module_obj.megahit_minlen = megahit_minlen
            module_obj.hasOptions = True
        if megahit_params: 
            module_obj.megahit_params = megahit_params
            module_obj.hasOptions = True
        if spades_params:
            module_obj.spades_params = spades_params
            module_obj.hasOptions = True
        if spades_memory:
            module_obj.spades_memory = spades_memory
            module_obj.hasOptions = True

        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("assembler", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'viral-identify',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Identify module'
)
@common_options
@click.option('--contig-minlen', required=False, default=None, help = 'Minimum contig length to filter BEFORE viral identification || default: 0 [INT]')
@click.option('--genomad-db', required=False, default=None, help = 'Path to geNomad databases || default: "workflow/database/genomad" [STR]')
@click.option('--genomad-minlen', required=False, default=None, help = 'Minimum viral contig length for the geNomad cutoff || default: 1500 [INT]')
@click.option('--genomad-params', required=False, default=None, help = 'Additional parameters to hand off to geNomad\'s analysis || default: "" [STR]')
@click.option('--genomad-cutoff', required=False, default=None, help = 'Parameters for hostile decontamination https://github.com/bede/hostile|| default: 0.7 [INT]')
@click.option('--checkv-original', required=False, default=None, help = 'Flag to use CheckV original instead of the much faster version in vOMIX-MEGA, CheckV-PyHMMER. || default: False [True or False]')
@click.option('--checkv-params', required=False, default=None, help = 'Additional parameters to pass on to CheckV. Read more at https://bitbucket.org/berkeleylab/CheckV/src || default: "" [STR]')
@click.option('--checkv-database', required=False, default=None, help = 'Path to CheckV\'s database || default: "workflow/database/checkv" [STR]')
@click.option('--clustering-fast', required=False, default=None, help = 'Flag to run fast clustering using CheckV\'s MEGABLAST approach. If set to False, CD-HIT will be used. Proceed with caution as it can be extremely slow at large sequence numbers. || default: True [True or False]')
@click.option('--cdhit-params', required=False, default=None, help = 'Additional parameters to pass on to CD-HIT if clustering-fast is set to False. Read more at https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki || default: "-c 0.95 -aS 0.85 -d 400 -M 0 -n 5" [STR]')
@click.option('--vOTU-ani', required=False, default=None, help = 'Minimum average nucleotide identity for fast clustering algorithm of viral contigs || default: 95 [INT]')
@click.option('--vOTU-targetcov', required=False, default=None, help = 'Minimum target coverage for fast clustering algorithm of viral contigs || default: 85 [NUM]')
@click.option('--vOTU-querycov', required=False, default=None, help = 'Minimum query coverage for fast clustering algorithm of viral contigs || default: 0 [NUM]')
@snakemake_options
def run_viral_identify(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, contig_minlen, genomad_db, genomad_minlen, genomad_params, genomad_cutoff, checkv_original, checkv_params, checkv_database, clustering_fast, cdhit_params, votu_ani, votu_targetcov, votu_querycov, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: viral-identify")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ViralIdentifyModule()
        module_obj.name = "viral-identify"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        if contig_minlen:
            module_obj.contig_minlen = contig_minlen
            module_obj.hasOptions = True
        if genomad_db:  
            module_obj.genomad_db = genomad_db
            module_obj.hasOptions = True
        if genomad_minlen:
            module_obj.genomad_minlen = genomad_minlen
            module_obj.hasOptions = True
        if genomad_params:  
            module_obj.genomad_params = genomad_params
            module_obj.hasOptions = True    
        if genomad_cutoff:
            module_obj.genomad_cutoff = genomad_cutoff
            module_obj.hasOptions = True
        if checkv_original:
            module_obj.checkv_original = checkv_original
            module_obj.hasOptions = True
        if checkv_params:
            module_obj.checkv_params = checkv_params
            module_obj.hasOptions = True
        if checkv_database: 
            module_obj.checkv_database = checkv_database
            module_obj.hasOptions = True
        if clustering_fast:
            module_obj.clustering_fast = clustering_fast
            module_obj.hasOptions = True
        if cdhit_params:
            module_obj.cdhit_params = cdhit_params
            module_obj.hasOptions = True
        if votu_ani:
            module_obj.vOTU_ani = votu_ani
            module_obj.hasOptions = True
        if votu_targetcov:
            module_obj.vOTU_targetcov = votu_targetcov
            module_obj.hasOptions = True
        if votu_querycov:
            module_obj.vOTU_querycov = votu_querycov
            module_obj.hasOptions = True

        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-identify", module_obj, snakemake_obj)
        logging.info(f"End module run")


@cli.command(
    'viral-taxonomy',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Taxonomy module'
)
@common_options
@click.option('--viphogs-hmmeval', required=False, default=None, help = 'Minimum e value for ViPhogs hmms to be considered a hit || default: 0.01 [NUM]')
@click.option('--viphogs-prop', required=False, default=None, help = 'Minimum proportion of annotated genes required for taxonomic assignment || default: 0.6 [NUM]')
@click.option('--PhaBox2-db', required=False, default=None, help = 'Path to phabox database directory || default: "workflow/database/phabox_db_v2"')
@click.option('--phagcn-minlen', required=False, default=None, help = 'Minimum contig length to filter before PaGCN taxonomy annotation || default: 1500 [INT]')
@click.option('--phagcn-params', required=False, default=None, help = 'Additional parameters to pass on to PhaGCN || default: "" [STR]')
@click.option('--diamond-params', required=False, default=None, help = 'Parameters for taxonomic classification using diamond || default: "--query-cover 50 --subject-cover 50 --evalue 1e-5 --max-target-seqs 1000" [INT] || WARNING: strongly recommend not changing this as it has been extensively tested by https://doi.org/10.1038/s41564-021-00928-6')
@click.option('--genomad-db', required=False, default=None, help = 'Path to geNomad database directory || default: "workflow/database/genomad" [STR]')
@click.option('--genomad-params', required=False, default=None, help = 'Additional parameters to pass on to geNomad || default: "--enable-score-calibration --relaxed" [INT]')
@snakemake_options
def run_viral_taxonomy(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, viphogs_hmmeval, viphogs_prop, phabox2_db, phagcn_minlen, phagcn_params, diamond_params, genomad_db, genomad_params, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: viral-taxonomy")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = ViralTaxonomyModule()
        module_obj.name = "viral-taxonomy"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        if viphogs_hmmeval:
            module_obj.viphogs_hmmeval = viphogs_hmmeval
            module_obj.hasOptions = True
        if viphogs_prop:
            module_obj.viphogs_prop = viphogs_prop
            module_obj.hasOptions = True
        if phabox2_db: 
            module_obj.PhaBox2_db = phabox2_db
            module_obj.hasOptions = True
        if phagcn_minlen:
            module_obj.phagcn_minlen = phagcn_minlen
            module_obj.hasOptions = True
        if phagcn_params:
            module_obj.phagcn_params = phagcn_params
            module_obj.hasOptions = True
        if diamond_params:  
            module_obj.diamond_params = diamond_params
            module_obj.hasOptions = True
        if genomad_db:
            module_obj.genomad_db = genomad_db
            module_obj.hasOptions = True
        if genomad_params:
            module_obj.genomad_params = genomad_params
            module_obj.hasOptions = True

        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-taxonomy", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'viral-host',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Host module'
)
@common_options
@click.option('--CHERRY-params', required=False, default=None, help = 'Parameters to pass on to CHERRY for host identification. Read more at https://phage.ee.cityu.edu.hk/wiki || default: "" [STR]')
@click.option('--PhaTYP-params', required=False, default=None, help = 'Parameters to pass on to MaxBin2 for lifestyle identification. Read more at https://phage.ee.cityu.edu.hk/wiki || default: "" [STR]')
@click.option('--iphop-cutoff', required=False, default=None, help = 'The number of correct host predictions was evaluated for 3 different score cutoffs corresponding to 20%, 10%, and 5% estimated FDR || default: 90 [NUM]')
@click.option('--iphop-params', required=False, default=None, help = 'Parameters to pass on to iPhOp for consensus host analysis. Read more at https://bitbucket.org/srouxjgi/iphop/src || default: "" [STR]')
@snakemake_options
def run_viral_host(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, cherry_params, phatyp_params, iphop_cutoff, iphop_params, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: viral-host")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = ViralHostModule()
        module_obj.name = "viral-host"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        if cherry_params:
            module_obj.CERRY_params = cherry_params
            module_obj.hasOptions = True
        if phatyp_params:
            module_obj.PhaTYP_params = phatyp_params
            module_obj.hasOptions = True
        if iphop_cutoff:
            module_obj.iphop_cutoff = iphop_cutoff
            module_obj.hasOptions = True
        if iphop_params:
            module_obj.iphop_params = iphop_params
            module_obj.hasOptions = True
        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-host", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'viral-community',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Community module'
)
@common_options
@click.option('--mpa-indexv', required=False, default=None, help = 'The version of the MetaPhlAn4 database to download || default: "mpa_vOct22_CHOCOPhlAnSGB_202212" [STR]')
@click.option('--mpa-params', required=False, default=None, help = 'Additional parameters to pass on to metaphlan function. See https://huttenhower.sph.harvard.edu/metaphlan/ for more. || default: "--ignore_eukaryotes" [STR]')
@snakemake_options
def run_viral_community(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, mpa_indexv, mpa_params, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: viral-community")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ViralCommunityModule()
        module_obj.name = "viral-community"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        if mpa_indexv:
            module_obj.mpa_indexv = mpa_indexv
            module_obj.hasOptions = True
        if mpa_params:
            module_obj.mpa_params = mpa_params
            module_obj.hasOptions = True
        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-community", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'viral-annotate',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Viral Annotate module'
)
@common_options
@click.option('--eggNOG-params', required=False, default=None, help = 'Parameters for running eggNOG-mapper v2. See more at https://github.com/eggnogdb/eggnog-mapper/wiki || default: "-m diamond --hmm_evalue 0.001 --hmm_score 60 --query-cover 20 --subject-cover 20 --tax_scope auto --target_orthologs all --go_evidence non-electronic --report_orthologs" [INT]')
@click.option('--PhaVIP-params', required=False, default=None, help = 'Minimum contig length to filter BEFORE viral identification || default: "" [STR]')
@snakemake_options
def run_viral_annotate(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, eggnog_params, phavip_params, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: viral-annotate")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ViralAnnotateModule()
        module_obj.name = "viral-annotate"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        if eggnog_params:
            module_obj.eggNOG_params = eggnog_params
            module_obj.hasOptions = True
        if phavip_params:
            module_obj.PhaVIP_params = phavip_params
            module_obj.hasOptions = True
        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("viral-annotate", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'prok-community',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Prokaryotic Community module'
)
@common_options
@click.option('--mpa-params', required=False, default=None, help = 'Parameters for metaphlan function. See more at https://huttenhower.sph.harvard.edu/metaphlan/ || default: "--ignore_eukaryotes" [STR]')
@click.option('--mpa-indexv', required=False, default=None, help = 'Database version for metaphlan to use. See more at https://huttenhower.sph.harvard.edu/metaphlan/ || default: "mpa_vOct22_CHOCOPhlAnSGB_202212" [STR]')
@snakemake_options
def run_prok_community(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, mpa_params, mpa_indexv, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: prok-community")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ProkaryoticCommunityModule()
        module_obj.name = "prok-community"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        if mpa_params:
            module_obj.mpa_params = mpa_params
            module_obj.hasOptions = True    
        if mpa_indexv:
            module_obj.mpa_indexv = mpa_indexv
            module_obj.hasOptions = True
        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("prok-community", module_obj, snakemake_obj)
        logging.info(f"End module run")

# TBD ProkaryoticBinningModule


@cli.command(
    'prok-annotate',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Prokaryotic Annotate module'
)
@common_options
@snakemake_options
def run_prok_annotate(workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: prok-annotate")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = ProkaryoticAnnotateModule()
        module_obj.name = "prok-annotate"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)
        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("prok-annotate", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'end-to-end',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the End-To-End module'
)
@common_options
@snakemake_options
def run_end_to_end(workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: end-to-end")
        logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")
        
        module_obj = EndToEndModule()
        module_obj.name = "end-to-end"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)
        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("end-to-end", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'cluster-fast',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Cluster Fast module'
)
@common_options
@click.option('--clustering-fast', required=False, default=None, help = 'Flag to run fast clustering using CheckV\'s MEGABLAST approach. If set to False, CD-HIT will be used. Proceed with caution as it can be extremely slow at large sequence numbers. || default: True [True or False]')
@click.option('--cdhit-params', required=False, default=None, help = 'Additional parameters to pass on to CD-HIT if clustering-fast is set to False. Read more at https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki || default: "-c 0.95 -aS 0.85 -d 400 -M 0 -n 5" [STR]')
@click.option('--vOTU-ani', required=False, default=None, help = 'Minimum average nucleotide identity for fast clustering algorithm of viral contigs || default: 95 [INT]')
@click.option('--vOTU-targetcov', required=False, default=None, help = 'Minimum target coverage for fast clustering algorithm of viral contigs || default: 85 [NUM]')
@click.option('--vOTU-querycov', required=False, default=None, help = 'Minimum query coverage for fast clustering algorithm of viral contigs || default: 0 [NUM]')
@snakemake_options
def run_cluster_fast(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, clustering_fast, cdhit_params, votu_ani, votu_targetcov, votu_querycov, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: cluster-fast")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = ClusterFastModule()
        module_obj.name = "cluster-fast"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        if clustering_fast:
            module_obj.clustering_fast = clustering_fast
            module_obj.hasOptions = True
        if cdhit_params:
            module_obj.cdhit_params = cdhit_params
            module_obj.hasOptions = True
        if votu_ani:
            module_obj.vOTU_ani = votu_ani
            module_obj.hasOptions = True
        if votu_targetcov:
            module_obj.votu_targetcov = votu_targetcov
            module_obj.hasOptions = True
        if votu_querycov:
            module_obj.votu_querycov = votu_querycov
            module_obj.hasOptions = True
        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("cluster-fast", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'checkv-pyhmmer',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the CheckV PyHMMER module'
)
@common_options
@click.option('--checkv-original', required=False, default=None, help = 'Flag to use CheckV original instead of the much faster version in vOMIX-MEGA, CheckV-PyHMMER. || default: False [True or False]')
@click.option('--checkv-params', required=False, default=None, help = 'Additional parameters to pass on to CheckV. Read more at https://bitbucket.org/berkeleylab/CheckV/src || default: "" [STR]')
@click.option('--checkv-database', required=False, default=None, help = 'Path to CheckV\'s database || default: "workflow/database/checkv" [STR]')
@snakemake_options
def run_checkv_pyhmmer(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, checkv_original, checkv_params, checkv_database, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: checkv-pyhmmer")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = CheckVPyHMMERModule()
        module_obj.name = "checkv-pyhmmer"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        if checkv_original:
            module_obj.checkv_original = checkv_original
            module_obj.hasOptions = True
        if checkv_params:
            module_obj.checkv_params = checkv_params
            module_obj.hasOptions = True
        if checkv_database:
            module_obj.checkv_database = checkv_database
            module_obj.hasOptions = True
        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("checkv-pyhmmer", module_obj, snakemake_obj)
        logging.info(f"End module run")

@cli.command(
    'setup-database',
    context_settings=dict(ignore_unknown_options=True),
    short_help='Run the Setup Database module'
)
@common_options
@click.option('--PhaBox2-db', required=False, default=None, help = 'Path to PhaBox2 database for download || default: "workflow/database/phabox_db_v2" [STR]')
@click.option('--genomad-db', required=False, default=None, help = 'Path to geNomad database for download || default: "workflow/database/genomad" [STR]')
@click.option('--checkv-db', required=False, default=None, help = 'Path to CheckV database for download || default: "workflow/database/phabox_db_v2" [STR]')
@click.option('--eggNOG-db', required=False, default=None, help = 'Path to eggNOG v2 database for download || default: "workflow/database/eggNOGv2" [STR]')
@click.option('--eggNOG-db-params', required=False, default=None, help = 'Parameters for downloading eggNOG v2 database || default: "" [STR]')
@click.option('--virsorter2-db', required=False, default=None, help = 'Path to VirSorter2 database for download || default: "workflow/database/virsorter2" [STR]')
@click.option('--iphop-db', required=False, default=None, help = 'Path to iPHoP database for download || default: "workflow/database/iphop/Aug_2023_pub_rw" [STR]')
@click.option('--humann-db', required=False, default=None, help = 'Path to HUMAnN3 databases for download || default: "workflow/database/humann" [STR]')
@snakemake_options
def run_setup_database(workdir, outdir, datadir, samplelist, custom_config, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, phabox2_db, genomad_db, checkv_db, eggnog_db, eggnog_db_params, virsorter2_db, iphop_db, humann_db, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args):
        logging.info(f"Running module: setup-database")
        logging.info(f"fasta: {fasta}, outdir: {outdir}")
        
        module_obj = SetupDatabaseModule()
        module_obj.name = "setup-pyhdatabasemmer"
        # Set the attributes of the module object
        module_obj = setOptions(module_obj, workdir, outdir, datadir, samplelist, fasta, fastadir, sample_name, assembly_ids, latest_run, splits, viral_binning, intermediate, setup_database, max_cores, email, ncbi_api_key, custom_config)

        if phabox2_db:  
            module_obj.PhaBox2_db = phabox2_db
            module_obj.hasOptions = True
        if genomad_db:
            module_obj.genomad_db = genomad_db
            module_obj.hasOptions = True
        if checkv_db:
            module_obj.checkv_db = checkv_db
            module_obj.hasOptions = True
        if eggnog_db:
            module_obj.eggNOG_db = eggnog_db
            module_obj.hasOptions = True
        if eggnog_db_params:
            module_obj.eggNOG_db_params = eggnog_db_params
            module_obj.hasOptions = True
        if virsorter2_db:   
            module_obj.virsorter2_db = virsorter2_db
            module_obj.hasOptions = True
        if iphop_db:
            module_obj.iphop_db = iphop_db
            module_obj.hasOptions = True
        if humann_db:
            module_obj.humann_db = humann_db
            module_obj.hasOptions = True

        snakemake_obj = SnakemakeFlags(dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, add_args)

        vomix_actions_instance = vomix_actions()
        vomix_actions_instance.run_module("setup-database", module_obj, snakemake_obj)
        logging.info(f"End module run")

