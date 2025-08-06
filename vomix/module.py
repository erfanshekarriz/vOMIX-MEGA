class Module:
    def __init__(self, module, workdir, outdir, datadir, samplelist, fasta="", fastadir="", sample_name="", assembly_ids="", latest_run="", splits=0, viral_binning=False, keep_intermediates=False, setup_database=True, max_cores=4, email="", NCBI_API_key="", custom_config=None):
        self.module = module
        self.workdir = workdir
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.fasta = fasta
        self.fastadir = fastadir
        self.sample_name = sample_name
        self.assembly_ids = assembly_ids
        self.latest_run = latest_run
        self.splits = splits
        self.viral_binning = viral_binning
        self.keep_intermediates = keep_intermediates
        self.setup_database = setup_database
        self.max_cores = max_cores
        self.email = email
        self.NCBI_API_key = NCBI_API_key
        self.custom_config = custom_config
