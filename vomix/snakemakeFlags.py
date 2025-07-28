class SnakemakeFlags:
    def __init__(self, dry_run, forceall, configfile, unlock, cores, jobs, latency_wait, rerun_incomplete, rerun_triggers, sdm, executor, quiet, snakemake_args):
        self.dry_run = dry_run
        self.forceall = forceall
        self.configfile = configfile
        self.unlock = unlock
        self.cores = cores
        self.jobs = jobs
        self.latency_wait = latency_wait
        self.rerun_incomplete = rerun_incomplete
        self.rerun_triggers = rerun_triggers
        self.sdm = sdm
        self.executor = executor
        self.quiet = quiet
        self.snakemake_args = snakemake_args
