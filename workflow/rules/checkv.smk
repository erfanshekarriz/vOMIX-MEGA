logdir=relpath("identify/viral/logs")
benchmarks=relpath("identify/viral/benchmarks")
tmpd=relpath("identify/viral/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

if isinstance(config['max-cores'], int):
 n_cores = config['max-cores']
else:
  console.print(Panel.fit(f"config['max-cores'] is not an integer: {config['max-cores']}, you can change the parameter in config/config.yml file", title="Error", subtitle="config['max-cores'] not integer"))
  sys.exit(1)


### Read single fasta file if input
if config['fasta'] != "":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
  samples, assemblies = parse_sample_list(config["workdir"] + config["samplelist"], datadir, outdir, email, nowstr)
  fastap = relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa")



### MASTER RULE

rule done_log:
  name: "checkv-pyhmmer.py Done. removing tmp files"
  localrule: True
  input:
    relpath("identify/viral/output/checkv/viruses.fna"),
    relpath("identify/viral/output/checkv/proviruses.fna"),
    relpath("identify/viral/output/checkv/quality_summary.tsv")
  output:
    os.path.join(logdir, "checkv-done.log")
  shell: "touch {output}"



### RULES

rule checkv:
  name: "checkv.smk CheckV dereplicated contigs"
  input:
    fna=fastap
  output:
    relpath("identify/viral/output/checkv/viruses.fna"),
    relpath("identify/viral/output/checkv/proviruses.fna"),
    relpath("identify/viral/output/checkv/quality_summary.tsv")
  params:
    checkvparams= config['checkv-params'],
    outdir=relpath("identify/viral/output/checkv"),
    tmpdir=os.path.join(tmpd, "checkv"),
    dbdir=config["checkv-database"]
  log: os.path.join(logdir, "checkv.log")
  benchmark: os.path.join(benchmarks, "checkv.log")
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt, input: attempt * 72 * 10**3 * 3
  conda: "../envs/checkv.yml"
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    checkv end_to_end \
        {input.fna} \
        {params.outdir} \
        -d {params.dbdir} \
        -t {threads} \
        {params.checkvparams} 2> {log}

    rm -rf {params.tmpdir}
    """
