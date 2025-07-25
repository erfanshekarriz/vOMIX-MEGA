logdir = relpath("identify/viral/logs")
tmpd = relpath("identify/viral/tmp")
benchmarks=relpath("identify/viral/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)

n_cores = config['cores']

### Read single fasta file if input
if config['fasta'] != "" and config["module"] == "clustering-sensitive":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
  samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, nowstr)
  fastap = relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa")



##################
# CD-HIT CLUSTER #
##################
rule cdhit_derep:
  name: "CD-HIT --clustering-sensitive [clustering-sensitive.smk]"
  input:
    fastap
  output:
    fa=relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"), 
    clstr=relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr")
  params:
    cdhitparams=config['cdhit-params'],
    outdir=relpath("identify/viral/output/derep"),
    tmpdir=os.path.join(tmpd, "cdhit")
  log: os.path.join(logdir, "clustering/cdhitderep.log")
  conda: "../envs/cd-hit.yml"
  benchmark: os.path.join(benchmarks, "identify/viral_cdhit.log")
  threads: 32
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}
    
    cd-hit -i {input} -o {params.tmpdir}/tmp.fa -T {threads} {params.cdhitparams} &> {log}

    mv {params.tmpdir}/tmp.fa {output.fa}
    mv {params.tmpdir}/tmp.fa.clstr {output.clstr}
    """
