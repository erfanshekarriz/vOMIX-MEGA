logdir=relpath("identify/viral/logs")
benchmarks=relpath("identify/viral/benchmarks")
tmpd=relpath("identify/viral/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

n_cores = config['max-cores']

### Read single fasta file if input
if config['fasta'] != "" and config["module"] == "checkv-pyhmmer":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
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

rule checkv_prodigalgv:
  name: "checkv-pyhmmer.smk CheckV run prodigal-gv"
  input: fastap
  output:
    relpath("identify/viral/output/checkv/tmp/proteins.faa")
  params:
    script="workflow/scripts/parallel_prodigal_gv.py",
    outdir=relpath("identify/viral/output/checkv/tmp"),
    tmpdir=os.path.join(tmpd, "checkv/prodigal-gv")
  log: os.path.join(logdir, "checkv_prodigal-gv.log")
  conda: "../envs/prodigal-gv.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input} \
        -a {params.tmpdir}/tmp.faa \
        -t {threads} &> {log}

    mv {params.tmpdir}/tmp.faa {output}
    rm -rf {params.tmpdir}/*
    """

rule checkv_pyhmmer:
  name: "checkv-pyhmmer.smk CheckV PyHMMER hmmsearch"
  localrule: False
  input:
    faa=relpath("identify/viral/output/checkv/tmp/proteins.faa"), 
    db=os.path.join(config["checkv-database"], "hmm_db/checkv_hmms/{index}.hmm"), 
    checkpoint=expand(os.path.join(config["checkv-database"], "hmm_db/checkv_hmms/{index}.hmm"), index=range(1, 81))
  output:
    relpath("identify/viral/output/checkv/tmp/hmmsearch/{index}.hmmout")
  params:
    script="workflow/scripts/pyhmmer_wrapper.py",
    outdir=relpath("identify/viral/output/checkv/tmp/hmmsearch"),
    tmpdir=os.path.join(tmpd, "checkv/hmmsearch/{index}"), 
    ecutoff=10.0
  log : os.path.join(logdir, "checkv_hmmsearch_{index}.log")
  benchmark: os.path.join(benchmarks, "checkv_hmmsearch_{index}.log")
  conda: "../envs/pyhmmer.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
         --proteins {input.faa} \
         --hmmdb {input.db} \
         --cores {threads} \
         --e_value {params.ecutoff} \
         --tblout {params.tmpdir}/tmp.hmmout &> {log}

    mv {params.tmpdir}/tmp.hmmout {output}
    rm -rf {params.tmpdir}
    """

rule checkv_hmm_merge:
  name: "checkv-pyhmmer.smk CheckV hmmsearch merge"
  localrule: True
  input:
    expand(relpath("identify/viral/output/checkv/tmp/hmmsearch/{index}.hmmout"), index = range(1, 81))
  output:
    relpath("identify/viral/output/checkv/tmp/hmmsearch.txt")
  shell:
    """
    cat {input} > {output}

    """

rule checkv_hmmer_checkpoint:
  name: "checkv-pyhmmer.smk CheckV hmmsearch checkpoint"
  localrule: True
  input:
    relpath("identify/viral/output/checkv/tmp/hmmsearch.txt")
  output:
    relpath("identify/viral/output/checkv/tmp/hmmsearch_checkpoint")
  shell:
    """
    touch {output}
    """

# This rule currently does not operate in tmpdir so that it matches the checkpoint
# Fix this later please thanks!

rule checkv:
  name: "checkv-pyhmmer.smk CheckV dereplicated contigs"
  input:
    checkpoint=relpath("identify/viral/output/checkv/tmp/hmmsearch_checkpoint"),
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
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: attempt * 72 * 10**3
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
