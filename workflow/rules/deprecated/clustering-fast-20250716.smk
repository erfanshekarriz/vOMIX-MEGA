logdir = relpath("identify/viral/logs")
tmpd = relpath("identify/viral/tmp")
benchmarks=relpath("identify/viral/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)


if isinstance(config['cores'], int):
 n_cores = config['cores']
else:
  console.print(Panel.fit(f"config['cores'] is not an integer: {config['cores']}, you can change the parameter in config/config.yml file", title="Error", subtitle="config['cores'] not integer"))
  sys.exit(1)


### Read single fasta file if input
if config['fasta'] != "" and config["module"] == "clustering-fast":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
  fastap = relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa")

### MASTER RULE

rule done_log:
  name: "clustering-fast.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("identify/viral/intermediate/derep/db.{suffix}"), suffix=["ntf", "ndb"]), 
    relpath("identify/viral/intermediate/derep/blast_out.csv"), 
    relpath("identify/viral/intermediate/derep/ani.tsv"), 
    relpath("identify/viral/output/derep/clusters.tsv"),
    relpath("identify/viral/output/derep/cluster_representatives.txt"),
    relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa")
  output:
    os.path.join(logdir, "clustering-fast-done.log")
  shell: "touch {output}"



### RULES

rule makeblastdb_derep:
  name: "clustering-fast.smk make blast db [--clustering-fast]"
  input: 
    fastap
  output: 
    expand(relpath("identify/viral/intermediate/derep/db.{suffix}"), suffix=["ntf", "ndb"])
  params:
    outdir=relpath("identify/viral/intermediate/derep/"), 
    dbtype='nucl', 
    tmpdir=tmpd
  log: os.path.join(logdir, "clustering/makeblastdb.log")
  benchmark: os.path.join(benchmarks, "makeblastdb.log")
  conda: "../envs/checkv.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    rm -rf {params.tmpdir}/* {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    makeblastdb -in {input} -dbtype {params.dbtype} -out {params.tmpdir}/db &> {log}

    mv {params.tmpdir}/* {params.outdir}
    rm -rf {params.tmpdir}/*

    """

rule megablast_derep:
  name: "clustering-fast.smk megablast [--clustering-fast]"
  input:
    fasta=fastap, 
    dbcheckpoints=expand(relpath("identify/viral/intermediate/derep/db.{suffix}"), suffix=["ntf", "ndb"])
  output:
    relpath("identify/viral/intermediate/derep/blast_out.csv")
  params:
    db=relpath("identify/viral/intermediate/derep/db"),
    outfmt="'6 std qlen slen'",
    maxtargetseqs=10000, 
    tmpdir=tmpd
  log: os.path.join(logdir, "megablastpairwise.log")
  benchmark: os.path.join(benchmarks, "megablastpairwise.log")
  conda: "../envs/checkv.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell: 
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    blastn -query {input.fasta} \
        -db {params.db} \
        -outfmt {params.outfmt} \
        -max_target_seqs {params.maxtargetseqs} \
        -out {params.tmpdir}/tmp.csv \
        -num_threads {threads} &> {log}
    
    mv {params.tmpdir}/tmp.csv {output}
    rm -rf {params.tmpdir}/*

    """
  
rule anicalc_derep:
  name : "clustering-fast.smk calculate ani [--clustering-fast]"
  input:
    relpath("identify/viral/intermediate/derep/blast_out.csv")
  output: 
    relpath("identify/viral/intermediate/derep/ani.tsv")
  params:
    script="workflow/scripts/clust_anicalc.py", 
    tmpdir=tmpd
  log: os.path.join(logdir, "anicalc.log")
  benchmark: os.path.join(benchmarks, "anicalc.log")
  conda: "../envs/checkv.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    rm -rf {params.tmpdir}/* 
    mkdir -p {params.tmpdir}

    python {params.script} \
        -i {input} \
        -o {params.tmpdir}/tmp.tsv &> {log}

    mv {params.tmpdir}/tmp.tsv {output}
    rm -rf {params.tmpdir}/*
    """


rule aniclust_derep:
  name : "clustering-fast.smk cluster [--clustering-fast]"
  input:
    fa=fastap, 
    ani=relpath("identify/viral/intermediate/derep/ani.tsv")
  output:
    tsv= relpath("identify/viral/output/derep/clusters.tsv"),
    reps=relpath("identify/viral/output/derep/cluster_representatives.txt")
  params:
    script="workflow/scripts/clust_ani.py",
    minani=config["vOTU-ani"],
    targetcov=config["vOTU-targetcov"],
    querycov =config["vOTU-querycov"], 
    tmpdir=tmpd
  log: os.path.join(logdir, "aniclust.log")
  benchmark: os.path.join(benchmarks, "aniclust.log")
  conda: "../envs/checkv.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    python {params.script} \
        --fna {input.fa} \
        --ani {input.ani} \
        --out {params.tmpdir}/tmp.tsv \
        --min_ani {params.minani} \
        --min_tcov {params.targetcov} \
        --min_qcov {params.querycov} &> {log}
    mv {params.tmpdir}/tmp.tsv {output.tsv}
    cut -f1 {output.tsv} > {output.reps}

    rm -rf {params.tmpdir}/*
    """


rule filtercontigs_derep:
  name: "clustering-fast.smk filter dereplicated viral contigs"
  input: 
    fna=fastap, 
    reps=relpath("identify/viral/output/derep/cluster_representatives.txt")
  output:
    relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa")
  params:
    outdir=relpath("identify/viral/checkv/output"),
    tmpdir=tmpd
  log: os.path.join(logdir, "filterderep.log")
  conda: "../envs/seqkit-biopython.yml" 
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    seqkit grep {input.fna} -f {input.reps} > {params.tmpdir}/tmp.fa 2> {log}
    mv {params.tmpdir}/tmp.fa {output}

    rm -rf {params.tmpdir}/*
    """
