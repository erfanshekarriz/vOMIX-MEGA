logdir = relpath("taxonomy/viral/logs")
tmpd = relpath("taxonomy/viral/tmp")
benchmarks = relpath("taxonomy/viral/benchmarks")
outdir_p = relpath("taxonomy/viral/output/")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)

n_cores = config['max-cores']


### Check if geNomad is run already 
if os.path.exists(relpath("identify/viral/output/classification_summary_vOTUs.csv")):
  genomad_out = relpath("identify/viral/output/classification_summary_vOTUs.csv")
  genomad_f = False
  console.print(Panel.fit(f"[dim] geNomad has already been run. Using its output in '{genomad_out}' for taxonomic annotation.", title = "Warning", subtitle="geNomad taxonomy"))
else:
  genomad_out = relpath("taxonomy/viral/intermediate/genomad/taxonomy.tsv")

### Read single fasta file if input
if config['fasta'] != "":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
  fastap = relpath("identify/viral/output/combined.final.vOTUs.fa")

### MASTER RULE 

rule done_log:
  name: "viral-taxonomy.smk Done. removing tmp files"
  localrule: True
  input:
    relpath("taxonomy/viral/output/merged_taxonomy.csv")
  output:
    os.path.join(logdir, "done.log")
  params:
    tmpdir=tmpd
  log: os.path.join(logdir, "done.log")
  shell:
    """
    rm -rf {params.tmpdir}/*
    touch {output}
    """

### RULES

if genomad_f:
  rule genomad_classify:
    name: "viral-taxonomy.smk geNomad classify"
    input:
      fna=fastap, 
      db=os.path.join(config['genomad-db'], "genomad_db.source")
    output:
      genomad=genomad_out
    params:
      genomadparams=config['genomad-params'],
      dbdir=config['genomad-db'],
      outdir=relpath("taxonomy/viral/intermediate/genomad"),
      tmpdir=os.path.join(tmpd, "genomad")
    log: os.path.join(logdir, "genomad_taxonomy.log")
    benchmark: os.path.join(benchmarks, "genomad_taxonomy.log")
    conda: "../envs/genomad.yml"
    threads: 64
    resources:
      mem_mb=lambda wildcards, attempt, input: 24 * 10**3 * attempt
    shell:
      """
      rm -rf {params.tmpdir} {params.outdir}
      mkdir -p {params.tmpdir} {params.outdir}

      genomad end-to-end \
          {input.fna} \
          {params.tmpdir} \
          {params.dbdir} \
          --threads {threads} \
          --cleanup \
          {params.genomadparams} &> {log}

      mv {params.tmpdir}/* {params.outdir}
      cp {params.outdir}/*_summary/*_virus_summary.tsv {output.genomad}
      rm -rf {params.tmpdir}
      """


rule genomad_taxonomy:
  name: "viral-taxonomy.smk geNomad parse taxonomy"
  localrule: True
  input:
    genomad_out
  output:
    relpath("taxonomy/viral/intermediate/genomad/taxonomy.csv")
  params:
    script="workflow/scripts/genomad_taxonomy.py",
    outdir=relpath("taxonomy/viral/intermediate/genomad"),
    tmpdir=os.path.join(tmpd, "genomad")
  conda: "../envs/ete3.yml"
  log: os.path.join(logdir, "genomad_parse.log")
  benchmark: os.path.join(benchmarks, "genomad_parse.log")
  shell:
    """
    rm -rf {params.tmpdir}/* 
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --input {input} \
        --output {params.tmpdir}/tmp.csv 2> {log}

    mv {params.tmpdir}/tmp.csv {output}
    rm -rf {params.tmpdir}/*

    """


rule phagcn_taxonomy:
  name: "viral-taxonomy.smk PhaGCN phage taxonomy"
  input:
    fna=fastap, 
    db=os.path.join(config['PhaBox2-db'], "genus2hostlineage.pkl")
  output:
    relpath("taxonomy/viral/intermediate/phagcn/taxonomy.tsv")
  params:
    parameters=config['phagcn-params'],
    dbdir=config['PhaBox2-db'],
    outdir=relpath("taxonomy/viral/intermediate/phagcn"),
    tmpdir=os.path.join(tmpd, "phagcn")
  log: os.path.join(logdir, "phagcn.log")
  benchmark: os.path.join(benchmarks, "phagcn.log")
  conda: "../envs/phabox2.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 128 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    phabox2 --task phagcn \
        --contigs {input.fna} \
        --threads {threads} \
        --outpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}/
    cp {params.outdir}/final_prediction/phagcn_prediction.tsv {output}

    rm -rf {params.tmpdir}
    """ 


rule merge_taxonomy:
  name: "viral-taxonomy.smk merge taxonomic classifications"
  localrule: True
  input:
    phagcn=relpath("taxonomy/viral/intermediate/phagcn/taxonomy.tsv"),
    genomad=relpath("taxonomy/viral/intermediate/genomad/taxonomy.csv"),
    contigs=fastap
  output:
    relpath("taxonomy/viral/output/merged_taxonomy.csv")
  params:
    script="workflow/scripts/taxonomy_merge.py",
    outdir=relpath("taxonomy/viral/output/"),
    tmpdir=os.path.join(tmpd, "merge")
  log: os.path.join(logdir, "merge_taxonomy.log")
  benchmark: os.path.join(logdir, "merge_taxonomy.log")
  conda: "../envs/ete3.yml"
  shell:
    """
    rm -rf {params.tmpdir}/* {output}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        --phagcnout {input.phagcn} \
        --genomadout {input.genomad} \
        --contigs {input.contigs} \
        --output {params.tmpdir}/tmp.csv &> {log}
    
    mv {params.tmpdir}/tmp.csv {output}
    """
        
