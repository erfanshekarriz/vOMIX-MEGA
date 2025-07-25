logdir = relpath("host/logs")
tmpd = relpath("host/tmp")
benchmarks = relpath("host/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)


### Read single fasta file if input
if config['fasta'] != "":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
  fastap = relpath("identify/viral/output/combined.final.vOTUs.fa")


### MASTER RULE 

rule done_log:
  name: "viral-host.smk Done. removing tmp files"
  localrule: True
  input:
    relpath("host/output/CHERRY/final_prediction/cherry_prediction.tsv"), 
    relpath("host/output/PhaTYP/final_prediction/phatyp_prediction.tsv"), 
    relpath("host/output/CHERRY/final_prediction/phavip_prediction.tsv"),
    relpath("host/output/merged_host.csv")
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


rule CHERRY:
  name: "viral-host.smk CHERRY host prediction"
  input:
    fna=fastap
  output:
    phavip=relpath("host/output/CHERRY/final_prediction/phavip_prediction.tsv"), 
    cherry=relpath("host/output/CHERRY/final_prediction/cherry_prediction.tsv"), 
    edges=relpath("host/output/CHERRY/final_prediction/cherry_supplementary/cherry_network_edges.tsv"),
    nodes=relpath("host/output/CHERRY/final_prediction/cherry_supplementary/cherry_network_nodes.tsv"),
  params:
    parameters=config['CHERRY-params'],
    dbdir=config['PhaBox2-db'],
    outdir=relpath("host/output/CHERRY"),
    tmpdir=os.path.join(tmpd, "CHERRY")
  conda: "../envs/phabox2.yml"
  log: os.path.join(logdir, "CHERRY.log")
  benchmark: os.path.join(benchmarks, "CHERRY.log")
  threads: 32
  resources: 
    mem_mb=lambda wildcards, attempt: attempt * 100 * 10**3
  shell:
    """
    rm -r {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    phabox2 --task cherry \
        --contigs {input.fna} \
        --threads {threads} \
        --outpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}/

    rm -rf {params.tmpdir}
    """

rule PhaVIP:
  name: "viral-host.smk PhaVIP results from CHERRY"
  localrule: True
  input:
    relpath("host/output/CHERRY/final_prediction/cherry_prediction.tsv")
  output:
    relpath("annotate/viral/output/PhaVIP/final_prediction/phavip_prediction.tsv")
  params:
    indir=relpath("host/output/CHERRY/final_prediction"),
    outdir=relpath("annotate/viral/output/PhaVIP/final_prediction"),
  log: os.path.join(logdir, "PhaVIP_copy.log")
  benchmark: os.path.join(benchmarks, "PhaVIP_copy.log")
  threads: 1
  shell:
    """
    rm -r {params.outdir}
    mkdir -p {params.outdir}/phavip_supplementary

    cp {params.indir}/phavip_prediction.tsv {params.outdir}
    cp {params.indir}/cherry_supplementary/a* {params.outdir}/phavip_supplementary
    cp {params.indir}/cherry_supplementary/gene_annotation.tsv {params.outdir}/phavip_supplementary
    """


rule PhaTYP:
  name: "viral-host.smk PhaTYP lifestyle prediction"
  input:
    fna=fastap, 
    db=os.path.join(config['PhaBox2-db'], "genus2hostlineage.pkl")
  output:
    relpath("host/output/PhaTYP/final_prediction/phatyp_prediction.tsv")
  params:
    parameters=config['PhaTYP-params'],
    dbdir=config['PhaBox2-db'],
    outdir=relpath("host/output/PhaTYP"),
    tmpdir=os.path.join(tmpd, "PhaTYP")
  conda: "../envs/phabox2.yml"
  log: os.path.join(logdir, "PhaTYP.log")
  benchmark: os.path.join(benchmarks, "PhaTYP.log")
  threads: 32
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -r {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    phabox2 --task phatyp \
        --contigs {input.fna} \
        --threads {threads} \
        --outpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}/

    rm -rf {params.tmpdir}
    """

rule merge_results:
  name: "viral-host.smk merge results"
  localrule: True
  input:
    contigs=fastap,
    cherry=relpath("host/output/CHERRY/final_prediction/cherry_prediction.tsv"),
    phatyp=relpath("host/output/PhaTYP/final_prediction/phatyp_prediction.tsv"),
    phavip=relpath("host/output/CHERRY/final_prediction/phavip_prediction.tsv"), 
    edges=relpath("host/output/CHERRY/final_prediction/cherry_supplementary/cherry_network_edges.tsv"),
    nodes=relpath("host/output/CHERRY/final_prediction/cherry_supplementary/cherry_network_nodes.tsv"),
  output:
    merged=relpath("host/output/merged_host.csv"), 
    edges=relpath("host/output/cherry_network_edges.tsv"), 
    nodes=relpath("host/output/cherry_network_nodes.tsv"),
  params: 
    script="workflow/scripts/host_merge.py",
    outdir=relpath("host/output/"), 
    tmpdir=os.path.join(tmpd, "merge")
  conda: "../envs/seqkit-biopython.yml"
  log: os.path.join(logdir, "merge_results.log")
  threads: 1
  shell:
    """
    mkdir -p {params.outdir} {params.tmpdir}

    python {params.script} \
        --cherryout {input.cherry} \
        --phatypout {input.phatyp} \
        --phavipout {input.phavip} \
        --contigs {input.contigs} \
        --output {params.tmpdir}/tmp.csv
    
    mv {params.tmpdir}/tmp.csv {output.merged}
    cp {input.edges} {output.edges}
    cp {input.nodes} {output.nodes}
    """
