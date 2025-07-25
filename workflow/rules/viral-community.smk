import os

logdir = relpath("community/viral/logs")
tmpd = relpath("community/viral/tmp")
benchmarks=relpath("community/viral/benchmarks")

os.makedirs(benchmarks, exist_ok=True)
os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

methodslist = config["coverm-methods"].split()
methods_c = ",".join(methodslist)

email=config["email"]
api_key=config["NCBI-API-key"]
nowstr=config["latest_run"]
outdir=config["outdir"]
datadir=config["datadir"]

samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, api_key, nowstr)

# MASTER RULE 

rule done:
  name: "viral-community.smk Done. removing tmp files"
  localrule: True
  input:
    pseudo=expand(relpath("community/viral/output/vOTU_table_{methods}.tsv"), methods = methodslist)
  output:
    os.path.join(logdir, "done.log")
  params:
    tmpdir=tmpd,
    interdir=relpath("community/viral/intermediate")
  log: os.path.join(logdir, "done.log")
  shell:
    """
    rm -rf {params.tmpdir}/*
    touch {output}
    """

# RULES

rule coverm_endtoend:
  name: "viral-community.smk CoverM calculate abundance"
  input:
    vOTUs=relpath("identify/viral/output/combined.final.vOTUs.fa"),
    R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz")
  output:
    tsv=relpath("community/viral/samples/{sample_id}/output/vOTU_table.tsv"),
    bam=relpath("community/viral/samples/{sample_id}/output/{sample_id}.bam")
  params:
    parameters=config["coverm-params"], 
    methods=config["coverm-methods"],
    outdir=relpath("community/viral/samples/{sample_id}/output"),
    bamdir=os.path.join(tmpd, "coverm/{sample_id}/bam"),
    tmpdir=os.path.join(tmpd, "coverm/{sample_id}")
  log: os.path.join(logdir, "coverm_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "coverm_{sample_id}.log")
  conda: "../envs/coverm.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 12 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}
  
    coverm contig \
        --coupled {input.R1} {input.R2} \
        --reference {input.vOTUs} \
        --methods {params.methods} \
         --bam-file-cache-directory {params.bamdir} \
        --output-file {params.tmpdir}/tmp.tsv \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/tmp.tsv {output.tsv}
    mv {params.bamdir}/*.bam {output.bam}

    rm -rf {params.tmpdir}
    """



rule coverm_merge:
  name:  "viral-community.smk merge vOTU tables"
  input:
    tables=expand(relpath("community/viral/samples/{sample_id}/output/vOTU_table.tsv"), sample_id = samples.keys())
  output:
    otu=expand(relpath("community/viral/output/vOTU_table_{methods}.tsv"), methods = methodslist)
  params:
    script="workflow/scripts/merged_vOTU_table.py",
    methods=methods_c,
    outdir=relpath("community/viral/output"),
    prefix="vOTU_table",
    tmpdir=os.path.join(tmpd, "merge")
  log: os.path.join(logdir, "coverm_merge.log")
  conda: "../envs/ete3.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: input.size_mb * 5 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}
    
    paste {input.tables} > {params.tmpdir}/tmp.tsv

    python {params.script} \
        --covermout {params.tmpdir}/tmp.tsv \
        --methods {params.methods} \
        --outdir {params.tmpdir} \
        --filename {params.prefix} &> {log}

    rm -f {params.tmpdir}/tmp.tsv
    mv {params.tmpdir}/* {params.outdir}

    rm -rf {params.tmpdir}
    """
