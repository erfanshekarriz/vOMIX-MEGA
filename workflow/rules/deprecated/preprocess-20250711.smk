logdir=relpath("preprocess/logs")
tmpd = relpath("preprocess/tmp")

email=config["email"]
api_key=config["NCBI-API-key"]
nowstr=config["latest-run"]
outdir=config["outdir"]
datadir=config["datadir"]

samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, api_key, nowstr)

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

def retrieve_accessions(wildcards):
  try:
    acc=samples[wildcards.sample_id]['accession']
  except KeyError:
    acc=wildcards.sample_id
  return acc

samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, nowstr)

# MASTER RULE

rule done:
    name: "preprocessing.py Done. deleting all tmp files"
    localrule: True
    input:
      expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"), sample_id=samples.keys()),
      expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"), sample_id=samples.keys()),
      expand(os.path.join(datadir, "{sample_id}_{i}.fastq.gz"), sample_id=samples.keys(), i=[1, 2]),
      expand(relpath("preprocess/samples/{sample_id}/output/{sample_id}_R{i}_cut.trim.filt.fastq.gz"), sample_id=samples.keys(), i=[1, 2]),
      relpath("preprocess/reports/preprocess_report.html"), 
      relpath("preprocess/reports/library_size_stats.csv")
    output:
      os.path.join(logdir, "done.log")
    shell:
      """
      touch {output}
      """

# RULES

rule download_fastq:
  name : "preprocessing.py download fastq from SRA"
  output:
    R1=os.path.join(datadir, "{sample_id}_1.fastq.gz"),
    R2=os.path.join(datadir, "{sample_id}_2.fastq.gz")
  params:
    download=config['dwnld-params'],
    pigz=config['pigz-params'],
    logdir=os.path.join(datadir, ".log"), 
    accessions= lambda wildcards: retrieve_accessions(wildcards),
    tmpdir=os.path.join(datadir, ".tmp/{sample_id}")
  log: os.path.join(datadir, ".log/{sample_id}.log")
  conda: "../envs/sratools-pigz.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 8 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.logdir}

    fasterq-dump {params.accessions} \
        {params.download} \
        --split-3 \
        --skip-technical \
        --outdir {params.tmpdir} \
        --temp {params.tmpdir} \
        --threads {threads} &> {log}

    pigz -p {threads} -c {params.pigz} {params.tmpdir}/*_1.fastq > {output.R1} 2> {log}
    pigz -p {threads} -c {params.pigz} {params.tmpdir}/*_2.fastq > {output.R2} 2> {log}

    rm -rf {params.tmpdir}

    """



rule fastp:
  name : "preprocessing.py fastp preprocess"
  input: 
    R1=os.path.join(datadir, '{sample_id}_1.fastq.gz'),
    R2=os.path.join(datadir, '{sample_id}_2.fastq.gz')
  output:
    R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"), 
    html=relpath("preprocess/samples/{sample_id}/report.fastp.html"),
    json=relpath("preprocess/samples/{sample_id}/report.fastp.json")
  params:
    fastp=config['fastp-params'],
    outdir=relpath("preprocess/samples/{sample_id}/output"),
    tmpdir=os.path.join(tmpd, "fastp/{sample_id}")
  log: os.path.join(logdir, "fastp_{sample_id}.log")
  threads: 12
  resources:
    mem_mb = lambda wildcards, input, attempt: attempt * max(5 * input.size_mb, 4000)
  conda: "../envs/fastp.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir} {params.outdir}

    fastp -i {input.R1} -I {input.R2} \
        -o {params.tmpdir}/R1.fastq.gz \
        -O {params.tmpdir}/R2.fastq.gz \
        --thread {threads} \
        --html {params.tmpdir}/tmp.html \
        --json {params.tmpdir}/tmp.json \
        {params.fastp} &> {log}

    mv {params.tmpdir}/R1.fastq.gz {output.R1}
    mv {params.tmpdir}/R2.fastq.gz {output.R2}
    mv {params.tmpdir}/tmp.html {output.html}
    mv {params.tmpdir}/tmp.json {output.json}

    rm -rf {params.tmpdir}
    """


rule aggregate_fastp:
  name: "preprocess.py summarise fastp stats"
  localrule: True
  input:
    jsons=expand(relpath("preprocess/samples/{sample_id}/report.fastp.json"), sample_id = samples.keys())
  output:
    relpath("preprocess/reports/library_size_stats.csv")
  params:
    script="workflow/scripts/fastp_parse.py",
    names=list(samples.keys()),
    outdir=relpath("preprocess/reports"),
    tmpdir=tmpd
  log: os.path.join(logdir, "fastp_summary_stats.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}
    
    echo "{params.names}" > {params.tmpdir}/tmp.names
    echo "{input.jsons}" > {params.tmpdir}/tmp.jsons

    python {params.script} \
        --names {params.tmpdir}/tmp.names \
        --jsons {params.tmpdir}/tmp.jsons > {params.tmpdir}/tmp.csv 2> {log}

    mv {params.tmpdir}/tmp.csv {output}
    rm -r {params.tmpdir}/*
    """ 

rule symlink:
  name: "preprocessing.py creating symbolic links"
  localrule: True
  input:
    R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz")
  output:
    R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz")
  shell:
    """
    ln -s $(pwd)/{input.R1} $(pwd)/{output.R1}
    ln -s $(pwd)/{input.R2} $(pwd)/{output.R2}
    """



rule multiqc:
  name: "preprocessing.py MultiQC preprocess report"
  input:
    R1s=expand(relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"), sample_id = samples.keys()),
    R2s=expand(relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"), sample_id = samples.keys()),
    logs=expand(relpath("preprocess/samples/{sample_id}/report.fastp.json"), sample_id = samples.keys())
  output:
    relpath("preprocess/reports/preprocess_report.html"),
    relpath("preprocess/reports/preprocess_report_data/multiqc.log")
  params:
    searchdir=relpath("preprocess/"),
    outdir=relpath("preprocess/reports"),
    tmpdir=os.path.join(tmpd, "multiqc")
  log: os.path.join(logdir, "multiqc.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, input, attempt: attempt * 8 * 10**3
  conda: "../envs/multiqc.yml"
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    multiqc {params.searchdir} -f -o {params.tmpdir} -n preprocess_report.html 2> {log}
    mv {params.tmpdir}/*.html {params.outdir}
    mv {params.tmpdir}/preprocess* {params.outdir}

    rm -rf {params.tmpdir}/*
    """
