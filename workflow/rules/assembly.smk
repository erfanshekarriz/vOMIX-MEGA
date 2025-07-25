assembler = config['assembler']
logdir = relpath(os.path.join("assembly", assembler, "logs"))
tmpd = relpath(os.path.join("assembly", assembler, "tmp"))
benchmarks = relpath(os.path.join("assembly", assembler, "benchmarks"))

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)

email=config["email"]
api_key=config["NCBI-API-key"]
nowstr=config["latest-run"]
outdir=config["outdir"]
datadir=config["datadir"]

samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, api_key, nowstr)

# Check if there are any co-assemblies and abort if using SPAdes and co-assembly
if (len(assemblies.keys()) != len(samples.keys())) and (assembler == "spades"):
  console.print(Panel.fit(f"[bold]Error[/bold]: [dim] SPAdes does not currently support co-assemblies, you may use assembler='megahit' instead. Please read more at https://ablab.github.io/spades/input.html", title="Error", subtitle="SPAdes Co-assembly Support"))
  sys.exit(1)


# MASTER RULE
rule done:
  name: "assembly.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa")),  assembly_id = assemblies.keys()),
    expand(relpath(os.path.join("assembly", assembler, "reports", "{summary_type}.tsv")), summary_type = ["assemblystats", "assembly_size_dist"])
  output:
    os.path.join(logdir, "done.log")
  shell:
    """
    touch {output}
    """

### RULES 

rule megahit:
  name : "assembly.smk MEGAHIT assembly"
  input:
    R1s=lambda wildcards: expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    R2s=lambda wildcards: expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"])
  output:
    fasta=relpath("assembly/megahit/samples/{assembly_id}/output/final.contigs.fa")
  params:
    parameters=config['megahit-params'],
    minlen=config["megahit-minlen"],
    outdir=relpath("assembly/megahit/samples/{assembly_id}/output"),
    tmpdir=os.path.join(tmpd, "megahit")
  log: os.path.join(logdir, "megahit_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "megahit_{assembly_id}.log")
  conda: "../envs/megahit.yml"
  threads: 24
  resources:
    mem_mb = lambda wildcards, attempt, threads, input: max(attempt * input.size_mb * 5, 8000)
  shell:
    """
    rm -rf {params.tmpdir}/{wildcards.assembly_id} {params.outdir}/*
    mkdir -p {params.outdir} {params.tmpdir}


    megahit \
        -1 $(echo "{input.R1s}" | tr ' ' ',') \
        -2 $(echo "{input.R2s}" | tr ' ' ',') \
        --min-contig-len {params.minlen} \
        -o {params.tmpdir}/{wildcards.assembly_id} \
        -t {threads} \
        {params.parameters} &> {log} 

    mv {params.tmpdir}/{wildcards.assembly_id}/final.contigs.fa {output.fasta}
    mv {params.tmpdir}/{wildcards.assembly_id}/* {params.outdir}
    """


rule spades:
  name : "assembly.smk SPAdes (--meta) assembly"
  input:
    R1s=lambda wildcards: expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    R2s=lambda wildcards: expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"])
  output:
    fasta=relpath("assembly/spades/samples/{assembly_id}/output/final.contigs.fa")
  params:
    parameters=config['spades-params'],
    memory=config['spades-memory'],
    outdir=relpath("assembly/spades/samples/{assembly_id}/output"),
    tmpdir=os.path.join(tmpd, "spades/{assembly_id}")
  log: os.path.join(logdir, "spades_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "spades_{assembly_id}.log")
  conda: "../envs/spades.yml"
  threads: 24
  resources:
    mem_mb = config['spades-memory'] * 1024
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}/*
    mkdir -p {params.outdir} {params.tmpdir}


    spades.py \
        -1 $(echo "{input.R1s}" | tr ' ' ',') \
        -2 $(echo "{input.R2s}" | tr ' ' ',') \
        -o {params.tmpdir} \
        -m {params.memory} \
        -t {threads} \
        {params.parameters} &> {log} 

    mv {params.tmpdir}/{wildcards.assembly_id}/scaffolds.fasta {output.fasta}
    mv {params.tmpdir}/{wildcards.assembly_id}/* {params.outdir}
    
    """

rule assembly_stats:
  name: "assembly.smk aggregate assembly statistics"
  localrule: False
  input:
    expand(relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa")), assembly_id = assemblies.keys())
  output:
    stats=relpath(os.path.join("assembly", assembler, "reports", "assemblystats.tsv")),
    sizedist=relpath(os.path.join("assembly", assembler, "reports", "assembly_size_dist.tsv"))
  params:
    script="workflow/scripts/assembly_stats.py",
    outdir=relpath(os.path.join("assembly", assembler, "reports")),
    tmpdir=os.path.join(tmpd, "reports")
  log: os.path.join(logdir, "stats.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt, threads, input: input.size_mb + 2000
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}/* 
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} -i {input} --size-dist-file {params.tmpdir}/tmp1.tsv > {params.tmpdir}/tmp2.tsv 2> {log} 
    
    mv {params.tmpdir}/tmp1.tsv {output.sizedist}
    mv {params.tmpdir}/tmp2.tsv {output.stats}
    """

