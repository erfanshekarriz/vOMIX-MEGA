import os 

logdir = relpath("binning/viral/logs")
benchmarks = relpath("binning/viral/benchmarks")
tmpd = relpath("binning/viral/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

email=config["email"]
api_key=config["NCBI-API-key"]
nowstr=config["latest_run"]
outdir=config["outdir"]
datadir=config["datadir"]

samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, api_key, nowstr)

rule prodigal_gv:
  name: "viral-binning.py prodigal-gv viral contigs"
  input:
    relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa")
  output:
    fna=relpath("binning/viral/intermediate/prodigal/proteins.vOTUs.fna"),
    faa=relpath("binning/viral/intermediate/prodigal/proteins.vOTUs.faa")
  params:
    script="workflow/software/prodigal-gv/parallel-prodigal-gv.py",
    outdir=relpath("binning/viral/intermediate/prodigal/"),
    tmpdir=os.path.join(tmpd, "prodigal")
  conda: "../envs/prodigal-gv.yml"
  log: os.path.join(logdir, "prodigal-gv.log")
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
        -d {params.tmpdir}/tmp.fna \
        -t {threads} &> {log}

    mv {params.tmpdir}/tmp.faa {output.faa}
    mv {params.tmpdir}/tmp.fna {output.fna}

    rm -rf {params.tmpdir}
    """



rule vRhyme:
  name: "viral-binning.py vRhyme run"
  input:
    contig=relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa"),
    fna=relpath("binning/viral/intermediate/prodigal/proteins.vOTUs.fna"), 
    faa=relpath("binning/viral/intermediate/prodigal/proteins.vOTUs.faa"), 
    fastq=lambda wildcards: expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R{i}.fastq.gz"), sample_id=assemblies[wildcards.assembly_id]["sample_id"], i=[1, 2])
  output:
    summary=relpath("binning/viral/samples/{assembly_id}/vrhyme_bin_summary.tsv"),
    membership=relpath("binning/viral/samples/{assembly_id}/vrhyme_bin_membership.tsv"),
    binned=relpath("binning/viral/samples/{assembly_id}/binned.list.txt")
  params:
    script="workflow/software/vRhyme/vRhyme/vRhyme",
    outdir=relpath("binning/viral/samples/{assembly_id}/vRhyme"), 
    parameters=config["vrhymeparams"], 
    minlen=config["vrhymeminlen"],
    tmpdir=os.path.join(tmpd, "vrhyme/{assembly_id}")
  log: os.path.join(logdir, "vRhyme_{assembly_id}.log")
  conda: "../envs/vrhyme.yml"
  threads: 16
  resources:
    mem_mb=lambda wildcards, attempt, threads: attempt * threads * 1/2 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.outdir}

    python {params.script} \
        -i {input.contig} \
        -r {input.fastq} \
        -g {input.fna} \
        -p {input.faa} \
        -l {params.minlen} \
        -t {threads} \
        -o {params.tmpdir} \
        --verbose \
        {params.parameters} &> {log}

    for f in {params.tmpdir}/vRhyme_best_bins_fasta/*fasta; do \
        newf=$(printf '%s\n' "${{f%.fasta}}_{wildcards.assembly_id}.fasta"); \
        mv $f $newf; done

    cp {params.tmpdir}/vRhyme_best_bins.*.summary.tsv {output.summary}
    cp {params.tmpdir}/vRhyme_best_bins.*.membership.tsv {output.membership}
    cut -f1 {output.membership} | tail -n+2 > {output.binned}
    

    mv {params.tmpdir}/* {params.outdir}
    """



rule merge_binned:
  name: "viral-binning.py add sample prefix"
  input: 
    binned=relpath("binning/viral/samples/{assembly_id}/binned.list.txt")
  output: 
    binned=relpath("binning/viral/samples/{assembly_id}/binned.list.txt")


rule drep_binds:
  name: "viral-binning.py dRep viral bins"



# Right now we are copying the fasta files instead of 
# moving them. Later implement a deletion of the second copy.

rule merge_results:
  name: "viral-binning.py merge binning results"
  localrule: True
  input:
    contig=relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa"),
    prot=relpath("binning/viral/intermediate/prodigal/proteins.vOTUs.faa"),
    binnedlist=expand(relpath("binning/viral/samples/{sample_id}/binned.list.txt"), sample_id=assemblies.keys())
  output:
    contig=relpath("binning/viral/output/combined.unbinned.fa"),
    binned=relpath("binning/viral/output/combined.binned.fa"), 
    binnedlist=relpath("binning/viral/output/binned.list.txt"),
    unbinnedlist=relpath("binning/viral/output/unbinned.list.txt")
  params:
    bindir=expand(relpath("binning/viral/samples/{sample_id}/vRhyme/vRhyme_best_bins_fasta"), sample_id=assemblies.keys()),
    outdir=relpath("binning/viral/output"),
    tmpdir=os.path.join(tmpd, "merge")
  log: os.path.join(logdir, "merge.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}/bins

    for dir in $(echo "{params.bindir}"); do ln -s $(pwd)/${{dir}}/*.fasta $(pwd)/{params.outdir}/bins/; done
    cat {input.binnedlist} | uniq -u > {output.binnedlist}

    seqkit grep --invert-match {input.contig} -f {input.binnedlist} > {params.tmpdir}/tmp.fna
    grep -oE '^>.+' {params.tmpdir}/tmp.fna | sed 's/^>//' > {output.unbinnedlist}
  
    mv {params.tmpdir}/tmp.fna {output.contig}
    """
  



rule checkv_input:
  name: "viral-binning.py prepare CheckV proteins input"


rule derep_bins:
  name: "viral-binning.py dRep viral bins"


