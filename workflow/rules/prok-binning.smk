import os 

logdir=relpath("binning/prok/logs")
benchmarks=relpath("binning/prok/benchmarks")
tmpd=relpath("binning/prok/tmp")

email=config["email"]
api_key=config["NCBI-API-key"]
nowstr=config["latest_run"]
outdir=config["outdir"]
datadir=config["datadir"]

samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, api_key, nowstr)

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

n_cores = config['max-cores'] 
assembler = config['assembler']


### MASTER RULE 
if config["binning-consensus"]:
  rule done_log:
    name: "prok-binning.smk Done. removing tmp files"
    localrule: True
    input:
      expand(relpath("binning/prok/samples/{sample_id}/strobealign/{sample_id}.sorted.{ext}"), sample_id=samples.keys(), ext = ["bam", "bai"]),
      expand(relpath("binning/prok/assemblies/{assembly_id}/{software}/depthfile.txt"), assembly_id=assemblies.keys(), software=["MetaBAT2", "MaxBin2"]),
      expand(relpath("binning/prok/assemblies/{assembly_id}/MaxBin2/bins/maxbin2.summary"), assembly_id=assemblies.keys()), 
      expand(relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/bins/metabat2.unbinned.fa"), assembly_id=assemblies.keys()), 
      expand(relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/concoct_clustering_merged.csv"), assembly_id=assemblies.keys()),
      expand(relpath("binning/prok/assemblies/{assembly_id}/finalbins_DASTool_summary.tsv"), assembly_id=assemblies.keys()),
      expand(relpath("binning/prok/output/no-drep/ids/{assembly_id}_MAGids.tsv"), assembly_id=assemblies.keys()),
      expand(relpath("binning/prok/output/no-drep/unbinned/{assembly_id}_unbinned.fasta"), assembly_id=assemblies.keys()),
      relpath("binning/prok/output/no-drep/checkm2/quality_report.tsv"), 
      relpath("binning/prok/output/clusters.tsv"), 
      relpath("binning/prok/output/taxonomy/gtdbtk/identify/gtdbtk.log")
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
else:
  rule done_log:
    name: "prok-binning.smk Done. removing tmp files"
    localrule: True
    input:
      expand(relpath("binning/prok/samples/{sample_id}/strobealign/{sample_id}.sorted.{ext}"), sample_id=samples.keys(), ext = ["bam", "bai"]),
      expand(relpath("binning/prok/assemblies/{assembly_id}/{software}/depthfile.txt"), assembly_id=assemblies.keys(), software=["MetaBAT2"]),
      expand(relpath("binning/prok/assemblies/{assembly_id}/VAMB/vae_clusters_metadata.tsv"), assembly_id=assemblies.keys())
      #expand(relpath("binning/prok/output/no-drep/ids/{assembly_id}_MAGids.tsv"), assembly_id=assemblies.keys()),
      #expand(relpath("binning/prok/output/no-drep/ids/unbinned/{assembly_id}_unbinned.fasta"), assembly_id=assemblies.keys()),
      #relpath("binning/prok/output/no-drep/checkm2/quality_report.tsv"),
      #relpath("binning/prok/output/clusters.tsv")
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

rule strobealign:
  name: "prok-binning.smk strobealign SR mapping"
  input:
    R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"),
    fasta=lambda wildcards: relpath(os.path.join("assembly", assembler, "samples", samples[wildcards.sample_id]["assembly"], "output/final.contigs.fa")),
  output:
    bam=relpath("binning/prok/samples/{sample_id}/strobealign/{sample_id}.sorted.bam")
  params:
    parameters=config["strobealign-params"],
    tmpdir=os.path.join(tmpd, "strobealign/{sample_id}"),
  log: os.path.join(logdir, "strobealign_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "strobealign_{sample_id}.log")
  conda: "../envs/strobealign.yml"
  threads: 4
  resources:
    mem_mb=lambda wildcards, attempt, input: 8 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}

    strobealign \
        -t {threads} \
        {input.fasta} \
        {input.R1} \
        {input.R2} \
        {params.parameters} 2>{log} | samtools sort - -o {params.tmpdir}/tmp.bam 2> {log}
    
    mv {params.tmpdir}/tmp.bam {output}
    """

rule indexbam:
  name: "prok-binning.smk index sorted BAM"
  localrule: True
  input:
    relpath("binning/prok/samples/{sample_id}/strobealign/{sample_id}.sorted.bam")
  output:
    relpath("binning/prok/samples/{sample_id}/strobealign/{sample_id}.sorted.bai")
  log: os.path.join(logdir, "indexbam_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "indexbam_{sample_id}.log")
  conda: "../envs/strobealign.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    samtools index {input} -o {output} 2> {log}
    """

rule binprep:
  name: "prok-binning.smk jgi summarize bam contig depths"
  input:
    bams=lambda wildcards: expand(relpath("binning/prok/samples/{sample_id}/strobealign/{sample_id}.sorted.bam"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    bais=lambda wildcards: expand(relpath("binning/prok/samples/{sample_id}/strobealign/{sample_id}.sorted.bai"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
  output:
    relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/depthfile.txt"), 
  params:
    parameters=config["jgi-summarize-params"],
    outdir=relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2"),
    tmpdir=os.path.join(tmpd, "MetaBAT2/{assembly_id}")
  log: os.path.join(logdir, "MetaBAT2_prep_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "MetaBAT2_prep_{assembly_id}.log")
  conda: "../envs/metabat2.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input: 20 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    jgi_summarize_bam_contig_depths \
        --outputDepth {params.tmpdir}/tmp.txt \
        {input.bams} \
        {params.parameters} 2>{log}
    
    mv {params.tmpdir}/tmp.txt {output}
    """


# VAMB BINNING 

rule vamb:
  name: "prok-binning.smk VAMB binning"
  input: 
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa")),
    txt=relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/depthfile.txt")
  output:
   relpath("binning/prok/assemblies/{assembly_id}/VAMB/vae_clusters_metadata.tsv")
  params:
    parameters=config["VAMB-params"],
    outdir=relpath("binning/prok/assemblies/{assembly_id}/VAMB"), 
    tmpdir=os.path.join(tmpd, "VAMB/{assembly_id}")
  log: os.path.join(logdir, "VAMB_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "VAMB_{assembly_id}.log")
  conda: "../envs/vamb.yml"
  threads: 16
  resources:
    mem_mb=lambda wildcards, attempt, input: 6 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.outdir}

    vamb \
        --outdir {params.tmpdir} \
        --fasta {input.fasta} \
        --jgi {input.txt} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}
    """


# CONSENSUS BINNING

rule metabat2maxbin:
  name: "prok-binning.smk MetaBAT2 2 MaxBin2 coverage file"
  input:
    relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/depthfile.txt")
  output:
    relpath("binning/prok/assemblies/{assembly_id}/MaxBin2/depthfile.txt")
  params:
    script="workflow/scripts/metabat2maxbin.py",
    outdir=relpath("binning/prok/assemblies/{assembly_id}/MaxBin2"),
  log: os.path.join(logdir, "MaxBin2_prep_{assembly_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: (input.size_mb + 2 * 10**3) * attempt
  shell:
    """
    mkdir -p {params.outdir}

    python {params.script} \
        -inputtxt {input} \
        -outputtxt {output} \
        -ending sorted 2> {log}
    """

rule metabat2:
  name: "prok-binning.smk MetaBAT2 binning"
  input:
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa")),
    txt=relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/depthfile.txt")
  output:
    relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/bins/metabat2.unbinned.fa")
  params:
    parameters=config["MetaBAT2-params"],
    outdir=relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2"), 
    tmpdir=os.path.join(tmpd, "MetaBAT2/{assembly_id}/bins")
  log: os.path.join(logdir, "MetaBAT2_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "MetaBAT2_{assembly_id}.log")
  conda: "../envs/metabat2.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input: 6 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}/bins

    metabat2 \
        -t {threads} \
        -i {input.fasta} \
        -a {input.txt} \
        -o {params.tmpdir}/metabat2 \
        --unbinned \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}/bins
    """

rule maxbin2:
  name: "prok-binning.smk MaxBin2 binning"
  input:
    jgi=relpath("binning/prok/assemblies/{assembly_id}/MaxBin2/depthfile.txt"),
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa"))
  output:
    relpath("binning/prok/assemblies/{assembly_id}/MaxBin2/bins/maxbin2.summary")
  params:
    parameters=config["MaxBin2-params"],
    outdir=relpath("binning/prok/assemblies/{assembly_id}/MaxBin2"),
    tmpdir=os.path.join(tmpd, "MaxBin2/{assembly_id}/bins")
  log: os.path.join(logdir, "MaxBin2_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "MaxBin2_{assembly_id}.log")
  conda: "../envs/maxbin2.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}/bins

    run_MaxBin.pl \
        -thread {threads} \
        -contig {input.fasta} \
        -abund {input.jgi} \
        -out {params.tmpdir}/maxbin2 \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}/bins
    """

rule concoctprep:
  name: "prok-binning.smk CONCOCT prepare"
  input:
    index=lambda wildcards: expand(relpath("binning/prok/samples/{sample_id}/strobealign/{sample_id}.sorted.bai"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    bams=lambda wildcards: expand(relpath("binning/prok/samples/{sample_id}/strobealign/{sample_id}.sorted.bam"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa"))
  output:
    bed=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/final.contigs.10k.bed"),
    fa=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/final.contigs.10k.fa"), 
    tsv=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/depth.tsv"),
  params:
    outdir=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT"), 
    tmpdir=os.path.join(tmpd, "CONCOCT/{assembly_id}")
  log: os.path.join(logdir, "CONCOCT_prep_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "CONCOCT_prep_{assembly_id}.log")
  conda: "../envs/concoct.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input: 2 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    cut_up_fasta.py {input.fasta} \
        -c 10000 \
        -o 0 \
        --merge_last \
        -b {params.tmpdir}/tmp.bed > {params.tmpdir}/tmp.fa 2> {log}

    concoct_coverage_table.py \
        {params.tmpdir}/tmp.bed \
        {input.bams} > {params.tmpdir}/tmp.tsv

    mv {params.tmpdir}/tmp.bed {output.bed}
    mv {params.tmpdir}/tmp.fa {output.fa}
    mv {params.tmpdir}/tmp.tsv {output.tsv}
    """


rule concoct:
  name: "prok-binning.smk CONCOCT binning"
  input:
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa")),
    bed=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/final.contigs.10k.bed"),
    fa=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/final.contigs.10k.fa"),       
    tsv=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/depth.tsv"),
  output:
    csv=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/concoct_clustering_merged.csv")
  params:
    parameters=config["CONCOCT-params"],
    outdir=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT"),
    tmpdir=os.path.join(tmpd, "CONCOCT/{assembly_id}")
  log: os.path.join(logdir, "CONCOCT_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "CONCOCT_{assembly_id}.log")
  conda: "../envs/concoct.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}/bins
    mkdir -p {params.tmpdir}/bins {params.outdir}

    concoct \
        -t {threads} \
        --composition_file {input.fa} \
        --coverage_file {input.tsv} \
        -b {params.tmpdir}/concoct \
        {params.parameters} 2> {log}

    merge_cutup_clustering.py \
        {params.tmpdir}/concoct_clustering_gt1000.csv > \
        {params.tmpdir}/tmp.csv 2> {log}
    
    extract_fasta_bins.py {input.fasta} \
        {params.tmpdir}/tmp.csv \
        --output_path {params.tmpdir}/bins
    
    mv {params.tmpdir}/bins {params.outdir}/
    mv {params.tmpdir}/tmp.csv {output.csv}
    mv {params.tmpdir}/* {params.outdir}
    """


rule contigs2bin:
  name: "prok-binning.smk Generate contigs2bin tables"
  input:
    MetaBAT2=relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/bins/metabat2.unbinned.fa"),
    MaxBin2=relpath("binning/prok/assemblies/{assembly_id}/MaxBin2/bins/maxbin2.summary"),
    CONCOCT=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/concoct_clustering_merged.csv")
  output:
    MetaBAT2=relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/contigs2bin.tsv"),
    MaxBin2=relpath("binning/prok/assemblies/{assembly_id}/MaxBin2/contigs2bin.tsv"),
    CONCOCT=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/contigs2bin.tsv")
  params:
    MetaBAT2=relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/bins"),
    MaxBin2=relpath("binning/prok/assemblies/{assembly_id}/MaxBin2/bins"),
  log: os.path.join(logdir, "contigs2bin_{assembly_id}.log")
  conda: "../envs/dastool.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    Fasta_to_Contig2Bin.sh -i {params.MetaBAT2}/ -e fa > {output.MetaBAT2}
    Fasta_to_Contig2Bin.sh -i {params.MaxBin2}/ -e fasta > {output.MaxBin2}
    tail -n+2 {input.CONCOCT} | tr ',' '\t' > {output.CONCOCT}
    """

rule dastool:
  name: "prok-binning.smk DAS Tool consensus binning"
  input:
    fasta=relpath(os.path.join("assembly", assembler, "samples/{assembly_id}/output/final.contigs.fa")),
    MetaBAT2=relpath("binning/prok/assemblies/{assembly_id}/MetaBAT2/contigs2bin.tsv"),
    MaxBin2=relpath("binning/prok/assemblies/{assembly_id}/MaxBin2/contigs2bin.tsv"),
    CONCOCT=relpath("binning/prok/assemblies/{assembly_id}/CONCOCT/contigs2bin.tsv")
  output:
    relpath("binning/prok/assemblies/{assembly_id}/finalbins_DASTool_summary.tsv")
  params:
    parameters=config["DASTool-params"], 
    outdir=relpath("binning/prok/assemblies/{assembly_id}"),
    tmpdir=os.path.join(tmpd, "DASTool/{assembly_id}")
  log: os.path.join(logdir, "DASTool_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "DASTool_{assembly_id}.log")
  conda: "../envs/dastool.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input: 8 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}/finalbins*
    mkdir -p {params.tmpdir}

    DAS_Tool \
        -i {input.MetaBAT2},{input.MaxBin2},{input.CONCOCT} \
        -l MetaBAT2,MaxBin2,CONCOCT \
        -c {input.fasta} \
        -o {params.tmpdir}/finalbins \
        --write_bin_evals \
        -t {threads} \
        --write_bins \
        --write_unbinned \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    """


rule movemags:
  name: "prok-binning.smk Combine MAGs from assemblies"
  localrule: True
  input:
    relpath("binning/prok/assemblies/{assembly_id}/finalbins_DASTool_summary.tsv")
  output:
    tsv=relpath("binning/prok/output/no-drep/ids/{assembly_id}_MAGids.tsv"), 
    unbinned=relpath("binning/prok/output/no-drep/unbinned/{assembly_id}_unbinned.fasta")
  params:
    indir=relpath("binning/prok/assemblies/{assembly_id}/finalbins_DASTool_bins"),
    outdir=relpath("binning/prok/output/no-drep")
  log: os.path.join(logdir, "move_mags_{assembly_id}.log")
  benchmark: os.path.join(logdir, "move_mags_{assembly_id}.log")
  threads: 1
  shell:
    """
    mkdir -p {params.outdir}/ids {params.outdir}/unbinned
    echo -e "ids\tDASTool_id" > {output.tsv}

    [ -f {params.indir}/unbinned* ] && mv {params.indir}/unbinned* {output.unbinned} || {{ touch {output.unbinned}; }}

    [ -f {params.indir}/unbinned*] && mv -f {params.indir}/unbinned* {output.unbinned}

    ls -v {params.indir} | grep -v "unbinned" | cat -n | \
        while read n f; do \
        echo -e "{wildcards.assembly_id}_${{n}}.fasta\t${{f}}" >> {output.tsv}; \
        cp "{params.indir}/$f" "{params.outdir}/{wildcards.assembly_id}_$n.fasta"; \
        done 2> {log}
    """


# REST OF ANALYSIS

rule checkm2:
  name: "prok-binning.smk Checkm2 predict"
  input:
    os.path.join(config["checkm2-db"], "uniref100.KO.1.dmnd"), 
    expand(relpath("binning/prok/output/no-drep/ids/{assembly_id}_MAGids.tsv"), assembly_id=assemblies.keys()),
    expand(relpath("binning/prok/output/no-drep/unbinned/{assembly_id}_unbinned.fasta"), assembly_id=assemblies.keys())
  output:
    relpath("binning/prok/output/no-drep/checkm2/quality_report.tsv"),
  params:
    parameters=config["checkm2-params"],
    indir=relpath("binning/prok/output/no-drep"),
    outdir=relpath("binning/prok/output/no-drep/checkm2"), 
    dbdir=config["checkm2-db"],
    tmpdir=os.path.join(tmpd, "CheckM2")
  log: os.path.join(logdir, "checkm2_predict.log")
  benchmark: os.path.join(logdir, "checkm2_predict.log")
  conda: "../envs/checkm2.yml"
  threads: 64
  resources: 
    mem_mb=lambda wildcards, attempt, input: 16 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.outdir} {params.tmpdir}/bins

    ln -s $(pwd)/{params.indir}/*.fasta $(pwd)/{params.tmpdir}/bins/

    checkm2 predict \
        --input {params.tmpdir}/bins \
        --output-directory {params.tmpdir}/results \
        --database_path {params.dbdir}/*.dmnd \
        --threads {threads} \
        --extension .fasta \
        {params.parameters} &> {log}

    mv {params.tmpdir}/results/* {params.outdir}
    """

rule galah:
  name: "prok-binning.smk Galah dereplicate bins"
  input:
    tsv=relpath("binning/prok/output/no-drep/checkm2/quality_report.tsv"),
  output:
    tsv=relpath("binning/prok/output/clusters.tsv")
  params:
    parameters=config["galah-params"], 
    outdir=relpath("binning/prok/output"),
    indir=relpath("binning/prok/output/no-drep"),
    tmpdir=os.path.join(tmpd, "galah")
  log: os.path.join(logdir, "galah_dereplicate.log")
  benchmark: os.path.join(benchmarks, "galah_dereplicate.log")
  conda: "../envs/galah.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt, input: 8 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}/drep
    mkdir -p {params.tmpdir} {params.outdir}/drep

    galah cluster \
        --genome-fasta-files {params.indir}/*fasta \
        --checkm2-quality-report {input.tsv} \
        --output-cluster-definition {params.tmpdir}/tmp.tsv \
        --output-representative-fasta-directory {params.outdir}/drep \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/tmp.tsv {output.tsv}
    """

rule GTDBTk_identify:
  name: "prok-binning.smk GTDB-Tk identify"
  input:
    tsv=relpath("binning/prok/output/clusters.tsv"), 
    db=os.path.join(config["GTDBTk-db"], "done.log")
  output:
    relpath("binning/prok/output/taxonomy/gtdbtk/identify/gtdbtk.log")
  params:
    parameters=config["GTDBTk-identify-params"], 
    outdir=relpath("binning/prok/output/taxonomy/gtdbtk/identify"), 
    indir=relpath("binning/prok/output/drep"), 
    tmpdir=os.path.join(tmpd, "gtdbtk/identify")
  log: os.path.join(logdir, "gtdbtk_identify.log")
  benchmark: os.path.join(benchmarks, "gtdbtk_identify.log")
  conda: "../envs/gtdbtk.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt, input: 8 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    gtdbtk identify \
        --genome_dir {params.indir} \
        --out_dir {params.tmpdir} \
        --extension fasta \
        --cpus {threads} \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    """

rule GTDBTk_align:
  name: "prok-binning.smk GTDB-Tk align"
  input:
    relpath("binning/prok/output/taxonomy/gtdbtk/identify/gtdbtk.log")
  output:
    relpath("binning/prok/output/taxonomy/gtdbtk/align/gtdbtk.log")
  params:
    parameters=config["GTDBTk-align-params"],
    outdir=relpath("binning/prok/output/taxonomy/gtdbtk/align"),
    indir=relpath("binning/prok/output/taxonomy/gtdbtk/identify"),
    tmpdir=os.path.join(tmpd, "gtdbtk/align")
  log: os.path.join(logdir, "gtdbtk_align.log")
  benchmark: os.path.join(benchmarks, "gtdbtk_align.log")
  conda: "../envs/gtdbtk.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt, input: 100 * 10**3 * attempt
  shell:
    """ 
    """

rule GTDBTk_classify:
  name: "prok-binning.smk GTDB-Tk classify"
  input:
    relpath("binning/prok/output/taxonomy/gtdbtk/align/gtdbtk.log")
  output:
    relpath("binning/prok/output/taxonomy/gtdbtk/classify/gtdbtk.log"),
  params:
    parameters=config["GTDBTk-classify-params"],
    outdir=relpath("binning/prok/output/taxonomy/gtdbtk/classify"),
    genomedir=relpath("binning/prok/output/drep"),
    aligndir=relpath("binning/prok/output/taxonomy/gtdbtk/align"),
    tmpdir=os.path.join(tmpd, "gtdbtk/classify")
  log: os.path.join(logdir, "gtdbtk_classify.log")
  benchmark: os.path.join(benchmarks, "gtdbtk_classify.log")
  conda: "../envs/gtdbtk.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt, input: 8 * 10**3 * attempt
  shell:
    """ 
    """
