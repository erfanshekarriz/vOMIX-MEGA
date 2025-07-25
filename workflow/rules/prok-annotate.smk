logdir=relpath("annotate/prok/logs")
benchmarks=relpath("annotate/prok/benchmarks")
tmpd = relpath("annotate/prok/tmp")

email=config["email"]
api_key=config["NCBI-API-key"]
nowstr=config["latest-run"]
outdir=config["outdir"]
datadir=config["datadir"]

samples, assemblies = parse_sample_list(config["workdir"] + config["samplelist"], datadir, outdir, email, api_key, nowstr)

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)


### MASTER RULE 
rule done_log:
  name: "prok-annotate.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("annotate/prok/samples/{sample_id}/input/{sample_id}_merged.fastq.gz"), sample_id = samples.keys()), 
    expand(relpath("annotate/prok/samples/{sample_id}/{sample_id}_{output_type}.tsv"), sample_id = samples.keys(), output_type = ["pathabundance", "pathcoverage", "genefamilies"]), 
    expand(relpath("annotate/prok/output/primary/{output_type}_merged.tsv"), output_type = ["pathabundance", "pathcoverage", "genefamilies"]),
    expand(relpath("annotate/prok/output/unnamed/genefamilies_merged-cpm-{database}.tsv"), database = ["ec", "eggnog", "go", "ko", "pfam"]),
    expand(relpath("annotate/prok/output/genefamilies_merged-cpm-{database}-named.tsv"), database = ["ec", "eggnog", "go", "ko", "pfam"]),
    expand(relpath("annotate/prok/output/{folder}/genefamilies_merged-cpm-{database}-named_{folder}.tsv"), folder = ["stratified", "unstratified"], database = ["ec", "eggnog", "go", "ko", "pfam"]),
    #os.path.join(benchmarks, "summary.tsv"),
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
rule pair_fastq:
  name: "prok-annotate.smk merge paired files"
  input: 
    R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"),
  output:
    relpath("annotate/prok/samples/{sample_id}/input/{sample_id}_merged.fastq.gz")
  params:
    outdir=relpath("annotate/prok/samples/{sample_id}/input"),
    tmpdir=os.path.join(tmpd, "cat/{sample_id}")
  log: os.path.join(logdir, "paired_concat_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "paired_concat_{sample_id}.log")
  threads: 1
  resources:
    mem_mb = lambda wildcards, input: max(4*input.size_mb, 2000)
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}

    cat {input.R1} {input.R2} > {params.tmpdir}/tmp.fastq.gz 2> {log}
    mv {params.tmpdir}/tmp.fastq.gz {output}
    """


rule humann3:
  name: "prok-annotate.smk HUMAnN3 run."
  input:
    fastq=relpath("annotate/prok/samples/{sample_id}/input/{sample_id}_merged.fastq.gz"),
    chocophlandb=os.path.join(config['humann-db'], "chocophlan/alaS.centroids.v201901_v31.ffn.gz"), 
    unirefdb=os.path.join(config['humann-db'], "uniref/uniref90_201901b_full.dmnd"), 
    utilitydb=os.path.join(config['humann-db'], "utility_mapping/map_ec_name.txt.gz")
  output:
    gene=relpath("annotate/prok/samples/{sample_id}/{sample_id}_genefamilies.tsv"),
    abund=relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathabundance.tsv"),
    cov=relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathcoverage.tsv")
  params:
    parameters=config["humann-params"],
    outdir=relpath("annotate/prok/samples/{sample_id}"), 
    dbdir=config['humann-db'], 
    name="{sample_id}",
    tmpdir=os.path.join(tmpd, "humann/{sample_id}"), 
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "HUMAnN3_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "HUMAnN3_{sample_id}.log")
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: max(42 * 10**3, 4000)
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    humann \
        --input {input.fastq} \
        --output {params.tmpdir} \
        --threads {threads} \
        --output-format tsv \
        --output-basename {params.name} \
        --protein-database {params.dbdir}/uniref \
        --nucleotide-database {params.dbdir}/chocophlan \
        --log-level INFO \
        --o-log {log} \
        {params.parameters} &>> {log}

    mv {params.tmpdir}/* {params.outdir}/
    rm -rf {params.tmpdir}
    """


rule humann_merge:
  name: "prok-annotate.smk HUMAnN3 merge results"
  input:
    gene=expand(relpath("annotate/prok/samples/{sample_id}/{sample_id}_genefamilies.tsv"), sample_id = assemblies.keys()),
    abund=expand(relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathabundance.tsv"), sample_id = assemblies.keys()),
    cov=expand(relpath("annotate/prok/samples/{sample_id}/{sample_id}_pathcoverage.tsv"), sample_id = assemblies.keys())
  output:
    gene=relpath("annotate/prok/output/primary/genefamilies_merged.tsv"),
    abund=relpath("annotate/prok/output/primary/pathabundance_merged.tsv"),
    cov=relpath("annotate/prok/output/primary/pathcoverage_merged.tsv")
  params:
    outdir=relpath("annotate/prok/output/primary"),
    tmpdir=os.path.join(tmpd, "humann_merge")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_merge.log")
  benchmark: os.path.join(benchmarks, "humann_merge.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: input.size_mb + 4000
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}
  
    for file in {input.gene}; do ln -s $(pwd)/$file $(pwd)/{params.tmpdir} &>> {log}; done
    humann_join_tables --input {params.tmpdir} --output {output.gene} &>> {log}
    rm {params.tmpdir}/* 

    for file in {input.abund}; do ln -s $(pwd)/$file $(pwd)/{params.tmpdir} &>> {log}; done
    humann_join_tables --input {params.tmpdir} --output {output.abund} &>> {log}
    rm {params.tmpdir}/*
    
    for file in {input.cov}; do ln -s $(pwd)/$file $(pwd)/{params.tmpdir} &>> {log}; done
    humann_join_tables --input {params.tmpdir} --output {output.cov} &>> {log}
    rm {params.tmpdir}/*
    """

rule humann_normalize:
  name: "prok-annotate.smk HUMAnN3 normalize results"
  input:
    relpath("annotate/prok/output/primary/genefamilies_merged.tsv")
  output:
    relpath("annotate/prok/output/primary/genefamilies_merged-cpm.tsv")
  params:
    outdir=relpath("annotate/prok/output/primary"),
    tmpdir=os.path.join(tmpd, "humann_normalize")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_normalize.log")
  benchmark: os.path.join(benchmarks, "humann_normalize.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: 400 * 10**3
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}
    
    humann_renorm_table --input {input} --output {params.tmpdir}/tmp.tsv --units cpm --mode community --special y --update-snames &> {log}

    mv {params.tmpdir}/tmp.tsv {output}
    """


rule humann_map_ec:
  name: "prok-annotate.smk HUMAnN3 map Level-4 Enzyme Comission"
  input:
    tsv=relpath("annotate/prok/output/primary/genefamilies_merged-cpm.tsv"),
    db=os.path.join(config['humann-db'], "utility_mapping/map_level4ec_uniref90.txt.gz"),
    dbname=os.path.join(config['humann-db'], "utility_mapping/map_ec_name.txt.gz")
  output:
    tsv=relpath("annotate/prok/output/unnamed/genefamilies_merged-cpm-ec.tsv"),
    tsvname=relpath("annotate/prok/output/genefamilies_merged-cpm-ec-named.tsv"),
  params:
    outdir=relpath("annotate/prok/output"),
    tmpdir=os.path.join(tmpd, "humann_map/ec")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_map_ec.log")
  benchmark: os.path.join(benchmarks, "humann_map_ec.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: input.size_mb + (200 * 10**3)
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.outdir}/unnamed {params.tmpdir}

    humann_regroup_table --input {input.tsv} --output {params.tmpdir}/tmp.tsv --custom {input.db} --function sum --ungrouped Y --protected Y &>> {log}
    humann_rename_table --input {params.tmpdir}/tmp.tsv --output {params.tmpdir}/tmpname.tsv --custom {input.dbname} &>> {log}

    mv {params.tmpdir}/tmp.tsv {output.tsv}
    mv {params.tmpdir}/tmpname.tsv {output.tsvname}
    """

rule humann_map_eggnog:
  name: "prok-annotate.smk HUMAnN3 map EggNOG including COGs"
  input:
    tsv=relpath("annotate/prok/output/primary/genefamilies_merged-cpm.tsv"),
    db=os.path.join(config['humann-db'], "utility_mapping/map_eggnog_uniref90.txt.gz"),
    dbname=os.path.join(config['humann-db'], "utility_mapping/map_eggnog_name.txt.gz"),
  output:
    tsv=relpath("annotate/prok/output/unnamed/genefamilies_merged-cpm-eggnog.tsv"),
    tsvname=relpath("annotate/prok/output/genefamilies_merged-cpm-eggnog-named.tsv"),
  params:
    outdir=relpath("annotate/prok/output"),
    tmpdir=os.path.join(tmpd, "humann_map/eggnog")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_map_eggnog.log")
  benchmark: os.path.join(benchmarks, "humann_map_eggnog.log")
  threads: 1 
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: input.size_mb + (200 * 10**3)
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.outdir}/unnamed {params.tmpdir}

    humann_regroup_table --input {input.tsv} --output {params.tmpdir}/tmp.tsv --custom {input.db} --function sum --ungrouped Y --protected Y &>> {log}
    humann_rename_table --input {params.tmpdir}/tmp.tsv --output {params.tmpdir}/tmpname.tsv --custom {input.dbname} &>> {log}

    mv {params.tmpdir}/tmp.tsv {output.tsv}
    mv {params.tmpdir}/tmpname.tsv {output.tsvname}
    """


rule humann_map_go:
  name: "prok-annotate.smk HUMAnN3 map Gene Ontology"
  input:
    tsv=relpath("annotate/prok/output/primary/genefamilies_merged-cpm.tsv"),
    db=os.path.join(config['humann-db'], "utility_mapping/map_go_uniref90.txt.gz"),
    dbname=os.path.join(config['humann-db'], "utility_mapping/map_go_name.txt.gz"),
  output:
    tsv=relpath("annotate/prok/output/unnamed/genefamilies_merged-cpm-go.tsv"),
    tsvname=relpath("annotate/prok/output/genefamilies_merged-cpm-go-named.tsv"),
  params:
    outdir=relpath("annotate/prok/output"),
    tmpdir=os.path.join(tmpd, "humann_map/go")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_map_go.log")
  benchmark: os.path.join(benchmarks, "humann_map_go.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: input.size_mb + (200 * 10**3)
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.outdir}/unnamed {params.tmpdir}

    humann_regroup_table --input {input.tsv} --output {params.tmpdir}/tmp.tsv --custom {input.db} --function sum --ungrouped Y --protected Y &>> {log}
    humann_rename_table --input {params.tmpdir}/tmp.tsv --output {params.tmpdir}/tmpname.tsv --custom {input.dbname} &>> {log}

    mv {params.tmpdir}/tmp.tsv {output.tsv}
    mv {params.tmpdir}/tmpname.tsv {output.tsvname}
    """

rule humann_map_ko:
  name: "prok-annotate.smk HUMAnN3 map KEGG Orthogroups"
  input:
    tsv=relpath("annotate/prok/output/primary/genefamilies_merged-cpm.tsv"),
    db=os.path.join(config['humann-db'], "utility_mapping/map_ko_uniref90.txt.gz"),
    dbname=os.path.join(config['humann-db'], "utility_mapping/map_ko_name.txt.gz"),
  output:
    tsv=relpath("annotate/prok/output/unnamed/genefamilies_merged-cpm-ko.tsv"),
    tsvname=relpath("annotate/prok/output/genefamilies_merged-cpm-ko-named.tsv"),
  params:
    outdir=relpath("annotate/prok/output"),
    tmpdir=os.path.join(tmpd, "humann_map/ko")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_map_ko.log")
  benchmark: os.path.join(benchmarks, "humann_map_ko.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: input.size_mb + (200 * 10**3)
  shell:
    """
    rm -rf {params.tmpdir} 
    mkdir -p {params.outdir}/unnamed {params.tmpdir}

    humann_regroup_table --input {input.tsv} --output {params.tmpdir}/tmp.tsv --custom {input.db} --function sum --ungrouped Y --protected Y &>> {log}
    humann_rename_table --input {params.tmpdir}/tmp.tsv --output {params.tmpdir}/tmpname.tsv --custom {input.dbname} &>> {log}

    mv {params.tmpdir}/tmp.tsv {output.tsv}
    mv {params.tmpdir}/tmpname.tsv {output.tsvname}
    """


rule humann_map_pfam:
  name: "prok-annotate.smk HUMAnN3 map Pfam domains"
  input:
    tsv=relpath("annotate/prok/output/primary/genefamilies_merged-cpm.tsv"),
    db=os.path.join(config['humann-db'], "utility_mapping/map_pfam_uniref90.txt.gz"),
    dbname=os.path.join(config['humann-db'], "utility_mapping/map_pfam_name.txt.gz"),
  output:
    tsv=relpath("annotate/prok/output/unnamed/genefamilies_merged-cpm-pfam.tsv"),
    tsvname=relpath("annotate/prok/output/genefamilies_merged-cpm-pfam-named.tsv"),
  params:
    outdir=relpath("annotate/prok/output"),
    tmpdir=os.path.join(tmpd, "humann_map/pfam")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_map_pfam.log")
  benchmark: os.path.join(benchmarks, "humann_map_pfam.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: input.size_mb + (200 * 10**3)
  shell:
    """
    rm -rf {params.tmpdir} 
    mkdir -p {params.outdir}/unnamed {params.tmpdir}

    humann_regroup_table --input {input.tsv} --output {params.tmpdir}/tmp.tsv --custom {input.db} --function sum --ungrouped Y --protected Y &>> {log}
    humann_rename_table --input {params.tmpdir}/tmp.tsv --output {params.tmpdir}/tmpname.tsv --custom {input.dbname} &>> {log}

    mv {params.tmpdir}/tmp.tsv {output.tsv}
    mv {params.tmpdir}/tmpname.tsv {output.tsvname}
    """

rule humann_map_metacyc:
  name: "prok-annotate.smk HUMAnN3 map MetaCyc reactions"
  input:
    tsv=relpath("annotate/prok/output/primary/genefamilies_merged-cpm.tsv"),
    db=os.path.join(config['humann-db'], "utility_mapping/map_pfam_uniref90.txt.gz"),
    dbname=os.path.join(config['humann-db'], "utility_mapping/map_pfam_name.txt.gz"),
  output:
    tsv=relpath("annotate/prok/output/unnamed/genefamilies_merged-cpm-metacycrxn.tsv"),
    tsvname=relpath("annotate/prok/output/genefamilies_merged-cpm-metacycrxn-named.tsv"),
  params:
    outdir=relpath("annotate/prok/output"),
    tmpdir=os.path.join(tmpd, "humann_map/metacycrxn")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_map_metacycrxn.log")
  benchmark: os.path.join(benchmarks, "humann_map_metacycrxn.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, threads, input: input.size_mb + (200 * 10**3)
  shell:
    """
    rm -rf {params.tmpdir} 
    mkdir -p {params.outdir}/unnamed {params.tmpdir}

    humann_regroup_table --input {input.tsv} --output {params.tmpdir}/tmp.tsv --groups uniref90_rxn --function sum --ungrouped Y --protected Y &>> {log}
    humann_rename_table --input {params.tmpdir}/tmp.tsv --output {params.tmpdir}/tmpname.tsv --names metacyc-rxn &>> {log}

    mv {params.tmpdir}/tmp.tsv {output.tsv}
    mv {params.tmpdir}/tmpname.tsv {output.tsvname}
    """


rule humann_split:
  name: "prok-annotate.smk HUMAnN3 split stratified table"
  localrule: True
  input:
    cpm=relpath("annotate/prok/output/primary/genefamilies_merged-cpm.tsv"), 
    ec=relpath("annotate/prok/output/genefamilies_merged-cpm-ec-named.tsv"), 
    eggnog=relpath("annotate/prok/output/genefamilies_merged-cpm-eggnog-named.tsv"), 
    ko=relpath("annotate/prok/output/genefamilies_merged-cpm-ko-named.tsv"), 
    go=relpath("annotate/prok/output/genefamilies_merged-cpm-go-named.tsv"), 
    pfam=relpath("annotate/prok/output/genefamilies_merged-cpm-pfam-named.tsv")
  output:
    expand(relpath("annotate/prok/output/{folder}/genefamilies_merged-cpm-{database}-named_{folder}.tsv"),
        folder = ["stratified", "unstratified"], database = ["ec", "eggnog", "go", "ko", "pfam"]),
  params:
    outdir=relpath("annotate/prok/output"),
    tmpdir=os.path.join(tmpd, "humann_split")
  conda: "../envs/biobakery3.yml"
  log: os.path.join(logdir, "humann_split.log")
  benchmark: os.path.join(benchmarks, "humann_split.log")
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.outdir}/stratified {params.outdir}/unstratified {params.tmpdir}
    
    humann_split_stratified_table --input {input.cpm} --output {params.tmpdir} &>> {log}
    humann_split_stratified_table --input {input.ec} --output {params.tmpdir} &>> {log}
    humann_split_stratified_table --input {input.eggnog} --output {params.tmpdir} &>> {log}
    humann_split_stratified_table --input {input.ko} --output {params.tmpdir} &>> {log}
    humann_split_stratified_table --input {input.go} --output {params.tmpdir} &>> {log}
    humann_split_stratified_table --input {input.pfam} --output {params.tmpdir} &>> {log}
    
    mv {params.tmpdir}/*_stratified.tsv {params.outdir}/stratified/
    mv {params.tmpdir}/*_unstratified.tsv {params.outdir}/unstratified/
    """

rule benchmark_summary:
  name: "prok-annotate.smk summaries performance benchmarks"
  localrule: True
  input:
    os.path.join(benchmarks, "humann_split.log")
  output:
    os.path.join(benchmarks, "summary.tsv")
  params:
    indir=benchmarks
  threads: 1
  shell:
    """
    cat {params.indir}/* | head -n 1 > {output}
    for file in $(find {params.indir}/*.log -type f); do echo -e "$(tail -n +2 $file)\t  $(basename $file)" >> {output}; done
    """
