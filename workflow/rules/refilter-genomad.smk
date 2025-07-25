logdir=relpath("identify/viral/logs")
benchmarks=relpath("identify/viral/benchmarks")
tmpd=relpath("identify/viral/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

if isinstance(config['cores'], int):
  n_cores = config['cores']
else:
   console.print(Panel.fit(f"config['cores'] is not an integer: {config['cores']}, you can change the parameter in config/config.yml file", title="Error", subtitle="config['cores'] not integer"))
   sys.exit(1)


############################
# Single-Sample Processing #
############################

# the "contigfile" is not the nested one

if config['inputdir']!="":

  indir = cleanpath(config['inputdir'])
  console.print(f"\nconfig['inputdir'] not empty, using '{indir}' as input for viral contig identification.")
  console.print("File names without the .fa extension will be used as sample IDs.")
  cwd = os.getcwd()
  indir_path = os.path.join(cwd, indir)

  if not os.path.exists(indir):
    console.print(Panel.fit(f"The input file path '{indir}' does not exist.", title="Error", subtitle="Contig Directory"))
    sys.exit(1)

  fasta_files = [f for f in os.listdir(indir_path) if f.endswith('.fa')]
  if len(fasta_files) == 0:
    console.print(Panel.fit(f"There are no files ending with .fa in '{indir_path}', other fasta extensions are not accepted (for now) :(", title="Error", subtitle="No .fa Files"))
    sys.exit(1)

  assembly_ids = [os.path.basename(fasta_file).rsplit(".", 1)[0] for fasta_file in fasta_files]
  wildcards_p = os.path.join(indir, "{sample_id}.fa")
  outdir_p = os.path.join(cwd, relpath("identify/viral/"))
  console.print(f"Creating output directory: '{outdir_p}'.\n")

else:
  wildcards_p = relpath("assembly/samples/{sample_id}/output/final.contigs.fa")
  assembly_ids = assemblies.keys()

### Delete all pre-existing files

if os.path.exists(os.path.join(logdir, "done.log")):
  os.remove(os.path.join(logdir, "done.log"))
  os.remove(os.path.join(logdir, "clustering-fast-done.log"))
  os.remove(os.path.join(logdir, "checkv-done.log"))
else:
  pass

files_to_delete = expand(relpath("identify/viral/samples/{sample_id}/output/merged_scores_filtered.csv"), sample_id=assembly_ids)
for file in files_to_delete:
    try:
        # Check if the file exists before attempting to delete it
        if os.path.exists(file):
            os.remove(file)
        else:
          pass
    except OSError as e:
        print(f"Error deleting file '{file}': {e}")

### MASTER RULE 

rule done_log:
  name: "genomad-refilter.smk Done. removing tmp files"
  localrule: True
  input:
    relpath("identify/viral/output/combined.final.vOTUs.fa"), 
    expand(relpath("identify/viral/samples/{sample_id}/output/merged_scores_filtered.csv"), sample_id=assembly_ids)
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
rule filter_contigs:
  name: "genomad-refilter.smk filter contigs [length]"
  localrule: True
  input:
    wildcards_p
  output:
    relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  params:
    minlen=config['contigminlen'],
    outdir=relpath("identify/viral/samples/{sample_id}/tmp"),
    tmpdir=os.path.join(tmpd, "contigs/{sample_id}")
  log: os.path.join(logdir, "filtercontig_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}/* {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}
    
    seqkit seq {input} --min-len {params.minlen} > {params.tmpdir}/tmp.fa

    mv {params.tmpdir}/tmp.fa {output}
    """

rule genomad_filter:
  name : "genomad-refilter.smk filter geNomad output [--genomad-only]"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
    tsv=relpath("identify/viral/samples/{sample_id}/intermediate/genomad/final.contigs.filtered_summary/final.contigs.filtered_virus_summary.tsv"),
  output:
    fna=relpath("identify/viral/samples/{sample_id}/output/viral.contigs.fa"),
    scrs=relpath("identify/viral/samples/{sample_id}/output/merged_scores_filtered.csv"),
    hits=relpath("identify/viral/samples/{sample_id}/output/viralhits_list")
  params:
    script="workflow/scripts/identify/viral/genomad_filter.py", 
    minlen=config['genomadminlen'],
    cutoff=config['genomadcutoff_p'],
    outdir=relpath("identify/viral/samples/{sample_id}/output/"),
    tmpdir=os.path.join(tmpd, "filter/{sample_id}")
  log: os.path.join(logdir, "genomad_filter_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir} {params.outdir}
    
    python {params.script} \
        --genomad_out {input.tsv} \
        --genomad_min_score {params.cutoff} \
        --genomad_min_len {params.minlen} \
        --output_path {params.tmpdir}/tmp.csv \
        --hitlist_path {params.tmpdir}/tmplist &> {log}

    mv {params.tmpdir}/tmp.csv {output.scrs}
    cat {params.tmpdir}/tmplist | uniq > {output.hits}
    
    seqkit grep {input.fna} -f {output.hits} | seqkit replace -p  "\s.*" -r "" | seqkit replace -p $ -r _{wildcards.sample_id}  > {params.tmpdir}/tmp.fa 2> {log}
    mv {params.tmpdir}/tmp.fa {output.fna}

    rm -rf {params.tmpdir}/*
    """
    

rule cat_contigs:
  name : "genomad-refilter.smk combine viral contigs"
  input:
    fna=expand(relpath("identify/viral/samples/{sample_id}/output/viral.contigs.fa"), sample_id=assembly_ids),
    scrs=expand(relpath("identify/viral/samples/{sample_id}/output/merged_scores_filtered.csv"), sample_id=assembly_ids)
  output: 
    fna=relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa"),
    scrs=relpath("identify/viral/intermediate/scores/combined_viral_scores.csv")
  params:
    script="workflow/scripts/identify/viral/mergeout_scores.py", 
    names=list(assembly_ids),
    tmpdir=tmpd
  log: os.path.join(logdir, "catcontigs.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}
    
    echo "{params.names}" > {params.tmpdir}/tmp.names
    echo "{input.scrs}" > {params.tmpdir}/tmp.csv.paths
    
    python {params.script} \
        --names {params.tmpdir}/tmp.names \
        --csvs {params.tmpdir}/tmp.csv.paths > {params.tmpdir}/tmp.csv 2> {log}
    cat {input.fna} > {params.tmpdir}/tmp.fa 2> {log}

    mv {params.tmpdir}/tmp.csv {output.scrs}
    mv {params.tmpdir}/tmp.fa {output.fna}

    rm -rf {params.tmpdir}
    """


rule combine_classifications:
  name: "genomad-refilter.smk combine derepped classification results"
  input:
    checkv_out=relpath("identify/viral/output/checkv/quality_summary.tsv"),
    classify_out=relpath("identify/viral/intermediate/scores/combined_viral_scores.csv")
  output:
    relpath("identify/viral/output/checkv/combined_classification_results.csv")
  params:
    script="workflow/scripts/identify/viral/combineclassify.py",
    tmpdir=tmpd
  log: os.path.join(logdir, "combine_classification.log")
  threads: 1
  conda : "../envs/seqkit-biopython.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    python {params.script} \
        --mergedclassify {input.classify_out} \
        --checkvsummary {input.checkv_out} \
        --output {params.tmpdir}/tmp.csv 2> {log}

    mv {params.tmpdir}/tmp.csv {output}
    """


# THEN GOES INTO 
# 1) VIRAL BINNING [optional]
# 2) CLUSTERING [sensitive or fast]
# 3) CHECKV-PYHMMER 
# THEN COMES BACK HERE TO GET vCONTIGS


rule consensus_filtering:
  name: "genomad-refilter.smk consensus vOTU filtering"
  input:
    relpath("identify/viral/output/checkv/combined_classification_results.csv")
  output:
    summary=relpath("identify/viral/output/classification_summary_vOTUs.csv"),
    proviruslist=relpath("identify/viral/output/provirus.list.txt"),
    viruslist=relpath("identify/viral/output/virus.list.txt")
  params:
    script="workflow/scripts/identify/viral/consensus_filtering_genomad.py",
    genomad=config['genomadcutoff_s'],
    tmpdir=tmpd
  log: os.path.join(logdir, "consensus_filtering.log")
  threads: 1
  conda: "../envs/seqkit-biopython.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    python {params.script} \
        --classification_results {input} \
        --genomad_min_score {params.genomad} \
        --summary_out {params.tmpdir}/tmp.csv \
        --provirus_list {params.tmpdir}/tmp.list.1 \
        --virus_list {params.tmpdir}/tmp.list.2  2> {log}

    mv {params.tmpdir}/tmp.csv {output.summary}
    mv {params.tmpdir}/tmp.list.1 {output.proviruslist}
    mv {params.tmpdir}/tmp.list.2 {output.viruslist}
    """


rule votu:
  name: "genomad-refilter.smk generate final vContigs"
  input:
    provirusfasta=relpath("identify/viral/output/checkv/proviruses.fna"),
    virusfasta=relpath("identify/viral/output/checkv/viruses.fna"), 
    provirushits=relpath("identify/viral/output/provirus.list.txt"),
    virushits=relpath("identify/viral/output/virus.list.txt")
  output:
    combined=relpath("identify/viral/output/combined.final.vOTUs.fa"),
    provirus=relpath("identify/viral/output/provirus.final.vOTUs.fa"),
    virus=relpath("identify/viral/output/virus.final.vOTUs.fa")
  params:
    outdir=relpath("identify/viral/output/"),
    tmpdir=tmpd
  log: os.path.join(logdir, "vOTUs.log")
  threads: 1
  conda: "../envs/seqkit-biopython.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir} {params.outdir}

    seqkit replace {input.provirusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} > {params.tmpdir}/tmp1.fa 2> {log}
    seqkit replace {input.virusfasta} --f-use-regexp -p "(.+)_\d\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} >> {params.tmpdir}/tmp1.fa 2> {log}

    seqkit grep {input.virusfasta} -f {input.virushits} > {params.tmpdir}/tmp2.fa 2> {log}
    seqkit grep {input.provirusfasta} -f {input.virushits} >> {params.tmpdir}/tmp2.fa 2> {log}

    cat {params.tmpdir}/tmp1.fa {params.tmpdir}/tmp2.fa  > {params.tmpdir}/tmp3.fa 2> {log}

    seqkit rmdup {params.tmpdir}/tmp1.fa > {output.provirus} 2> {log}
    seqkit rmdup {params.tmpdir}/tmp2.fa > {output.virus} 2> {log}
    seqkit rmdup {params.tmpdir}/tmp3.fa > {output.combined} 2> {log}

    rm -rf {params.tmpdir}/*
    """


