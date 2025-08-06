import os

logdir=os.path.join(config['basedir'], "workflow/database/.logs")
benchmarks=os.path.join(config['basedir'], "workflow/database/.benchmarks")
tmpd=os.path.join(config['basedir'], "workflow/database/.tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

n_cores = config['max-cores']


rule done:
  name: "setup-database.smk Done."
  localrule: True
  input:
    os.path.join(config['genomad-db'], "genomad_db.source"), 
    expand("workflow/database/checkv/hmm_db/checkv_hmms/{index}.hmm", index=range(1, 81)), 
    os.path.join(config['PhaBox2-db'], "genus2hostlineage.pkl"), 
    os.path.join(config["GTDBTk-db"], ("gtdbtk_r" + config["GTDBTk-db-version"] + "_data.tar.gz"))
  output:
    os.path.join(benchmarks, "done.log")
  shell:
    """
    touch {output}
    """


rule genomad_db:
  name: "setup-database.smk geNomad database (1.3 G)"
  localrule: True
  output: os.path.join(config['genomad-db'], "genomad_db.source")
  params:
    outdir=config['genomad-db'],
    tmpdir=os.path.join(tmpd, "genomad/db")
  log: os.path.join(logdir, "genomad_db.log")
  benchmark: os.path.join(benchmarks, "genomad_db.log")
  conda: "../envs/genomad.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    genomad download-database {params.tmpdir} &> {log}

    mv {params.tmpdir}/genomad_db/* {params.outdir}
    """


rule checkv_db:
  name: "setup-database.smk CheckV database (7.3 G)"
  localrule: True
  output: expand("workflow/database/checkv/hmm_db/checkv_hmms/{index}.hmm", index=range(1, 81))
  params:
    outdir=config['checkv-db'],
    tmpdir=os.path.join(tmpd, "checkv/db")
  log: os.path.join(logdir, "checkv_db.log")
  benchmark: os.path.join(benchmarks, "checkv_db.log")
  conda: "../envs/checkv.yml"
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    checkv download_database {params.tmpdir} &> {log}

    mv {params.tmpdir}/checkv-db*/* {params.outdir}
    """


rule phabox2_db:
  name: "setup-database.smk PhaBox2 database (1.6 G)"
  localrule: True
  output: os.path.join(config['PhaBox2-db'], "genus2hostlineage.pkl")
  params:
    outdir=config['PhaBox2-db'],
    tmpdir=os.path.join(tmpd, "PhaBox2/db")
  log: os.path.join(logdir, "PhaBox2_db.log")
  benchmark: os.path.join(benchmarks, "PhaBox2_db.log")
  conda: "../envs/phabox2.yml"
  threads: 1
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}/phabox_db_v2
    mkdir -p {params.tmpdir} {params.outdir}

    wget -O {params.tmpdir}/phabox_db_v2.zip https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2.zip &> {log}
    unzip -o {params.tmpdir}/phabox_db_v2.zip -d {params.tmpdir} > {log}

    mv {params.tmpdir}/phabox_db_v2/* {params.outdir}
    """


rule virsorter2_db:
  name: "setup-database.smk VirSorter2 setup database (X.X G)"
  output: os.path.join(config['virsorter2-db'], "db.tgz")
  params:
    outdir=config['virsorter2-db'],
    tmpdir=os.path.join(tmpd, "virsorter2/db")
  log: os.path.join(logdir, "virsorter2_database.log")
  benchmark: os.path.join(benchmarks, "virsorter2_db.log")
  conda: "../envs/virsorter2.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 8000
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir}

    virsorter setup \
        -d {params.tmpdir} \
        -j {threads} \
        -s 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    """

rule vibrant_db:
  name: "setup-database.smk VirSorter2 setup database (11.1 G)"
  output: os.path.join(config['vibrant-db'], "files/VIBRANT_machine_model.sav")
  params:
    outdir=config['vibrant-db'],
    tpmdir=os.path.join(tmpd, "vibrant/db")
  log: os.path.join(logdir, "vibrant_database.log")
  benchmark: os.path.join(benchmarks, "vibrant_database.log")
  conda: "../envs/vibrant.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 8000
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir}

    download-db.sh {params.tmpdir}
    python3 {params.tmpdir}/database/VIBRANT_setup.py

    mv {params.tmpdir}/* {params.output}
    """


rule chocophlan_db:
  name: "setup-database.smk HUMAnN3 chocophlan database (16.4 G)"
  output: os.path.join(config['humann-db'], "chocophlan/alaS.centroids.v201901_v31.ffn.gz")
  params:
    outdir=config['humann-db'],
    tmpdir=os.path.join(tmpd, "humann/db")
  log: os.path.join(logdir, "humann_chocophlan_db.log")
  benchmark: os.path.join(benchmarks, "humann_chocophlan_db.log")
  conda: "../envs/biobakery3.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 4000
  shell:
    """
    rm -rf {params.outdir}/chocophlan {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}

    humann_databases --download chocophlan full {params.tmpdir}/ --update-config yes &> {log}

    mv {params.tmpdir}/chocophlan {params.outdir}/
    """


rule uniref_db:
  name: "setup-database.smk HUMAnN3 uniref database (19.7 G)"
  output: os.path.join(config['humann-db'], "uniref/uniref90_201901b_full.dmnd")
  params:
    outdir=config['humann-db'],
    tmpdir=os.path.join(tmpd, "humann/db")
  log: os.path.join(logdir, "humann_uniref_db.log")
  benchmark: os.path.join(benchmarks, "humann_uniref_db.log")
  conda: "../envs/biobakery3.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 4000
  shell:
    """
    rm -rf {params.outdir}/uniref {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}
 
    humann_databases --download uniref uniref90_diamond {params.tmpdir}/ --update-config yes &> {log}

    mv {params.tmpdir}/uniref {params.outdir}/
    """


rule utilitymap_db:
  name: "setup-database.smk HUMAnN3 utility mapping database (2.5 G)"
  output: 
    ecdb=os.path.join(config['humann-db'], "utility_mapping/map_level4ec_uniref90.txt.gz"),
    eggnogdb=os.path.join(config['humann-db'], "utility_mapping/map_eggnog_uniref90.txt.gz"),
    godb=os.path.join(config['humann-db'], "utility_mapping/map_go_uniref90.txt.gz"),
    kodb=os.path.join(config['humann-db'], "utility_mapping/map_ko_uniref90.txt.gz"),
    pfamdb=os.path.join(config['humann-db'], "utility_mapping/map_pfam_uniref90.txt.gz"),
    ecnamedb=os.path.join(config['humann-db'], "utility_mapping/map_ec_name.txt.gz"),
    eggnognamedb=os.path.join(config['humann-db'], "utility_mapping/map_eggnog_name.txt.gz"),
    gonamedb=os.path.join(config['humann-db'], "utility_mapping/map_go_name.txt.gz"),
    konamedb=os.path.join(config['humann-db'], "utility_mapping/map_ko_name.txt.gz"),
    pfamnamedb=os.path.join(config['humann-db'], "utility_mapping/map_pfam_name.txt.gz")
  params:
    outdir=config['humann-db'],
    tmpdir=os.path.join(tmpd, "humann/db")
  log: os.path.join(logdir, "humann_utility_db.log")
  benchmark: os.path.join(benchmarks, "humann_utility_db.log")
  conda: "../envs/biobakery3.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 4000
  shell:
    """
    rm -rf {params.outdir}/utility_mapping {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}
 
    humann_databases --download utility_mapping full {params.tmpdir}/ --update-config yes &> {log}

    mv {params.tmpdir}/utility_mapping {params.outdir}/
    """


rule eggnog_download:
  name: "setup-database.smk eggNOG-mapper v2 Database (43.0 G)"
  output:
    os.path.join(config['eggNOG-db'], "eggnog.db"),
    os.path.join(config['eggNOG-db'], "eggnog_proteins.dmnd")
  params:
    dbdir=config['eggNOG-db'],
    outdir="workflow/database/eggNOGv2",
    tmpdir=os.path.join(tmpd, "eggNOGdb")
  conda: "../envs/eggnog-mapper.yml"
  log: os.path.join(logdir, "eggNOGv2_db.log")
  benchmark: os.path.join(benchmarks, "eggNOGv2_db.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 8 * 10**3
  shell:
    """
    rm -r {params.outdir}
    mkdir -p {params.tmpdir} {params.dbdir}
    
    download_eggnog_data.py -y --data_dir {params.tmpdir} 2> {log}
    mv {params.tmpdir}/* {params.outdir}/

    rm -rf {params.tmpdir}
    """

rule checkm2_download:
  name: "setup-database.smk CheckM2 Database (2.7 G)"
  output: 
    os.path.join(config["checkm2-db"], "uniref100.KO.1.dmnd")
  params:
    outdir=config['checkm2-db'],
    tmpdir=os.path.join(tmpd, "checkm2db")
  conda: "../envs/checkm2.yml"
  log: os.path.join(logdir, "checkm2_db.log")
  benchmark: os.path.join(benchmarks, "checkm2_db.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4 * 10**3
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} 
    
    checkm2 database --download --path {params.tmpdir} 2> {log}

    mv {params.tmpdir}/CheckM2_database/* {params.outdir}/
    """


gtdbtk_db_base_url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release"
gtdbtk_db_version = config["GTDBTk-db-version"]
gtdbtk_dwnld_lnk =f"{gtdbtk_db_base_url}{gtdbtk_db_version}/{gtdbtk_db_version}.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r{gtdbtk_db_version}_data.tar.gz"

rule GTDBTk_download:
  name: "setup-database.smk GTDB-Tk Database (63.3 G)"
  output:
    os.path.join(config["GTDBTk-db"], ("gtdbtk_r" + config["GTDBTk-db-version"] + "_data.tar.gz"))
  params:
    outdir=config["GTDBTk-db"], 
    dbdir=os.path.join(config["GTDBTk-db"], "release" + config["GTDBTk-db-version"]),
    dwnldlink=gtdbtk_dwnld_lnk,
    tmpdir=os.path.join(tmpd, "gtdbtkdb")
  conda: "../envs/gtdbtk.yml"
  log: os.path.join(logdir, "gtdbtk_db.log")
  benchmark: os.path.join(benchmarks, "gtdbtk_db.log")
  threads: 1
  resources: 
    mem_mb=lambda wildcards, attempt: attempt * 4 * 10**3
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} 

    wget -P {params.tmpdir} {params.dwnldlink} 2> {log}
    #download-db.sh {params.tmpdir} 2> {log}
    #touch {output}

    mv {params.tmpdir}/* {params.outdir}/
    tar -xzvf {output} -C {params.outdir}
    conda env config vars set GTDBTK_DATA_PATH="{params.dbdir}"
    """


rule iPHoP_download:
  name: "setup-database.smk iPHoP Database (318.1 G)"
  output: 
    os.path.join(config["iphop-db"], config["iphop-db-basename"], "md5checkfile.txt")
  params:
    outdir=os.path.abspath(config["iphop-db"]),
    dbversion=config["iphop-db-version"],
    tmpdir=os.path.join(tmpd, "iphop-db")
  conda: "../envs/iphop.yml"
  log: os.path.join(logdir, "iphop_db.log")
  benchmark: os.path.join(benchmarks, "iphop_db.log")
  threads: 1
  resources: 
    mem_mb=lambda wildcards, attempt: attempt * 4 * 10**3
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}

    iphop download \
        --db_dir {params.tmpdir} \
        --db_version {params.dbversion} \
        --no_prompt

    mv {params.tmpdir}/* {params.outdir}
    """
