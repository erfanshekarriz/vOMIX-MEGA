os.makedirs(relpath(".vomix/log"), exist_ok=True)

email=config["email"]
api_key=config["NCBI-API-key"]
nowstr=config["latest-run"]
outdir=config["outdir"]
datadir=config["datadir"]

samples, assemblies = parse_sample_list(config["workdir"] + config["samplelist"], datadir, outdir, email, api_key, nowstr)

###################
# DELETE SYMLINKS #
###################

symlinklist = []
symlinklist += [relpath(".vomix/log/symlink_done.log")]
symlinklist += expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R{i}.fastq.gz"), 
    sample_id=samples.keys(), i=[1,2])

for symlink in symlinklist:
  try:
    os.remove(symlink)
  except OSError:
    pass



###################
# REMAKE SYMLINKS #
###################

rule done:
  name: "symlink.py Done. fixed all broken symlinks"
  localrule: True
  input: 
    expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R{i}.fastq.gz"), sample_id=samples.keys(), i=[1,2])
  output:
    relpath(".vomix/log/symlink_done.log")
  shell:
    """
    touch {output}
    """


rule symlink:
  name: "preprocessing.py creating symbolic links"
  localrule: True
  input:
    R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz")
  output:
    R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"), 
  shell:
    """
    ln -s $(pwd)/{input.R1} $(pwd)/{output.R1}
    ln -s $(pwd)/{input.R2} $(pwd)/{output.R2}
    """

