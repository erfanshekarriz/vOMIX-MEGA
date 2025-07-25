import os
import re
import argparse

import pandas as pd


# input
#covermout = "results/abundance/intermediate/vOTU_table.tsv"
#methods = "tpm rpkm"

# output
#outdir = "results/abundance/output"
#filename = "vOTU_table"

# function 
def parse_OTU(coverm_out, methods_l, out_dir, file_name):

	df = pd.read_csv(coverm_out, index_col=0, delimiter="\t")
	df.index.name = "vOTU_id"

	pattern = r'\/(.+)_R[12]_'

	for method in methods_l:
		cols = [col for col in df.columns if method.lower() in col.lower()]
		df_f = df[cols]
		df_fr = df_f.rename(columns=lambda col: re.search(pattern, col).group(1))
		outfile = file_name + "_" + method + ".tsv"
		outpath = os.path.join(out_dir, outfile)
		df_fr.to_csv(outpath, sep="\t")


# argument parsing
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.description = "Takes coverm output and a list of abundance methods used, parses as TSV tables and saves multiple OTU table with unique suffixes according to the method. It also cleans column names to contain sample_id name before _R1_ or _R2_ in sequencing read names"

	parser.add_argument("--covermout", required=True, type=str, help="Full path to coverM contig output")
	parser.add_argument("--methods", required=True, type=str, help="Comma separated list of methods used for coverm output (e.g. --methods tpm,rpkm,mean")
	parser.add_argument("--outdir", required=True, type=str, help="Output directory to save files")
	parser.add_argument("--filename", required=True, type=str, help="Output file name without extension")

	args = parser.parse_args()
	methods_l = args.methods.split(",")
	
	parse_OTU(args.covermout, methods_l, args.outdir, args.filename)


