import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np

def merge_taxonomy(cherry_out, phatyp_out, phavip_out,  contigs_fa, out_csv):
	""" Takes the output of different host annotation tools and returns a merged dataframe in csv format with column prefixes."""
	cherryout = cherry_out
	phatypout = phatyp_out
	phavipout = phavip_out
	
	contigsfa = contigs_fa
	
	# Function
	headers = []
	with open(contigsfa, "r") as f:
		for record in SeqIO.parse(f, "fasta"):
			headers.append(record.description)
	
	cherrydf = pd.read_csv(cherryout, delimiter="\t", index_col=0)
	phatypdf = pd.read_csv(phatypout, delimiter="\t", index_col=0)
	phavipdf = pd.read_csv(phavipout, delimiter="\t", index_col=0)
	
	cherrydf = cherrydf.add_prefix("cherry_")
	phatypdf = phatypdf.add_prefix("phatyp_")
	phavipdf = phavipdf.add_prefix("phavip_")
	
	mergedf = pd.DataFrame(index=headers)
	
	mergedf = mergedf.merge(cherrydf, left_index=True, right_index=True, how='outer')
	mergedf = mergedf.merge(phatypdf, left_index=True, right_index=True, how='outer')
	mergedf = mergedf.merge(phavipdf, left_index=True, right_index=True, how='outer')
	
	mergedf = mergedf.replace(["Na", "Unclassified"], np.nan)
	mergedf = mergedf.drop(columns=mergedf.columns[mergedf.columns.str.contains('subfamily', case=False)])
	
	
	mergedf.to_csv(out_csv)
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.description = "Takes the output of different host annotation tools and returns a merged dataframe in csv format with column prefixes."
	parser.add_argument("--cherryout", required=True, type=str, help="CHERRY final taxonomic annotation tsv file")
	parser.add_argument("--phatypout", required=True, type=str, help="PhaTYP final taxonomic annotation tsv file")
	parser.add_argument("--phavipout", required=True, type=str, help="PhaVIP final taxonomic annotation tsv file")
	parser.add_argument("--contigs", required=True, type=str, help="Contigs fasta file used as input for taxonomic annotation")
	parser.add_argument("--output", required=True, type=str, help="Path to merged output file csv")

	args = parser.parse_args()

	merge_taxonomy(args.cherryout, args.phatypout, args.phavipout, args.contigs, args.output)

