#!/usr/bin/env python3

# this code was directly retrieved from https://github.com/EBI-Metagenomics/emg-viral-pipeline/blob/master/bin/viral_contigs_annotation.py
# it was modified to take output file name and path instead of just directory 

import argparse
import operator
import re
from pathlib import Path

from Bio import SeqIO

import pandas as pd


def prot_annot_tbl(protein_file, ratio_evalue_file):
	"""This function takes a fasta file containing the proteins predicted in a
	set of putative viral contigs and a dataframe that collates the results
	obtained with hmmscan against the ViPhOG database for the same proteins
	NOTE: if the protein_file is empty it will create an empty tsv
	"""
	annotation_list = []
	ratio_evalue_df = pd.read_csv(ratio_evalue_file, sep="\t")
	if Path(protein_file).exists():
		for protein in SeqIO.parse(protein_file, "fasta"):
			contig_id = re.split(r"_\d+$", protein.id)[0]
			protein_prop = protein.description.split(" # ")[:-1]
			if protein_prop[0] in ratio_evalue_df["query"].values:
				filtered_df = ratio_evalue_df[
					ratio_evalue_df["query"] == protein_prop[0]
				]
				if len(filtered_df) > 1:
					best_value_index = max(
					    filtered_df["Abs_Evalue_exp"].items(),
					    key=operator.itemgetter(1),
					)[0]
					protein_prop.extend(
					    list(
					        filtered_df.loc[
					            best_value_index, ["ViPhOG", "Abs_Evalue_exp", "Taxon"]
					        ]
					    )
					)
				else:
					protein_prop.extend(
					    list(
					        filtered_df.loc[
					            filtered_df.index[0],
					            ["ViPhOG", "Abs_Evalue_exp", "Taxon"],
					        ]
					    )
					)
			else:
				protein_prop.extend(["No hit", "NA", ""])
			annotation_list.append([contig_id] + protein_prop)
	protein_annot_df = pd.DataFrame(
		annotation_list,
		columns=[
			"Contig",
			"CDS_ID",
			"Start",
			"End",
			"Direction",
			"Best_hit",
			"Abs_Evalue_exp",
			"Label",
		],
	)
	return protein_annot_df


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Generate tabular file with ViPhOG annotation results for proteins predicted in viral contigs"
	)
	parser.add_argument(
		"-p",
		"--prot",
		dest="prot_file",
		help="Relative or absolute path to protein file of predicted viral contigs",
		required=True,
	)
	parser.add_argument(
		"-t",
		"--table",
		dest="ratio_file",
		help="Relative or absolute path to ratio_evalue tabular file generated for predicted viral contigs",
		required=True,
	)
	parser.add_argument(
		"-o",
		"--outpath",
		dest="output_path",
		help="Default name with full relative or absolute path that you'd like to save the results with"
	)
	args = parser.parse_args()

	output_name = Path(args.prot_file).stem
	csv_output = Path(args.output_path)

	final_df = prot_annot_tbl(args.prot_file, args.ratio_file)
	final_df.to_csv(csv_output, sep="\t", index=False)

