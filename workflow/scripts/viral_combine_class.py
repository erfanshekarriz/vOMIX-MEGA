import sys
import argparse
import os
import pandas as pd


def merge_classifications(classify_csv, checkv_tsv, output_file):

	# read files 
	classifydf = pd.read_csv(classify_csv, index_col=0)
	checkvdf = pd.read_csv(checkv_tsv, delimiter='\t', index_col=0)
	
	# add checkv to the column names 
	checkvdf.columns = ['checkv_' + col if not col.startswith("checkv") else col for col in checkvdf.columns]
	
	# merge the two results
	merged = checkvdf.merge(classifydf, how='left', left_index=True, right_index=True)

	merged.to_csv(output_file, index=True)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Merge outputs of viral contig identification')

	parser.add_argument('--mergedclassify', type=str, help='combined classification output path')
	parser.add_argument('--checkvsummary', type=str, help='checkv summary_quality.tsv file path')
	parser.add_argument('--output', type=str, help='output csv file path and name')

	args = parser.parse_args()

	merge_classifications(args.mergedclassify, args.checkvsummary, args.output)


	


