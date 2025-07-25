#!/usr/bin/env python

import argparse
import pandas as pd
import os

def filter_csv(genomad_out, genomad_score_cutoff, genomad_min_len, output_path, hitlist_path):
	
	# Read the CSV file into a DataFrame
	df = pd.read_csv(genomad_out, sep="\t")
	df = df[(df['length'] >= genomad_min_len) & (df['virus_score'] >= genomad_score_cutoff)]

	df['seq_name'] = df['seq_name'].str.replace(r'\|provir.+$', '', regex=True)
	df.columns = ['genomad_' + col for col in df.columns]
	df['sequence_id'] = df['genomad_seq_name']


	# Save filtered data to CSV
	output_dir = os.path.dirname(output_path)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	df.to_csv(output_path, index=False)
	df['sequence_id'].to_csv(hitlist_path, index=False, header=False)

if __name__ == '__main__':
	# Parse command-line arguments
	parser = argparse.ArgumentParser(description='viral contig identification CSV filtering script based on software cutoffs')
	parser.add_argument('--genomad_out', help='Path to the CSV file of merged output')
	parser.add_argument('--genomad_min_score', type=float, default=0.7, help='genomad_min_score cutoff')
	parser.add_argument('--genomad_min_len', type=float, default=0, help='genomad_min_length cutoff')
	parser.add_argument('--output_path', help='Path to save the filtered CSV file')
	parser.add_argument('--hitlist_path', help='Path to save positive hits sequence id as a TXT file')
	args = parser.parse_args()
	
	# Filter the CSV file based on cutoff scores and save filtered data
	filter_csv(args.genomad_out, args.genomad_min_score, args.genomad_min_len, args.output_path, args.hitlist_path)



