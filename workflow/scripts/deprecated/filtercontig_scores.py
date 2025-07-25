#!/usr/bin/env python

import argparse
import pandas as pd

def filter_csv(csv_path, genomad_score_cutoff, dvf_score_cutoff, dvf_pvalue_cutoff, phamer_pred_value, phamer_score_cutoff, output_path, hitlist_path):
	# Read the CSV file into a DataFrame
	df = pd.read_csv(csv_path)
	# Apply filters based on cutoff scores
	filtered_df = df[
			((df['genomad_virus_score'] >= genomad_score_cutoff) & (df['dvf_score'] >= dvf_score_cutoff) 
				& (df['dvf_pvalue'] <= dvf_pvalue_cutoff)) |
			((df['genomad_virus_score'] >= 0.9) | (df['dvf_score'] >= 0.9)) |
			((df['phamer_Pred'] == phamer_pred_value) & (df['phamer_Score'] >= phamer_score_cutoff))
			]

	# Save filtered data to CSV
	filtered_df.to_csv(output_path, index=False)
	filtered_df['sequence_id'].to_csv(hitlist_path, index=False, header=False)

if __name__ == '__main__':
	# Parse command-line arguments
	parser = argparse.ArgumentParser(description='viral contig identification CSV filtering script based on software cutoffs')
	parser.add_argument('--csv_path', help='Path to the CSV file of merged output')
	parser.add_argument('--genomad_min_score', type=float, default=0.7, help='genomad_min_score cutoff')
	parser.add_argument('--dvf_min_score', type=float, default=0.9, help='dvf_min_score cutoff')
	parser.add_argument('--dvf_max_pval', type=float, default=0.05, help='dvf_max_pvalue cutoff')
	parser.add_argument('--phamer_pred', default='phage', help='phamer_Pred value')
	parser.add_argument('--phamer_min_score', type=float, default=0.9, help='phamer_min_score cutoff')
	parser.add_argument('--output_path', help='Path to save the filtered CSV file')
	parser.add_argument('--hitlist_path', help='Path to save positive hits sequence id as a TXT file')
	args = parser.parse_args()
	
	# Filter the CSV file based on cutoff scores and save filtered data
	filter_csv(args.csv_path, args.genomad_min_score, args.dvf_min_score, args.dvf_max_pval, args.phamer_pred, args.phamer_min_score, args.output_path, args.hitlist_path)



