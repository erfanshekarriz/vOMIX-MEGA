#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np

def merge_files(genomad_path, dvfout_path, phamerout_path, genomad_min_len, dvf_min_len, phamer_min_len, output_file):
	
	
	# Read the input files
	try:
		genomadout = pd.read_csv(genomad_path, delimiter='\t')
	except FileNotFoundError:
		columns = ['seq_name', 'length', 'topology', 'coordinates', 'n_genes', 'genetic_code', 'virus_score', 'fdr', 'n_hallmarks', 'marker_enrichment', 'taxonomy']
		genomadout = pd.DataFrame({}, columns=columns)
		
	#vs2out = pd.read_csv(vs2out_path, delimiter='\t')
	try:
		dvfout = pd.read_csv(dvfout_path, delimiter='\t')
	except FileNotFoundError:
		columns = ['name', 'len', 'score', 'pvalue']
		dvfout = pd.DataFrame({}, columns=columns)
	#virbotout = pd.read_csv(virbotout_path)
	try:
		phamerout = pd.read_csv(phamerout_path)
	except FileNotFoundError:
		columns = ['Accession', 'Length', 'Pred', 'Score']
		phamerout = pd.DataFrame({}, columns=columns)


	# length filtering 
	genomadout = genomadout[genomadout['length'] >= genomad_min_len]
	dvfout = dvfout[dvfout['len'] >= dvf_min_len]
	
	# in some instances if all contigs are filtered for phamer, it will 
	# not produce a length column, so we must consider that. 

	if 'Length' in phamerout.columns:
		phamerout = phamerout[phamerout['Length'] >= phamer_min_len]

	# Add prefixes to column title
	#vs2out['seqname'] = vs2out['seqname'].str.replace(r'\|\|.*', '')
	genomadout['seq_name'] = genomadout['seq_name'].str.replace(r'\|provir.+$', '', regex=True)
	dvfout['name'] = dvfout['name'].str.replace(r'\sflag.+$', '', regex=True)

	genomadout.columns = ['genomad_' + col for col in genomadout.columns]
	#vs2out.columns = ['vs2_' + col for col in vs2out.columns]
	dvfout.columns = ['dvf_' + col for col in dvfout.columns]
	#virbotout.columns = ['virbot_' + col for col in virbotout.columns]
	phamerout.columns = ['phamer_' + col for col in phamerout.columns]

	# Create the 'sequence_id' column
	genomadout['sequence_id'] = genomadout['genomad_seq_name']
	#vs2out['sequence_id'] = vs2out['vs2_seqname']
	dvfout['sequence_id'] = dvfout['dvf_name']
	#virbotout['sequence_id'] = virbotout['virbot_Contig_acc']
	phamerout['sequence_id'] = phamerout['phamer_Accession']


	# Filter based on length 


	# Merge the files based on the 'sequence_id' column using outer join
	merged = pd.merge(genomadout, dvfout, on='sequence_id', how='outer')
	merged = pd.merge(merged, phamerout, on='sequence_id', how='outer')

	# Replace missing values with "NA"
	merged = merged.replace(np.nan, 'NA')
	merged.set_index('sequence_id', inplace=True) # Set sequence_id as rownames
	merged.sort_index(inplace=True)

	# Save the merged dataframe as a CSV file
	output_dir = os.path.dirname(output_file)
	if not os.path.exists(output_dir):
		    os.makedirs(output_dir)

	merged.to_csv(output_file, index=True)
	print("Merged dataframe saved to", output_file)

if __name__ == '__main__':
	# Create the argument parser
	parser = argparse.ArgumentParser(description='Merge outputs of viral contig identification')
	parser.add_argument('--genomadout', type=str, help='Path to the vs2out TSV file')
	parser.add_argument('--dvfout', type=str, help='Path to the dvfout TXT file')
	parser.add_argument('--phamerout', type=str, help='Path to the phamerout CSV file')
	parser.add_argument('--genomadminlen', type=float, help='Min length for genomad annotation')
	parser.add_argument('--dvfminlen', type=float, help='Min length for dvf annotation')
	parser.add_argument('--phamerminlen', type=float, help='Min length for phamer annotation')
	parser.add_argument('--output', type=str, default='merged_output.csv', help='Output file name (default: merged_output.csv)')

	# Parse the command-line arguments
	args = parser.parse_args()

	# Merge the files
	merge_files(args.genomadout, args.dvfout, args.phamerout, args.genomadminlen, args.dvfminlen, args.phamerminlen, args.output)
