import argparse
import os
import pandas as pd
import json


def parse_json_files(names, csvs):
	merged_data = pd.DataFrame()
	df = pd.DataFrame(columns = ["total_reads_before_filtering",
		"total_reads_after_filtering", 
		"total_bases_before_filtering", 
		"total_bases_after_filtering", 
		"gc_content_before_filtering", "gc_content_after_filtering"], index = names)
	
	for sample_id, csv in zip(names, csvs):
		
		with open(csv, 'r') as jsonf:
			data = json.load(jsonf)
			databefore = data['summary']['before_filtering']
			dataafter = data['summary']['after_filtering']

			df.at[sample_id, 'total_reads_before_filtering'] = databefore['total_reads']
			df.at[sample_id, 'total_bases_before_filtering'] = databefore['total_bases']
			df.at[sample_id, 'gc_content_before_filtering'] = databefore['gc_content']
			df.at[sample_id, 'total_reads_after_filtering'] = dataafter['total_reads']
			df.at[sample_id, 'total_bases_after_filtering'] = dataafter['total_bases']
			df.at[sample_id, 'gc_content_after_filtering'] = dataafter['gc_content']
	
	# Print merged table to STDOUT
	print(df.to_csv(index=True))



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Merge outputs of viral contig identification tools output')
	parser.add_argument('--names', type=str, help='Names of all sample IDs')
	parser.add_argument('--jsons', type=str, help='CSV output files of merged viral contig identification tools')
	
	args = parser.parse_args()

	with open(args.names, "r") as namesfile:
		namestring = namesfile.read()

	with open(args.jsons, "r") as pathsfile:
		jsonstring = pathsfile.read()

	parse_json_files(namestring.split(), jsonstring.split())
