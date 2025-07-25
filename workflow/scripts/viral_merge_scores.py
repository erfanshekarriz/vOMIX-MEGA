import argparse
import os
import pandas as pd


def merge_csv_files(names, csvs):
	merged_data = pd.DataFrame()
	data_frames = []
	for sample_id, csv in zip(names, csvs):
		df = pd.read_csv(csv, index_col=0)
		df['sample_id'] = sample_id  # Add a new column with the file name
		df.index = df.index.astype(str) + "_" + str(sample_id)
		data_frames.append(df)

	merged_data = pd.concat(data_frames, ignore_index=False)
	# Print merged table to STDOUT
	print(merged_data.to_csv(index=True))



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Merge outputs of viral contig identification tools output')
	parser.add_argument('--names', type=str, help='Names of all sample IDs')
	parser.add_argument('--csvs', type=str, help='CSV output files of merged viral contig identification tools')
	
	args = parser.parse_args()

	with open(args.names, "r") as namesfile:
		namestring = namesfile.read()
	
	with open(args.csvs, "r") as pathsfile:
		csvstring = pathsfile.read()

	merge_csv_files(namestring.split(), csvstring.split())

