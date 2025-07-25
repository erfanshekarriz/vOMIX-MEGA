import os 
import argparse
import pandas as pd


def extract_taxonomy(input_f, output_f):

	_, extension = os.path.splitext(input_f)
	if extension.lower() == ".tsv":
		delimiter = "\t"
	elif extension.lower() == ".csv":
		delimiter = ","
	else:
		raise ValueError(f"Unsupported file extension: {extension}")


	df = pd.read_csv(input_f, index_col=0, delimiter=delimiter)
	try:
		taxdf = df.loc[:, ['genomad_taxonomy', 'genomad_virus_score']]
	except KeyError:
		try:
			#print(df)
			taxdf = df.loc[:, ['taxonomy', 'virus_score']]
			taxdf.rename(columns={'taxonomy': 'genomad_taxonomy', 'virus_score': 'genomad_virus_score'}, inplace=True)
		except KeyError:
			raise ValueError("Neither 'genomad_taxonomy' nor 'taxonomy' column exists in the DataFrame.")
	taxdf['genomad_taxonomy'] = taxdf['genomad_taxonomy'].replace('Unclassified', 'Unassigned')

	# Assuming your DataFrame is called df
	taxdf[['virus', 'realm', 'kingdom', 'phylum', 'class', 'order', 'family']] = taxdf['genomad_taxonomy'].str.split(';', expand=True)
	taxdf = taxdf.drop('genomad_taxonomy', axis=1)

	taxdf.to_csv(output_f, index=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.description = "Takes the classification summary and the genomad taxonomy column and creates a taxonomy table"

	parser.add_argument("--input", dest="infile", help="input classification summary file", required=True)
	parser.add_argument("--output", dest="outfile", help="ouput csv file path", required=True)
	
	args = parser.parse_args()
	extract_taxonomy(args.infile, args.outfile)



