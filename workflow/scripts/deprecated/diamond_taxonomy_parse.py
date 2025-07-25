import argparse 
import pandas as pd 

def parse_diamond(dmndout, ncbitaxcsv, thresh, taxcols, output_csv):
	dmndtsvout = dmndout
	ncbitaxcsv= ncbitaxcsv
	thresh = thresh
	taxcols = taxcols

	metadf = pd.read_csv(ncbitaxcsv, delimiter=",")
	metadf['Accession_nv'] = metadf['Accession'].str.rsplit('.', n=1).apply(lambda x: x[0])
	metadf = metadf.replace('^.*unclassified.*$', 'Unclassified', regex=True)
	
	
	taxdf = metadf[taxcols].drop_duplicates()
	
	dmndheaders = [
			"query_id", 
			"subject_id", 
			"percentage_match", 
			"ali_length", 
			"n_mistmatch", 
			"n_gaps", 
			"start_query", 
			"end_query", 
			"start_subject", 
			"end_subject", 
			"eval", 
			"bit_score" 
			]
	df = pd.read_csv(dmndtsvout, delimiter="\t", comment='#', names=dmndheaders)
	
	# add extra columns
	df['database_id'] = df['subject_id'].str.rsplit('.', n=1).apply(lambda x: x[0])
	df['contig_id'] = df['query_id'].str.rsplit('_', n=1).apply(lambda x: x[0])
	
	mergedf=df.merge(metadf, left_on="database_id", right_on="Accession_nv")
	
	# group by query id and sort by bit-score
	groupdf = mergedf.sort_values(['query_id', 'bit_score'], ascending=False).groupby('query_id')
	
	# 1) FILTER TOP HITS PER PROTEIN
	# per protein diamond hit, pick the highest score hit 
	# if it does not have taxonomy annotation, pick the second if it is within 
	# 25% of the top hit and has taxonomy annotation 
	
	filterlist = []
	for group_n, group_df in groupdf:
		# check if highest bit-score has taxonomy
		firstrow = group_df.iloc[0]
		tax = firstrow.loc[taxcols]
		notax = tax.isnull().all()
		nhits = group_df.shape[0]
		# if there is taxonomy or only one top hit
		if (not notax) or (nhits == 1):
			filterlist.append(firstrow)
		else:
			secondrow = group_df.iloc[1]
			tax = secondrow.loc[taxcols]
			notax = tax.isnull().all()
			firsthitbit = firstrow['bit_score']
			secondhitbit = secondrow['bit_score']
			if (secondhitbit * 1.25 >= firsthitbit) and (not notax):
				filterlist.append(secondrow)
			else:
				filterlist.append(firstrow)
	
	tophitsdf = pd.DataFrame(filterlist)
	
	# 2) AGGREGATE PER CONTIG & BIT-SCORE PER TAX RANK
	# if the sum of the bit score for a species is 70% of the total 
	# accumulated bitscore [since bitscores are independent of query length]
	# the contig will be assigned that value
	
	taxdict = {}
	# Iterate through each contig 
	for contig in tophitsdf['contig_id'].unique():
		taxdict[contig] = {}
		# Iterate through each tax level
		col_len = len(taxcols)
		for i, taxlev in enumerate(taxcols):
			contigdf = tophitsdf[tophitsdf['contig_id'].values == contig]
			groupdf = contigdf.groupby(taxlev)
			bitsumdf = groupdf['bit_score'].sum()
			bitsumcontig = bitsumdf.sum()
			bitpercdf = pd.DataFrame(bitsumdf / bitsumcontig * 100).reset_index()
			bitpercdf  = bitpercdf.rename(columns={'bit_score':'bit_perc'})
			# if the function fails it mean there is not tax annotation at that level
			try:
				maxhit = bitpercdf.iloc[bitpercdf['bit_perc'].idxmax()]
			except ValueError:
				tax = "Na"
				taxdict[contig][taxlev] = tax
				continue
			score = maxhit['bit_perc']
			if score < (thresh * 100):
				tax = "Unclassified"
				taxdict[contig][taxlev] = tax
				continue
			else:
				tax = maxhit[taxlev]
				# retrieve entire lineage
				taxdfslice = taxdf.iloc[:, i:(col_len+1)].drop_duplicates()
				lineage = taxdfslice[taxdfslice[taxlev] == tax].iloc[0].rename_axis(None).to_dict()
				taxdict[contig].update(lineage)
				break
	
	outdf = pd.DataFrame(taxdict).T
	outdf.to_csv(output_csv)



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.description = "Takes in diamond output and table with taxonomic labels for query IDs and give them taxonomic assignments according to https://doi.org/10.1038/s41564-021-00928-6 filtering thresholds"
	
	parser.add_argument(
			"--diamondout", 
			required=True, 
			dest="diamond_out", 
			help="Diamond tsv file output file", 
			type=str)

	parser.add_argument("--taxcsv", 
			required=True, 
			dest="taxcsv", 
			help="CSV file with 'Accession' column that matches diamond query id", 
			type=str)

	parser.add_argument(
			"--threshold", 
			required=True, 
			dest="threshold_v", 
			help="The minimum ratio of bit scores that should match a taxonomic level for it to be classified as that before moving on to higher taxonomic level (default 0.7)",
			default=0.7, 
			type=float)

	parser.add_argument("--taxcolumns", 
			required=True, 
			dest="taxcolstr", 
			help="A ',' separated string indicating the taxonomic column names in --taxcsv in ORDER to be assesed (e.g. --taxcolumns genus,family,order,class)", 
			type=str)

	parser.add_argument("--outputcsv", 
			required=True, 
			dest="out_path", 
			help="The full file path of csv file to store final classifications", 
			type=str)

	args = parser.parse_args()
	taxcolslist = args.taxcolstr.split(",")

	parse_diamond(args.diamond_out, args.taxcsv, args.threshold_v, taxcolslist, args.out_path)
