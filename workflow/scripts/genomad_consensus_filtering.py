import os 
import argparse
import pandas as pd 


def consensus_filtering(classification_summary, genomad_min, summary_out, provirus_out, virus_out):

	df = pd.read_csv(classification_summary, index_col=0)

	# add a new column called hit-score which shows how many softwares it has a hit from 
	# according to the config file thresholds set

	df['hit-score'] = 0
	df.loc[(df['genomad_virus_score'] > genomad_min), 'hit-score'] += 1
	
	df['provirus-hit-score'] = 0
	df.loc[(df['checkv_provirus'] == 'Yes'), 'provirus-hit-score'] += 1
	df.loc[(df['genomad_topology'] == 'Provirus'), 'provirus-hit-score'] += 1
	
	# set sequences identified as proviruses by both checkV and genoad as true pro-viruses
	# if only one tool says it's a provirus, it will label as lower confidence
	df['type'] = 'Virus'
	df.loc[(df['provirus-hit-score'] == 2), 'type'] = 'Provirus (geNomad+CheckV)'
	df.loc[(df['provirus-hit-score'] == 1) & (df['genomad_topology'] == 'Provirus'), 'type'] = 'Provirus (geNomad)'
	df.loc[(df['provirus-hit-score'] == 1) & (df['checkv_provirus'] == 'Yes'), 'type']  = 'Provirus (CheckV)'

	# filter the dataframe if not 'type' is 'Virus' but contamination is higher than 10%
	df = df[(df['type'] != 'Virus') | (df['checkv_contamination'] <= 10)]

	exclusion_warnings = [
		#"low-confidence DTR",
		#"low-confidence ITR",
		#"no viral genes detected",
		"no viral genes detected; low-confidence DTR",
		"no viral genes detected; low-confidence ITR",
		">1 viral region detected; contig >1.5x longer than expected genome length",
		">1 viral region detected",
		"high kmer_freq may indicate large duplication; contig >1.5x longer than expected genome length; low-confidence DTR",
		"high kmer_freq may indicate large duplication; contig >1.5x longer than expected genome length; low-confidence ITR",
		"no viral genes detected; contig >1.5x longer than expected genome length",
		"high kmer_freq may indicate large duplication; contig >1.5x longer than expected genome length",
		"no viral genes detected; high kmer_freq may indicate large duplication; low-confidence ITR",
		"high kmer_freq may indicate large duplication; low-confidence ITR",
		"contig >1.5x longer than expected genome length",
		">1 viral region detected; low-confidence ITR",
		">1 viral region detected; low-confidence DTR",
		"contig >1.5x longer than expected genome length; low-confidence Provirus",
		"no viral genes detected; high kmer_freq may indicate large duplication",
		"low-confidence Provirus"]

	df_exclusions = df[df["checkv_warnings"].isin(exclusion_warnings)]
	df_filt = df[~df["checkv_warnings"].isin(exclusion_warnings)]

	df_highqual = df_filt[~df_filt["checkv_quality"].isin(["Not-determined", "Low-quality"])]
	df_lowqual = df_filt[df_filt["checkv_quality"].isin(["Not-determined", "Low-quality"])]
	
	df_lowqual_pass = df_lowqual[
			((df_lowqual["genomad_n_hallmarks"] == 0) & (df_lowqual["genomad_virus_score"] >= 0.99))
			| ((df_lowqual["genomad_n_hallmarks"] >= 1) & (df_lowqual["genomad_n_hallmarks"] <= 5) & (df_lowqual["genomad_virus_score"] >= 0.95))
			| ((df_lowqual["genomad_n_hallmarks"] > 5) & (df_lowqual["genomad_virus_score"] >= 0.90))
			| ((df_lowqual["genomad_marker_enrichment"] >= df["genomad_marker_enrichment"].quantile(0.90)))
			]

	df_pass = pd.concat([df_highqual, df_lowqual_pass])

	# drop any sequence where more than one proviruse region has been identified by geNomad
	# but checkV does not see a provirus
	#filtered_df = filtered_df[~(filtered_df.index.duplicated()) & (filtered_df['provirus-hit-score'] != 2)]
	
	# drop duplicate indices 
	# these regions represent the ones where geNomad found two proviruses 
	# we can find a cleaner way to handle this later
	df_pass = df_pass[~df_pass.index.duplicated(keep='first')]

	# make all indecies unique 
	provirusdf = df_pass[(df_pass['type'] != 'Virus')]
	virusdf = df_pass[df_pass['type'] == 'Virus']


	# save outputs
	df_pass.to_csv(summary_out, index=True, sep=",")
	provirusdf.index.to_frame().to_csv(provirus_out, header=False, index=False, sep='\t')
	virusdf.index.to_frame().to_csv(virus_out, header=False, index=False, sep='\t')


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Perform consensus viral classification on the outputs to find vOTUs')

	parser.add_argument('--classification_results', type=str, help='Path to final combined classification summary')
	parser.add_argument('--genomad_min_score', type=float, help='Minimum threshold for geNomad to be considered as a hit')
	parser.add_argument('--summary_out', type=str, help='Output file path to final filtered summary of vOTUs')
	parser.add_argument('--provirus_list', type=str, help='Output file path to be list of proviral sequences')
	parser.add_argument('--virus_list', type=str, help='Output file path to be list of viral (non-proviral) sequences')

	args = parser.parse_args()

	consensus_filtering(args.classification_results, args.genomad_min_score, 
			args.summary_out, args.provirus_list, args.virus_list)

