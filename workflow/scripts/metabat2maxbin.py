import pandas as pd
import argparse

def main():
	parser = argparse.ArgumentParser(description="Process a text file with pandas")
	parser.add_argument("-inputtxt", type=str, help="Path to MetaBAT2 jgi_summarize_bam_contig_depths txt file")
	parser.add_argument("-outputtxt", type=str, help="Path to MaxBIN2 coverage output txt file")
	parser.add_argument("-ending", type=str, default="sorted", help="The ending of the column names to extract, usually .bam or .sorted")
	args = parser.parse_args()
	
	try:
		df = pd.read_csv(args.inputtxt, sep='\t', index_col=0)
		
		# Print any column ending with the specified ending
		dffilt = df.loc[:, df.columns.str.endswith(".sorted")]
		dffilt.to_csv(args.outputtxt, sep='\t', index=True)

	except FileNotFoundError:
		print(f"Error: File not found at {args.file_path}")
	except pd.errors.ParserError:
		print(f"Error: Could not parse the file at {args.file_path}")

if __name__ == "__main__":
	main()
