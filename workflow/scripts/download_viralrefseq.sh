#!/bin/bash

if [ "$#" -ne 4 ]; then
	echo "Usage: $0 \ninput=ncbi_file.csv (IDs in first column)\noutput=proteins_file.faa\nmax_batch_size (max=90000)\ncpus (int)"
	exit 1
fi

input_file="$1"
output_file="$2"
max_batch_size="$3"
n_cpus="$4"
tmpdir="tmp"

if [ -e "$output_file" ]; then
	echo "Output file '$output_file' already exists. Cannot overwrite."
	exit 1
fi

total_lines=$(wc -l < "$input_file")
batch_size=$((total_lines / n_cpus))
if ((batch_size > max_batch_size)); then
	batch_size=$max_batch_size
fi

num_batches=$((total_lines / batch_size))
if ((total_lines % batch_size != 0)); then
	  num_batches=$((num_batches + 1))
  fi

echo "Processing $total_lines sequences with $n_cpus CPUs..."
echo "The files will be downloaded in $num_batches batches..."
echo "Number of sequences per run: $batch_size"

download_batch() {
	start_line=$1
	end_line=$2
	batch_no=$3
	input_file=$4
	tmpdir=$5
	
	echo "Processing batch #$3 ..."
	echo -e "Downloading sequences $start_line-$end_line \n"
	
	ids=$(cut -f1 -d "," ${input_file} | sed -n "${start_line},${end_line}p")
	efetch -id ${ids} -db protein -format fasta > ${tmpdir}/${batch_no}.faa


}


export -f download_batch

# make a tmp directory to store all parallel tasks 
if [ ! -d "$tmpdir" ]; then
        mkdir "$tmpdir"
        echo -e "A $tmpdir directory has been created.\n"
else
	{
		rmdir $tmpdir
		mkdir "$tmpdir"
		echo -e "A $tmpdir directory has been created.\n"
	} || {
	echo "Directory $tmpdir already exists and is not empty. Cannot overwite."
        exit 1
}
fi


# Run parallel job
starts=$(seq 1 $batch_size $total_lines)
ends=$(seq $batch_size $batch_size $total_lines && echo $total_lines)
batchno=$(seq $num_batches)
parallel --link --jobs $4 --delay 0 download_batch {1} {2} {3} {4} {5} ::: $starts ::: $ends ::: $batchno ::: "$input_file" ::: "$tmpdir"

# Concatenate all parallel job .faa outputs into one final output
echo "All files done! Merging files into $output_file"
cat ${tmpdir}/*.faa > $output_file
echo "Done!"
