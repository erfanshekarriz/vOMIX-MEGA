taxcols = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']

scoredict = {}
scoredict = {key: {} for key in headers}

for key in scoredict:
	for taxlev in taxcols:
		scoredict[taxlev] = set(mergedf.filter(like=taxlev).values.flatten())
