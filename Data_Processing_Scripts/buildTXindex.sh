input="/scratch/alice/h/hp329/Steered_Project/hg38/ncbi_dataset/data/GCF_000001405.40/genomic.gff"
# Script to extract transcript to gene mappings from a GFF file

grep "ID=rna-" "$input" | \ # Pull every row labelled as RNA
	awk '{
		match($0, /ID=([^;]+)/, a); #Look for matches followed by non ; character
		match($0, /Parent=([^;]+)/, b);
		gsub(/^rna-/, "", a[1]); # Remove rna prefix from transcripts
		if (a[1] != "" && b[1] != "") print a[1] "\t" b[1]; # If both IDs are valid, print into tab separated rows
	}' > tx2gene.tsv


