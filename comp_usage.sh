function comp_usage(){
    echo "fasta2comp.sh calculates the fractions of synonymous codon usage and %G+C content of protein-coding sequences"
    echo
    echo "Options:"
    echo "  --fasta_file -> File containing CDS in fasta format (mandatory)"
    echo "  --out_file   -> Output file for the results (mandatory)"
    echo "  --threads    -> Number of CPU threads to use (mandatory)"
    echo
    echo "Example:"
    echo "  fasta2comp.sh --fasta_file /path/to/my/sequences.fasta --out_file results.tsv --threads 16"
	}
export -f comp_usage