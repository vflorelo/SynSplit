# SynSplit

## Overview

Working with high-throughput data from non-model organisms is challenging in many levels, even more so when there might be mixed populations of genes in the analysed data such as contaminants, genes of different evolutionary origins, genes from symbionts, etc.

Differences in %G+C content has been shown to be indicative of potential contamination or presence of mixed populations in a given dataset. %G+C content directly influences the coding capacity of an organism, which ultimately might reflect availability of tRNAs, translation efficiency, or regulatory mechanisms.

We designed SynSplit as a way to separate populations using machine learning applied to %G+C and synonymous codon usage. SynSplit workflow is straightforward:

- %G+C and codon usage is calculated for a given dataset (65 metrics).
- The user provides two gene lists as the training sets, genes of known origin.
- Two distributions are built for each metric, one for genes in set A and one for genes in set B
- A 2-sample Kolmogorov-Smirnov test is performed to select metrics that are significantly different between the two sets
- A Random forest classifier is built using the selected metrics.
- An internal evaluation of the model is performed, calculating precision, recall, accuracy, true/false positives and true/false negatives
- The rest of genes is classified as part of set A or B.

SynSplit produces tsv and svg files to visualise the classified genes

The program was successfully applied to separate genes from the hybrid genome of [*Perkinsus chesapeaki*](https://www.ebi.ac.uk/ena/browser/view/GCA_978019755.1).

## Dependencies

- [EMBOSS](https://emboss.sourceforge.net/)
- [seqtk](https://github.com/lh3/seqtk)
- [gnu-parallel](https://www.gnu.org/software/parallel/)
- [matplotlib](https://pypi.org/project/matplotlib/)
- [numpy](https://pypi.org/project/numpy/)
- [pandas](https://pypi.org/project/pandas/)
- [seaborn](https://pypi.org/project/seaborn/)
- [scipy](https://pypi.org/project/scipy/)
- [scikit-learn](https://pypi.org/project/scikit-learn/)

## Usage

The first step is to calculate the fractions of synonymous codons as well as the %G+C from the coding sequences. The coding sequences can be transcriptome- or genome-derived CDS. To remove unnecessary bias it is a good idea to remove alternative isoforms of the same gene, if the CDS come from genome annotations just look for the main isoforms in the gff files; if the CDS come from Trinity transcriptomes, look for the identifiers that end with "_i1" [[source](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Output-of-Trinity-Assembly#explaining-the-identifiers-genes-vs-transcripts)]

Once the CDS fasta file is ready, simply run fasta2comp.sh to get the composition metrics from each CDS.

```
fasta2comp.sh --fasta_file cds.fasta --out_file cds.comp.tsv --threads 16
```

Since the splitting process uses a supervised machine learning approach (Random Forest classifier), there should be *a priori* information on at least 50-100 genes per group. This information might be hard to obtain, BLAST searches and/or phylogenetics can be useful to assign labels to at least some genes in the dataset. For SynSplit to work two files are needed, each containing a list of genes of known origin.

Once the files are ready, run SynSplit as follows.
```
SynSplit.py --comp_file cds.comp.tsv --train_list_a group_a.txt --train_list_b group_b.txt --test_list unknown.txt --prefix cds_classified
```

## Output

- prefix.ks.json -> A JSON file with the results of Kolmogorov-Smirnov test for all the evaluated metrics when comparing genes coming from group A and group B in the training step
- prefix.model.json -> A JSON file with the internal quality control of the Random Forest classifier including the accuracy, precision and recall of the RF-classifier.
- prefix.model.eval.json -> A JSON file with the fractions of true/false positives and true/false negatives
- prefix.svg -> Violin-plots in svg format showing the distribution of the different metrics *after* the unknown genes were classified.
- prefix.out.tsv -> A tsv file with the metrics of the complete dataset (training + testing genes), the corresponding metrics and the classifications.

## Example

See the contents of the [example](example) directory.