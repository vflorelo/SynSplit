#!/usr/bin/python3
codon_dict = {"AAA":"Lys","AAC":"Asn","AAG":"Lys","AAT":"Asn","ACA":"Thr","ACC":"Thr","ACG":"Thr","ACT":"Thr",
              "AGA":"Arg","AGC":"Ser","AGG":"Arg","AGT":"Ser","ATA":"Ile","ATC":"Ile","ATT":"Ile","CAA":"Gln",
              "CAC":"His","CAG":"Gln","CAT":"His","CCA":"Pro","CCC":"Pro","CCG":"Pro","CCT":"Pro","CGA":"Arg",
              "CGC":"Arg","CGG":"Arg","CGT":"Arg","CTA":"Leu","CTC":"Leu","CTG":"Leu","CTT":"Leu","GAA":"Glu",
              "GAC":"Asp","GAG":"Glu","GAT":"Asp","GCA":"Ala","GCC":"Ala","GCG":"Ala","GCT":"Ala","GGA":"Gly",
              "GGC":"Gly","GGG":"Gly","GGT":"Gly","GTA":"Val","GTC":"Val","GTG":"Val","GTT":"Val","TAC":"Tyr",
              "TAT":"Tyr","TCA":"Ser","TCC":"Ser","TCG":"Ser","TCT":"Ser","TGC":"Cys","TGT":"Cys","TTA":"Leu",
              "TTC":"Phe","TTG":"Leu","TTT":"Phe"}
amino_list = list(set(codon_dict.values()))
amino_list.sort()
reverse_codon_dict = {}
for amino_acid in amino_list:
    codon_list = []
    for codon in codon_dict.keys():
        amino = codon_dict[codon]
        if (amino == amino_acid):
            codon_list.append(codon)
    reverse_codon_dict[amino_acid] = codon_list
codon_list = list(set(codon_dict.keys()))
codon_list.sort()