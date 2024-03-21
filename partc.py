from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import Bio.Data.CodonTable
import seaborn as sns

# 1. בחלק ג' בפרויקט הסופי אפשר להשתמש בפונקציה המובנת עבור החישוב של dnds.


# 1. Get the standard_dna_table
standard_table = Bio.Data.CodonTable.standard_dna_table
print(standard_table)

# 2. Get the codons
codons = standard_table.forward_table
print(codons)


# 3. group by amino acid
aa_codons = {}
for codon, aa in codons.items():
    if aa not in aa_codons:
        aa_codons[aa] = []
    aa_codons[aa].append(codon)

print(aa_codons)

synonimous_codons = {}

for key, strings in aa_codons.items():
    for string in strings:
        if key not in synonimous_codons:
            synonimous_codons[key] = []
        synonimous_codons[key].append(string)
