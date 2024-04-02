from Bio import Entrez, SeqIO, codonalign
import pandas as pd
import matplotlib.pyplot as plt
import Bio.Data.CodonTable
import seaborn as sns


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

synonimous_codons_first_char = {}
synonimous_codons_second_char = {}
synonimous_codons_third_char = {}

for key, strings in aa_codons.items():
    # remove first char from each string
    strings_minus_first_char = [x[1:] for x in strings]
    strings_minus_second_char = [x[0] + x[2] for x in strings]
    strings_minus_third_char = [x[:2] for x in strings]

    # check if cell have duplicates in the list
    duplicates_first_char = [
        x for x in strings_minus_first_char if strings_minus_first_char.count(x) > 1
    ]
    synonimous_codons_first_char[key] = duplicates_first_char

    duplicates_second_char = [
        x for x in strings_minus_second_char if strings_minus_second_char.count(x) > 1
    ]
    synonimous_codons_second_char[key] = duplicates_second_char

    duplicates_third_char = [
        x for x in strings_minus_third_char if strings_minus_third_char.count(x) > 1
    ]
    synonimous_codons_third_char[key] = duplicates_third_char


print("==========")
for key, value in synonimous_codons_first_char.items():
    print(key, value)
print("=====2nd=====")
for key, value in synonimous_codons_second_char.items():
    print(key, value)

print("=====3rd=====")

for key, value in synonimous_codons_third_char.items():
    print(key, value)

count_first_char = 0
count_second_char = 0
count_third_char = 0

for key, value in synonimous_codons_first_char.items():
    count_first_char += sum([value.count(x) - 1 for x in value])

for key, value in synonimous_codons_second_char.items():
    count_second_char += sum([value.count(x) - 1 for x in value])

for key, value in synonimous_codons_third_char.items():
    count_third_char += sum([value.count(x) - 1 for x in value])

print(count_first_char)
print(count_second_char)
print(count_third_char)
total = count_first_char + count_second_char + count_third_char
print(total)


# שאלה 2
# if file is not found, download it from genbank
# if file is found, read it from the file
april_2021_covid_file_path = "./Data/MZ054892.gb"
feb_2024_covid_file_path = "./Data/PP348372.gb"
try:
    april_2021_covid_gb_record = SeqIO.read(april_2021_covid_file_path, "genbank")
    print(f"File exists: {april_2021_covid_file_path}")

except FileNotFoundError:
    Entrez.email = "kobieshka@gmail.com"
    with Entrez.efetch(
        db="nucleotide", id="MZ054892.1", rettype="gb", retmode="text"
    ) as in_handle:
        with open(april_2021_covid_file_path, "w") as out_handle:
            out_handle.write(in_handle.read())
        print(f"Saved file {april_2021_covid_file_path}")
    april_2021_covid_gb_record = SeqIO.read(april_2021_covid_file_path, "genbank")

try:
    feb_2024_covid_gb_record = SeqIO.read(feb_2024_covid_file_path, "genbank")
    print(f"File exists: {feb_2024_covid_file_path}")
except FileNotFoundError:
    Entrez.email = "kobieshka@gmail.com"
    with Entrez.efetch(
        db="nucleotide", id="PP348372.1", rettype="gb", retmode="text"
    ) as in_handle:
        with open(feb_2024_covid_file_path, "w") as out_handle:
            out_handle.write(in_handle.read())
        print(f"Saved file {feb_2024_covid_file_path}")
    feb_2024_covid_gb_record = SeqIO.read(feb_2024_covid_file_path, "genbank")

# a. genome length
april_2021_covid_genome_length = len(april_2021_covid_gb_record.seq)
print(f"April 2021 covid genome length: {april_2021_covid_genome_length}")

feb_2024_covid_genome_length = len(feb_2024_covid_gb_record.seq)
print(f"Feb 2024 covid genome length: {feb_2024_covid_genome_length}")


# b. number of genes
april_2021_genes = set()

for feature in april_2021_covid_gb_record.features:
    if feature.type == "gene":
        april_2021_genes.add(feature.qualifiers["gene"][0])

print(f"April 2021 covid number of genes: {len(april_2021_genes)}")

feb_2024_genes = set()
for feature in feb_2024_covid_gb_record.features:
    if feature.type == "gene":
        feb_2024_genes.add(feature.qualifiers["gene"][0])

print(f"Feb 2024 covid number of genes: {len(feb_2024_genes)}")

# c. number of proteins
april_2021_num_of_proteins = 0
for feature in april_2021_covid_gb_record.features:
    if feature.type == "CDS":
        april_2021_num_of_proteins += 1
print(f"April 2021 covid number of proteins: {april_2021_num_of_proteins}")

feb_2024_num_of_proteins = 0
for feature in feb_2024_covid_gb_record.features:
    if feature.type == "CDS":
        feb_2024_num_of_proteins += 1
print(f"Feb 2024 covid number of proteins: {feb_2024_num_of_proteins}")

# d. shared genes
print(april_2021_genes)
print(feb_2024_genes)

shared_genes = april_2021_genes.intersection(feb_2024_genes)
print(f"Shared genes: {shared_genes}")


def some_other_intersting():
    # get translations
    april_2021_translations = []
    for feature in april_2021_covid_gb_record.features:
        if feature.type == "CDS":
            april_2021_translations.append(feature.qualifiers["translation"][0])

    feb_2024_translations = []
    for feature in feb_2024_covid_gb_record.features:
        if feature.type == "CDS":
            feb_2024_translations.append(feature.qualifiers["translation"][0])

    # april_2021_translations md5
    april_2021_translations_md5 = []
    for translation in april_2021_translations:
        april_2021_translations_md5.append(hash(translation))

    feb_2024_translations_md5 = []
    for translation in feb_2024_translations:
        feb_2024_translations_md5.append(hash(translation))

    # print(april_2021_translations_md5)
    # print(feb_2024_translations_md5)


# dnds
