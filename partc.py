from Bio import Entrez, SeqIO
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
    gb_record = SeqIO.read(april_2021_covid_file_path, "genbank")
except FileNotFoundError:
    Entrez.email = "kobieshka@gmail.com"
    with Entrez.efetch(
        db="nucleotide", id="MZ054892.1", rettype="gb", retmode="text"
    ) as in_handle:
        with open(april_2021_covid_file_path, "w") as out_handle:
            out_handle.write(in_handle.read())
        print(f"Saved file {april_2021_covid_file_path}")
    gb_record = SeqIO.read(april_2021_covid_file_path, "genbank")

try:
    gb_record = SeqIO.read(feb_2024_covid_file_path, "genbank")
except FileNotFoundError:
    Entrez.email = "kobieshka@gmail.com"
    with Entrez.efetch(
        db="nucleotide", id="PP348372.1", rettype="gb", retmode="text"
    ) as in_handle:
        with open(feb_2024_covid_file_path, "w") as out_handle:
            out_handle.write(in_handle.read())
        print(f"Saved file {feb_2024_covid_file_path}")
    gb_record = SeqIO.read(feb_2024_covid_file_path, "genbank")
