import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

hydrophopic_amino_acids = ["A", "F", "I", "L", "M", "P", "V", "W"]

uniprot_data = pd.read_csv("./Data/uniprotkb_Bacillus_clausii_2024_03_16.tsv", sep="\t")

genbank_data = [rec for rec in SeqIO.parse("./Data/Bacillus clausii.gb", "genbank")]


# get locus_tags from genbank file
genbank_locus_tags = []
for record in genbank_data:
    for feature in record.features:
        if feature.type == "CDS":
            genbank_locus_tags.append(feature.qualifiers["locus_tag"][0])

genbank_locus_tags = list(set(genbank_locus_tags))

# get locus_tags from uniprot file
uniprot_locus_tags = uniprot_data["Gene Names (ordered locus)"].tolist()
uniprot_locus_tags = list(set(uniprot_locus_tags))

common_locus_tags = list(set(genbank_locus_tags) & set(uniprot_locus_tags))
unique_locus_tags_genbank = list(set(genbank_locus_tags) - set(uniprot_locus_tags))
unique_locus_tags_uniprot = list(set(uniprot_locus_tags) - set(genbank_locus_tags))

print(f"Number of common locus tags: {len(common_locus_tags)}")
print(f"Total number of locus tags in genbank file: {len(genbank_locus_tags)}")
print(f"Total number of locus tags in uniprot file: {len(uniprot_locus_tags)}")
print(f"Number of unique locus tags in genbank file: {len(unique_locus_tags_genbank)}")
print(f"Number of unique locus tags in uniprot file: {len(unique_locus_tags_uniprot)}")


#  difference between the number of locus tags in genbank and uniprot file in percentage
fig, ax = plt.subplots()
ax.bar("Uniprot", len(unique_locus_tags_uniprot), label="Unique in Uniprot")
ax.bar("Genbank", len(unique_locus_tags_genbank), label="Unique in Genbank")
ax.set_ylabel("Number of locus tags")
ax.set_title("Number of unique locus tags in Genbank and Uniprot files")
ax.legend()
plt.show()

similarity_genbank = (len(common_locus_tags) / len(genbank_locus_tags)) * 100
print(
    f"Percentage similarity between the number of locus tags in genbank and uniprot file: {similarity_genbank}%"
)

similarity_uniprot = (len(common_locus_tags) / len(uniprot_locus_tags)) * 100
print(
    f"Percentage similarity between the number of locus tags in genbank and uniprot file: {similarity_uniprot}%"
)


# matplotlib plot
fig, ax = plt.subplots()
ax.bar("similarity_genbank", similarity_genbank, label="Similarity")
ax.bar("similarity_uniprot", similarity_uniprot, label="Similarity")
ax.set_ylabel("Percentage")

ax.text(
    "similarity_genbank",
    similarity_genbank,
    f"{similarity_genbank}%",
    ha="center",
    va="bottom",
)
ax.text(
    "similarity_uniprot",
    similarity_uniprot,
    f"{similarity_uniprot}%",
    ha="center",
    va="bottom",
)

ax.set_title("Similarity between the number of locus tags (%)")
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
# set dimmensions of the plot
plt.subplots_adjust(right=0.7)

plt.show()


transmembrane_sequences = uniprot_data[
    uniprot_data["Transmembrane"].notna() & uniprot_data["Sequence"].notna()
]
# print len
print(f"Number of transmembrane sequences: {len(transmembrane_sequences)}")
print(transmembrane_sequences["Transmembrane"].head())

transmembrane_sequences_list = []
# loop through the transmembrane sequences dataframe
# add data in the format of {"entry": row["Entry"], "sequences": [(1, 1), (2, 2)]}
for index, row in transmembrane_sequences.iterrows():
    transmembrane_sequences_list.append(
        {
            "entry": row["Entry"],
            "sequences": [
                (int(i), int(j))
                for i, j in row.str.extractall(r"(\d+)\.\.(\d+)").values
            ],
        }
    )

lengths = []
count = 0
for sequence in transmembrane_sequences_list:
    for start, end in sequence["sequences"]:
        lengths.append(end - start)
        count += 1

print(f"Number of transmembrane sequences: {count}")
print(f"Average length: {sum(lengths) / len(lengths)}")
print(f"Minimum length: {min(lengths)}")
print(f"Maximum length: {max(lengths)}")

# matplotlib histogram, Length and their occurences
fig, ax = plt.subplots()
ax.hist(lengths, edgecolor="black", bins=20, alpha=0.7, rwidth=0.85, align="left")
ax.set_xlabel("Length")
ax.set_ylabel("Occurences")
ax.set_title(
    f"Length and their occurences from \ntotal of {count} transmembrane sequences"
)

# show count of each length
for i in range(min(lengths), max(lengths) + 1):
    ax.text(
        i,
        lengths.count(i),
        lengths.count(i),
        ha="center",
        va="bottom",
        fontsize=8,
        color="black",
    )

plt.show()


hydrophobic_amino_acids_count = dict.fromkeys(hydrophopic_amino_acids, 0)

total_amino_acids = 0
for sequence in transmembrane_sequences_list:
    seq = uniprot_data[uniprot_data["Entry"] == sequence["entry"]]["Sequence"].values[0]

    for start, end in sequence["sequences"]:
        segment = seq[start - 1 : end]
        for aa in hydrophopic_amino_acids:
            hydrophobic_amino_acids_count[aa] += segment.count(aa)
        total_amino_acids += len(segment)

print(hydrophobic_amino_acids_count)
print(f"Total amino acids: {total_amino_acids}")

# sum of hydrophobic_amino_acids_count from total
total = sum(hydrophobic_amino_acids_count.values())
print(f"Total hydrophobic amino acids: {total}")

print(
    f"Percentage of hydrophobic amino acids in transmembrane sequences: {(total / total_amino_acids) * 100}%"
)
