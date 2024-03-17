from Bio import SeqIO

genbank_data = [rec for rec in SeqIO.parse("Bacillus clausii.gb", "genbank")]

# count feature.type
# 1 get all feature types
feature_type = [feature.type for record in genbank_data for feature in record.features]
feature_type_count = {i: feature_type.count(i) for i in feature_type}
print(feature_type_count)
# {'source': 1, 'gene': 4203, 'CDS': 4108, 'rRNA': 22, 'tRNA': 73, 'misc_feature': 27}


# 2.a genes length
gene_length = [
    len(feature)
    for record in genbank_data
    for feature in record.features
    if feature.type == "gene"
]
print(f"Total number of genes: {len(gene_length)}")
