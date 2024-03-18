from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from os import path
from pathlib import Path



DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Data')

HYDRO_AMINO_ACIDS = ["A", "F", "I", "L", "M", "P", "V", "W"]
DNA = ['A', 'C', 'T', 'G']

TRANS_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


def read_and_parse_genbank(gb_file_path):
    """
    Parses a GenBank file and returns the first record.

    Parameters:
    - bs_path (str): Path to the GenBank file.

    Returns:
    - Bio.SeqRecord.SeqRecord: The first GenBank record found in the file.

    Raises:
    - FileNotFoundError: If the GenBank file does not exist.
    - ValueError: If no records are found in the file.
    """

    gb_path_obj = Path(gb_file_path)
    if not gb_path_obj.exists():
        raise FileNotFoundError(f"Bacillus Clausii GenBank file not found at {gb_file_path}")
    
    with gb_path_obj.open('r') as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        try:
            gb_record = next(gen)
        except StopIteration:
            raise ValueError(f"No records found in {gb_file_path}")
    
    return gb_record


def parse_genbank_to_dataframe(gb_file_path, id_header='locus_tag'):
    """
    Parses a GenBank file to extract features into a pandas DataFrame.
    
    This function reads a GenBank file, extracts various features and compiles them 
    into a DataFrame for further analysis. 

    Parameters:
    - gb_file_path (str): The file path to the GenBank file.
    - id_header (str, optional): The header under which the gene identifier is stored,
      defaulting to 'locus_tag'.

    Returns:
    - pd.DataFrame: The DataFrame.
    - str: The uppercased genomic sequence.
    - list: A list of gene names extracted from the GenBank file.
    """
    gb_record = read_and_parse_genbank(gb_file_path)
    genom_seq = gb_record.seq.upper()  
    
    data = {
        'id': [], 'start': [], 'end': [], 'strand': [], 'type': [],
        'table': [], 'translation': [], 'codon_start': [],
        'cell_wall': [], 'product': [], 'protein_id': []
    }
    
    gene_names = []

    for feature in gb_record.features[1:]:  # Skip the source feature
        qualifiers = feature.qualifiers
        
        g_name = qualifiers.get(id_header, [None])[0]
        
        if feature.type == 'gene':
            gene_names.append(g_name)
            continue  
        
        data['id'].append(g_name)
        data['type'].append(feature.type)
        data['start'].append(int(feature.location.start))
        data['end'].append(int(feature.location.end))
        data['strand'].append(feature.location.strand)
        data['cell_wall'].append('Yes' if 'cell wall' in qualifiers.get('product', [''])[0] else 'No')
        
        # For CDS features, extract additional details
        if feature.type == 'CDS':
            data['table'].append(qualifiers.get('transl_table', [1])[0])
            data['translation'].append(qualifiers.get('translation', [None])[0])
            data['codon_start'].append(qualifiers.get('codon_start', [None])[0])
            data['product'].append(qualifiers.get('product', [None])[0])
            data['protein_id'].append(qualifiers.get('protein_id', [None])[0])
        else:
            for key in ['table', 'translation', 'codon_start', 'product', 'protein_id']:
                data[key].append(None)

    # Construct DataFrame object from dictionary
    df = pd.DataFrame(data)
    # Add sub_sequence and check columns
    df['sub_sequence'] = df.apply(lambda row: genom_seq[row['start']:row['end']], axis=1)
    
    # Assuming check_translation is a predefined function
    df['check'] = [validate_translation(row['type'], row['translation'], row['sub_sequence'], row['table'], row['strand'], row['codon_start']) for index, row in df.iterrows()]

    return df, genom_seq, gene_names

def validate_translation(gene_type, org_trans, seq, trans_table_id, strand, codon_start_index):
    """
    Validates the translation of a given sequence against an original translation,
    considering the specified codon start position, translation table, and strand.

    Parameters:
    - gene_type (str): Type of gene (only processes "CDS" types).
    - org_trans (str): Original translation string to validate against.
    - seq (str): DNA sequence to be translated.
    - trans_table_id (int): ID of the translation table to use.
    - strand (int): Strand of the DNA sequence (1 for forward, -1 for reverse).
    - codon_start_index (int): Starting index of the codon (1-based).

    Returns:
    - str: "OK" if the translation matches the original, otherwise returns the error message.
    - None: If gene_type is not "CDS".
    """
    if gene_type != "CDS":
        return None

    codon_start_index = int(codon_start_index) - 1
    sequence = Seq(seq[(codon_start_index - 1):])

    if strand == -1:
        sequence = sequence.reverse_complement()

    try:
        # Ensure trans_table_id is an integer
        trans_table_id = int(trans_table_id)
        trans = sequence.translate(table=trans_table_id, to_stop=True)
        if trans == org_trans:
            return "OK"
    except Exception as err:  
        return str(err)