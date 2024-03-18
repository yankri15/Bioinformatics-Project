from template import *

BC_FILE_PATH = os.path.join(DATA_PATH, 'Bacillus clausii.gb')

class BacillusClausiiGB:
    def __init__(self, gb_path, id_header="locus_tag"):
        dataframe, sequence, gene_names = parse_genbank_to_dataframe(gb_path, id_header)
        self.dataframe = dataframe
        self.sequence = sequence
        self.gene_names = gene_names
        self.cds = None
        self.other = None

    def print_dictionary(self):
        """
        Prints a summary of element types including genes and other types from the DataFrame.
        """
        elements_summary = {'gene': len(self.gene_names), **self.dataframe['type'].value_counts().to_dict()}

        print("Elements Summary:")
        for element_type, count in elements_summary.items():
            print(f"{element_type}: {count}")



if __name__ == "__main__":
    genbank = BacillusClausiiGB(BC_FILE_PATH)

    genbank.print_dictionary()
