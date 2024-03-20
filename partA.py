from template import *

BC_FILE_PATH = os.path.join(DATA_PATH, 'Bacillus clausii.gb')

class BacillusClausiiGB:
    def __init__(self, gb_path, id_header="locus_tag"):
        dataframe, sequence, gene_names = parse_genbank_to_dataframe(gb_path, id_header)
        self.dataframe = dataframe
        self.sequence = sequence
        self.gene_names = gene_names
        self.cds = None
        # self.other = None
        self.plus_strand_all_genes = None
        self.minus_strand_all_genes = None
        self.plus_strand_cds = None
        self.minus_strand_cds = None


    def print_dictionary(self):
        """
        Prints a summary of element types including genes and other types from the DataFrame.
        """
        elements_summary = {'gene': len(self.gene_names), **self.dataframe['type'].value_counts().to_dict()}

        print("Elements Summary:")
        for element_type, count in elements_summary.items():
            print(f"{element_type}: {count}")
        
        print(elements_summary)

    def find_cds(self):
        """
        Find all the genes that are proteins (Type = 'cds') and save them in the arrtibute.
        """
        self.cds = self.dataframe[self.dataframe['type'] == 'CDS']
        # self.other = self.dataframe[self.dataframe['type'] != 'CDS']
    
    def add_gene_length_to_df(self):
        """
        Adds a 'length' column to the dataframe attribute of the instance, representing
        the length of each gene. The gene length is calculated as the absolute difference
        between the 'start' and 'end' positions of the gene.
        """
        self.dataframe['length'] = (self.dataframe['end'] - self.dataframe['start']).abs()
        
    def plot_gene_length_histogram(self):
        """
        Plots histograms of gene lengths for all genes, coding sequences (CDS),
        and non-coding sequences (other), including statistical summaries.
        """
        # Check if the necessary data is present
        if self.dataframe.empty:
            print("Data missing in one or more attributes ('dataframe', 'cds', 'other').")
            return

        # Use direct pandas operations for efficiency
        total_len = self.dataframe['length']

        x_max = total_len.max()
        y_max = len(self.gene_names)
        x_label = 'length'
        y_label = 'count'
        plot_histograma("All Genes", total_len, x_label, y_label, x_max, y_max)

    def set_all_strands(self):
        if self.cds is None:
            self.find_cds()
            
        self.plus_strand_all_genes = self.dataframe[self.dataframe['strand'] == 1]
        self.minus_strand_all_genes = self.dataframe[self.dataframe['strand'] == -1]
        self.plus_strand_cds = self.cds[self.cds['strand'] == 1]
        self.minus_strand_cds = self.cds[self.cds['strand'] == -1]


    def count_genes_by_strand_and_cds(self):
        """
        Divides genes into two groups based on their strand (plus or minus),
        counts how many genes there are in each group, and then counts how many genes in each
        group encode a protein (type = 'CDS').
        """
        if self.cds is None:
            self.find_cds()
        
        plus_strand_all_genes_count = self.plus_strand_all_genes.shape[0]
        plus_strand_cds_count = self.plus_strand_cds.shape[0]
        
        minus_strand_all_genes_count = self.minus_strand_all_genes.shape[0]
        minus_strand_cds_count = self.minus_strand_cds.shape[0]
        
        print(f"Total genes on the plus strand: {plus_strand_all_genes_count}, of which {plus_strand_cds_count} are encoded as proteins.")
        print(f"Total genes on the minus strand: {minus_strand_all_genes_count}, of which {minus_strand_cds_count} are encoded as proteins.")

    def plot_cds_length_histogram_by_strand(self):
        """
        For genes that encode proteins, draws two histograms of gene lengths,
        one for genes on the plus strand and another for genes on the minus strand.
        """
        if self.cds is None:
            self.find_cds()
                
        plus_lengths_cds = self.plus_strand_cds['length']
        minus_lengths_cds = self.minus_strand_cds['length']
        
        if not plus_lengths_cds.empty:
            plot_histograma("CDS Length Distribution on Plus Strand", plus_lengths_cds, 'Gene Length', 'Frequency', color='blue')
        else:
            print("No CDS data on plus strand to plot.")
        
        if not minus_lengths_cds.empty:
            plot_histograma("CDS Length Distribution on Minus Strand", minus_lengths_cds, 'Gene Length', 'Frequency', color='green')
        else:
            print("No CDS data on minus strand to plot.")

    def report_cds_stats_by_strand(self):
        """
        Generates and prints statistical reports for lengths of genes on the plus and minus strands that are CDS.

        """
        if self.cds is None:
            self.find_cds()
        
        plus_strand_lengths_cds = self.plus_strand_cds['length']
        minus_strand_lengths_cds = self.minus_strand_cds['length']
        
        plus_strand_stats = get_stats_table(plus_strand_lengths_cds, 'Plus Strand CDS Lengths Group')
        minus_strand_stats = get_stats_table(minus_strand_lengths_cds, 'Minus Strand CDS Lengths Group')
        
        combined_stats = pd.concat([plus_strand_stats, minus_strand_stats])
        print(combined_stats)


if __name__ == "__main__":
    genbank = BacillusClausiiGB(BC_FILE_PATH)

    # Question 1
    genbank.print_dictionary()

    # Question 2
    genbank.add_gene_length_to_df()
    genbank.plot_gene_length_histogram()
    genbank.set_all_strands()
    genbank.count_genes_by_strand_and_cds()
    # Need to determine if the distributions are similar in the Word file.
    genbank.plot_cds_length_histogram_by_strand()
    genbank.report_cds_stats_by_strand()
    