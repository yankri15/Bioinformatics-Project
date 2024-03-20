from template import *

BC_FILE_PATH = os.path.join(DATA_PATH, 'Bacillus clausii.gb')

class BacillusClausiiGB:
    def __init__(self, gb_path, id_header="locus_tag"):
        dataframe, sequence, gene_names = parse_genbank_to_dataframe(gb_path, id_header)
        self.dataframe = dataframe
        self.sequence = sequence
        self.gene_names = gene_names
        self.cds = None
        self.no_cds = None
        self.plus_strand_all_genes = None
        self.minus_strand_all_genes = None
        self.plus_strand_cds = None
        self.minus_strand_cds = None
        self.plus_strand_no_cds = None
        self.minus_strand_no_cds = None

    # General functions for all questions
    def export_dataframe_to_csv(self):
        """
        Exports the dataframe attribute to a CSV file named 'part_A.csv'.
        """
        csv_file_path = os.path.join(DATA_PATH, 'part_A.csv')
        self.dataframe.to_csv(csv_file_path, index=False)
        print(f'Dataframe exported to CSV at {csv_file_path}')

    def find_cds(self):
        """
        Find all the genes that are proteins (Type = 'cds') and save them in the arrtibute.
        """
        self.cds = self.dataframe[self.dataframe['type'] == 'CDS']

    def find_no_cds(self):
        """
        Find all the genes that are not proteins (Type != 'cds') and save them in the arrtibute.
        """
        self.no_cds = self.dataframe[self.dataframe['type'] != 'CDS']

    def set_all_strands(self):
        if self.cds is None:
            self.find_cds()

        if self.no_cds is None:
            self.find_no_cds()
            
        self.plus_strand_all_genes = self.dataframe[self.dataframe['strand'] == 1]
        self.minus_strand_all_genes = self.dataframe[self.dataframe['strand'] == -1]
        self.plus_strand_cds = self.cds[self.cds['strand'] == 1]
        self.minus_strand_cds = self.cds[self.cds['strand'] == -1]
        self.plus_strand_no_cds = self.no_cds[self.no_cds['strand'] == 1]
        self.minus_strand_no_cds = self.no_cds[self.no_cds['strand'] == -1]

    # Question 1    
    def print_dictionary(self):
        """
        Prints a summary of element types including genes and other types from the DataFrame.
        """
        elements_summary = {'gene': len(self.gene_names), **self.dataframe['type'].value_counts().to_dict()}

        print("Elements Summary:")
        for element_type, count in elements_summary.items():
            print(f"{element_type}: {count}")
        
        print(elements_summary)

    # Question 2 section 1
    def add_gene_length_to_df(self):
        """
        Adds a 'length' column to the dataframe attribute of the instance, representing
        the length of each gene. The gene length is calculated as the absolute difference
        between the 'start' and 'end' positions of the gene.
        """
        self.dataframe['length'] = (self.dataframe['end'] - self.dataframe['start']).abs()

    def plot_gene_length_histogram(self):
        """
        Plots histograms of gene lengths for all genes.
        """
        if self.dataframe.empty:
            print("Data missing in one or more attributes ('dataframe', 'cds', 'other').")
            return

        total_len = self.dataframe['length']

        x_max = total_len.max()
        y_max = len(self.gene_names)
        x_label = 'length'
        y_label = 'count'
        plot_single_histograma("All Genes", total_len, x_label, y_label, x_max, y_max)

    # Question 2 section 2
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
        
    # Question 2 section 3,5
    def plot_cds_and_no_cds_length_histogram_by_strand(self):
        """
        For genes on both strands, draws histograms of gene lengths,
        one for genes that encode proteins (CDS) and another for genes that do not encode proteins (non-CDS).
        Plots are displayed in subplots with specified y-axis limits.
        """
        if self.cds is None:
            self.find_cds()
        if self.no_cds is None:
            self.find_no_cds()
            
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # Adjust the figure size here
        plot_multiple_histograma("CDS Length Distribution on Plus Strand", self.plus_strand_cds['length'],
                        'Gene Length', 'Frequency', color='blue', ax=axs[0, 0])
        axs[0, 0].set_ylim([0, 1300])  
        plot_multiple_histograma("CDS Length Distribution on Minus Strand", self.minus_strand_cds['length'],
                        'Gene Length', 'Frequency', color='green', ax=axs[0, 1])
        axs[0, 1].set_ylim([0, 1300]) 
        plot_multiple_histograma("None CDS Length Distribution on Plus Strand", self.plus_strand_no_cds['length'],
                        'Gene Length', 'Frequency', color='blue', ax=axs[1, 0])
        axs[1, 0].set_ylim([0, 80])
        plot_multiple_histograma("None CDS Length Distribution on Minus Strand", self.minus_strand_no_cds['length'],
                        'Gene Length', 'Frequency', color='green', ax=axs[1, 1])
        axs[1, 1].set_ylim([0, 80])  
        plt.tight_layout(pad=3.0)
        plt.show()

        # Save the figure as JPEG
        jpeg_file_path = os.path.join(DATA_PATH, 'gene_length_distributions.jpeg')
        fig.savefig(jpeg_file_path, format='jpeg')
        print(f"Histograms saved as JPEG at {jpeg_file_path}")

    # Question 2 section 4
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

        csv_file_path = os.path.join(DATA_PATH, 'cds_stats_report.csv')    
        combined_stats.to_csv(csv_file_path, index=True)
        print(f"Statistical reports written to: {csv_file_path}")

    def report_no_cds_stats_by_strand(self):
        """
        Generates and prints statistical reports for lengths of genes on the plus and minus strands that are none CDS.
        """
        if self.no_cds is None:
            self.find_no_cds()
        
        plus_strand_lengths_no_cds = self.plus_strand_no_cds['length']
        minus_strand_lengths_no_cds = self.minus_strand_no_cds['length']
        
        plus_strand_stats = get_stats_table(plus_strand_lengths_no_cds, 'Plus Strand none CDS Lengths Group')
        minus_strand_stats = get_stats_table(minus_strand_lengths_no_cds, 'Minus Strand none CDS Lengths Group')
        
        combined_stats = pd.concat([plus_strand_stats, minus_strand_stats])
        print(combined_stats)

        csv_file_path = os.path.join(DATA_PATH, 'none_cds_stats_report.csv')    
        combined_stats.to_csv(csv_file_path, index=True)
        print(f"Statistical reports written to: {csv_file_path}")

if __name__ == "__main__":
    genbank = BacillusClausiiGB(BC_FILE_PATH)

    # Question 1
    genbank.print_dictionary()

    # Question 2
    genbank.add_gene_length_to_df()
    genbank.plot_gene_length_histogram()
    genbank.set_all_strands()
    genbank.count_genes_by_strand_and_cds()
    # An explanation about the distributions is located in Q2S3.txt under Data folder for now.
    genbank.report_cds_stats_by_strand()
    genbank.report_no_cds_stats_by_strand()
    genbank.plot_cds_and_no_cds_length_histogram_by_strand()

    # At the end export the DataFrame holding all the information into part_A.csv file inder Data folder.
    genbank.export_dataframe_to_csv()