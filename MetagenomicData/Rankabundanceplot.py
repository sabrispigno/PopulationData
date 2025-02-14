import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import math
class RankAbundancePlotter:
    def __init__(self, data, metadata, column_name, bins=None, pdf_filename='histograms.pdf'):
        """
        Initialize the HistogramPlotter with data and metadata.

        :param data: 2D numpy array, each row corresponds to data for a histogram
        :param metadata: Metadata DataFrame containing additional information
        :param column_name: Column in metadata to group histograms by
        :param bins: Bin edges for the histograms, default is logspace(-5, 1, 10)
        :param pdf_filename: Name of the PDF file to save the plots
        """
        self.data = data
        self.metadata = metadata
        self.column_name = column_name
        self.bins = bins if bins is not None else np.logspace(-5, 1, 10)
        self.pdf_filename = pdf_filename

    def plot_histograms(self):
        """
        Generate and save histograms grouped by the specified column to a PDF.
        """
        # Extract unique values from the specified column
        unique_values = self.metadata[self.column_name].unique()
        
        # Calculate optimal grid size
        num_plots = len(unique_values)
        self.rows = math.ceil(math.sqrt(num_plots))
        self.cols = math.ceil(num_plots / self.rows)
        
        with PdfPages(self.pdf_filename) as pdf:
            fig, axes = plt.subplots(self.rows, self.cols, figsize=(18, 6), squeeze=False)
            axes = axes.flatten()

            for idx, value in enumerate(unique_values):
                # Filter data by the selected column value
                value_indices = self.metadata[self.column_name] == value
                value_data = self.data[value_indices, :]

                ax = axes[idx]

                for i in range(value_data.shape[0]):
                    # Remove NaN values
                    masked_data = value_data[i, :][~np.isnan(value_data[i, :])]
                    vec = np.sort(masked_data)[::-1]  # Sort in descending order
                    ranks = np.arange(1, len(vec) + 1)  # Calculate ranks
                    # Scatter plot for the current sample
                    ax.scatter(ranks, vec, edgecolors="blue", facecolors='none')      
            

                # Set the x and y axes to logarithmic scale
                ax.set_xscale('log')
                ax.set_yscale('log')


                ax.set_title(f'{value}')

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)