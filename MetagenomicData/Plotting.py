import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import math

class HistogramPlotter:
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


    def plot_histograms(self, row=None, col=None):
        """
        Generate and save histograms grouped by the specified column to a PDF,
        with axes set dynamically according to data range.
        """
        # Extract unique values from the specified column
        unique_values = self.metadata[self.column_name].unique()
        
        # Calculate optimal grid size
        num_plots = len(unique_values)
        self.row = row  
        self.col = col 
        
        with PdfPages(self.pdf_filename) as pdf:
            fig, axes = plt.subplots(self.row, self.col, figsize=(18, 6))
            axes = axes.flatten()

            for idx, value in enumerate(unique_values):
                # Filter data by the selected column value
                value_indices = self.metadata[self.column_name] == value
                value_data = self.data[value_indices, :]

                ax = axes[idx]
                all_data=[]
                for i in range(value_data.shape[0]):
                    # Remove NaN values
                    masked_data = value_data[i, :][~np.isnan(value_data[i, :])]
                    all_data.extend(masked_data)
                    # Create histogram with specified bins
                    hist, __ = np.histogram(masked_data, bins=self.bins, density=True)
                    bins_mid = (self.bins[1:] + self.bins[:-1]) / 2
                    ax.scatter(bins_mid, hist, facecolors='none', edgecolors='r', alpha=0.7)
                if all_data:
                    ax.set_xlim([min(all_data), max(all_data)])
                    #ax.set_ylim([min(hist) , max(hist) ])

                
                ax.set_xscale('log')
                ax.set_yscale('log')
                ax.set_title(f'{value}')

                # Hide unused subplots
                for j in range(len(value_data), len(axes)):
                    axes[j].axis('off')

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)