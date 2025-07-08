import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

def load_degree_distribution(filename):
    df = pd.read_csv(filename, sep='\t')
    return df['degree'].values, df['count'].values

def plot_degree_distribution(degrees, counts, output_pdf):
    with PdfPages(output_pdf) as pdf:
        # Figure 1: Regular degree distribution
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.bar(degrees, counts)
        ax.set_xlim(0, 800)
        ax.set_xlabel('Degree (k)')
        ax.set_ylabel('Number of nodes')
        ax.set_title('Degree Distribution (linear)')
        ax.grid(True, alpha=0.3)
        pdf.savefig(fig)
        plt.close()
        
        # Figure 2: Log-log plot for scale-free analysis
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Filter out zeros for log plot
        log_degrees = np.log2(degrees, out=np.zeros_like(degrees, np.float64), where=(degrees != 0))
        log_counts = np.log2(counts, out=np.zeros_like(counts, np.float64), where=(counts != 0))
        
        ax.scatter(log_degrees, log_counts, alpha=0.6, s=30)
        
        # Fit a line to estimate gamma (power law exponent)
        # For scale-free: P(k) ~ k^(-gamma)
        # log(P(k)) ~ -gamma * log(k)
        slope, intercept, rvalue, p, se = linregress(log_degrees, log_counts)
        
        fit_line = slope * log_degrees + intercept
        ax.plot(log_degrees, fit_line, 'r-', linewidth=2, 
                label=f'Fit: γ = {-slope:.2f}, R² = {rvalue**2:.3f}') # type: ignore
        
        ax.set_xlabel('log₂(Degree)')
        ax.set_ylabel('log₂(Count)')
        ax.set_title('Degree Distribution (log-log)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        pdf.savefig(fig)
        plt.close()

        print(f"Wrote plots to '{output_pdf}'")


def main():
    input_data_path = "degree_distribution.txt"
    plot_output_path = "degree_distribution.pdf"

    degrees, counts = load_degree_distribution(input_data_path)
    
    plot_degree_distribution(degrees, counts, plot_output_path)

if __name__ == "__main__":
    main()