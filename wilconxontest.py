# Test whether values at one loci (eg. first exon) is significantly different from values at another loci (eg. 4th exon)
import pandas as pd
from scipy.stats import ranksums

# Load the datasets from the files
file1 = 'ORF Net GC chimp exon1.txt'
file2 = 'ORF Net GC chimp exon4.txt'

data1 = pd.read_csv(file1)
data2 = pd.read_csv(file2)

# Perform the Wilcoxon rank-sum test
statistic, p_value = ranksums(data1, data2)

# Print the results
print(f"Wilcoxon rank-sum test statistic: {statistic}")
print(f"P-value: {p_value}")