'''
Use the below codes to find the distribution of GC content of the middle 100
of a particular set of sequences of interest
Get a set of intergenic sequences that have the same distribution of GC content
'''
import matplotlib.pyplot as plt
import numpy as np


## For counting 100 in the middle ##
def calculate_gc_percentage(sequence):
    gc_count = 0
    total_count = 0
    start_position = 450
    end_position = 550

    for index, base in enumerate(sequence):
        if start_position <= index + 1 <= end_position:
            if base != "_":
                total_count += 1
                if base.upper() == "G" or base.upper() == "C":
                    gc_count += 1

    if total_count > 0:
        gc_percentage = (gc_count / total_count) * 100
    else:
        gc_percentage = 0

    return gc_percentage

'''
## For counting all ##
def calculate_gc_percentage(sequence):
    gc_count = 0
    total_count = 0

    for base in sequence:
        if base != "_":
            total_count += 1
            if base.upper() == "G" or base.upper() == "C":
                gc_count += 1

    if total_count > 0:
        gc_percentage = (gc_count / total_count) * 100
    else:
        gc_percentage = 0

    return gc_percentage
'''

def read_fasta_file(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    sequences = []
    current_sequence = ""
    current_title = ""

    for line in lines:
        if line.startswith(">"):
            if current_sequence:
                if len(current_sequence) >= 100:
                    sequences.append((current_title, current_sequence))
                current_sequence = ""
            current_title = line.strip()[1:]
        else:
            current_sequence += line.strip()

    if current_sequence:
        if len(current_sequence) >= 100:
            sequences.append((current_title, current_sequence))

    return sequences


def write_gc_percentage(filename, sequences):
    with open(filename, "w") as file:
        for title, sequence in sequences:
            gc_percentage = calculate_gc_percentage(sequence)
            file.write(f">{title}\n")
            file.write(f"GC percentage: {gc_percentage:.2f}%\n\n")


input_filename = "Homologous Rodent Genes All Alignment 60.0% Cutoff_mod.txt"  # Replace with the path to your input FASTA file
output_filename = "Homologous Rodent Genes All Alignment 60.0% Cutoff_mod_GC_percent_100.txt"  # Replace with the desired output file path

sequences = read_fasta_file(input_filename)
write_gc_percentage(output_filename, sequences)

# Filter GC percentages for histogram and normal distribution plot
filtered_gc_percentages = [calculate_gc_percentage(sequence) for _, sequence in sequences]

# Plotting the histogram
plt.subplot(1, 2, 1)
plt.hist(filtered_gc_percentages, bins=50, edgecolor='black')
plt.xlabel('GC Percentage')
plt.ylabel('Frequency')
plt.title('Histogram of GC Percentage')

# Plotting the normal distribution curve
plt.subplot(1, 2, 2)
mean = np.mean(filtered_gc_percentages)
std = np.std(filtered_gc_percentages)
x = np.linspace(mean - 3*std, mean + 3*std, 100)
y = 1/(std * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mean)/std)**2)
plt.plot(x, y, color='blue')
plt.xlabel('GC Percentage')
plt.ylabel('Probability Density')
plt.title('Normal Distribution of GC Percentage')

plt.tight_layout()
plt.show()


#### 

import matplotlib.pyplot as plt
import numpy as np

# Function to parse GC percentages from a file
def parse_gc_percentages(filename):
    file = open(filename, 'r')
    gc_percentages = []
    for line in file:
        if line.startswith("GC percentage:"):
            percentage = float(line.split(":")[1].strip("%\n"))
            gc_percentages.append(percentage)
    file.close()
    return gc_percentages

# Define the filenames
file1_name = 'mousehotspot_mm39_modified_GC_percent.txt'
file2_name = 'Homologous Rodent Genes All Alignment 60.0% Cutoff_mod_GC_percent_100.txt'


# Parse GC percentages from each file
file1_gc_percentages = parse_gc_percentages(file1_name)
file2_gc_percentages = parse_gc_percentages(file2_name)

# Create histograms for each file
bins = np.arange(0, 110, 10)  # Define bin ranges

# Histogram for File 1
file1_frequencies, _ = np.histogram(file1_gc_percentages, bins=bins)
file1_total = len(file1_gc_percentages)
file1_percentages = file1_frequencies / file1_total * 100
file1_table = zip(bins[:-1], bins[1:], file1_frequencies, file1_percentages)
file1_header = ("Bin Start", "Bin End", "Frequency", "Percentage")

# Histogram for File 2 (Old)
file2_frequencies_old, _ = np.histogram(file2_gc_percentages, bins=bins)
file2_total = len(file2_gc_percentages)
file2_percentages_old = file2_frequencies_old / file2_total * 100
file2_table_old = zip(bins[:-1], bins[1:], file2_frequencies_old, file2_percentages_old)
file2_header_old = ("Bin Start", "Bin End", "Frequency", "Percentage")

# Calculate new frequencies in File 2 based on File 1 percentages and File 2 total
file2_percentages_new = file1_percentages
file2_frequencies_new = np.round(file2_percentages_new / 100 * file2_total).astype(int)
file2_table_new = zip(bins[:-1], bins[1:], file2_frequencies_new, file2_percentages_new)
file2_header_new = ("Bin Start", "Bin End", "Frequency", "Percentage")

# Calculate ratio of frequencies (updated/old)
frequency_ratio = file2_frequencies_new / file2_frequencies_old

# Create updated table with the frequency ratio
file2_table_updated = zip(bins[:-1], bins[1:], file2_frequencies_new, file2_percentages_new, frequency_ratio)
file2_header_updated = ("Bin Start", "Bin End", "Frequency", "Percentage", "Frequency Ratio")

# Write tables to separate output files
output_file1 = 'recombination_gc_match_file1.txt'
output_file2_old = 'recombination_gc_match_file2_old.txt'
output_file2_new = 'recombination_gc_match_file2_new.txt'

with open(output_file1, 'w') as file:
    file.write("Table for " + file1_name + "\n")
    file.write("{}\t{}\t{}\t{}\n".format(*file1_header))
    for row in file1_table:
        file.write("{}\t{}\t{}\t{:.2f}\n".format(*row))

with open(output_file2_old, 'w') as file:
    file.write("Table for " + file2_name + " (Old)\n")
    file.write("{}\t{}\t{}\t{}\n".format(*file2_header_old))
    for row in file2_table_old:
        file.write("{}\t{}\t{}\t{:.2f}\n".format(*row))

with open(output_file2_new, 'w') as file:
    file.write("Table for " + file2_name + " (Updated)\n")
    file.write("{}\t{}\t{}\t{}\n".format(*file2_header_updated))
    for row in file2_table_updated:
        file.write("{}\t{}\t{}\t{:.2f}\n".format(*row))

# Plot histograms
plt.figure(figsize=(12, 4))

# Histogram for File 1
plt.subplot(1, 3, 1)
plt.hist(file1_gc_percentages, bins=100, edgecolor='black')
plt.xlabel('GC Percentage')
plt.ylabel('Frequency')
plt.title('Mouse Recombination')

# Histogram for File 2 (Old)
plt.subplot(1, 3, 2)
plt.hist(file2_gc_percentages, bins=100, edgecolor='black')
plt.xlabel('GC Percentage')
plt.ylabel('Frequency')
plt.title('Mouse Intergenic')

# Histogram for File 2 (Updated)
plt.subplot(1, 3, 3)
plt.hist(file2_gc_percentages, bins=100, edgecolor='black')
plt.xlabel('GC Percentage')
plt.ylabel('Frequency')
plt.title('Mouse Intergenic (Updated)')

plt.tight_layout()
plt.show()

### intergenic frequency match
def parse_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            header = lines[i].strip()
            sequence = ""
            gc_percentage = None

            # Find the line with GC percentage
            while i + 1 < len(lines) and not lines[i + 1].startswith(">"):
                if "GC percentage:" in lines[i + 1]:
                    percentage_line = lines[i + 1].split(":")[1].strip()
                    gc_percentage = float(percentage_line.rstrip('%'))
                else:
                    sequence += lines[i + 1].strip()
                i += 1

            entry = {
                "header": header,
                "sequence": sequence,
                "gc_percentage": gc_percentage
            }
            sequences.append(entry)

        i += 1

    return sequences



def select_entries(entries, gc_ranges, num_entries):
    selected_entries = []

    for gc_range, num in zip(gc_ranges, num_entries):
        start, end = gc_range
        count = 0

        for entry in entries:
            gc_percentage = entry['gc_percentage']

            if start <= gc_percentage <= end:
                selected_entries.append(entry)
                count += 1

            if count >= num:
                break

    return selected_entries


def write_selected_entries_to_fasta(selected_entries, output_file):
    with open(output_file, 'w') as file:
        for entry in selected_entries:
            file.write(entry['header'] + '\n')
            file.write(entry['sequence'])

# Usage example
fasta_file = "Homologous Rodent Genes All Alignment 60.0% Cutoff_mod_GC_percent_100.txt"
output_file = "Homologous Rodent Genes All Alignment 60.0% Cutoff_mod_GC_percent_100_match_to_recombination_mm39.txt"

entries = parse_fasta(fasta_file)


gc_ranges = [(0,10), (10,20), (20,30), (30,40), (40,50), (50,60), (60,70), (70,80), (80,90), (90,100)]  # Example GC percentage ranges
num_entries = [2, 5, 183, 1111, 2361, 1280, 123, 25, 2, 0] # Replace with actual values (see below)

'''
The above values will need to be adjusted for different datasets
the values are based on the output_file2_new file from above
'''

selected_entries = select_entries(entries, gc_ranges, num_entries)

write_selected_entries_to_fasta(selected_entries, output_file)


###extract the gc_percent_matched lines from the original aligned file, and perform SNP analysis again


def extract_lines(file_a, file_b, output_file):
    lines_to_extract = set()
    with open(file_a, 'r') as file_a:
        for line in file_a:
            lines_to_extract.add(line.strip())

    with open(file_b, 'r') as file_b, open(output_file, 'w') as output:
        extract = False
        extracted_lines = 0
        for line in file_b:
            if line.strip() in lines_to_extract:
                extract = True

            if extract:
                output.write(line)
                extracted_lines += 1

                if extracted_lines == 6:
                    extract = False
                    extracted_lines = 0

            if extract and line == '\n':
                extract = False

# Usage example
file_a = 'Homologous Rodent Genes All Alignment 60.0% Cutoff_mod_GC_percent_100_match_to_recombination_mm39.txt'
file_b = 'Homologous Rodent Genes All Alignment 60.0% Cutoff_mod.txt'
output_file = 'Alignment_GC_percent_matched_to_recombination_mm39_mid100.txt'
extract_lines(file_a, file_b, output_file)

###now use the find_substituions.py file to get the substitutions from these gc-matched alignments