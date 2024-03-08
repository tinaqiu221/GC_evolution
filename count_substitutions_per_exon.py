'''
##########################
Find ORF sequences that only have 4 or more exons, 
then use that to filter the previous substitutions results for genes that have 4 or more exons
'''
def process_file(input_file, output_file):
    with open(input_file, 'r') as file_in, open(output_file, 'w') as file_out:
        lines = file_in.readlines()
        for i in range(len(lines)):
            if lines[i].startswith(">"):
                # Splitting column 7 by ","
                values = lines[i].split("\t")[6].split(",")
                if len(values) >= 4:
                    # Writing the current line and the following line to the output file
                    file_out.write(lines[i])
                    if i + 1 < len(lines):
                        file_out.write(lines[i + 1])

# Usage example
input_file = "human_orf.fa"
output_file = "human_orf_4exons.txt"
process_file(input_file, output_file)

def process_files(file1, file2, output_file1, output_file2):
    # Open output files for writing
    with open(output_file1, 'w') as out_file1, open(output_file2, 'w') as out_file2:
        with open(file2, 'r') as f2:
            match_found = False
            for line2 in f2:
                if line2.startswith(">"):
                    match_found = False
                    search_term = line2.strip().split("\t")[1]
                    with open(file1, 'r') as f1:
                        for line1 in f1:
                            if search_term in line1:
                                out_file1.write(line1)
                                out_file1.write(next(f1))  # Write the following line
                                match_found = True
                                break
                if match_found:
                    out_file2.write(line2)
                    out_file2.write(next(f2))  # Write the following line

# Example usage:
file2 = "human_orf_4exons.txt"
file1 = "Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP Human GC4 Final.txt"
output_file2 = "human_orf_4exons_final.txt"
output_file1 = "Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP human 4exons.txt"
process_files(file1, file2, output_file1, output_file2)

'''
######################
Find the relative start and end coordinates for exon 1 or exon 4
'''

# Open the input file
with open("human_orf_4exons_final.txt", "r") as file:
    # Open output files for writing
    with open("human_orf_4exons_final_exon1_coordinates.txt", "w") as outfile1, open("human_orf_4exons_final_exon4_coordinates.txt", "w") as outfile2:
        # Iterate through each line in the file
        for line in file:
            # Check if the line starts with ">"
            if line.startswith(">"):
                # Split the line by tabs
                parts = line.strip().split("\t")
                # Assign gene_name from the second column
                gene_name = parts[1]
                # Split exon coordinates by comma
                exon_coordinates = parts[6].split(",")
                # Initialize variables to store exon start and end values
                exon1_start, exon1_end = "", ""
                exon2_start, exon2_end = "", ""
                exon3_start, exon3_end = "", ""
                exon4_start, exon4_end = "", ""
                # Check the strand
                if parts[3] == "+":
                    # Assign values for + strand
                    exon1_start, exon1_end = map(int, exon_coordinates[0].split("-"))
                    exon2_start, exon2_end = map(int, exon_coordinates[1].split("-"))
                    exon3_start, exon3_end = map(int, exon_coordinates[2].split("-"))
                    exon4_start, exon4_end = map(int, exon_coordinates[3].split("-"))
                elif parts[3] == "-":
                    # Assign values for - strand
                    exon1_end, exon1_start = map(int, exon_coordinates[-1].split("-"))
                    #print (exon1_start, exon1_end)
                    exon2_end, exon2_start = map(int, exon_coordinates[-2].split("-"))
                    exon3_end, exon3_start = map(int, exon_coordinates[-3].split("-"))
                    exon4_end, exon4_start = map(int, exon_coordinates[-4].split("-"))

                # Calculate differences for output file 1
                diff_exon1 = abs(exon1_end - exon1_start) + 1
                # Write to output file 1
                outfile1.write(f"{gene_name}\t1\t{diff_exon1}\n")

                # Calculate differences for output file 2
                diff_start = abs(exon1_end - exon1_start) + abs(exon2_end - exon2_start) + abs(exon3_end - exon3_start) + 3
                diff_end = abs(exon4_end - exon4_start) + diff_start + 1
                # Write to output file 2
                outfile2.write(f"{gene_name}\t{diff_start}\t{diff_end}\n")

'''
####################
Write all exon1 or exon4 start/end coordinate ranges into positions
'''

input_filename = "chimp_orf_4exons_final_exon1_coordinates.txt"
output_filename = "chimp_orf_4exons_final_exon1_coordinates_all_positions.txt"

with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
    for line in input_file:
        parts = line.strip().split("\t")
        gene_id = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        
        output_file.write(">"+gene_id + "\n")
        output_file.write(",".join(str(i) for i in range(start, end + 1)) + "\n")
'''
################
Add * after substitutions that match the previous positions
'''
import re

# Function to parse numbers inside square brackets
def extract_numbers(text):
    return [int(match.group()) for match in re.finditer(r'\b(\d+)\b', text)]

# Define the output file
output_file = 'Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP chimp 4exons exon1.txt'

# Open the output file for writing
with open(output_file, 'w') as output:
    # Read and parse both files simultaneously
    with open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP chimp 4exons.txt', 'r') as file1, open('chimp_orf_4exons_final_exon1_coordinates_all_positions.txt', 'r') as file2:
        for line1, line2 in zip(file1, file2):
            if line1.startswith('>') and line2.startswith('>'):
                current_entry1 = line1.strip()[1:]
                current_entry2 = line2.strip()[1:]
                numbers_file2 = set(map(int, next(file2).strip().split(',')))

                # Extract original column data
                original_data = next(file1).strip()

                # Find matching numbers and apply the conditions
                matching_numbers = [num for num in extract_numbers(original_data) if num in numbers_file2]

                # Replace the original numbers with the matching numbers in the output
                updated_data = re.sub(r'\b(\d+)\b', lambda x: str(x.group()) if int(x.group()) not in matching_numbers else f"{x.group()}*", original_data)

                # Write the result to the output file
                output.write(f">{current_entry1}\n{updated_data}\n")

'''
#### 
Delete the ones with * added
'''
import re

# Read the input file
with open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP chimp 4exons exon1.txt', 'r') as file:
    lines = file.readlines()

# Function to process SNP information
def process_snp_info(snp_info):
    # Remove numbers not followed by '*'
    snp_info = re.sub(r'\b\d+(?!\*)\b', '', snp_info)
    # Replace "[," with "[" and ",]" with "]"
    snp_info = snp_info.replace("[,", "[").replace(",]", "]")
    return snp_info

# Output file path
output_file_path = 'Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP chimp 4exons exon1 mod.txt'

# Open the output file for writing
with open(output_file_path, 'w') as output_file:
    # Iterate through lines and process SNP information
    for i in range(len(lines)):
        if lines[i].startswith('>'):
            output_file.write(lines[i])
            # Make sure we don't go beyond the end of the file
            if i + 1 < len(lines):
                snp_info_line = process_snp_info(lines[i + 1].strip().replace(",", ""))
                # Write the processed line to the output file
                output_file.write(snp_info_line + '\n')

def process_file(input_file, output_file):
    with open(input_file, 'r') as f_in:
        lines = f_in.readlines()

    # Remove triple space bars
    lines = [line.replace("   ", "") for line in lines]
    # Remove double space bars
    lines = [line.replace("  ", "") for line in lines]
    # Remove single space bars
    lines = [line.replace(" ", "") for line in lines]
    # Replace "*]" with "]"
    lines = [line.replace("*]", "]") for line in lines]
    # Replace * with ", "
    lines = [line.replace("*", ", ") for line in lines]
    # Replace "SNPsin" with "SNPs in "
    lines = [line.replace("SNPsin", "SNPs in ") for line in lines]
    # Remove lines that start with ">"
    lines = [line for line in lines if not line.startswith(">")]

    # Add the header as the first line
    header = "Species\tA to T\tA to G\tA to C\tT to A\tT to G\tT to C\tG to A\tG to T\tG to C\tC to A\tC to T\tC to G\tCG to TG\tCG to CA\tHomologous Genes\n"
    lines.insert(0, header)

    with open(output_file, 'w') as f_out:
        f_out.write("".join(lines))

# Example usage
input_filename = "Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP chimp 4exons exon1 mod.txt"  # Change this to your input file name
output_filename = "Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP chimp 4exons exon1 mod Final.txt"  # Change this to your desired output file name

process_file(input_filename, output_filename)

'''
#######
Count the total number of substituions for exon 1 or exon 4
'''

import re

# Initialize a line counter
total_lines = 0

# Read the file
with open("Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP human 4exons exon1 mod Final.txt", "r") as file:
    lines = file.readlines()

# Initialize a counter for column 3 integer values
AtoG_count = 0
AtoC_count = 0
TtoG_count = 0
TtoC_count = 0
GtoA_count = 0
GtoT_count = 0
CtoA_count = 0
CtoT_count = 0
        
for line in lines:
    total_lines += 1
    # Split the line by tab delimiter
    columns = line.strip().split("\t")
    # Extract the value in the third column
    AtoG_value = columns[2]
    # Use regular expression to extract integers within square brackets
    integers_in_brackets = re.findall(r'\[(\d+)\]', AtoG_value)
    # Iterate through extracted integers and count those that are integers
    for integer in integers_in_brackets:
        if integer.isdigit():
            AtoG_count += 1

    AtoC_value = columns[3]
    integers_in_brackets = re.findall(r'\[(\d+)\]', AtoC_value)
    for integer in integers_in_brackets:
        if integer.isdigit():
            AtoC_count += 1

    TtoG_value = columns[5]
    integers_in_brackets = re.findall(r'\[(\d+)\]', TtoG_value)
    for integer in integers_in_brackets:
        if integer.isdigit():
            TtoG_count += 1

    TtoC_value = columns[6]
    integers_in_brackets = re.findall(r'\[(\d+)\]', TtoC_value)
    for integer in integers_in_brackets:
        if integer.isdigit():
            TtoC_count += 1

    GtoA_value = columns[7]
    integers_in_brackets = re.findall(r'\[(\d+)\]', GtoA_value)
    for integer in integers_in_brackets:
        if integer.isdigit():
            GtoA_count += 1

    GtoT_value = columns[8]
    integers_in_brackets = re.findall(r'\[(\d+)\]', GtoT_value)
    for integer in integers_in_brackets:
        if integer.isdigit():
            GtoT_count += 1

    CtoA_value = columns[10]
    integers_in_brackets = re.findall(r'\[(\d+)\]', CtoA_value)
    for integer in integers_in_brackets:
        if integer.isdigit():
            CtoA_count += 1

    CtoT_value = columns[11]
    integers_in_brackets = re.findall(r'\[(\d+)\]', CtoT_value)
    for integer in integers_in_brackets:
        if integer.isdigit():
            CtoT_count += 1

total_gain = AtoG_count + AtoC_count + TtoG_count + TtoC_count
total_loss = GtoA_count + GtoT_count + CtoA_count + CtoT_count
# Print the total count of integer values in column 3
print("Total number of gain of GC substitutions:", total_gain)
print("Total number of loss of GC substitutions:", total_loss)
print("Net Gain-Loss:", total_gain-total_loss)
print("Total number of genes:", total_lines)
print("Net Gain-Loss per gene analyzed:", (total_gain-total_loss)/total_lines)

'''
#######
Count the total number of substitutions for exon 1 or exon 4 taking into account the size of the exons
'''

# Open the first input file for reading
with open('human_orf_4exons_final_exon1_coordinates.txt', 'r') as file1:
    # Read lines from the first input file
    lines_file1 = file1.readlines()

# Open the second input file for reading
with open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP human 4exons exon1 mod Final.txt', 'r') as file2:
    # Skip the first line of the second input file
    next(file2)
    # Read lines from the second input file
    lines_file2 = file2.readlines()

# Open a new output file for writing
with open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP human 4exons exon1 mod Final size.txt', 'w') as output_file:
    output_file.write("Species\tA to T\tA to G\tA to C\tT to A\tT to G\tT to C\tG to A\tG to T\tG to C\tC to A\tC to T\tC to G\tCG to TG\tCG to CA\tHomologous Genes\tExon Size\n")
    # Iterate over the lines of both input files simultaneously
    for line_file1, line_file2 in zip(lines_file1, lines_file2):
        # Split the lines by tabs
        parts_file1 = line_file1.strip().split('\t')
        parts_file2 = line_file2.strip().split('\t')
        
        # Calculate the difference between the third and second columns of the first file
        difference = int(parts_file1[2]) - int(parts_file1[1])
        
        # Append the difference as a new column to the line of the second file
        output_line = '\t'.join(parts_file2 + [str(difference)]) + '\n'
        
        # Write the updated line to the output file
        output_file.write(output_line)

# Define a function to parse each line and perform the required calculations

def process_line(line):

    # Split the line by tabs
    parts = line.strip().split('\t')
    
    # Extract values from columns 3, 4, 6, and 7 and count the number of values
    count_values1 = sum(1 for x in parts[2].split(',') if x.strip('[]')) + \
                   sum(1 for x in parts[3].split(',') if x.strip('[]')) + \
                   sum(1 for x in parts[5].split(',') if x.strip('[]')) + \
                   sum(1 for x in parts[6].split(',') if x.strip('[]'))
    count_values2 = sum(1 for x in parts[7].split(',') if x.strip('[]')) + \
                   sum(1 for x in parts[8].split(',') if x.strip('[]')) + \
                   sum(1 for x in parts[10].split(',') if x.strip('[]')) + \
                   sum(1 for x in parts[11].split(',') if x.strip('[]'))
    
    difference = count_values1 - count_values2
    #print (count_values)
    # Extract the divisor from column 17 and convert it to float
    divisor = float(parts[16])
    
    # Check if divisor is zero to avoid division by zero
    if divisor == 0:
        return 0  # or handle it according to your use case
    
    # Divide the count of values by the divisor and return the result
    return difference / divisor


# Initialize a variable to store the sum of results
total_sum = 0
total_lines = 0

# Open the input file for reading and the output file for writing
with open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod Primate SNP human 4exons exon1 mod Final size.txt', 'r') as input_file, open('ORF Net GC human exon1.txt', 'w') as output_file:
    # Skip the header line
    next(input_file)
    
    # Iterate over each line in the input file
    for line in input_file:
        total_lines += 1
        
        # Process the line and add the result to the total sum
        result = process_line(line)
        total_sum += result
        
        # Write the result to the output file
        output_file.write(f'{result}\n')

# Print the total sum
print("Gain of GC minus loss of GC per exon size:", total_sum)
print("Gain of GC minus loss of GC per exon size per gene:", total_sum/total_lines)
print (total_lines)