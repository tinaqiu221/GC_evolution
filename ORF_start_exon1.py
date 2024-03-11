'''
##########
Previously in the code "get_homologou_sequences.py", we generated annotation files 
that contain the best TSS for each gene and its associated transcript information.
Using this file, we will generate a file contains the coordinates for the TSS and the first exon/intron boundary
'''
file_path = 'Chimpanzee_Protein_Coding_Genes_GTF_Best_TSS.txt'

# This function will process the file according to the described logic
def process_file(file_path):
    # Dictionary to hold exon lines for each gene ID, with strand consideration
    gene_exons = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into columns
            columns = line.strip().split('\t')
            # Check if the line contains an exon feature and has the required number of columns
            if len(columns) >= 9 and columns[2] == 'exon':
                gene_id = columns[8]
                #strand = columns[6]
                # For genes on the + strand, keep the first exon encountered
                if gene_id not in gene_exons:
                    gene_exons[gene_id] = line
                # For genes on the - strand, keep replacing with the latest exon encountered
                #elif strand == '-':
                    #gene_exons[gene_id] = line

    # Output file path
    output_file_path = 'Chimpanzee_Protein_Coding_Genes_GTF_Best_TSS_EIB.txt'
    # Write the selected exon lines to the output file
    with open(output_file_path, 'w') as output_file:
        for exon_line in gene_exons.values():
            output_file.write(exon_line)
    # Return the path to the output file for download or further processing
    return output_file_path

# Process the file and get the path to the output file
output_file_path = process_file(file_path)
output_file_path

'''
###################
Now we just want to check if the beginning of the ORF falls in between the TSS and the exon/intron boundary as outlined above.
The ORF sequences are generated based on Github Repository https://github.com/glarue/cdseq
'''
# Define a function to read and process the input files
def process_files(file1_path, file2_path, output_path):
    # Read file 2 and store relevant info in a dictionary {Gene ID: [chromosome, start, end, strand]}
    gene_info = {}
    with open(file2_path, 'r') as file2:
        for line in file2:
            parts = line.strip().split('\t')
            chromosome, _, _, start, end, _, strand, _, gene_id, transcript_id = parts
            gene_info[gene_id] = [chromosome, int(start), int(end), strand]

    # Open file 1 and the output file
    with open(file1_path, 'r') as file1, open(output_path, 'w') as output:
        for line in file1:
            if line.startswith('>'):  # Header line in file 1
                parts = line.strip().split('\t')
                gene_id = parts[1]
                strand = parts[3]
                genomic_range = parts[4]
                # Correctly parse the range start and end from the genomic range
                range_start, range_end = map(int, genomic_range.split(':'))

                # Check if Gene ID from file 1 exists in file 2
                if gene_id in gene_info:
                    _, file2_start_pos, file2_end_pos, _ = gene_info[gene_id]

                    # Apply the logic based on the strand from file 1
                    if strand == "+":
                        # Check if the start of the genomic range falls within file 2's positions
                        if file2_start_pos <= range_start <= file2_end_pos:
                            output.write(line)  # Write the header line
                            output.write(next(file1))  # Write the nucleotide sequence
                    elif strand == "-":
                        # Check if the end of the genomic range falls within file 2's positions
                        if file2_start_pos <= range_end <= file2_end_pos:
                            output.write(line)  # Write the header line
                            output.write(next(file1))  # Write the nucleotide sequence


# Specify the paths to your input files and the output file
file1_path = 'dog_orf.fa'
file2_path = 'dog_Protein_Coding_Genes_GTF_Best_TSS_EIB.txt'
output_path = 'dog_orf_firstexoninorf.txt'

# Call the function with the file paths
process_files(file1_path, file2_path, output_path)
