#### Use codes from: https://github.com/glarue/cdseq to get sequences of the open reading frame given a gtf file

python3 orf_seq.py Homo_sapiens.GRCh38.dna.toplevel.fa Homo_sapiens.GRCh38.108.chr.gtf >> human_orf.fa
python3 orf_seq.py Gorilla_gorilla.gorGor4.dna.toplevel.fa Gorilla_gorilla.gorGor4.108.chr.gtf >> gorilla_orf.fa
python3 orf_seq.py Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa Pan_troglodytes.Pan_tro_3.0.108.chr.gtf >> chimp_orf.fa

'''
##### 
Filter for ORF sequences that start with ATG and the sequence length is divisible by 3
'''

def is_valid_sequence(sequence):
    return sequence.startswith("ATG") and len(sequence) % 3 == 0

def process_entries(input_filename, valid_output_filename, invalid_output_filename):
    with open(input_filename, 'r') as input_file, \
         open(valid_output_filename, 'w') as valid_output_file, \
         open(invalid_output_filename, 'w') as invalid_output_file:

        current_entry = ''
        current_sequence = ''

        for line in input_file:
            if line.startswith('>'):
                # Process the previous entry, if any
                if current_entry and is_valid_sequence(current_sequence):
                    valid_output_file.write(current_entry + current_sequence + '\n')
                elif current_entry:
                    invalid_output_file.write(current_entry + current_sequence + '\n')

                # Start a new entry
                current_entry = line
                current_sequence = ''
            else:
                current_sequence += line.strip()

        # Process the last entry in the file
        if current_entry and is_valid_sequence(current_sequence):
            valid_output_file.write(current_entry + current_sequence + '\n')
        elif current_entry:
            invalid_output_file.write(current_entry + current_sequence + '\n')

# Usage
process_entries('human_orf.fa', 'human_orf_atg.txt', 'human_orf_nonatg.txt')

'''
##### From the previously aligned files (see needleman_alignment.py), 
filter for the sequences that meet the above requirements
'''
def extract_lines(input_file1, input_file2, matching_output_file, non_matching_output_file):
    with open(input_file1, 'r') as infile1:
        search_values = set()
        for line in infile1:
            if line.startswith('>') and len(line.split('\t')) >= 2:
                ensg_value = line.split('\t')[1].split()[0]
                if ensg_value.startswith("ENSG"):
                    search_values.add(ensg_value)

    with open(input_file2, 'r') as infile2, open(matching_output_file, 'w') as match_outfile, open(non_matching_output_file, 'w') as non_match_outfile:
        lines = infile2.readlines()

        for i, line in enumerate(lines):
            if line.startswith('>'):
                columns = line.strip().split('\t')
                if len(columns) >= 2 and any(value in columns[0] for value in search_values):
                    # Write the matching line and the next 5 lines to the matching output file
                    match_outfile.write(line)
                    for j in range(1, 6):
                        match_outfile.write(lines[i + j])
                else:
                    # Write the non-matching line to the non-matching output file
                    non_match_outfile.write(line)

if __name__ == "__main__":
    # Specify file names here
    input_file1 = "human_orf_atg.txt"
    input_file2 = "Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod.txt"
    matching_output_file = "Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg.txt"
    non_matching_output_file = "Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_noatg.txt"

    extract_lines(input_file1, input_file2, matching_output_file, non_matching_output_file)

##### Filter for aligned sequences that start aligning at position 1
input_file = "Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg.txt"
output_file = "Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1.txt"

with open(input_file, "r") as file:
    entries = file.read().split(">")[1:]

with open(output_file, "w") as output:
    for entry in entries:
        lines = entry.strip().split("\n")
        identifier_line = lines[0]
        gene_info = identifier_line.split(",")
        sequence_lines = lines[1:5]
        coordinate_line = lines[5]

        if coordinate_line.startswith("1:"):
            # Write the entire identifier line to the output file
            output.write(">" + identifier_line + "\n")
            output.write("\n".join(sequence_lines) + "\n")
            output.write(coordinate_line + "\n")

### Find the positions (within the aligned sequence) of the 3rd position synonymous codons (GC4)
# Define the codons of interest
codons_of_interest = ["GCT", "GCC", "GCA", "GCG", "GGT", "GGC", "GGA", "GGG",
                      "CCT", "CCC", "CCA", "CCG", "CTT", "CTC", "CTA", "CTG",
                      "TCT", "TCC", "TCA", "TCG", "ACT", "ACC", "ACA", "ACG",
                      "GTT", "GTC", "GTA", "GTG", "CGT", "CGC", "CGA", "CGG"]

# Open input and output files
with open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1.txt', 'r') as input_file, open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 GC4 Positions Chimp.txt', 'w') as output_file:
    # Iterate through lines in the input file
    for line in input_file:
        # Check if the line starts with ">"
        if line.startswith(">"):
            # Write the identifier to the output file
            output_file.write(line)

            # Update current sequence and position
            current_sequence = ""
            current_position = 1
            positions = []

            #*** Use this to skip the first sequence line (reference sequence) after the identifier to get the cousin sequence ***
            next(input_file)

            # Read the second sequence line after the identifier
            sequence_line = next(input_file).strip()
            #print (sequence_line)

            # Remove underscores from the sequence
            sequence_line = sequence_line.replace('_', '')

            # Check if the length of the sequence line is at least 3
            while len(sequence_line) >= 3:
                # Get the codon from the sequence line
                codon = sequence_line[:3]

                # Check if the codon is in the list of interest
                if codon in codons_of_interest:
                    # Collect the 3rd position of the codon
                    positions.append(current_position + 2)

                # Move to the next codon
                sequence_line = sequence_line[3:]
                current_position += 3

            # Write the positions to the output file if there are any
            if positions:
                output_file.write(','.join(map(str, positions)) + '\n')

'''
### 
Perform the substitution position analysis for the above filtered sequences
'''
def InDel_Count (n, x):
	count = 0
	for a in range (x):
		if n[a] == "_":
			count = count + 1
	return (count)
	
def nth_nonInDel_Character_Mutations (n, x):
	count = 0
	position = 0
	while count < x:
		if n[position] == "_":
			position = position + 1
		else:
			count = count + 1
			position = position + 1
	return (position)
	
# InDel Location:InDel Size   
def mutations (n):
	input_file = open(n+".txt", 'r')
	output_file = open(n+" Primate InDels.txt", 'w')
	output_file2 = open(n+" Primate SNP.txt", 'w')
	hAT = []
	hAG = []
	hAC = []
	hTA = []
	hTG = []
	hTC = []
	hGA = []
	hCGCA = []
	hGT = []
	hGC = []
	hCA = []
	hCT = []
	hCGTG = []
	hCG = []
	cAT = []
	cAG = []
	cAC = []
	cTA = []
	cTG = []
	cTC = []
	cGA = []
	cCGCA = []
	cGT = []
	cGC = []
	cCA = []
	cCT = []
	cCGTG = []
	cCG = []
	human_Deletion = []
	human_Insertion = []
	chimp_Deletion = []
	chimp_Insertion = []
	output_file.write ("Species" + '\t' + "Insertions" + '\t' + "Deletions" + '\t' + "Homologous Genes" + '\r')
	output_file2.write ("Species" + '\t' + "A to T" + '\t' + "A to G" + '\t' + "A to C" + '\t' + "T to A" + '\t' + "T to G" + '\t' + "T to C" + '\t' + "G to A" + '\t' + "G to T" + '\t' + "G to C" + '\t' + "C to A" + '\t' + "C to T" + '\t' + "C to G" + '\t' + "CG to TG" + '\t' + "CG to CA" + '\t' + "Homologous Genes" + '\r')
	name = input_file.readline()
	while name != "*":
		indel = 0  # size of the InDel
		tracker = 0  # location of the start of InDel
		human1 = input_file.readline()
		chimp = input_file.readline()
		human2 = input_file.readline()
		gorilla = input_file.readline()
		crop_record = input_file.readline()
		'''
		print (name)
		print (human1)
		print (chimp)
		print (human2)
		print (gorilla)
		print (crop_record)
		'''
		crop_parts = crop_record.split(':')
		crop_start = int(crop_parts[0])
		#print (crop_start)
		for a in range (len(human1)):
			if (human1[a] != chimp[a]):
				'''
				print (human1[a])
				print (chimp[a])
				print (a)
				'''
				if (human1[a] == "_" or chimp[a] == "_"):
					if indel == 0:
						tracker = a
					indel = indel + 1
				elif indel > 0:
					human1_Count = InDel_Count(human1, tracker)
					chimp_Count = InDel_Count(chimp, tracker)
					human_InDel = human1[tracker:a]
					chimp_InDel = chimp[tracker:a]
					human2_Count = nth_nonInDel_Character_Mutations (human2, (tracker-human1_Count))
					gorilla_InDel = gorilla [human2_Count:(human2_Count+indel+1)]
					if (gorilla_InDel == chimp_InDel) and (gorilla_InDel != human2_InDel):
						if gorilla_InDel.count('_') == 0:
							human_Deletion.append (str(tracker-human1_Count+crop_start) + ":" + str(indel))
						else:
							human_Insertion.append (str(tracker-human1_Count+crop_start) + ":" + str(indel))
					elif (gorilla_InDel == human2_InDel):
						if human_InDel.count('_') == 0:
							chimp_Deletion.append (str(tracker-human1_Count+crop_start) + ":" + str(indel))
						else:
							chimp_Insertion.append (str(tracker-human1_Count+crop_start) + ":" + str(indel))
					indel = 0
				else:
					human1_Count = InDel_Count(human1, a)
					chimp_Count = InDel_Count(chimp, a)
					human2_Count = nth_nonInDel_Character_Mutations (human2, (a-human1_Count))
					if gorilla[human2_Count] == chimp [a] and gorilla[human2_Count] != human1[a]:
						if (gorilla[human2_Count] == "A") and (human1[a] == "T"):
							hAT.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "A") and (human1[a] == "G"):
							hAG.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "A") and (human1[a] == "C"):
							hAC.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "T") and (human1[a] == "A"):
							hTA.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "T") and (human1[a] == "G"):
							hTG.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "T") and (human1[a] == "C"):
							hTC.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "G") and (human1[a] == "A"):
							if gorilla[human2_Count-1] == "C" and human1[a-1] == "C" and chimp[a-1] == "C":
								hCGCA.append (a-human1_Count+crop_start)
							else:
								hGA.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "G") and (human1[a] == "T"):
							hGT.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "G") and (human1[a] == "C"):
							hGC.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "C") and (human1[a] == "A"):
							hCA.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "C") and (human1[a] == "T"):
							if gorilla[human2_Count+1] == "G" and human1[a+1] == "G" and chimp[a+1] == "G":
								hCGTG.append(a-human1_Count+crop_start)
							else:
								hCT.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "C") and (human1[a] == "G"):
							hCG.append (a-human1_Count+crop_start)
					elif gorilla[human2_Count] == human1 [a] and gorilla[human2_Count] != chimp[a]:
						if (gorilla[human2_Count] == "A") and (chimp[a] == "T"):
							cAT.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "A") and (chimp[a] == "G"):
							cAG.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "A") and (chimp[a] == "C"):
							cAC.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "T") and (chimp[a] == "A"):
							cTA.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "T") and (chimp[a] == "G"):
							cTG.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "T") and (chimp[a] == "C"):
							cTC.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "G") and (chimp[a] == "A"):
							if gorilla[human2_Count-1] == "C" and human1[a-1] == "C" and chimp[a-1] == "C":
								cCGCA.append (a-human1_Count+crop_start)
							else:
								cGA.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "G") and (chimp[a] == "T"):
							cGT.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "G") and (chimp[a] == "C"):
							cGC.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "C") and (chimp[a] == "A"):
							cCA.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "C") and (chimp[a] == "T"):
							if gorilla[human2_Count+1] == "G" and human1[a+1] == "G" and chimp[a+1] == "G":
								cCGTG.append(a-human1_Count+crop_start)
							else:
								cCT.append (a-human1_Count+crop_start)
						elif (gorilla[human2_Count] == "C") and (chimp[a] == "G"):
							cCG.append (a-human1_Count+crop_start)
			else:
				if indel > 0:
					human1_Count = InDel_Count(human1, tracker)
					chimp_Count = InDel_Count(chimp, tracker)
					human_InDel = human1[tracker:a]
					chimp_InDel = chimp[tracker:a]
					#print (tracker)
					#print (human1_Count)
					human2_Count = nth_nonInDel_Character_Mutations (human2, (tracker-human1_Count))
					gorilla_InDel = gorilla [human2_Count:(human2_Count+indel)]
					human2_InDel = human2 [human2_Count:(human2_Count+indel)]
					if (gorilla_InDel == chimp_InDel) and (gorilla_InDel != human2_InDel):
						if gorilla_InDel.count('_') == 0:
							human_Deletion.append (str(tracker-human1_Count+crop_start) + ":" + str(indel))
						else:
							human_Insertion.append (str(tracker-human1_Count+crop_start) + ":" + str(indel))
					elif (gorilla_InDel == human2_InDel):
						if human_InDel.count('_') == 0:
							chimp_Deletion.append (str(tracker-human1_Count+crop_start) + ":" + str(indel))
						else:
							chimp_Insertion.append (str(tracker-human1_Count+crop_start) + ":" + str(indel))
					indel = 0
		output_file.write (name[0:-1] + '\t' + crop_record)
		output_file.write ("InDels in Humans"+ '\t' + str(human_Insertion) + '\t' + str(human_Deletion) + '\t' + "." + '\r')
		output_file.write ("InDels in Chimps"+ '\t' + str(chimp_Insertion) + '\t' + str(chimp_Deletion) + '\t' + "." + '\r')
		output_file2.write (name[0:-1] + '\t' + crop_record)
		output_file2.write ("SNPs in Humans" + '\t' + str(hAT) + '\t' + str(hAG) + '\t' + str(hAC) + '\t' + str(hTA) + '\t' + str(hTG) + '\t' + str(hTC) + '\t' + str(hGA) + '\t' + str(hGT) + '\t' + str(hGC) + '\t' + str(hCA) + '\t' + str(hCT) + '\t' + str(hCG) + '\t' + str(hCGTG) + '\t' + str(hCGCA) + '\t' + "END" + '\r')
		output_file2.write ("SNPs in Chimps" + '\t' + str(cAT) + '\t' + str(cAG) + '\t' + str(cAC) + '\t' + str(cTA) + '\t' + str(cTG) + '\t' + str(cTC) + '\t' + str(cGA) + '\t' + str(cGT) + '\t' + str(cGC) + '\t' + str(cCA) + '\t' + str(cCT) + '\t' + str(cCG) + '\t' + str(cCGTG) + '\t' + str(cCGCA) + '\t' + "END" + '\r')
		hAT = []
		hAG = []
		hAC = []
		hTA = []
		hTG = []
		hTC = []
		hGA = []
		hCGCA = []
		hGT = []
		hGC = []
		hCA = []
		hCT = []
		hCGTG = []
		hCG = []
		cAT = []
		cAG = []
		cAC = []
		cTA = []
		cTG = []
		cTC = []
		cGA = []
		cCGCA = []
		cGT = []
		cGC = []
		cCA = []
		cCT = []
		cCGTG = []
		cCG = []
		human_Deletion = []
		human_Insertion = []
		chimp_Deletion = []
		chimp_Insertion = []
		name = input_file.readline()
	output_file.write("*")
	output_file2.write("*")
	input_file.close()
	output_file.close()

#mutations ("Canid Genes Aligned 60.0% Cutoff")
mutations ("Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1")

'''
##### 
Separate the triple aligned files into the reference and the cousin
'''
def process_input_file(input_filename, output_filename1, output_filename2):
    with open(input_filename, 'r') as input_file, \
         open(output_filename1, 'w') as output_file1, \
         open(output_filename2, 'w') as output_file2:

        # Read the first line from the input file
        first_line = input_file.readline()

        # Write the first line to both output files
        output_file1.write(first_line)
        output_file2.write(first_line)

        # Read and process the input file in groups of 3 lines
        while True:
            line1 = input_file.readline().strip()
            line2 = input_file.readline().strip()
            line3 = input_file.readline().strip()

            # Break the loop if there are no more lines to read
            if not line1 or not line2 or not line3:
                break

            # Write the first and second lines to output file 1
            output_file1.write(f"{line1}\n{line2}\n")

            # Write the first and third lines to output file 2
            output_file2.write(f"{line1}\n{line3}\n")

# Example usage
input_filename = 'Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP.txt'
output_filename1 = 'Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP Human.txt'
output_filename2 = 'Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP Chimp.txt'

process_input_file(input_filename, output_filename1, output_filename2)

'''
##### 
In the substitutions file, identify the substitutions that fall on a GC4 position
'''
import re

# Function to parse numbers inside square brackets
def extract_numbers(text):
    return [int(match.group()) for match in re.finditer(r'\b(\d+)\b', text)]

# Define the output file
output_file = 'Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP Chimp GC4.txt'

# Open the output file for writing
with open(output_file, 'w') as output:
    # Read and parse both files simultaneously
    with open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP Chimp.txt', 'r') as file1, open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 GC4 Positions Chimp.txt', 'r') as file2:
        for line1, line2 in zip(file1, file2):
            if line1.startswith('>') and line2.startswith('>'):
                current_entry1 = line1.strip()[1:]
                current_entry2 = line2.strip()[1:]
                numbers_file2 = set(map(int, next(file2).strip().split(',')))

                # Extract original column data
                original_data = next(file1).strip()

                # Find matching numbers and apply the conditions
                matching_numbers = [num for num in extract_numbers(original_data) if num in numbers_file2 and num - 1 not in numbers_file2 and num - 2 not in numbers_file2]

                # Replace the original numbers with the matching numbers in the output
                updated_data = re.sub(r'\b(\d+)\b', lambda x: str(x.group()) if int(x.group()) not in matching_numbers else f"{x.group()}*", original_data)

                # Write the result to the output file
                output.write(f">{current_entry1}\n{updated_data}\n")


# Read the input file
with open('Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP Human GC4.txt', 'r') as file:
    lines = file.readlines()

# Function to process SNP information
def process_snp_info(snp_info):
    # Remove numbers not followed by '*'
    snp_info = re.sub(r'\b\d+(?!\*)\b', '', snp_info)
    # Replace "[," with "[" and ",]" with "]"
    snp_info = snp_info.replace("[,", "[").replace(",]", "]")
    return snp_info

# Output file path
output_file_path = 'Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP Human GC4 Final.txt'

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
'''
##### Manually delete all spaces and *

##### Since the ORF has different length, need to split all ORFs into 40 equal-sized bins
'''
import re
from pathlib import Path
import sys

# Read data from a file
file_path = Path("Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP Chimp GC4 Final.txt")

with file_path.open("r") as file:
    data = file.read()

# Split the data into individual gene entries
gene_entries = data.split('>')

# Output to a file
output_file_path = Path("ORF SNP 40Bins Chimp GC4.txt")
with output_file_path.open("w") as output_file:
    # Redirect standard output to the file
    sys.stdout = output_file

    # Iterate through each gene entry
    for gene_entry in gene_entries:
        if gene_entry.strip():
            # Extract second column values from Gene Information lines
            gene_info_lines = [line for line in gene_entry.split('\n') if line.startswith("ENS")]

            # Extract "SNPs in Humans" lines
            snps_lines = [line for line in gene_entry.split('\n') if line.startswith("SNPs in Chimps")]

            # Check if there is at least one "SNPs in Humans" line
            if snps_lines:
                snps_line = snps_lines[0]  # Assuming there's only one "SNPs in Humans" line

                # Extract values between square brackets using regular expression
                snp_columns = re.findall(r'\[([^]]*)\]', snps_line)

                # Check if there are any SNP values
                if not snp_columns or not any(snp_columns):
                    # Handle the case when there are no SNP values
                    continue  # Skip to the next gene

                # Split values in each column
                snp_columns_values = [column.split(', ') for column in snp_columns]

                # Extract numbers from the second column of gene_info_lines
                second_column_values = [int(line.split('\t')[1]) for line in gene_info_lines]

                # Use the maximum value from the second column as the basis for bin ranges
                max_value = max(second_column_values)
                bin_width = max_value / 40
                bin_ranges = [(int(i * bin_width), int((i + 1) * bin_width)) for i in range(40)]

                # Count the number of values in each bin for each column
                bin_counts_per_column = [
                    [sum(start <= int(value) < end for value in column if value) for start, end in bin_ranges]
                    for column in snp_columns_values
                ]

                # Display bin counts for each column
                for bin_count in zip(*bin_counts_per_column):
                    print("\t".join(map(str, bin_count)))
            else:
                continue

    # Reset standard output
    sys.stdout = sys.__stdout__

print(f"Output written to: {output_file_path}")


##### Split the data such that each file contains the values for one ORF divided into 40 bins
import os

# Create a new folder if it doesn't exist
output_folder = 'ORF SNP 40Bins Human GC4'
os.makedirs(output_folder, exist_ok=True)

with open('ORF SNP 40Bins Human GC4.txt', 'r') as input_file:
    lines = input_file.readlines()

# Define the number of lines per file
lines_per_file = 40

# Calculate the number of files needed
num_files = len(lines) // lines_per_file + (len(lines) % lines_per_file > 0)

# Create and write to each file in the output folder
for i in range(num_files):
    start_idx = i * lines_per_file
    end_idx = (i + 1) * lines_per_file
    output_filename = os.path.join(output_folder, f'ORF SNP 40Bins Human GC4_{i + 1}.txt')

    with open(output_filename, 'w') as output_file:
        output_file.writelines(lines[start_idx:end_idx])

'''
######### 
Count the number of times each coordinate position is accounted for in the alignment in the bin
'''
import re
from pathlib import Path
import sys

# Read data from a file
file_path = Path("Homologous Primate ORF Genes All Alignment 60.0% Cutoff_mod_atg_start1 Primate SNP Human GC4 Final.txt")

with file_path.open("r") as file:
    data = file.read()

# Split the data into individual gene entries
gene_entries = data.split('>')

# Output to a file for SNP counts
output_file_path = Path("ORF_SNP_Counts_per_Gene GC4.txt")
with output_file_path.open("w") as output_file:
    # Iterate through each gene entry
    for gene_entry in gene_entries:
        if gene_entry.strip():
            # Extract second column values from Gene Information lines
            gene_info_lines = [line for line in gene_entry.split('\n') if line.startswith("ENSG")]

            # Extract "SNPs in Humans" lines
            snps_lines = [line for line in gene_entry.split('\n') if line.startswith("SNPs in Humans")]

            # Check if there is at least one "SNPs in Humans" line
            if snps_lines and gene_info_lines:
                snps_line = snps_lines[0]  # Assuming there's only one "SNPs in Humans" line

                # Extract values between square brackets using regular expression
                snp_columns = re.findall(r'\[([^]]*)\]', snps_line)

                # Check if there are any SNP values
                if not snp_columns or not any(snp_columns):
                    # Handle the case when there are no SNP values
                    continue  # Skip to the next gene

                # Split values in each column
                snp_columns_values = [column.split(', ') for column in snp_columns]

                # Extract numbers from the second column of gene_info_lines
                second_column_values = [int(line.split('\t')[1]) for line in gene_info_lines]

                # Extract values from the third column in the format of "start:end"
                third_column_values = [line.split('\t')[2] for line in gene_info_lines]

                # Use the maximum value from the second column as the basis for bin ranges
                max_value = max(second_column_values)
                bin_width = max_value / 40
                bin_ranges = [(int(i * bin_width), int((i + 1) * bin_width)) for i in range(40)]

                # Count the number of times values in the third column fall into each bin
                bin_counts_per_gene = [0] * 40
                for start_end in third_column_values:
                    start, end = map(int, start_end.split(':'))
                    bin_flags = [any(bin_start <= x <= bin_end for x in range(start, end + 1)) for bin_start, bin_end in bin_ranges]
                    bin_counts_per_gene = [count + flag for count, flag in zip(bin_counts_per_gene, bin_flags)]

                # Display bin counts for each gene
                #output_file.write(f"Gene: {gene_info_lines[0].split()[0]}\n")
                output_file.write("\t".join(map(str, bin_counts_per_gene)) + "\n")
print(f"SNP Counts per Gene written to: {output_file_path}")