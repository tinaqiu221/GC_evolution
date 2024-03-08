'''
##########################################
After doing the Needleman Alignment and applying the cutoffs,
we will look at the aligned sequences and find the postisitions for all the substitutions 
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
mutations ("Homologous Canid Genes Aligned 60.0% Cutoff")

'''
######################################
Separate the cousin and reference substitutions into 2 files
'''
def separate (n):
	input_file = open(n+".txt", 'r')
	output_file1 = open(n+" Dog.txt", 'w')
	output_file2 = open(n+" Fox.txt", 'w')
	name = input_file.readline()
	output_file1.write (name)
	output_file2.write (name)
	name = input_file.readline()
	while name != "*":
		dogSNP = input_file.readline()
		dogSNP = dogSNP.replace ("[", "")
		dogSNP = dogSNP.replace ("]", "")
		foxSNP = input_file.readline()
		foxSNP = foxSNP.replace ("[", "")
		foxSNP = foxSNP.replace ("]", "")
		output_file1.write (dogSNP)
		output_file2.write (foxSNP)
		name = input_file.readline()
	input_file.close()
	output_file1.close()
	output_file2.close()

separate ("Homologous Canid Genes Aligned 60.0% Cutoff Primate SNP")

'''
###############################################
Count the number of times each coordinate position is accounted for in the alignment
'''

# Open the input file in read mode
input_file_path = "Homologous Canid Genes Aligned 60.0% Cutoff.txt"  # Update with your input file path
output_file_path = "Homologous Canid Genes Aligned 60.0% Cutoff_Coor.txt"  # Update with your output file path

with open(input_file_path, "r") as input_file:
    # Read all lines from the input file
    lines = input_file.readlines()

    # Extract every 6th line
    extracted_lines = lines[5::6]

    # Open the output file in write mode
    with open(output_file_path, "w") as output_file:
        # Write the extracted lines to the output file
        for line in extracted_lines:
            output_file.write(line)

# Function to count occurrences of numbers in a range
def count_occurrences(range_str, num_list):
    """
    Given a range string in the format "start:end", this function
    counts how many times the numbers from start to end (inclusive)
    appear in the provided num_list.
    """
    start, end = map(int, range_str.split(':'))
    count = 0
    for num in num_list:
        if start <= num <= end:
            count += 1
    return count

# Read file and count occurrences of each number from 1 to 1000
filename = "Homologous Canid Genes Aligned 60.0% Cutoff_Coor.txt"  # Replace with your input file name
num_list = list(range(1, 5001))  # Numbers from 1 to 1000
count_dict = {num: 0 for num in num_list}  # Dictionary to store counts
with open(filename, "r") as file:
    for line in file:
        range_str = line.split()[0]  # Extract range string from first column
        for num in num_list:
            count_dict[num] += count_occurrences(range_str, [num])

# Write counts to a new file
output_filename = "Homologous Canid Genes Aligned 60.0% Cutoff_Coor_Count.txt"  # Replace with your output file name
with open(output_filename, "w") as outfile:
    for num, count in count_dict.items():
        outfile.write(f"{num}\t {count}\n")


'''
##################################################
Calculate the net change in GC per gene per nucleotide and output a final file with desired format
'''

input_file = 'Homologous Canid Genes Aligned 60.0% Cutoff Primate SNP Fox.txt'
output_file = 'Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_GainofGC.txt'

columns_to_extract = [3, 4, 6, 7]  # Adjusted for 0-based indexing

# Read input file and extract specified columns
data = []
with open(input_file, 'r') as f:
    for line in f:
        columns = line.strip().split('\t')
        extracted_data = [columns[i - 1] if i <= len(columns) else "" for i in columns_to_extract]
        data.append(extracted_data)

# Write extracted columns to output file as a single column
with open(output_file, 'w') as f:
    for row in data:
        f.write('\n'.join(row) + '\n')

# Replace ", " with a new line character '\n' in the output file
with open(output_file, 'r') as f:
    content = f.read()
    content = content.replace(', ', '\n')

with open(output_file, 'w') as f:
    f.write(content)


input_file = 'Homologous Canid Genes Aligned 60.0% Cutoff Primate SNP Fox.txt'
output_file = 'Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_LossofGC.txt'

columns_to_extract = [8, 9, 11, 12]  # Adjusted for 0-based indexing

# Read input file and extract specified columns
data = []
with open(input_file, 'r') as f:
    for line in f:
        columns = line.strip().split('\t')
        extracted_data = [columns[i - 1] if i <= len(columns) else "" for i in columns_to_extract]
        data.append(extracted_data)

# Write extracted columns to output file as a single column
with open(output_file, 'w') as f:
    for row in data:
        f.write('\n'.join(row) + '\n')

# Replace ", " with a new line character '\n' in the output file
with open(output_file, 'r') as f:
    content = f.read()
    content = content.replace(', ', '\n')

with open(output_file, 'w') as f:
    f.write(content)

def count_numbers(filename):
    # Initialize a dictionary to store counts for numbers 1 to 5000
    counts = {i: 0 for i in range(1, 5001)}

    # Open the file for reading
    with open(filename, 'r') as file:
        for line in file:
            # Split the line into words
            words = line.split()

            # Loop through the words and update the counts
            for word in words:
                try:
                    num = int(word)
                    if 1 <= num <= 5000:
                        counts[num] += 1
                except ValueError:
                    pass

    return counts

def write_table_to_file(counts, output_filename):
    with open(output_filename, 'w') as file:
        file.write("Gain of GC\n")
        for num, count in counts.items():
            file.write(f"{num}\t{count}\n")

if __name__ == "__main__":
    input_filename = "Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_GainofGC.txt"  # Replace with your input file name
    output_filename = "Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_GainofGC_counts.txt"   # Replace with your desired output file name
    number_counts = count_numbers(input_filename)
    write_table_to_file(number_counts, output_filename)

def write_table_to_file(counts, output_filename):
    with open(output_filename, 'w') as file:
        file.write("Loss of GC\n")
        for num, count in counts.items():
            file.write(f"{num}\t{count}\n")

if __name__ == "__main__":
    input_filename = "Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_LossofGC.txt"  # Replace with your input file name
    output_filename = "Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_LossofGC_counts.txt"   # Replace with your desired output file name
    number_counts = count_numbers(input_filename)
    write_table_to_file(number_counts, output_filename)



# Read data from File A
with open('Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_GainofGC_counts.txt', 'r') as file_a:
    data_a = [line.strip().split('\t') for line in file_a]

# Read data from File B
with open('Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_LossofGC_counts.txt', 'r') as file_b:
    data_b = [line.strip().split('\t') for line in file_b]

# Skip headers in both data lists
header_a = data_a[0]
header_b = data_b[0]
data_a = data_a[1:]
data_b = data_b[1:]



merged_data = []

for row_a, row_b in zip(data_a, data_b):
    idx_a, count_a = int(row_a[0]), int(row_a[1])
    idx_b, count_b = int(row_b[0]), int(row_b[1])

    difference = count_a - count_b
    merged_data.append([idx_a, count_a, count_b, difference])

# Write merged data to output file
with open('Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_NetGC.txt', 'w') as output_file:
    # Write header
    header_out = ["Nucleotide position", "Gain of GC", "Loss of GC"]
    output_file.write('\t'.join(header_out) + '\tDifference\n')
    
    # Write data rows
    for row in merged_data:
        output_file.write('\t'.join(map(str, row)) + '\n')

# Read the contents of the first file into a list of lines
with open('Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_NetGC.txt', 'r') as f1:
    lines1 = f1.readlines()

# Read the contents of the second file into a list of lines
with open('Homologous Canid Genes Aligned 60.0% Cutoff_Coor_Count.txt', 'r') as f2:
    lines2 = f2.readlines()

# Create a dictionary to store nucleotide positions and corresponding values from the second file
position_value_dict = {}
for line in lines2:
    parts = line.strip().split('\t')
    if len(parts) == 2:
        position, value = parts
        position_value_dict[position] = value

# Open a new output file for writing
with open('Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_NetGC_Coor_Count.txt', 'w') as output_file:
    # Write the header line (assuming the header is present in the first file)
    output_file.write(lines1[0].strip() + '\tCoor_Count\tRelative_Position\n')

    # IteFoxe over the lines in the first file (skipping the header)
    for line in lines1[1:]:
        parts = line.strip().split('\t')
        position = int(parts[0])

        # Get the corresponding value from the second file's dictionary
        appended_value = position_value_dict.get(str(position), '')

        # Calculate the relative position based on your description
        relative_position = (position - 1) - 2500

        # Write the original line from the first file, the appended value, and the relative position to the output file
        output_file.write(line.strip() + f'\t{appended_value}\t{relative_position}\n')

# Read the data from the input file
with open('Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_NetGC_Coor_Count.txt', 'r') as f:
    lines = f.readlines()

# Extract and write the header to the output file
header = lines[0].strip() + '\tDivide\n'

# Process the data and calculate the new values
new_lines = [header]
for line in lines[1:]:
    columns = line.strip().split('\t')
    difference = int(columns[3])
    coor_count = int(columns[4])
    relative_position = int(columns[5])
    
    new_value = difference / coor_count
    new_relative_position = f"{relative_position}\t{new_value:.10f}"
    
    new_line = '\t'.join(columns[:5]) + '\t' + new_relative_position + '\n'
    new_lines.append(new_line)

# Write the new data to the output file
with open('Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_NetGC_Coor_Count_Divide.txt', 'w') as f:
    f.writelines(new_lines)

import pandas as pd

# Read input data from file
input_file = 'Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_NetGC_Coor_Count_Divide.txt'
output_file = 'Homologous Canid Genes All Alignment 60.0% Cutoff_mod Primate SNP Fox_NetGC_Coor_Count_Divide_Final.txt'
window_size = 100

# Read the input data with tab sepaFoxion
df = pd.read_csv(input_file, sep='\t')

# Calculate half of the window size
half_window = window_size // 2

# Calculate sliding window average centered around the middle
df['Sliding_Window_Avg'] = df['Divide'].rolling(window=window_size, center=True).mean()

# Calculate the average of the Sliding_Window_Avg column
average_sliding_window_avg = df['Sliding_Window_Avg'].mean()

# Add the average value to a new column
df['Average_Sliding_Window_Avg'] = average_sliding_window_avg

# Write the output data to a new tab-sepaFoxed text file
df.to_csv(output_file, sep='\t', index=False)
