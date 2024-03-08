'''
######################################
Get the sequences of homologous genes 
download a triplet set of homologous genes in related species from Ensembl > make sure to click in this order reference, cousin, outgroup
However some of the entries in the list don't have homologous counterpart
Use the below codes to get rid of these

'''

def Index (n):
    i_file = open(n+".txt", 'r')
    o_file = open(n+" Ordered.txt", 'w')
    input_line = i_file.readline()
    o_file.write(input_line)
    input_line = i_file.readline()
    while input_line != "*":
        species1 = input_line.find ("ENSCAFG")
        species2 = input_line.find ("ENSUAMG")
        species3 = input_line.find ("ENSVVUG")
        if species1 != -1 and species2 != -1 and species3 != -1:
            o_file.write(input_line)
        input_line = i_file.readline()
    o_file.write ("*")
    i_file.close()
    o_file.close()

Index ("Index")

'''
#################################
Download GTF file from Ensembl
Remember to add star at the end of file and resave as txt
Use below codes to get rid of extrac unnecessary info to reduce file size
'''

def EnsemblFormatting (n, g, t):
    i_file = open(n+".fa", 'r')
    o_file = open(n+" Formatted.fa", 'w')
    input_line = i_file.readline()
    input_line = i_file.readline()
    input_line = i_file.readline()
    input_line = i_file.readline()
    input_line = i_file.readline()
    input_line = i_file.readline()
    geneLocation = 0
    transcriptLocation = 0
    while input_line != "*":
        input_parts = input_line.split ('\t')
        if input_parts[2]=="gene" or input_parts[2]=="transcript" or input_parts[2]=="exon":
        # if input_parts[2]=="gene" or input_parts[2]=="transcript" or input_parts[2]=="exon" or input_parts[2]=="CDS" or input_parts[2]=="five_prime_utr" or input_parts[2]=="start_codon":
            geneLocation = input_parts[8].find(g)
            transcriptLocation = input_parts[8].find(t)
            o_file.write (input_parts[0] + '\t' + input_parts[1] + '\t' + input_parts[2] + '\t' + input_parts[3] + '\t' + input_parts[4] + '\t' + input_parts[5] + '\t' + input_parts[6] + '\t' + input_parts[7] + '\t' + input_parts[8][geneLocation:(geneLocation+len(g)+11)] + '\t' + input_parts[8][transcriptLocation:(transcriptLocation+len(t)+11)] + '\r')
        input_line = i_file.readline()
    o_file.write ("*")
    i_file.close()
    o_file.close()

EnsemblFormatting ("Canis_lupus_familiaris.ROS_Cfam_1.0.108.chr", "ENSCAFG", "ENSCAFT")
EnsemblFormatting ("Vulpes_vulpes.VulVul2.2.108", "ENSVVUG", "ENSVVUT")
EnsemblFormatting ("Ursus_americanus.ASM334442v1.108", "ENSUAMG", "ENSUAMT")

'''
##################################
Use the below codes to get protein-coding genes only
'''
def PCGs_Only (GTF, PCG, species, g):
    GTF_file = open(GTF+".fa", 'r')
    PCG_file = open(PCG+".fa", 'r')
    output_file = open(species + "_Protein_Coding_Genes_GTF.txt", 'w')
    dimension = 3000000000
    array = [0]*dimension
    save = ""
    input_line = PCG_file.readline()
    while input_line != "*":
        if input_line[0] == ">":
            ENSG_find = input_line.find(g)
            ID = int(input_line[ENSG_find+len(g):(ENSG_find+len(g)+11)])
            array[ID] = 1
        input_line = PCG_file.readline()
    PCG_file.close()
    input_line = GTF_file.readline()
    while input_line != "*":
        input_parts = input_line.split ('\t')
        if input_parts[2] == "gene":
            ENSG_find = input_line.find(g)
            ID = int(input_line[ENSG_find+len(g):(ENSG_find+len(g)+11)])
            if array[ID] == 1:
                output_file.write (input_line)
                input_line = GTF_file.readline()
                input_parts = input_line.split ('\t')
                while input_parts[2] != "gene":
                    output_file.write (input_line)
                    input_line = GTF_file.readline()
                    if input_line == "*":
                        break
                    input_parts = input_line.split ('\t')
            else:
                input_line = GTF_file.readline()
                input_parts = input_line.split ('\t')
                while input_parts[2] != "gene":
                    input_line = GTF_file.readline()
                    if input_line == "*":
                        break
                    input_parts = input_line.split ('\t')
    output_file.write ("*")
    GTF_file.close()
    output_file.close()
    
PCGs_Only ("Canis_lupus_familiaris.ROS_Cfam_1.0.108.chr Formatted", "Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all", "Dog", "ENSCAFG")
PCGs_Only ("Ursus_americanus.ASM334442v1.108 Formatted", "Ursus_americanus.ASM334442v1.pep.all", "Bear", "ENSUAMG")
PCGs_Only ("Vulpes_vulpes.VulVul2.2.108 Formatted", "Vulpes_vulpes.VulVul2.2.pep.all", "Fox", "ENSVVUG")

'''
#####################################3
Some genes have multiple transcripts
this can be due to alternative splicing or alternative transcription start sites (TSS)
We define the best TSS ourselves as:
	-the most commonly used TSS
	-if TSS use is tied select one at random
Use the below codes to go through the GTF file and select the first alternatively spliced transcript that uses the TSS that we just selected
'''
import random

def mostFrequent(arr, n): # randomize when tied
	maxcount = 0
	element_having_max_freq = 0
	for i in range(0, n):
		count = 0
		for j in range(0, n):
			if(arr[i] == arr[j]):
				count += 1
		if(count > maxcount):
			maxcount = count
			element_having_max_freq = arr[i]
		elif count == maxcount:
			randNum = random.randint (1,2)
			if randNum == 1:
				element_having_max_freq = arr[i]
	return element_having_max_freq

def Best_TSS (n, code):
	input_file = open(n+".txt", 'r')
	o_file = open(n+"_Best_TSS.txt", 'w')
	strand = ""
	array = [""]
	ID = ""
	dimension = 1000000000
	most_Common = [""]*dimension
	input_line = input_file.readline()
	while input_line != "*":
		input_parts = input_line.split ('\t')
		if input_parts[2] == "gene":
			ID = int (input_line[input_line.find(code)+len(code):input_line.find(code)+len(code)+11])
			if input_parts[6] == "+":
				strand = "+"
			elif input_parts[6] == "-":
				strand = "-"
			input_line = input_file.readline()
			input_parts = input_line.split ('\t')
			while input_parts[2] != "gene":
				if input_parts[2] == "transcript":
					if strand == "+":
						array.append(input_parts[3])
					elif strand == "-":
						array.append(input_parts[4])
				input_line = input_file.readline()
				if input_line == "*":
					break
				input_parts = input_line.split ('\t')
			array.remove("")
			x = len(array)
			most_Common[ID] = mostFrequent (array, x)
			# print (str(array) + " " + most_Common[ID])
			strand = ""
			array = [""]
			ID = ""
	input_file.close()
	input_file = open(n+".txt", 'r')
	checker = "no"
	input_line = input_file.readline()
	while input_line != "*":
		input_parts = input_line.split ('\t')
		if input_parts[2] == "gene":
			o_file.write (input_line)
			ID = int (input_line[input_line.find(code)+len(code):input_line.find(code)+len(code)+11])
			if input_parts[6] == "+":
				strand = "+"
			elif input_parts[6] == "-":
				strand = "-"
			input_line = input_file.readline()
			input_parts = input_line.split ('\t')
			while input_parts[2] != "gene":
				if checker == "no":
					if input_parts[2] == "transcript":
						if strand == "+":
							if input_parts[3] == most_Common[ID]:
								o_file.write (input_line)
								input_line = input_file.readline()
								input_parts = input_line.split ('\t')
								while input_parts[2] == "exon":
									o_file.write (input_line)
									input_line = input_file.readline()
									input_parts = input_line.split ('\t')
									if input_line == "*":
										break
								checker = "yes"
							else:
								input_line = input_file.readline()
								input_parts = input_line.split ('\t')
								while input_parts[2] == "exon":
									input_line = input_file.readline()
									input_parts = input_line.split ('\t')
									if input_line == "*":
										break
						elif strand == "-":
							if input_parts[4] == most_Common[ID]:
								o_file.write (input_line)
								input_line = input_file.readline()
								input_parts = input_line.split ('\t')
								while input_parts[2] == "exon":
									o_file.write (input_line)
									input_line = input_file.readline()
									input_parts = input_line.split ('\t')
									if input_line == "*":
										break
								checker = "yes"
							else:
								input_line = input_file.readline()
								input_parts = input_line.split ('\t')
								while input_parts[2] == "exon":
									input_line = input_file.readline()
									input_parts = input_line.split ('\t')
									if input_line == "*":
										break
				else:
					input_line = input_file.readline()
					input_parts = input_line.split ('\t')
				if input_line == "*":
					break
			checker = "no"
	o_file.write ("*")
	input_file.close()
	o_file.close()

Best_TSS ("Dog_Protein_Coding_Genes_GTF", "ENSCAFG")
Best_TSS ("Bear_Protein_Coding_Genes_GTF", "ENSUAMG")
Best_TSS ("Fox_Protein_Coding_Genes_GTF", "ENSVVUG")

'''
#############################################
Identify the genomic coordinates of sequences we want to download
	we want sequences surrounding the TSS
	the best studied species will be used as reference 
		download a smaller track for this
	for the other two species download a larger track
		the best annotated species will scan the other two larger track until alignment looks good and locks on
	simply add +/- x Kb surrounding the TSS
'''

def KbDownload (n, x, code):
    input_file = open(n+".txt", 'r')
    o_file = open(n+"_" + str(int(x/1000)) + "kb_Seq.gff", 'w')
    #p_file.write ("track name=Plus" + '\r')
    #m_file.write ("track name=Minus" + '\r')
    input_line = input_file.readline()
    while input_line != "*":
        input_parts = input_line.split ('\t')
        if input_parts[2] == "transcript": 
            ENSG_find = input_line.find (code)
            if input_parts [6] == "+":
                o_file.write (input_parts[0] + '\t' + str(int(input_parts[3])-int(x/2)) + '\t' + str(int(input_parts[3])+int(x/2)) + '\t' + input_line[ENSG_find:ENSG_find+len(code)+11] + '\t' + "1" + '\t' + "+" + '\n')
            elif input_parts [6] == "-":
                o_file.write (input_parts[0] + '\t' + str(int(input_parts[4])-int(x/2)) + '\t' + str(int(input_parts[4])+int(x/2)) + '\t' + input_line[ENSG_find:ENSG_find+len(code)+11] + '\t' + "1" + '\t' + "-" + '\n')
            while input_parts[2] != "gene":
                input_line = input_file.readline()
                if input_line == "*":
                    break
                input_parts = input_line.split ('\t')
        else:
            input_line = input_file.readline()
    input_file.close()
    o_file.close()
    
KbDownload ("Bear_Protein_Coding_Genes_GTF_Best_TSS", 10000, "ENSUAMG")                 # Chimpanzee +/- 5kb
KbDownload ("Fox_Protein_Coding_Genes_GTF_Best_TSS", 10000, "ENSVVUG")                     # Gorilla +/- 5kb
KbDownload ("Dog_Protein_Coding_Genes_GTF_Best_TSS", 5000, "ENSCAFG")                          # Human +/- 2.5kb

'''
############################################
Go to Ensembl to download the genome (toplevel)
write the lines below into the command line (don't run it as a python code)
Now we have sequences surrounding the TSS
'''
bedtools getfasta -fi Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa -bed Dog_Protein_Coding_Genes_GTF_Best_TSS_5kb_Seq.gff -fo Dog_TSS.txt -s -name

bedtools getfasta -fi Ursus_americanus.ASM334442v1.dna.toplevel.fa -bed Bear_Protein_Coding_Genes_GTF_Best_TSS_10kb_Seq.gff -fo Bear_TSS.txt -s -name

bedtools getfasta -fi Vulpes_vulpes.VulVul2.2.dna.toplevel.fa -bed Fox_Protein_Coding_Genes_GTF_Best_TSS_10kb_Seq.gff -fo Fox_TSS.txt -s -name

'''
###########################################
Combine the homolgous gene sequences in sets of three 
	ordered like this: gene name
			   outgroup sequences
			   cousin sequences
			   reference sequences
'''
def Index (human, chimp, gorilla, annotated, cousin, outgroup, index, family):
    human_file = open(human+".txt", 'r')
    chimp_file = open(chimp+".txt", 'r')
    gorilla_file = open(gorilla+".txt", 'r')
    index_file = open(index+".txt", 'r') # index file tells you which genes are homologous
    output_file = open("Homologous " + family + " Genes.txt", 'w')
    human_array = [""]*1000000000
    chimp_array = [""]*1000000000
    gorilla_array = [""]*1000000000
    ID = ""
    seq = ""
    count = 0
    # store trio sequences into appropriate arrays
    human_line = human_file.readline()
    while human_line != "*":
        if human_line[0] == ">":
            geneLocation = human_line.find(annotated)
            ID = human_line[(geneLocation+len(annotated)):(geneLocation+len(annotated)+11)]
            human_line = human_file.readline()
            while human_line[0] != ">":
                if human_line == "*":
                    break
                seq = seq + human_line[0:-1]
                human_line = human_file.readline()
            if int(ID) <= 1000000000 and human_array[int(ID)] == "":
                human_array[int(ID)] = (seq)
            seq = ""
    human_file.close()
    chimp_line = chimp_file.readline()
    while chimp_line != "*":
        if chimp_line[0] == ">":
            geneLocation = chimp_line.find(cousin)
            ID = chimp_line[(geneLocation+len(cousin)):(geneLocation+len(cousin)+11)]
            chimp_line = chimp_file.readline()
            while chimp_line[0] != ">":
                if chimp_line == "*":
                    break
                seq = seq + chimp_line[0:-1]
                chimp_line = chimp_file.readline()
            if int(ID) <= 1000000000 and chimp_array[int(ID)] == "":
                chimp_array[int(ID)] = (seq)
            seq = ""
    chimp_file.close()
    gorilla_line = gorilla_file.readline()
    while gorilla_line != "*":
        if gorilla_line[0] == ">":
            geneLocation = gorilla_line.find(outgroup)
            ID = gorilla_line[(geneLocation+len(outgroup)):(geneLocation+len(outgroup)+11)]
            gorilla_line = gorilla_file.readline()
            while gorilla_line[0] != ">":
                if gorilla_line == "*":
                    break
                seq = seq + gorilla_line[0:-1]
                gorilla_line = gorilla_file.readline()
            if int(ID) <= 1000000000 and gorilla_array[int(ID)] == "":
                gorilla_array[int(ID)] = (seq)
            seq = ""
    gorilla_file.close()
    index_line = index_file.readline()
    index_line = index_file.readline()
    while index_line != "*":
        index_parts = index_line.split ('\t')
        if len(index_parts[0]) > 2 and len(index_parts[1]) > 2 and len(index_parts[2]) > 2:
            human_Location = index_line.find(annotated)
            human_ID = int(index_line[(human_Location+len(annotated)):(human_Location+len(annotated)+11)])
            chimp_Location = index_line.find(cousin)
            chimp_ID = int(index_line[(chimp_Location+len(cousin)):(chimp_Location+len(cousin)+11)])
            gorilla_Location = index_line.find(outgroup)
            gorilla_ID = int(index_line[(gorilla_Location+len(outgroup)):(gorilla_Location+len(outgroup)+11)])
            if human_array[human_ID] != "" and chimp_array[chimp_ID] != "" and gorilla_array[gorilla_ID] != "":
                output_file.write (">" + index_parts[2][0:-1] + "," + index_parts[1] + "," + index_parts[0] + '\r')
                output_file.write (gorilla_array[gorilla_ID] + '\r')
                output_file.write (chimp_array[chimp_ID] + '\r')
                output_file.write (human_array[human_ID] + '\r')
        index_line = index_file.readline()
        count = count + 1
        print (count)
    output_file.write ("*")
    index_file.close()
    output_file.close()
    
Index ("Dog_TSS", "Fox_TSS", "Bear_TSS", "ENSCAFG", "ENSVVUG", "ENSUAMG", "Index Ordered", "Canid")
