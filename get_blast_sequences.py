'''
# convert coordinates to sequnces in UCSC:
https://genome.ucsc.edu/cgi-bin/hgCustom
1. Head to the Custom track upload page: https://genome.ucsc.edu/cgi-bin/hgCustom, and select your organism of interest (Mouse mm10) 
2. Use the upload file button or paste your BED input file directly into the text box, the click submit. 
3. On the resulting "Manage Custom Tracks" page, select "Table Browser" from the "view in" dropdown and click go. 
Have to chunk into 1000 lines at a time
4. The mm10 assembly and your uploaded custom track should be pre-selected, so you should be able to just change the "output format" dropdown to "sequence", and click "get output". Optionally, if you would like the resulting output in a file you can enter a file name into the "output file" text box and the next steps will download to a file. 
5. On the resulting page select any extra formatting options you would like and the click "get sequence", which will result in FASTA format output, for example:
>mm10_ct_UserTrack_3545_0 range=chrX:150001-200000 5'pad=0 3'pad=0 strand=+ repeatMasking=none NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

Then need to download and install local blast https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
'''

# make fasta into blast database (Hamster example below)
makeblastdb -in Hamster_genome.fasta -dbtype nucl

# perform local blast for the 3 specoes (Rat example below)

blastn -query Mouse_Recombination_Hotspots.fasta -db ratgenome.fasta -num_threads 32 -evalue 1e-5 -outfmt "6 qseqid sseqid qstart qend sstart send" -max_target_seqs 1 -out rat_blast.txt

#merge mouse, rat and hasmter blast results together
def merge (Cousin, Outgroup, Output):
    file_a = open (Cousin, 'r')
    file_b = open (Outgroup, 'r')
    output_file = open (Output, 'w')
    file_b_dict = {}
    for line_b in file_b:
        line_b = line_b.strip().split()
        file_b_dict[line_b[0]] = line_b[1:]
    for line_a in file_a:
        line_a = line_a.strip().split()
        if line_a[0] in file_b_dict:
            output_line = line_a[:]
            output_line.extend(file_b_dict[line_a[0]][:])
            output_file.write('\t'.join(output_line) + '\n')
        else:
            pass
    file_a.close()
    file_b.close()
    output_file.close()

merge ('chimp_blast.txt', 'gorilla_blast.txt', 'primate_blast.txt')


# now with trio blast aligned sequences, we will perform needleman alignment again to standardize our alignment methods. See codes in find_substitutions.py