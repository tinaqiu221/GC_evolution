# weighted GC calculation into x bins without the window function

def nuc_count(n, nuc):
  """ take in sequence n, and determine GC percentage """
  Gs = n.count('G') + n.count('g')
  Cs = n.count('C') + n.count('c')
  As = n.count('A') + n.count('a')
  Ts = n.count('T') + n.count('t')
  CpGs = n.count ('CG') + n.count ("cg")
  if (Gs+Cs+As+Ts) > 0:
    if nuc == "GC":
        return (Gs+Cs)/(Gs+Cs+As+Ts)
    elif nuc == "G":
        return (Gs)/(Gs+Cs+As+Ts)
    elif nuc == "C":
        return (Cs)/(Gs+Cs+As+Ts)
    elif nuc == "A":
        return (As)/(Gs+Cs+As+Ts)
    elif nuc == "T":
        return (Ts)/(Gs+Cs+As+Ts)
    elif nuc == "CpG":
        return (CpGs)/(2*(Gs+Cs+As+Ts))
  else:
    return (0)

def weighted_Nucleotide_analysis_NoWindow (n, x, g, nuc):
  """ take in a file name, n, and produce an output file with nucleotide count for x bins"""
  input_file = open(n+".txt", 'r')
  output_file = open(n+" Weighted "+ nuc +"_"+ str(x)+"Bins.txt", 'w')
  gene_name = input_file.readline()
  seq_line = input_file.readline()
  seq = ''
  beginning = 0
  dimension = 10000000
  array = [[0 for i in range(x+2)] for j in range (dimension+1)]
  count_array = [0]*(dimension+1)
  while gene_name[0] != '*':
    while (seq_line[0] != '>') and (seq_line[0] != '*'):
      seq = seq+seq_line[0:-1]
      seq_line = input_file.readline()
    ENSG_find = gene_name.find(g)
    ID = int(gene_name[ENSG_find+len(g):(ENSG_find+len(g)+11)])
    if seq.count('N') == 0:
        bin_size = len(seq)/x
        if bin_size>2:
            if gene_name[-1] == '\n' or gene_name[-1] == '\r':
                array[ID][0] = gene_name[0:-1]
            else:
                array[ID][0] = gene_name	     
            for a in range(x):
                end = int(round(((a+1)*bin_size),0))
                next_beginning = end
                bin_nucleotide = nuc_count(seq[beginning:end], nuc)
                array[ID][a+1] = array[ID][a+1] + bin_nucleotide
                beginning = next_beginning
            array[ID][x+1] = len(seq)
            count_array[ID] = count_array[ID] + 1
            beginning = 0
    seq = ''
    gene_name = seq_line
    seq_line = input_file.readline()
  for a in range (dimension+1):
    if count_array[a] != 0:
        output_file.write (array[a][0] + '\t')
        for b in range (1, x+1):
            output_file.write (str(array[a][b]/count_array[a]) + '\t')
        output_file.write (str(array[a][x+1]) + '\r')
  output_file.write ("*")  
  input_file.close()
  output_file.close()

def Nucleotide_Analysis(n, bins, code):
    #weighted_Nucleotide_analysis_NoWindow (n, bins, code, "CpG")
    weighted_Nucleotide_analysis_NoWindow (n, bins, code, "GC")
    #weighted_Nucleotide_analysis_NoWindow (n, bins, code, "G")
    #weighted_Nucleotide_analysis_NoWindow (n, bins, code, "C")
    #weighted_Nucleotide_analysis_NoWindow (n, bins, code, "A")
    #weighted_Nucleotide_analysis_NoWindow (n, bins, code, "T")
     
Nucleotide_Analysis ("Human PCG TSS Exon1 Normalized Highest GC", 41, "ENSG")