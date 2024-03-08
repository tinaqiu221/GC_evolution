'''
! Now we are ready to begin alignment !

****************************************************

Needleman Alignment.py: uses the Needleman Wunsch algorithm to perform 2 sets of pairwise alignment
	reference organism vs cousin
	reference organism vs outgroup
	does "end mapping" where the code scans where the reference sequence starts to resemble the other homologous sequence from both ends
		could also play around with cutoffs to increase the number of alignments
'''
def delta (x,y):        # Penalty for misalignment
    if x==y:
        return (0)
    else:
        return (-1)
        
def evaluate (x, y):    # Calculate how similar 2 sequences are with each other
    score = 0
    length = len (x)
    for a in range (length):
        if x[a] == y[a]:
            score = score + 1
    return (score/length)
        
def fwd_window(a, b, window, cutoff):   # Map when the alignment starts to look good
  window = int(window)
  cutoff = float(cutoff)
  alignment_score = 0
  start = 0
  for x in range(1000):
  # for x in range(len(a)-window+1):
    while alignment_score < cutoff and start <= (len(b)-window):
        alignment_score = evaluate(a[x:x+window], b[start:start+window])
        # seq_align_parts = seq_align.split ('\t')
        '''
        print (seq_align_parts[0])
        print (seq_align_parts[1])
        '''
        # alignment_score = evaluate(seq_align_parts[0], seq_align_parts[1])
        start = start + 1
    if alignment_score >= cutoff:
        '''
        print (a[x:x+window])
        print (b[start-1:start+window-1])
        print (alignment_score)
        '''
        # return (start)
        return (x+1, start)
    start = 0
  return (-1, -1)
  
#fwd_window ("HHHHHHHATGCACAGTC", "BBAATGCACAGTC", 5, 0.8)
        
def bwd_window(a, b, window, cutoff):   # Map where the alignment finishes being good
    window = int(window)
    cutoff = float(cutoff)
    alignment_score = 0
    start = len(b)
    for x in range (1000):
    #   for x in range(len(a)-window+1):
        while alignment_score < cutoff and start >= window:
            alignment_score = evaluate(a[(len(a)-window-x):len(a)-x], b[(start-window):start])
            # seq_align_parts = seq_align.split ('\t')
            '''
            print (seq_align_parts[0])
            print (seq_align_parts[1])
            '''
            # alignment_score = evaluate(seq_align_parts[0], seq_align_parts[1])
            start = start - 1
        if alignment_score >= cutoff:
            '''
            print (a[(len(a)-window-1):len(a)-1])
            print (b[(start-window+1):start+1])
            print (alignment_score)
            '''
            # return (start+1)
            return (len(a)-x, start+1)
        start = len(b)  
    return (-1, -1)
    
#bwd_window ("ATGCACAGTCHHHHHHA", "AATGCACAGTCBB", 5, 0.8)
        
# 2D Alignment with Gap Penalty
def alignment (x, y):       # Alignment via Needleman-Wunsch Algorithm
    string1 = ""
    string2 = ""
    n = len(y)
    m = len(x)
    OPT = [[0 for i in range(n+1)] for j in range (m+1)]
    #print (OPT)
    for i in range(1,m+1):
        OPT[i][0] = -(i)
    for j in range(1,n+1):
        OPT[0][j] = -(j)
    for i in range(1, m+1):
        for j in range(1,n+1):
            OPT[i][j] = max(OPT[i-1][j-1]+ delta(x[i-1],y[j-1]), OPT[i-1][j] - 1, OPT[i][j-1] - 1) # align, delete, insert     
    i = m
    j = n
    #print (i)
    #print (j)
    prevDel = "No"
    prevIns = "No"
    while i != 0 or j != 0:
        if prevIns == "No": 
            insert = OPT[i][j-1] - 1.1
        else:
            insert = OPT[i][j-1] - 1
        if prevDel == "No":
            delete = OPT[i-1][j] - 1.1
        else:
            delete = OPT[i-1][j] - 1
        align = OPT[i-1][j-1] + delta(x[i-1],y[j-1])
        best_choice = max(align, insert, delete)
        if best_choice == insert:
            string1 = "_" + string1
            string2 = y[j-1]  + string2
            j = j-1
            prevDel = "No"
            prevIns = "Yes"
        elif best_choice == delete:
            string1 = x[i-1] + string1
            string2 = "_" + string2
            i = i-1
            prevDel = "Yes"
            prevIns = "No"
        elif best_choice == align:
            string1 = x[i-1] + string1
            string2 = y[j-1] + string2
            i = i-1
            j = j-1
            prevDel = "No"
            prevIns = "No"
    return (string1 + '\t' + string2)
    
#alignment ("ATTTGACGACATGATCGACTTTACGGCA", "ATGACGACATGGGATCGCCTATACGGCAAA")
    
'''(Human, Chimp, Human, Gorilla)'''
def triple_Align (n):       # 2 sets of pair-wise alignment (annotated, cousin, outgroup)
    input_file = open (n+".txt", 'r')
    output_file = open(n+" Alignment.txt", 'w')
    output_file2 = open(n+" Alignment Error.txt", 'w')
    output_file3 = open(n+" No Alignment.txt", 'w')
    name = input_file.readline()
    count = 0
    while name != "*":
        name_Parts = name.split(',')
        gorilla = input_file.readline()
        chimp = input_file.readline()
        human = input_file.readline()
        start1 = fwd_window(human, chimp, 20, 0.9)        # May want to change the stringency
        end1 = bwd_window (human, chimp, 20, 0.9)
        start2 = fwd_window(human, gorilla, 20, 0.9)
        end2 = bwd_window(human, gorilla, 20, 0.9)
        '''
        if start1[0] == -1:
            print ("msg1")
            start1 = fwd_window(human, chimp, 40, 0.80)
        if end1[0] == -1:
            print ("msg2")
            end1 = bwd_window (human, chimp, 40, 0.80)
        if start2[0] == -1:
            print ("msg3")
            start2 = fwd_window(human, gorilla, 40, 0.80)
        if  end2[0] == -1:
            print ("msg4")
            end2 = bwd_window(human, gorilla, 40, 0.80)
        if end1[1] <= start1[1] and end1[1] != -1:
            print ("msg5")
            start1 = fwd_window(human, chimp, 40, 0.80)
            end1 = bwd_window (human, chimp, 40, 0.80)
        if end2[1] <= start2[1] and end2[1] != -1:
            print ("msg6")
            start2 = fwd_window(human, gorilla, 40, 0.80)
            end2 = bwd_window(human, gorilla, 40, 0.80)
        print (start1)
        print (end1)
        print (start2)
        print (end2)
        print (start1[0])
        print (end1[0])
        print (start2[0])
        print (end2[0])
        print (start1[1])
        print (end1[1])
        print (start2[1])
        print (end2[1])
        '''

        if start1[0] != -1 and end1[0] != -1 and start2[0] != -1 and end2[0] != -1:
            #if end1 >= (start1+4000) and end2 >= (start2+4000):
            if end1[1] - start1[1] > 20 and end2[1] - start2[1] > 20:
                Edit_Seq = alignment(human[start1[0]-1:end1[0]-1],chimp[start1[1]-1:end1[1]-1])
                Edit_Parts = Edit_Seq.split('\t')
                output_file.write (">" + name_Parts[2][0:-1] + "," + name_Parts[1] + "," + name_Parts[2][0:-1] + "," + name_Parts[0][1:] + '\r' + Edit_Parts[0] + '\r' + Edit_Parts[1] + '\r')
                Edit_Seq = alignment(human[start2[0]-1:end2[0]-1],gorilla[start2[1]-1:end2[1]-1])
                Edit_Parts = Edit_Seq.split('\t')
                output_file.write (Edit_Parts[0] + '\r' + Edit_Parts[1] + '\r')
            else:
                print ("Alignment error")
                output_file2.write (name)
                output_file2.write (gorilla)
                output_file2.write (chimp)
                output_file2.write (human)
        else:
            print ("No alignment")
            output_file3.write (name)
            output_file3.write (gorilla)
            output_file3.write (chimp)
            output_file3.write (human)
        count = count + 1
        print (count)
        name = input_file.readline()
    output_file.write ("*")
    input_file.close()
    output_file.close()
    
triple_Align ("Homologous Canid Genes2")

'''
###########################################
Below codes are a quality control step where alignments are compared and only returns alignments that are x% identical
	also back calculates how much of the ends were cropped off due to end mapping in the previous step
	returns the cropped window
'''
def evaluate (x, y):	# Calculate how similar 2 sequences are with each other
	score = 0
	length = len (x)
	for a in range (length):
		if x[a] == y[a]:
			score = score + 1
	return (score/length)

def noIndel (n):
	count = 0
	string = ""
	for a in range (len(n)):
		if n[a] != "_":
			string = string + n[a]
	return (string)

def nth_nonInDel_Character (n, x):
	count = 1
	position = 0
	while count < x:
		if n[position] == "_":
			position = position + 1
		else:
			count = count + 1
			position = position + 1
	return (position)

def Alignment_Cutoff (n, unaligned, cutoff, code):
	input_file = open(n+".txt", 'r')
	unaligned_file = open(unaligned+".txt", 'r')
	output_file = open(n+" " + str(cutoff*100) + "% Cutoff.txt", 'w')
	array = ["" for i in range(1000000000)]
	name = unaligned_file.readline()
	while name != "*":
		ENSG_loc = name.find(code)
		ENSG = int(name[ENSG_loc+7:ENSG_loc+18])
		outgroup = unaligned_file.readline()
		cousin = unaligned_file.readline()
		annotated = unaligned_file.readline()
		array[ENSG] = annotated
		name = unaligned_file.readline()
	unaligned_file.close()
	name = input_file.readline()
	while name != "*":
		ENSG_loc = name.find(code)
		ENSG = int(name[ENSG_loc+7:ENSG_loc+18])
		fwd_crop_record = 0
		bwd_crop_record = 5000
		human1 = input_file.readline()
		human1 = human1.strip()
		chimp = input_file.readline()
		chimp = chimp.strip()
		human2 = input_file.readline()
		human2 = human2.strip()
		gorilla = input_file.readline()
		gorilla = gorilla.strip()
		check1 = evaluate (human1, chimp)
		check2 = evaluate (human2, gorilla)
		if check1 > cutoff and check2 > cutoff:
			#print (name[0:-1])
			H1 = noIndel(human1)
			H2 = noIndel(human2)
			H1 = H1.strip()
			H2 = H2.strip()
			#print (H1)
			#print (H2)
			H1_start = array[ENSG].find (H1)
			H2_start = array[ENSG].find (H2)
			#print (H1_start)
			#print (H2_start)
			if H1_start > H2_start:
				HG_start = nth_nonInDel_Character (human2, H1_start-H2_start)
				human2 = human2[HG_start+1:]
				gorilla = gorilla[HG_start+1:]
				fwd_crop_record = H1_start
			elif H2_start > H1_start:
				HC_start = nth_nonInDel_Character (human1, H2_start-H1_start)
				human1 = human1[HC_start+1:]
				chimp = chimp[HC_start+1:]
				fwd_crop_record = H2_start
			else:
				fwd_crop_record = H1_start
			H1_end = H1_start + len(H1)
			H2_end = H2_start + len(H2)
			#print (H1_end)
			#print (H2_end)
			#print (human1)
			#print (human2)
			if H1_end < H2_end:
				HG_end = nth_nonInDel_Character (human2, (H1_end-fwd_crop_record))
				human2 = human2[0:HG_end+1]
				gorilla = gorilla[0:HG_end+1]
				bwd_crop_record = H1_end
			elif H2_end < H1_end:
				HC_end = nth_nonInDel_Character (human1, (H2_end-fwd_crop_record))
				human1 = human1[0:HC_end+1]
				chimp = chimp[0:HC_end+1]
				bwd_crop_record = H2_end
			else:
				bwd_crop_record = H1_end
			output_file.write (name)
			output_file.write (human1 + '\r')
			output_file.write (chimp + '\r')
			output_file.write (human2 + '\r')
			output_file.write (gorilla + '\r')
			output_file.write (str(fwd_crop_record+1) + ":" + str(bwd_crop_record) + '\r')
		name = input_file.readline()
	output_file.write("*")
	input_file.close()
	output_file.close()

Alignment_Cutoff ("Homologous Canid Genes Aligned", "Homologous Canid Genes", 0.6, "ENSCAFG")
Alignment_Cutoff ("Homologous Canid Genes Aligned", "Homologous Canid Genes", 0.7, "ENSCAFG")
#Alignment_Cutoff ("sample", "Homologous Canid Genes", 0.6, "ENSCAFG")
