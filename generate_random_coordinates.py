'''
The is an example where we input all protein coding gene coordinates for chromosome 1 for humans, 
'''

import random

# open the input file for reading
with open("chr1.txt", "r") as f:
    # create an empty set to store all values within the ranges
    values = set()
    # read each line in the file
    for line in f:
        # split the line into columns
        columns = line.strip().split()
        # extract the start and end values from columns 2 and 3
        start = int(columns[1])
        end = int(columns[2])
        # add all values within the range to the set
        values.update(range(start, end+1))

# create an empty list to store the new ranges
new_ranges = []

# generate random ranges until we have 10000 in total
while len(new_ranges) < 10000:
    # generate a random start value for the range
    start = random.randint(1, 248951422)
    end = start + 5000
    # check if any value within the range is already in the set of stored values
    if not any(i in values for i in range(start, end+1)):
        # check if the new range overlaps with any existing range
        if not any(start <= r_end and end >= r_start for r_start, r_end in new_ranges):
            # if not, add the range to the new list
            new_ranges.append((start, end))

# open a new file for writing the new ranges
with open("chr1_random.txt", "w") as f:
    # write each range to a new line in the file
    for start, end in new_ranges:
        f.write(f"chr1\t{start}\t{end}\n")
