'''
##############################
Example of doing a permuation for mouse recombination hotspots vs. mouse intergenic gc-matched to recombination hotspot
'''

import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def split_and_analyze(input_file, analysis_file, num_iterations, value_A, value_B, observed_delta):
    net_frequency_differences = []
    for iteration in range(num_iterations):
        with open(input_file, 'r') as f:
            lines = f.readlines()

        random.shuffle(lines)

        lines_A = lines[:value_A]
        #print (lines_A)
        lines_B = lines[value_A:(value_A+value_B)]
        #print (lines_B)

        num_columns = len(lines_A[0].split('\t'))
        count_A = [0] * num_columns
        count_B = [0] * num_columns

        for line in lines_A:
            #print (line)
            values = line.split('\t')
            for i, value in enumerate(values):
                sub_values = value.split(',')  # Split by comma
                if len(sub_values) == 1:
                    sub_value = sub_values[0].strip()
                    if sub_value and sub_value.isdigit() and 450 <= int(sub_value) <= 550:
                        count_A[i] += 1
                else:
                    for sub_value in sub_values:
                        sub_value = sub_value.strip()
                        if sub_value and sub_value.isdigit() and 450 <= int(sub_value) <= 550:
                            count_A[i] += 1
        #print (count_A)

        for line in lines_B:
            #print (line)
            values = line.split('\t')
            for i, value in enumerate(values):
                sub_values = value.split(',')  # Split by comma
                if len(sub_values) == 1:
                    sub_value = sub_values[0].strip()
                    if sub_value and sub_value.isdigit() and 450 <= int(sub_value) <= 550:
                        count_B[i] += 1
                else:
                    for sub_value in sub_values:
                        sub_value = sub_value.strip()
                        if sub_value and sub_value.isdigit() and 450 <= int(sub_value) <= 550:
                            count_B[i] += 1
        #print (count_B)

        expression_A = (count_A[2] + count_A[3] + count_A[5] + count_A[6]) - (count_A[7] + count_A[8] + count_A[10] + count_A[11])
        #print (expression_A)
        expression_B = (count_B[2] + count_B[3] + count_B[5] + count_B[6]) - (count_B[7] + count_B[8] + count_B[10] + count_B[11])
        #print (expression_B)
        frequency_expression_A = expression_A / (value_A * 100)
        #print (frequency_expression_A)
        frequency_expression_B = expression_B / (value_B * 100)
        #print (frequency_expression_B)
        net_frequency_difference = frequency_expression_A - frequency_expression_B
        #print (net_frequency_difference)
        net_frequency_differences.append(net_frequency_difference)
    #output_file = open("Brandon Test.txt", 'w')
    #output_file.write (str(net_frequency_differences))
    #print (net_frequency_differences)
    with open(analysis_file, 'w') as f:
        f.write("\n".join(str(diff) for diff in net_frequency_differences) + "\n")

    # Plotting the normal distribution
    plt.hist(net_frequency_differences, bins=100, density=True)

    mu = np.mean(net_frequency_differences)
    sigma = np.std(net_frequency_differences)

    x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
    plt.plot(x, stats.norm.pdf(x, mu, sigma), color='r')
    plt.axvline(x=observed_delta, color='b', linestyle='--')

    plt.xlabel('Mouse Recombinatio vs. Mouse Intergenic GC-matched')
    plt.ylabel('Density')
    plt.title('Mouse Recombinatio vs. Mouse Intergenic GC-matched')
    plt.grid(True)

    # Calculate the area under the curve before x=0.00226
    area = stats.norm.cdf(observed_delta, mu, sigma)
    area_text = f"P-value: {1-area:.10f}"
    plt.text(0.5, 0.95, area_text, transform=plt.gca().transAxes, ha='center')

    plt.savefig('Mouse Recombinatio vs. Mouse Intergenic GC-matched.png')  # Save the figure as an image file

    plt.show()


input_file = "mouse_recombination_gcmatched.txt" #make this file manually by combining raw substitution data from the two groups we are testing
analysis_file = "mouse_recombination_gcmatched_frequency_analysis.txt"
num_iterations = 1000
value_A = 9248  # Customize this value
value_B = 20209  # Customize this value
observed_delta = -0.0007981

split_and_analyze(input_file, analysis_file, num_iterations, value_A, value_B, observed_delta)

import matplotlib.pyplot as plt

# Read data from a file
with open('mouse_recombination_gcmatched_frequency_analysis.txt', 'r') as file:
    data = [float(line) for line in file]

# Plot histogram
plt.hist(data, bins=100, edgecolor='black')

# Add labels and title
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram')

# Show the plot
plt.show()

'''
###########
Based on the above frequency_analysis permutation results, use the codes below to graph the distribution
then calculate the p-value based on the area under the curve before the actual delta
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import quad

# Your list of values obtained from the previous frequency_analysis output file
data = """
-4.82E-06
3.07E-06
6.36E-06
1.10E-06
-2.16E-07
1.76E-06
-9.43E-06
4.39E-06
-2.85E-06
5.05E-06
-6.80E-06
-2.16E-07
3.07E-06
1.76E-06
8.34E-06
2.42E-06
1.76E-06
-2.16E-07
-2.19E-06
-8.77E-06
-8.74E-07
3.07E-06
-2.85E-06
-8.74E-07
3.73E-06
4.42E-07
-8.11E-06
5.71E-06
1.10E-06
7.02E-06
-2.85E-06
-4.16E-06
7.68E-06
-6.14E-06
-4.16E-06
7.02E-06
-2.16E-07
-5.48E-06
-8.11E-06
-8.11E-06
3.73E-06
3.07E-06
-2.85E-06
1.76E-06
-2.19E-06
-6.14E-06
3.73E-06
1.76E-06
4.42E-07
-3.51E-06
-4.16E-06
-4.82E-06
-2.19E-06
-8.74E-07
-4.16E-06
-7.45E-06
4.42E-07
-4.82E-06
-1.53E-06
-2.16E-07
5.05E-06
2.42E-06
-2.16E-07
1.10E-06
4.39E-06
1.10E-06
-3.51E-06
4.39E-06
1.76E-06
1.10E-06
8.34E-06
-2.16E-07
-3.51E-06
1.76E-06
4.39E-06
-2.19E-06
-8.74E-07
9.65E-06
1.76E-06
-4.16E-06
4.42E-07
6.36E-06
3.73E-06
7.02E-06
-2.16E-07
9.00E-06
-7.45E-06
-6.14E-06
7.02E-06
-7.45E-06
4.39E-06
-9.43E-06
-2.16E-07
-1.14E-05
-8.74E-07
3.73E-06
-7.45E-06
9.00E-06
3.07E-06
1.10E-06
1.10E-06
-4.16E-06
3.07E-06
-6.14E-06
1.76E-06
1.76E-06
-8.74E-07
-1.53E-06
-1.53E-06
-4.82E-06
7.68E-06
5.71E-06
3.07E-06
-1.53E-06
2.42E-06
-3.51E-06
-5.48E-06
1.76E-06
7.02E-06
3.07E-06
-2.85E-06
-2.85E-06
-2.85E-06
1.23E-05
-8.74E-07
1.76E-06
2.42E-06
2.42E-06
5.05E-06
5.71E-06
3.07E-06
-1.53E-06
-8.74E-07
-4.16E-06
1.76E-06
-2.19E-06
1.76E-06
9.00E-06
2.42E-06
-5.48E-06
-2.85E-06
-4.82E-06
4.42E-07
-8.74E-07
3.73E-06
1.76E-06
-2.16E-07
7.02E-06
3.07E-06
1.76E-06
7.02E-06
-2.85E-06
3.07E-06
2.42E-06
-2.85E-06
4.42E-07
5.05E-06
-1.53E-06
4.42E-07
-4.82E-06
-1.53E-06
2.42E-06
-2.16E-07
-4.82E-06
-4.16E-06
-6.14E-06
-2.19E-06
6.36E-06
2.42E-06
-4.82E-06
-2.19E-06
4.39E-06
-2.85E-06
4.39E-06
1.10E-06
4.42E-07
1.10E-06
5.71E-06
-2.19E-06
-2.19E-06
-2.85E-06
-2.85E-06
1.10E-06
-9.43E-06
-4.16E-06
-3.51E-06
-2.19E-06
3.07E-06
-2.85E-06
7.68E-06
1.76E-06
-1.53E-06
-6.14E-06
1.76E-06
4.39E-06
-1.53E-06
4.42E-07
-2.19E-06
-1.53E-06
-2.16E-07
-1.53E-06
3.73E-06
1.10E-06
4.39E-06
-8.74E-07
-7.45E-06
-1.53E-06
3.73E-06
-6.80E-06
-2.19E-06
-4.82E-06
-1.53E-06
-3.51E-06
1.10E-06
1.10E-06
5.05E-06
3.73E-06
3.07E-06
1.10E-06
5.71E-06
3.73E-06
6.36E-06
3.07E-06
1.10E-06
-2.16E-07
9.65E-06
9.00E-06
1.10E-06
3.07E-06
1.76E-06
9.65E-06
-2.16E-07
2.42E-06
1.76E-06
3.73E-06
-2.85E-06
3.07E-06
-8.74E-07
3.73E-06
-4.16E-06
2.42E-06
-2.85E-06
4.42E-07
-4.16E-06
7.68E-06
-1.07E-05
7.68E-06
4.42E-07
-2.16E-07
-2.16E-07
5.05E-06
-1.53E-06
-1.53E-06
2.42E-06
-3.51E-06
-2.16E-07
4.39E-06
3.73E-06
-2.19E-06
5.05E-06
8.34E-06
1.76E-06
3.07E-06
-2.85E-06
-1.53E-06
2.42E-06
-8.74E-07
-4.82E-06
4.39E-06
2.42E-06
8.34E-06
4.42E-07
4.42E-07
-2.19E-06
1.10E-06
-2.16E-07
-8.74E-07
4.42E-07
-2.16E-07
4.39E-06
3.07E-06
4.42E-07
-2.19E-06
-4.82E-06
2.42E-06
-4.16E-06
1.10E-06
4.42E-07
-5.48E-06
-3.51E-06
-4.16E-06
5.05E-06
-2.85E-06
-2.16E-07
9.00E-06
3.07E-06
7.02E-06
-2.85E-06
1.76E-06
-2.85E-06
1.76E-06
-2.19E-06
-2.85E-06
-2.85E-06
6.36E-06
4.39E-06
-7.45E-06
3.07E-06
-8.11E-06
-2.85E-06
3.07E-06
5.05E-06
-1.07E-05
-4.16E-06
-2.16E-07
6.36E-06
-6.80E-06
3.73E-06
3.73E-06
-4.82E-06
-1.53E-06
-2.16E-07
-4.82E-06
7.02E-06
-2.19E-06
4.42E-07
-4.16E-06
-2.85E-06
-2.16E-07
2.42E-06
-5.48E-06
-8.74E-07
-2.85E-06
-8.74E-07
-2.16E-07
-4.82E-06
-4.82E-06
-1.53E-06
-3.51E-06
5.71E-06
5.05E-06
-7.45E-06
-2.19E-06
-5.48E-06
6.36E-06
4.42E-07
-4.82E-06
5.05E-06
-2.16E-07
-2.16E-07
-1.53E-06
-3.51E-06
-4.16E-06
4.42E-07
-8.74E-07
1.76E-06
-2.85E-06
-2.16E-07
1.10E-06
5.05E-06
4.39E-06
-2.85E-06
-4.82E-06
1.10E-06
4.42E-07
-4.16E-06
2.42E-06
5.71E-06
1.10E-06
4.39E-06
-2.85E-06
-3.51E-06
-8.74E-07
-3.51E-06
4.39E-06
-1.53E-06
3.07E-06
-2.85E-06
-1.53E-06
-2.85E-06
-3.51E-06
-4.82E-06
4.39E-06
5.71E-06
-4.82E-06
-7.45E-06
-5.48E-06
1.76E-06
4.42E-07
-8.74E-07
1.76E-06
-4.16E-06
-4.82E-06
-8.74E-07
-6.80E-06
-2.16E-07
-5.48E-06
5.05E-06
3.73E-06
1.10E-06
3.73E-06
-5.48E-06
5.05E-06
-8.74E-07
2.42E-06
-4.16E-06
5.71E-06
-2.85E-06
-4.82E-06
-8.77E-06
5.71E-06
4.39E-06
3.73E-06
-2.16E-07
3.07E-06
6.36E-06
2.42E-06
1.76E-06
5.71E-06
-7.45E-06
5.71E-06
4.42E-07
5.05E-06
4.42E-07
3.73E-06
-1.53E-06
5.05E-06
5.71E-06
-1.01E-05
-8.74E-07
3.07E-06
-2.16E-07
-2.85E-06
5.05E-06
3.73E-06
1.76E-06
-3.51E-06
-2.85E-06
7.02E-06
3.07E-06
-2.19E-06
-5.48E-06
1.76E-06
1.76E-06
-3.51E-06
1.10E-06
-2.16E-07
1.16E-05
-2.16E-07
-2.85E-06
4.42E-07
5.71E-06
1.76E-06
1.10E-06
-1.01E-05
6.36E-06
3.07E-06
1.76E-06
-6.80E-06
2.42E-06
-8.74E-07
3.07E-06
2.42E-06
1.76E-06
-2.19E-06
-1.53E-06
7.02E-06
1.76E-06
-3.51E-06
-8.77E-06
-2.16E-07
1.76E-06
-5.48E-06
7.68E-06
3.73E-06
1.10E-06
3.07E-06
1.10E-06
1.76E-06
5.71E-06
5.71E-06
-7.45E-06
4.39E-06
-4.16E-06
2.42E-06
1.10E-06
-6.14E-06
1.76E-06
-2.85E-06
-8.77E-06
-1.53E-06
-8.11E-06
-6.14E-06
3.73E-06
-1.53E-06
-2.85E-06
-8.74E-07
3.73E-06
4.42E-07
-3.51E-06
-2.85E-06
-2.16E-07
-8.74E-07
-2.19E-06
-2.85E-06
-2.85E-06
1.76E-06
-3.51E-06
7.68E-06
-1.53E-06
3.73E-06
-8.74E-07
-8.74E-07
-1.53E-06
1.03E-05
-6.80E-06
-9.43E-06
-8.77E-06
3.07E-06
6.36E-06
-1.07E-05
-5.48E-06
3.07E-06
3.73E-06
-4.16E-06
1.10E-06
1.76E-06
-2.85E-06
-8.74E-07
3.73E-06
1.76E-06
5.71E-06
3.73E-06
-3.51E-06
4.42E-07
3.07E-06
-2.16E-07
-2.85E-06
-4.16E-06
-2.85E-06
7.02E-06
1.76E-06
1.10E-06
-3.51E-06
2.42E-06
-8.74E-07
-1.53E-06
-1.53E-06
-5.48E-06
-3.51E-06
1.76E-06
7.02E-06
-5.48E-06
-2.16E-07
5.05E-06
-3.51E-06
-3.51E-06
4.39E-06
7.02E-06
-8.74E-07
1.10E-06
-2.16E-07
-2.16E-07
1.10E-06
-8.74E-07
-2.16E-07
3.07E-06
-2.19E-06
7.68E-06
4.42E-07
4.42E-07
4.42E-07
-3.51E-06
-3.51E-06
-3.51E-06
-3.51E-06
4.42E-07
1.76E-06
-2.16E-07
-2.85E-06
3.73E-06
2.42E-06
-2.19E-06
-7.45E-06
-6.14E-06
1.10E-06
-4.16E-06
5.05E-06
8.34E-06
5.05E-06
6.36E-06
-2.85E-06
-2.85E-06
3.07E-06
3.07E-06
-9.43E-06
5.05E-06
-4.82E-06
1.10E-06
-3.51E-06
1.76E-06
-2.16E-07
6.36E-06
8.34E-06
1.76E-06
3.73E-06
5.05E-06
1.16E-05
7.68E-06
-8.74E-07
-6.14E-06
2.42E-06
-3.51E-06
3.07E-06
2.42E-06
6.36E-06
4.42E-07
-1.40E-05
4.39E-06
-8.74E-07
-2.16E-07
2.42E-06
3.07E-06
1.76E-06
2.42E-06
4.42E-07
1.10E-06
-2.19E-06
2.42E-06
-1.53E-06
-2.85E-06
5.05E-06
-8.77E-06
6.36E-06
-8.74E-07
1.10E-06
-8.74E-07
1.76E-06
3.07E-06
5.05E-06
4.39E-06
1.10E-06
-1.53E-06
5.05E-06
-2.16E-07
-8.74E-07
-6.14E-06
4.42E-07
-8.74E-07
-2.85E-06
-3.51E-06
7.02E-06
-5.48E-06
-3.51E-06
7.02E-06
6.36E-06
5.05E-06
-2.85E-06
7.68E-06
3.73E-06
-3.51E-06
-8.74E-07
1.76E-06
-2.19E-06
2.42E-06
-2.85E-06
5.05E-06
5.71E-06
1.10E-06
1.10E-06
1.10E-06
2.42E-06
-1.53E-06
3.07E-06
-3.51E-06
-8.74E-07
5.05E-06
-2.85E-06
-8.74E-07
-2.19E-06
1.10E-06
5.71E-06
-8.11E-06
2.42E-06
-2.16E-07
-8.74E-07
-6.14E-06
3.73E-06
-3.51E-06
-6.14E-06
-3.51E-06
5.05E-06
-3.51E-06
-8.74E-07
4.39E-06
-2.19E-06
-4.82E-06
4.42E-07
-2.19E-06
2.42E-06
2.42E-06
9.00E-06
2.42E-06
-2.85E-06
-5.48E-06
-2.16E-07
-6.14E-06
4.39E-06
-4.16E-06
-5.48E-06
-1.01E-05
-8.74E-07
1.10E-06
9.00E-06
3.07E-06
3.07E-06
-8.74E-07
-1.53E-06
-3.51E-06
-7.45E-06
1.10E-06
2.42E-06
-3.51E-06
3.73E-06
1.10E-06
1.76E-06
4.42E-07
3.07E-06
-2.16E-07
-2.85E-06
5.05E-06
3.73E-06
1.10E-06
-1.53E-06
7.68E-06
4.39E-06
3.07E-06
-1.53E-06
-8.74E-07
6.36E-06
-7.45E-06
7.02E-06
3.07E-06
2.42E-06
8.34E-06
1.10E-05
1.76E-06
1.10E-06
-1.53E-06
-8.74E-07
-8.74E-07
3.07E-06
-4.82E-06
-2.19E-06
1.10E-06
-2.85E-06
-4.82E-06
3.73E-06
-1.53E-06
7.68E-06
-1.01E-05
-2.85E-06
-3.51E-06
1.76E-06
-2.19E-06
7.02E-06
5.71E-06
-8.74E-07
1.76E-06
3.73E-06
4.42E-07
-2.85E-06
-8.74E-07
-2.16E-07
5.05E-06
-3.51E-06
-3.51E-06
-2.85E-06
1.76E-06
2.42E-06
-1.53E-06
3.73E-06
4.39E-06
-8.74E-07
1.10E-06
1.10E-06
4.42E-07
-2.16E-07
4.42E-07
1.76E-06
3.73E-06
3.07E-06
3.73E-06
1.10E-05
3.07E-06
-2.16E-07
-3.51E-06
-4.82E-06
1.10E-06
2.42E-06
-2.19E-06
-2.16E-07
1.10E-06
2.42E-06
6.36E-06
1.10E-06
4.39E-06
2.42E-06
-4.16E-06
-2.85E-06
4.42E-07
-3.51E-06
1.76E-06
1.10E-06
-2.19E-06
1.10E-06
-4.16E-06
-2.16E-07
-2.16E-07
-1.14E-05
-3.51E-06
2.42E-06
-4.16E-06
3.73E-06
3.07E-06
-1.53E-06
-6.14E-06
3.73E-06
1.10E-06
-1.53E-06
1.10E-06
-2.19E-06
3.73E-06
3.73E-06
-3.51E-06
-4.82E-06
-7.45E-06
-4.16E-06
4.42E-07
-3.51E-06
1.10E-06
6.36E-06
-2.85E-06
-5.48E-06
-4.16E-06
-2.85E-06
-3.51E-06
3.73E-06
-2.16E-07
1.76E-06
-8.74E-07
-3.51E-06
1.76E-06
5.05E-06
2.42E-06
6.36E-06
-1.53E-06
-8.77E-06
4.42E-07
-2.85E-06
-1.07E-05
-1.53E-06
4.39E-06
-1.07E-05
6.36E-06
2.42E-06
1.10E-06
-4.82E-06
-8.74E-07
-2.85E-06
-2.16E-07
2.42E-06
-2.19E-06
-4.16E-06
2.42E-06
-6.80E-06
-2.85E-06
-4.16E-06
-8.74E-07
-4.16E-06
-2.19E-06
-2.85E-06
7.02E-06
-2.19E-06
-8.74E-07
2.42E-06
-2.19E-06
-2.19E-06
-4.82E-06
5.05E-06
-4.82E-06
4.42E-07
-1.07E-05
3.73E-06
-4.16E-06
-4.16E-06
4.39E-06
-2.19E-06
-2.85E-06
3.73E-06
5.05E-06
5.71E-06
-2.16E-07
-6.80E-06
5.05E-06
-4.82E-06
5.05E-06
-4.16E-06
-6.14E-06
-4.16E-06
-9.43E-06
-2.85E-06
3.07E-06
3.73E-06
1.10E-06
5.71E-06
5.05E-06
1.76E-06
-1.53E-06
-2.85E-06
1.10E-06
-6.80E-06
-5.48E-06
-2.19E-06
3.07E-06
-2.16E-07
-3.51E-06
-2.19E-06
-2.19E-06
2.42E-06
4.42E-07
1.10E-06
-2.85E-06
4.42E-07
-2.16E-07
-2.16E-07
-8.74E-07
-8.74E-07
4.42E-07
4.42E-07
1.76E-06
-8.74E-07
4.42E-07
1.10E-06
-3.51E-06
-4.16E-06
6.36E-06
-2.16E-07
4.42E-07
2.42E-06
-1.53E-06
1.76E-06
5.71E-06
-8.74E-07
-8.11E-06
-2.85E-06
-6.14E-06
-2.19E-06
3.73E-06
-9.43E-06
-2.16E-07
4.42E-07
-2.16E-07
-2.85E-06
-2.16E-07
-1.53E-06
-2.19E-06
3.07E-06
-2.19E-06
-8.74E-07
-2.85E-06
3.07E-06
-4.16E-06
1.76E-06
-3.51E-06
5.05E-06
5.71E-06
-2.85E-06
1.10E-06
2.42E-06
3.07E-06
-2.16E-07
-4.16E-06
1.76E-06
3.73E-06
2.42E-06
1.10E-06
-1.53E-06
-9.43E-06
2.42E-06
5.05E-06
-2.19E-06
4.42E-07
-8.74E-07
6.36E-06
4.42E-07
-5.48E-06
1.10E-05
-6.14E-06
2.42E-06
"""

# Convert the data to a NumPy array
data_array = np.fromstring(data, sep='\n')

# Calculate mean and standard deviation
mean_value = np.mean(data_array)
std_dev = np.std(data_array)

# Create a normal distribution with mean and standard deviation
x = np.linspace(min(data_array), max(data_array), 1000)
y = norm.pdf(x, mean_value, std_dev)

# Plot the normal distribution curve
plt.plot(x, y, label='Normal Distribution Curve')

# Set up the vertical line
vertical_line_x = -2.25212E-05  # Replace with the actuual calculated delta between the 2 groups
plt.axvline(x=vertical_line_x, color='r', linestyle='--', label=f'Vertical Line at x={vertical_line_x}')

# Calculate the area under the curve to the left of the vertical line
def integrand(x):
    return norm.pdf(x, mean_value, std_dev)

area, _ = quad(integrand, min(data_array), vertical_line_x)

# Display the result
print(f'Area under the curve to the left of x={vertical_line_x}: {area}')

# Show the plot
plt.legend()
plt.title('Normal Distribution Curve and Vertical Line')
plt.xlabel('X-axis')
plt.ylabel('Probability Density')
plt.show()
