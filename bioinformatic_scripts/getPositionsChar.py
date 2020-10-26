#
import pandas as pd
import numpy as np
import argparse
#import csv

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required = True, help = 'Input file')
parser.add_argument('-o', '--output', required = True, help = 'Output file')
parser.add_argument('-c', '--character', required = True, help = 'Special char to count')

args = parser.parse_args()

headers = []
lines = []
results = []
iupac = ['A','C','G','T','U','R','Y','S','W','K','M','B','D','H','V','N','-']

#Parse alignment file to memory
with open(args.input, 'r') as f:
    for line in f:
        if line.startswith(">"):
            headers.append(line.split(" ")[0])
        if not line.startswith(">") and not line.startswith('\n'):
            lines.append(line.replace(' ',''))


#combine into 2d-list
for i in lines:
    results.append(list(i))

#convert to pandas dataframe and tidy column headers
y = np.array([np.array(xi) for xi in results], dtype='object')

#define new matrix for counts with first column just the range of number of sequences
#and the other columns populated with 0'
matLen = len(results)
t = np.zeros((matLen,2))

#define pandas dataframe
o = pd.DataFrame(t, columns = ['start','end'])
o['headers'] = headers
o.set_index('headers',inplace=True)

#Do search
counter = 0
for i in y:
    temp = np.in1d(i, iupac)
    temp = temp*1
    o.loc[headers[counter],'start'] = (temp!=0).argmax(0) + 1
    o.loc[headers[counter],'end'] = np.max(np.nonzero(temp)) + 1
    counter+=1

# Write to file
pd.DataFrame(o).to_csv(args.output)
#write a stripped alignment file
# with open(args.output + '.csv', 'w', newline = '\n') as f:
#     writer = csv.writer(f)
#     writer.writerows(lines[1:5])
