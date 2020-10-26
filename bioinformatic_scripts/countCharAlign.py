#
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required = True, help = 'Input file')
parser.add_argument('-o', '--output', required = True, help = 'Output file')
parser.add_argument('-c', '--character', required = True, help = 'Special char to count')

args = parser.parse_args()

lines = []
results = []

#Parse alignment file to memory
with open(args.input, 'r') as f:
    for line in f:
        if not line.startswith(">"):
            lines.append(line)
#combine into 2d-list
for i in lines:
    results.append(list(i))

#convert to pandas dataframe and tidy column headers
y = np.array([np.array(xi) for xi in results], dtype='object')

t = np.ndarray.transpose(np.array([list(range(1,len(results[1])+1)),[0]*len(results[1])]))

o = pd.DataFrame(t,columns = ['position','count'])
o.set_index('position', inplace=True)

#Do search for character
for i in y:
    temp = np.where(i == args.character)
    o.loc[o.index[temp],'count'] += 1

#Tidy (remove carriage returns) and write to file
o.drop(o.tail(1).index,inplace=True)

pd.DataFrame(o).to_csv(args.output)
