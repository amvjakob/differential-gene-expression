import numpy as np
import pickle
import re

def write_dataset(X, y, metadata, gene_names):
    # X is too large to be stored with the other values
    np.save('data_X.npy', X)
    # Store other values
    data = { 'y': y, 'metadata': metadata, 'gene_names': gene_names }
    with open('data_y.pkl', 'wb') as f:
        pickle.dump(data, f)

filename = "../dataX.txt"
filenameMeta = "../dataY.txt"
samples = []
meta = []
genes = []
expression_values = []
y = []
reading = False

firstline = True
with open(filename, encoding="utf-8") as f:
    for line in f:
        # split at tabs
        vals = re.split(r'\t+', line.rstrip('\t'))
        if firstline:
            samples = [val.strip('"') for val in vals]
            firstline = False
        else:
            genes.append(vals[0].strip('"'))
            expression_values.append(np.array(vals[1:]).astype(float))

firstline = True
with open(filenameMeta, encoding="utf-8") as f:
    for line in f:
        if firstline:
            firstline = False
        else:
            vals = re.split(r'\t+', line.rstrip('\t'))
            vals = [val.strip('"') for val in vals]
            meta.append(vals[1:])
            y.append(np.sum([fact == "AUTISM" for fact in vals]))

X = np.transpose(expression_values)
y = np.array(y)
meta = np.array(meta)

print("Number of samples (y): %i" % len(y))
print("Number of samples: %i" % len(samples))
print("Number of genes: %i" % len(genes))

write_dataset(X, y, meta, genes)

