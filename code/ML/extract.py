import numpy as np
import pickle

"""
Anthony Jakob, 11.02.2019

Extract data from time series
"""

def write_dataset(X, y, metadata, gene_names):
    # X is too large to be stored with the other values
    np.save('data_X.npy', X)
    # Store other values
    data = { 'y': y, 'metadata': metadata, 'gene_names': gene_names }
    with open('data_y.pkl', 'wb') as f:
        pickle.dump(data, f)

filename = "GSE18123-GPL570_series_matrix.txt"
samples = []
meta = []
genes = []
expression_values = []
reading = False

with open(filename, encoding="utf-8") as f:
    for line in f:
        if not reading:
            if line.startswith("!Sample_characteristics_ch1"):
                components = [comp.strip('"') for comp in line.split("\t")]
                meta.append(components[1:])
            
            reading = line.startswith("!series_matrix_table_begin")
            
        elif not line.startswith("!series_matrix_table_end"):
            # split elements per line and strip quotes
            components = [comp.strip('"') for comp in line.split("\t")]

            if not samples:
                samples = components[1:]
            else:
                genes.append(components[0])
                expression_values.append(np.array(components[1:]).astype(float))

X = np.transpose(expression_values)
meta.append(samples)
metadata = np.transpose(meta)
y = [np.sum([characteristic.startswith("diagnosis:") and not characteristic.startswith("diagnosis: CONTROL") for characteristic in data]) for data in metadata]
y = np.array(y)

print("Number of samples (y): %i" % len(y))
print("Number of samples: %i" % len(samples))
print("Number of genes: %i" % len(genes))

write_dataset(X, y, metadata, genes)
