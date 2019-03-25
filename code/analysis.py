import numpy as np
import pickle
from sklearn.linear_model import LogisticRegression, ElasticNet, SGDClassifier

def load_dataset_np(filename):
    return np.load(filename)

def load_dataset_pickle(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

def classification_error(y, yhat):
    return np.mean(y!=yhat)

def classification_accuracy(y, yhat):
    return 1 - np.mean(y!=yhat)


X = load_dataset_np('data_X_pa.npy')
dataset = load_dataset_pickle('data_y_pa.pkl')
y = dataset["y"]
metadata = np.array(dataset["metadata"])
genes = dataset["gene_names"]

modelGenes = load_dataset_pickle("bestSGD.pkl")
indices = []


bestLogReg = load_dataset_pickle("bestLogReg.pkl")
bestSGD = load_dataset_pickle("bestSGD.pkl")
best = []
for gene in bestLogReg:
    if gene in bestSGD:
        best.append(gene)
for g in best: print('|' + g.strip("'") + '|')


index = 0
for gene in genes:
    if gene in modelGenes:
        indices.append(index)
    index += 1

# Reduce X
X = np.array([[example[idx] for idx in indices] for example in X])
split = int(len(y)*7/10)
n = 100
acc = 0.0
nnz = 0.0

# format for markdown
# for g in modelGenes: print('|' + g.strip("'") + '|')

for i in range(n):
    # Shuffle data
    randomize = np.arange(len(y))
    np.random.shuffle(randomize)
    X, y, metadata = X[randomize], y[randomize], metadata[randomize]

    Xtrain, ytrain = X[:split], y[:split]
    Xtest, ytest = X[split:], y[split:]
            
    model = LogisticRegression()
    model.fit(Xtrain, ytrain)

    acc += classification_accuracy(model.predict(Xtest), ytest)
    nnz += (model.coef_ != 0).sum()

print("avg validation accuracy: %.3f with %d non-zeros" % (acc/n, nnz/n))    
