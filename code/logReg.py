import numpy as np
import pickle
from sklearn.linear_model import LogisticRegression, ElasticNet, SGDClassifier
from sklearn import svm
from sklearn.preprocessing import scale

import multiprocessing as mp

"""
Anthony Jakob, 11.02.2019

Perform analysis
"""

def load_dataset_np(filename):
    return np.load(filename)

def load_dataset_pickle(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

def classification_error(y, yhat):
    return np.mean(y!=yhat)

def classification_accuracy(y, yhat):
    return 1 - np.mean(y!=yhat)

def cross_validate_async(X, y, C):
    n = len(y)
    n_blocks = 4
    len_blocks = int(n * 1.0 / n_blocks)
    # print("cross_validate_async with C = %.3f" % C)

    # do crossvalidation leaving out one example at the time
    cumulative_error = 0
    for i in range(0, n, len_blocks):
        indices_val = np.array([j >= i and j <= i + len_blocks - 1 for j in range(n)])
        indices_train = np.invert(indices_val)

        model = LogisticRegression(penalty='l1', C=C, max_iter=100)
        model.fit(X[indices_train, :], y[indices_train])
        cumulative_error += classification_error(model.predict(X[indices_val,:]), y[indices_val])

    cumulative_error /= 1.0 * n_blocks

    return { 'error': cumulative_error, 'C': C }

if __name__ == "__main__":  # always guard your multiprocessing code

    hasControl = True
    verbose = True

    X = load_dataset_np('data_X.npy')
    dataset = load_dataset_pickle('data_y.pkl')
    y = dataset["y"]
    metadata = dataset["metadata"]
    genes = dataset["gene_names"]  

    """
    # caution: gene names might not match
    if hasControl:
        X_ctrl = load_dataset_np('data_X_control.npy')
        dataset_ctrl = load_dataset_pickle('data_y_control.pkl')
        genes_ctrl = dataset_ctrl["gene_names"] 

        print(genes[100:110])
        print(genes_ctrl[1000:1100])
        
        X_new = []
        genes_new = []
        for i in range(X_ctrl.shape[1]):
            try:
                idx = genes.index(genes_ctrl[i])
                X_new.append(X[:,idx])
                genes_new.append(genes[idx])
            except:
                # skip
                continue

        print(np.array(X_new).transpose().shape)
        print(X_ctrl.shape)
        X = np.append(np.array(X_new).transpose(), X_ctrl, axis=0)
        y = y + dataset_ctrl["y"]
        genes = genes_new
    """

    # Remove lowly expressed genes
    X_new = []
    genes_new = []
    for i in range(X.shape[1]):
        if all(gene >= 150 for gene in X[:,i]):
            X_new.append(X[:,i])
            genes_new.append(genes[i])

    X = np.array(X_new).transpose()
    genes = genes_new

    # Center and scale data
    X = scale(X)

    if verbose: print("Shuffling")
    randomize = np.arange(len(y))
    np.random.seed(12)
    np.random.shuffle(randomize)
    X, y, metadata = X[randomize], y[randomize], metadata[randomize]

    cores = max(mp.cpu_count() - 1, 1)  # ensure at least one process
    if verbose: print("Running on %i cores" % cores)
    C = 0.01
    with mp.Pool(processes=cores) as pool:
        Cs = [pool.apply_async(cross_validate_async, (X, y, 10 ** c,)) for c in range(-3, 3)]
        results = [result.get() for result in Cs]

        best_err = np.inf
        for result in results:
            if result['error'] < best_err:
                best_err = result['error']
                C = result['C']


    if verbose: print("Best C: %.3f" % C)
    

    split = int(len(y)*7/10)
    Xtrain, ytrain = X[:split], y[:split]
    Xtest, ytest = X[split:], y[split:]

    nModels = 100
    best_acc = 0
    best_acc_nnz = 10000
    best_acc_genes = []
    best_nnz = 10000
    best_nnz_acc = 0
    best_nnz_genes = []
    for i in range(nModels):
        if verbose: print("Iteration %d" % i)
		
		# Shuffle data
        randomize = np.arange(len(y))
        # 12 np.random.seed(i)
        np.random.shuffle(randomize)
        X, y, metadata = X[randomize], y[randomize], metadata[randomize]

        Xtrain, ytrain = X[:split], y[:split]
        Xtest, ytest = X[split:], y[split:]
        
        model = LogisticRegression(penalty='l1', C=C, max_iter=100)
        model.fit(Xtrain, ytrain)

        acc = classification_accuracy(model.predict(Xtest), ytest)
        nnz = (model.coef_ != 0).sum()
        if acc > best_acc:
            best_acc = acc
            best_acc_nnz = nnz
            best_acc_genes = (model.coef_ != 0)[0]
        if best_nnz > nnz:
            best_nnz = nnz
            best_nnz_acc = acc
            best_nnz_genes = (model.coef_ != 0)[0]

    print("Best validation accuracy: %.3f with %d non-zeros" % (best_acc, best_acc_nnz))
    print("Best non-zeros: %d with validation accuracy %.3f" % (best_nnz, best_nnz_acc))
    
    if best_acc_nnz < 100:
        dt = np.array(genes)[best_acc_genes]
        with open('best_genes_logReg_' + str(best_acc_nnz) + '_' + str(best_acc) + '.pkl', 'wb') as f:
            pickle.dump(dt, f)
