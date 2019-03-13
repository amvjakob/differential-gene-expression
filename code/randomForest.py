import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
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

        model = SGDClassifier(loss='hinge', penalty='elasticnet', alpha=C, l1_ratio=0.3, max_iter=100)
        model.fit(X[indices_train, :], y[indices_train])
        cumulative_error += classification_error(model.predict(X[indices_val,:]), y[indices_val])

    cumulative_error /= 1.0 * n_blocks

    return { 'error': cumulative_error, 'C': C }

if __name__ == "__main__":  # always guard your multiprocessing code

    X = load_dataset_np('data_X.npy')
    dataset = load_dataset_pickle('data_y.pkl')
    y, metadata, gene_names = dataset["y"], dataset["metadata"], dataset["gene_names"]
    verbose = True

    # Remove lowly expressed genes
    X_new = []
    for i in range(X.shape[1]):
        if all(gene >= 150 for gene in X[:,i]):
            X_new.append(X[:,i])

    X = np.array(X_new).transpose()

    # Center and scale data
    X = scale(X)

    if verbose: print("Shuffling... ", end='')
    randomize = np.arange(len(y))
    np.random.seed(12)
    np.random.shuffle(randomize)
    X, y, metadata = X[randomize], y[randomize], metadata[randomize]
    if verbose: print("Done")

    """
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
    """

    split = 120
    Xtrain, ytrain = X[:split], y[:split]
    Xtest, ytest = X[split:], y[split:]

    nModels = 100
    best_acc = 0
    best_acc_nnz = 10000
    best_nnz = 10000
    best_nnz_acc = 0
    for i in range(nModels):
        if verbose: print("Iteration %d" % i)

        model = RandomForestClassifier(n_estimators=100, max_depth=20)
        model.fit(Xtrain, ytrain)

        acc = classification_accuracy(model.predict(Xtest), ytest)
        nnz = model.n_features_ #(model.coef_ != 0).sum()
        if acc > best_acc:
            best_acc = acc
            best_acc_nnz = nnz
        if best_nnz > nnz:
            best_nnz = nnz
            best_nnz_acc = acc

    print("Best validation accuracy: %.3f with %d non-zeros" % (best_acc, best_acc_nnz))
    print("Best non-zeros: %d with validation accuracy %.3f" % (best_nnz, best_nnz_acc))
