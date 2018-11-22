import numpy as np

from sklearn.metrics.pairwise import pairwise_kernels

from scipy.stats import pearsonr

from netanalytics.thresholding import thresholding

def pearson_correlation(x, y):
    c, _ = pearsonr(x, y)
    return c

def mutual_rank_graph(X, threshold='1'):
    """
    Parameters
    ----------
    X: array-like, shape=(n, d)
        The input data matrix with n samples and d variables.
    threshold: string, optional default='1'
        The type of threshold to use for eliminate the weak connections.
    Returns
    -------
    array-like:
        The adjacency matrix of the mutual rank graph.
    """

    pc = pairwise_kernels(X.T, metric=pearson_correlation)
    pc -= np.diag(np.diag(pc))
    mr_r = np.zeros_like(pc)
    mr_c = np.zeros_like(pc)
    ix = np.arange(1, pc.shape[0]+1)
    for i in range(pc.shape[0]):
        mr_r[i,:] = ix[np.argsort(pc[i,:])]
        mr_c[:,i] = ix[np.argsort(pc[:,i])]
    mr = (mr_r + mr_c)/2
    mr /= np.max(mr)
    res = thresholding(mr, mode='1')
    return res + res.T - np.diag(np.diag(res))
