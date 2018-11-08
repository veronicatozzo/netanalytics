import warnings
import numpy as np


def _check_axis(axis):
    if axis<0 or axis>1:
        warnings.warn("Found a value for axis different than 0 and 1. Setting"
                       "default value 0.")
    axis=0
    return axis


def _normalize_degree_distribution(degrees, count):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        S = count/degrees
        S[np.where(degrees==0)[0]] = count[np.where(degrees==0)[0]]
        S[np.where(count==0)[0]] = 0
        T = np.sum(S)
        N = S/T
        return N


def _remove_loops(edges_list):
    aux = edges_list.copy()
    for e in edges_list:
        if e[0]==e[1]:
            aux.remove(e)
    return aux


def _from_edges_list_to_adjacency(edges_list, values=None):
    edges = np.array(edges_list)
    if values == None:
        values = np.ones(edges_list.shape[0])
    M = coo_matrix((values, (edges[:,0],edges[:,1])),
                    shape=(len(edges), len(edges))).todense()
    return M


def _from_adjacency_to_edges_list(A, undirected=False):
    edges_list = []
    values = []
    A = A.copy()
    if undirected:
        ix = np.triu_indices(A.shape[0])
        A[ix] = 0
    for i in range(A.shape[0]):
        for j in np.where(A[i,:]!=0)[0]:
            edges_list.append([i,j])
            values.append(A[i,j])
    return edges_list, values
