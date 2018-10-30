import warnings

import numpy as np

from itertools import product

from netanalytics.utils import _check_axis


def local_clustering_coefficient(A, i, directed=False, axis=0):
    axis = _check_axis(axis)
    A = A.copy()
    A -= np.diag(np.diag(A))
    neighbours = np.nonzero(A[:, i])[0] if axis else np.nonzero(A[i, :])[0]
    k = neighbours.shape[0]
    if directed:
        how_many_links = 0
        for ixs in product(neighbours, neighbours):
            print(ixs)
            i,j = ixs
            if axis:
                how_many_links += int(A[j,i] != 0)
            else:
                how_many_links += int(A[i,j] != 0)
    else:
        sub_net = A[neighbours, :]
        sub_net = sub_net[:, neighbours]
        how_many_links = np.sum((sub_net!=0).astype(int))
    if k == 1 or how_many_links == 0:
        return 0
    return how_many_links/float(k*(k-1))


def graph_clustering_coefficient(A, directed=False, axis=0):
    cc = [local_clustering_coefficient(A, i, directed, axis)
            for i in range(A.shape[0])]
    return np.mean(cc)


#def betweenness_centrality():


#def closeness_centrality():
