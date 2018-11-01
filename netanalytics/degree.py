import warnings

import numpy as np
from collections import Counter

from netanalytics.utils import _normalize_degree_distribution


def laplacian_matrix(A):
    A = A.copy()
    A[np.where(A!=0)] = 1
    D = np.diag(np.sum(A, axis=1))
    return D - A


def number_of_edges(A):
    """
    Parameters
    ------

    A: array_like, shape=(n,n)
        Adjacency matrix of the graph.

    Returns
    -------
    int :
        The number of edges.
    """
    aux = A.copy()
    aux -= np.diag(np.diag(aux))
    return np.sum((aux!=0).astype(int))


def node_degree(A, pos, axis=0):
    """
    Parameters
    ------

    A: array_like, shape=(n,n)
        Adjacency matrix of the graph.
    pos: int
        Number of node.
    axis: int, default=0
        The axis to consider to compute the degree. To use in case of directed
        networks.

    Returns
    -------
    float :
        The weighted degree of the node. If A is binary is the same as the
        unweighted degree.
    int :
        The unweighted degree of the node.
    """
    if axis==0:
        edges = A[pos, :].copy()
    elif axis==1:
        edges = A[:, pos].copy()
    else:
        raise ValueError("Not allowed axis value.")
    if pos < 0:
        raise ValueError("Negative node position.")
    if pos >= A.shape[axis]:
        raise ValueError("Position too big, given position %d, number of "
                         "nodes %d" %(pos, A.shape[axis]))
    edges = np.reshape(edges, A.shape[0])
    edges[0, pos] = 0
    return np.sum(edges), np.sum((edges!=0).astype(int))


def degree_distribution(A, axis=0):
    """
    Parameters
    ------

    A: array_like, shape=(n,n)
        Adjacency matrix of the graph.
    axis: int, default=0
        The axis to consider to compute the degree. To use in case of directed
        networks.

    Returns
    -------
    array-like: shape=(n,)
        The vector containing for each node of the network its degree.
    dict :
        A dictionary that has as keys the degree and as values the occurence of
        that degree in the network.
    """

    degrees = np.zeros(A.shape[axis])
    for i in range(A.shape[axis]):
        degrees[i] = node_degree(A, i, axis)[1]
    return degrees, Counter(degrees)





def degree_distribition_distance(G1, G2):
    _, d_G1 = degree_distribution(G1)
    _, d_G2 = degree_distribution(G2)
    degrees_G1 = list(d_G1.keys())
    degrees_G2 = list(d_G2.keys())

    degrees = np.array(sorted(list(set(degrees_G1).union(degrees_G2))))
    count_G1 = [d_G1.get(d, 0) for d in degrees]
    count_G2 = [d_G2.get(d, 0) for d in degrees]
    N_G1 = _normalize_degree_distribution(degrees, np.array(count_G1))
    N_G2 = _normalize_degree_distribution(degrees, np.array(count_G2))
    distance = np.sqrt(np.sum(np.square(N_G1 - N_G2)))/np.sqrt(2)
    return distance


#def degree_distribution_distance_kolmogorov_smirnov(G1, G2):
