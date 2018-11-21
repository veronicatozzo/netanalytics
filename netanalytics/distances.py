import warnings

import numpy as np
import pandas as pd

from scipy.linalg import eigh

from netanalytics.graph import Graph
from netanalytics.degree import laplacian_matrix, degree_distribition_distance
from netanalytics.graphlets import GDD_agreement, GCD
from netanalytics.subnetworks import common_subgraph
from netanalytics.utils import jaccard_index


def jaccard_similarity(G1, G2):
    return jaccard_index(G1.edges, G2.edges)

def spectral_distance(A, B):
    La = laplacian_matrix(A)
    Lb = laplacian_matrix(B)
    ev_a = list(eigh(La, eigvals_only=True))
    ev_b = list(eigh(Lb, eigvals_only=True))

    _sum = 0
    for i in range(max(len(ev_a), len(ev_b))):
        if i >= len(ev_a):
            _sum += ev_b[-i]**2
        elif i>= len(ev_b):
            _sum += ev_a[-i]**2
        else:
            _sum += (ev_a[-i] - ev_b[-i])**2
    return np.sqrt(_sum)


def all_distances(G1, G2, verbose=0, label1=None, label2=None):
    """
    # TODO
    """
    if label1 is None:
        label1 = '1'
    if label2 is None:
        label2 = '2'
    if not isinstance(G1, Graph):
        raise ValueError("The graphs must be instances of the class Graph")

    if not isinstance(G2, Graph):
        raise ValueError("The graphs must be instances of the class Graph")

    if not G1.is_fitted:
        warning.warn("The graph analysis was not fit, doing it now.")
        G1.fit()
    if not G2.is_fitted:
        warning.warn("The graph analysis was not fit, doing it now.")
        G2.fit()

    spectral = spectral_distance(G1.adjacency, G2.adjacency)
    ddd = degree_distribition_distance(G1.adjacency, G2.adjacency)
    GDDAa, GDDAgeo= GDD_agreement(G1.GDV, G2.GDV)
    gcd73, gdc11 = GCD(G1.GDV, G2.GDV)
    _, _, n_edges = common_subgraph(G1.adjacency, G2.adjacency,
                                    G1.nodes, G2.nodes)
    if verbose:
        print("Computed sitances between graphs:")
        print("Spectral distance=%.4f" %spectral)
        print("Degree distribution distance= %.4f"%ddd)
        print("GDDA=%.4f related distance=%.4f"%(GDDAa, 1-GDDAa))
        print("GCD=%.4f"%gcd11)
    df = pd.DataFrame([[G1.clustering_coefficient, G2.clustering_coefficient,
                       spectral, ddd, GDDAa, GDDAgeo, gdc11, gcd73, n_edges]],
                      columns=['cl_coeff_1', 'cl_coeff_2', 'SD', 'DDD',
                      'GDDA (arithmetic)',
                      'GDDA (geometric)', 'GCD-11', 'GDC-73', 'No commmon edges'],
                      index=[label1+'-'+label2])
    return df
