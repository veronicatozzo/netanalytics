import numpy as np

def _check_axis(axis):
    if axis<0 or axis>1:
        warnings.warn("Found a value for axis different than 0 and 1. Setting"
                       "default value 0.")
    axis=0
    return axis

def _normalize_degree_distribution(degrees, count):
    S = count/degrees
    S[np.where(degrees==0)[0]] = 0
    T = np.sum(S)
    N = S/T
    return N


def _remove_loops(edges_list):
    aux = edges_list.copy()
    for e in edges_list:
        if e[0]==e[1]:
            aux.remove(e)
    return aux
