import sys
import os
import warnings

import numpy as np
import pandas as pd

from scipy.stats import spearmanr
from collections import Counter

from netanalytics.io import to_ORCA
from netanalytics.utils import _normalize_degree_distribution


def graphlet_degree_vectors(nodes_list, edges_list, graphlet_size=5):
    if not (graphlet_size == 4 or graphlet_size == 5):
        raise ValueError("The maximum graphlet size must be either 4 or 5.")

    to_ORCA(nodes_list, edges_list, '_orca.txt')

    cmd = './../_orca/orca_mac '+ str(graphlet_size)  + ' _orca.txt _temp.ndump2'
    os.system(cmd)
    res = []
    with open('_temp.ndump2', 'rb') as _from:
        i = 0
        for line in _from:
            res.append(str(line).split("'")[-2].split('\\')[0].split(' '))
            i +=1
    os.remove("_orca.txt")
    os.remove("_temp.ndump2")
    df = pd.DataFrame(res, index=nodes_list)
    return df


#def graphlet_frequency_count(A):


#def relative_graphlet_frequency_distance():


def _graphlet_degree_distribution(GDV):
    _max = np.max(np.max(GDV))
    degrees = []
    dicts = []
    for orbit in range(GDV.shape[1]):
        dicts.append(Counter(GDV.iloc[orbit,:]))
        degrees.append(list(dicts[-1].keys()))
   # print(dicts)
    total_degrees = set()
    for dg in degrees:
        total_degrees = total_degrees.union(set(dg))
    total_degrees = np.array(sorted(list(total_degrees)))
    Ns = []
    for orbit in range(GDV.shape[1]):
        aux = np.array([dicts[orbit].get(d, 0) for d in total_degrees])
        Ns.append(_normalize_degree_distribution(total_degrees.astype(int), aux))
    res = pd.DataFrame(np.array(Ns).T, index=total_degrees)
    res.sort_index(inplace=True)
    return res


def GDD(nodes_list, edges_list):
    """
    Graphlet Degree Distribution.

    It computes the Graphlet Degree Vector and then computes the distributions
    for all the 73 graphlets orbits.
    """
    counts = graphlet_degree_vectors(nodes_list, edges_list, 5)
    return _graphlet_degree_distribution(counts)


def _graphlet_distribution_distance(GDD_1, GDD_2):
    if GDD_1 is None or GDD_2 is None:
        warnings.warn("Empty graphlet degree vector")
        return np.nan
    if (GDD_1.shape[1] != 73 and GDD_2.shape[1] != 73 and
       GDD_1.shape[1] != 15 and GDD_2.shape[1] !=15):
        raise ValueError("The number of orbits must be either 73 pr 15, "
                         "found %d for the first graph and %d for the second"
                         %(GDD_1.shape[1], GDD_2.shape[1]))

    indices = list(set(GDD_1.index.values).union(set(GDD_2.index.values)))
    distance = np.zeros(GDD_1.shape[1])
    v1 = np.zeros(len(indices))
    v2 = v1.copy()
    for orbit in range(GDD_1.shape[1]):
        for i, ix in enumerate(indices):
            try:
                v1[i] = GDD_1.loc[ix].iloc[orbit]
            except KeyError:
                v1[i] = 0
            try:
                v2[i] = GDD_2.loc[ix].iloc[orbit]
            except KeyError:
                v2[i] = 0
        distance[orbit] = np.sqrt(np.sum(np.square(v1-v2)))
    return distance


def GDD_distance(GDV1, GDV2):
    GDD_1 = _graphlet_degree_distribution(GDV1)
    GDD_2 = _graphlet_degree_distribution(GDV2)
    return _graphlet_distribution_distance(GDD_1, GDD_2)


def GDD_agreement(GDV1, GDV2):
    """
    Graphlet Degree Distribution Agreement.

    This measures uses the Graphlet Degree Distribution to compare two networks.
    """
    distances = GDD_distance(GDV1, GDV2)
    agreements = 1 - distances
    arithmetic = np.mean(agreements)
    geometric = np.product(agreements)
    return arithmetic, geometric


def graphlet_correlation_matrix(GDV):
    if GDV.shape[1] == 73:
        GDV = GDV.values.astype(int)
        GDV = np.vstack((GDV, np.ones((1,GDV.shape[1]))))
        spearm_corr = spearmanr(GDV)
        GCM_73 = spearm_corr[0]
        to_consider = [0,1,2,4,5,6,7,8,9,10,11]
        GCM_11 = GCM_73[to_consider, ]
        GCM_11 = GCM_11[:,to_consider]
        return GCM_73, GCM_11
    else:
        GDV = GDV.values.astype(int)
        GDV = np.vstack((GDV, np.ones((1,GDV.shape[1]))))
        spearm_corr = spearmanr(GDV)
        GCM_11 = spearm_corr[0]
        return None, GCM_11



def _graphlet_correlation_distance(GCM1, GCM2):
    _sum = 0
    if GCM1 is None or GCM2 is None:
        warnings.warn("Empty correlation matrix")
        return 0
    if GCM1.shape != GCM2.shape:
        raise ValueError("bla bla bla")
    for i in range(GCM1.shape[0]):
        for j in range(i, GCM1.shape[0]):
            _sum += (GCM1[i,j] - GCM2[i,j])**2
    return np.sqrt(_sum)


def GCD(GDV1, GDV2):
    GCM1_73, GCM1_11 =  graphlet_correlation_matrix(GDV1)
    GCM2_73, GCM2_11 =  graphlet_correlation_matrix(GDV2)
    return (_graphlet_correlation_distance(GCM1_73, GCM2_73),
            _graphlet_correlation_distance(GCM1_11, GCM2_11))
