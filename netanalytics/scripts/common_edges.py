import numpy as np
import pickle as pkl

from joblib import Parallel, delayed
import multiprocessing

from os import listdir
from os.path import isfile, join
from netanalytics.subnetworks import multi_common_subgraph
from netanalytics.io import from_edge_list
from netanalytics.degree import number_of_edges


def count(path):
    files = [join(path, f) for f in listdir(path) if isfile(join(path, f))]
    As = []
    nodes_lists= []
    n_edges_single = []
    for f in files:
        A, nodes, _= from_edge_list(f, sep=' ', only_adjacency=False)
        n_edges_single.append(number_of_edges(A))
        As.append(A)
        nodes_lists.append(nodes)
    common, nodes, n_edges = multi_common_subgraph(As, nodes_lists)
    return [n_edges_single, common, nodes, n_edges]


#TODO: adjust with parameters inputing
if __name__ == "__main__":
    path = "/cs/research/bioinf/bionet1/Coexpression_Study/BCStages-SubtypesNetworks/"
    folders = ['HER2Networks', 'LuminalANetworks', 'LuminalBNetworks', 'TripleNegativeNetworks',
           'Stage1Networks', 'Stage2Networks', 'Stage3Networks', 'Stage4Networks']

    results = Parallel(n_jobs=8)(delayed(count)(path+fol) for fol in folders)
    with open("results_common_edges.pkl", 'wb') as f:
        pkl.dump(results, f)
