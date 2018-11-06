import numpy as np
import pickle as pkl

from joblib import Parallel, delayed
import multiprocessing

from os import listdir
from os.path import isfile, join
from netanalytics.subnetworks import common_subgraph
from netanalytics.io import from_edge_list
from netanalytics.degree import number_of_edges

path = "/cs/research/bioinf/bionet1/Coexpression_Study/BCStages-SubtypesNetworks/"
folders = ['HER2Networks', 'LuminalANetworks', 'LuminalBNetworks', 'TripleNegativeNetworks',
           'Stage1Networks', 'Stage2Networks', 'Stage3Networks', 'Stage4Networks']

glob_files = []
glob_labels = []
for fol in folders:
    files = [f for f in listdir(path+fol) if
             (isfile(join(path+fol, f)) and f.endswith('txt'))]
    labels = []
    for f in files:
        if str(f).startswith("Adap"):
            labels.append(fol+'adaptive_lasso')
        elif str(f).startswith('Lasso'):
            labels.append(fol+'lasso')
        elif str(f).startswith('aracnea'):
            labels.append(fol+'aracnea')
        elif str(f).startswith('aracnem'):
            labels.append(fol+'aracnem')
        elif str(f).startswith('c3'):
            labels.append(fol+'c3net')
        elif str(f).startswith('gene'):
            labels.append(fol+'genenet')
        elif str(f).startswith('clr'):
            labels.append(fol+'clr')
        elif str(f).startswith('genie'):
            labels.append(fol+'genie3')
        elif str(f).startswith('mrnetb'):
            labels.append(fol+'mrnetb')
        elif str(f).startswith('mrnet') and str(f)[5]!='b':
            labels.append(fol+'mrnet')
        elif str(f).startswith('wgcna'):
            labels.append(fol+'wgcna')
    glob_files += [join(path+fol, f) for f in files]
    glob_labels += labels


def count(i):
    res = []
    A, nodesA, _ = from_edge_list(glob_files[i], sep=' ', only_adjacency=False)

    for j in range(i+1, len(glob_files)):
        B, nodesB, _ = from_edge_list(glob_files[i], sep=' ', only_adjacency=False)
        common, nodes, n_edges = common_subgraph(A, B, nodesA, nodesB)
        res.append([glob_labels[i]+'-'+glob_labels[j], number_of_edges(A),
                    number_of_edges(B), common, nodes, n_edges])
    return res


#TODO: adjust with parameters inputing
if __name__ == "__main__":
    results = Parallel(n_jobs=4)(delayed(count)(i) for i in range(len(glob_files)-1))
    with open("results_common_edges2.pkl", 'wb') as f:
        pkl.dump(results, f)
