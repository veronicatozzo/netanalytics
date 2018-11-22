
import numpy as np
import pickle as pkl
import pandas as pd

from joblib import Parallel, delayed
import multiprocessing

from os import listdir
from os.path import isfile, join
from netanalytics.subnetworks import common_subgraph
from netanalytics.io import from_edge_list, get_adjacency_csv
from netanalytics.degree import number_of_edges

path = "/cs/research/bioinf/bionet1/Veronica/NMF/network_inference/networks_coex"

files = [join(path, f) for f in listdir(path) if (isfile(join(path, f))
         and f.endswith('csv'))]

def count(i):
    A = get_adjacency_csv(files[i])
    name = files[i].split('/')[-1].split('.')[-2]
    no_edges = number_of_edges(A)
    return name, no_edges


#TODO: adjust with parameters inputing
if __name__ == "__main__":
    results = Parallel(n_jobs=4)(delayed(count)(i) for i in range(len(files)))
    with open("number_of_edges.pkl", 'wb') as f:
        pkl.dump(results, f)
