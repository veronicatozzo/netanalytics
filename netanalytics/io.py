import warnings
import numpy as np
import pandas as pd

from scipy.sparse import coo_matrix

from netanalytics.utils import _check_axis


def get_adjacency_csv(file):
    data = pd.read_csv(file, index_col=0)
    return data.values

    
def _params_check(G, filename, labels, axis):
    n1, n2 = G.shape
    if n1 != n2:
        raise ValueError("The adjacency matrix of the graph must be a square "
                         "matrix. Found dimensions (%d, %d)" %(n1, n2))

    if(str(filename).split('.')[-1]!='csv'):
        warnings.warn("The filename is not in csv format. Changing it.")
        filename = "".join(str(filename).split('.')[:-1]+['csv'])
    n = len(labels)
    if n1 != n:
        raise ValueError("The number of labels of the nodes must be the same "
                         "as the number of nodes. Found %d labels for %d "
                         "nodes" %(n, n1))
    axis = _check_axis(axis)
    return filename, axis


def to_LEDA(G, filename, labels=None, keep_loops=False, directed=False, axis=0):
    """
    Write a graph G to a file in LEDA formatself.

    Parameters
    ----------

    G: array-like, shape=(n,n)
        The graph.
    filename: string,
        The filename in which to write the network.
    labels: list of string, optional default=None
        The labels of the nodes to insert in the file.
    keep_loops: boolean, optional default=False
        If True self loops are stored.
    directed: boolean, default=False
        If the graph has to be considered directed or undirected.
    axis: int, default=0
        If undirected is not used. If directed the directions in which to
        consider the edges (rows->columns axis=0, columns-> rows axis=1)
    """
    # TODO
    filename, axis = _params_check(G, filename, labels, axis)


def to_edge_list(G, filename, labels=None, keep_loops=False,
                 directed=False, axis=0):
    """
    Write a graph G to a file as an edge list.

    Parameters
    ----------

    G: array-like, shape=(n,n)
        The graph.
    filename: string,
        The filename in which to write the network.
    labels: list of string, optional default=None
        The labels of the nodes to insert in the file.
    keep_loops: boolean, optional default=False
        If True self loops are stored.
    directed: boolean, default=False
        If the graph has to be considered directed or undirected.
    axis: int, default=0
        If undirected is not used. If directed the directions in which to
        consider the edges (rows->columns axis=0, columns-> rows axis=1)
    """
    filename, axis = _params_check(G, filename, labels, axis)
    res = []
    if not directed:
        A = A.copy()
        indices = np.triu_indices(A.shape[0])
        A[indices] = 0
    if labels==None:
        labels = np.arange(1, G.shape[0]+1)
    for i, l in enumerate(labels):
        non_zero = np.where(A[i,:]!= 0)[0]
        label_nz = label[non_zero]
        values_nz = A[i, non_zero]
        for lab, val in zip(label_nz, values_nz):
            res.append([l, lab, val])
    df = pd.DataFrame(res, columns=['Source', 'Target', 'Value'])
    df.to_csv(filename)


def from_edge_list(filename, sep=',', only_adjacency=True):
    """
    From a tabular file reads the graph.

    parameters
    ----------
    filename: string,
        The file from which the network is retrieved. It must be a tabular file.
    sep: string, optional
	The type of separator used in the file to read.

    only_adjacency: boolean, optional, default=True
        If False returns also the edge list as a list of indices.

    returns
    -------
    array-like:
        The adjacency matrix of the graph.
    list:
        The labels of the nodes.
    list:
        The list of edges as indices of the list of labels. Only if
        only_adjacency=False.
    """
    data = pd.read_table(filename,  sep=sep)
    nodes = data.iloc[:, 0].tolist() + data.iloc[:, 1].tolist()
    nodes_labels = sorted(list(set(nodes)))
    nodes = [(i,nodes_labels[i]) for i in range(len(nodes_labels))]
    for i in range(len(nodes)):
        data = data.replace(nodes[i][1], nodes[i][0])
    M = coo_matrix((data.iloc[:,2], (data.iloc[:,0],data.iloc[:,1])),
                    shape=(len(nodes), len(nodes))).todense()
    edge_list = []
    if not only_adjacency:
        for _from, _to in zip(data.iloc[:,0], data.iloc[:,1]):
            edge_list.append([_from, _to])
        return M, nodes_labels, edge_list
    return M, nodes_labels


def to_ORCA(nodes_list, edges_list, filename):
    with open(filename, 'w') as f:
        f.write(str(len(nodes_list)) + ' ' + str(len(edges_list)) + '\n')
        for e in edges_list:
            f.write(str(e[0]) + ' ' + str(e[1]) + '\n')
