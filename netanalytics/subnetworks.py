
import numpy as np


def common_subgraph(A, B, nodes_list_A, nodes_list_B, directed=False):
    """
    Retrievs the common sub graph of two networks.

    Given the adjacency matrices of two networks and their nodes lists it finds
    the nodes in common and consequentely the edges.

    Parameters
    ----------
    A: array-like, shape=(n1, n1)
        Adjacency matrix of the first network.
    B: array-like, shape=(n2, n2)
        Adjacency matrix of the second network.
    nodes_list_A: list, length=n1
        List of nodes identities for the network A.
    nodes_list_B: list, length=n2
        List of nodes identities for the network A.
    directed: boolean, optional defualt=False
        If the network has to be considered as directed or undirected. Relevant
        only for the computation of the number of common edges.
    Returns
    -------
    array-like:
        The adjacency matrix of the common sub graph.
    list:
        The list of common nodes.
    int:
        The number of common edges in the graph.
    """
    nodes_list_A = [str(n).lower() for n in nodes_list_A]
    nodes_list_B = [str(n).lower() for n in nodes_list_B]
    common_nodes = sorted(list(set(nodes_list_A).intersection(set(nodes_list_B))))
    ix_a = [nodes_list_A.index(n) for n in common_nodes]
    ix_b = [nodes_list_B.index(n) for n in common_nodes]
    A_sub = A[ix_a, :]
    A_sub = A_sub[:, ix_a]
    B_sub = B[ix_b, :]
    B_sub = B_sub[:, ix_b]
    assert A_sub.shape == B_sub.shape
    A_sub = (A_sub != 0).astype(int)
    B_sub = (B_sub != 0).astype(int)
    common = A_sub + B_sub
    common -= np.diag(np.diag(common))
    ix = sorted(list(set(np.where(common==2)[0]).union(set(np.where(common==2)[1]))))

    nodes = np.array(nodes_list_A)[ix_a][ix]
    common_sub = common[ix,:]
    common_sub = common_sub[:, ix]
    common_sub[np.where(common_sub==1)] = 0
    common_sub = common_sub / 2
    n_edges = np.sum(common_sub) if directed else np.sum(common_sub)/2
    return common_sub, nodes, n_edges


def multi_common_subgraph(Gs, nodes_lists, directed=False):
    common, nodes, _ = common_subgraph(Gs[0], Gs[1],
                                       nodes_lists[0], nodes_lists[1],
                                       directed)
    for i in range(2, len(Gs)):
        common, nodes, n_edges = common_subgraph(Gs[i], common,
                                                 nodes_lists[i], nodes,
                                                 directed)
    return common, nodes, n_edges
