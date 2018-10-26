
import numpy as np

def ER(n, m):
    """
    Erdos-Renyi random graph.
    n: int,
        Number of nodes
    m: int,
        Number of edges
    """
    comb = np.array(list(combinations(np.arange(0, n), 2)))
    np.random.shuffle(comb)
    m = int(m)
    selected_comb = comb[:m]
    x = [c[0] for c in selected_comb]
    y = [c[1] for c in selected_comb]
    network = np.zeros((n,n))
    network[x, y] = 1
    network[y, x] = 1

    return network


def ER_DD():
    """
    Generalized random model
    """


def GEO():
    """
    Geometric model
    """


def GEO_GD():
    """
    Geometric models with gene duplication.
    """

def SF_BA():
    """
    Barabasi-Albet Scale free model
    """

def SF_GD():
    """
    Scale free model with gene duplication and divergence
    """

def STICKY():
    """
    Stickiness-index based model
    """
