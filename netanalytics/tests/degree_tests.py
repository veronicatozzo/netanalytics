import numpy as np

from netanalytics.degree import node_degree, number_of_edges, \
                                degree_distribution


def node_degree_test():
    A = np.array([[1,0,0], [1,1,1], [0,1,1]])
    assert(node_degree(A, 0)[1]==0)
    assert(node_degree(A, 1)[1]==2)
    assert(node_degree(A, 1)[0]==2)

    A = A * 0.5
    assert(node_degree(A, 0)[1]==0)
    assert(node_degree(A, 1)[1]==2)
    assert(node_degree(A, 1)[0]==1.)


def number_of_edges_test():
    A = np.array([[1,0,0], [1,1,1], [0,1,1]])
    assert(number_of_edges(A)==3)
    A = A * 0.5
    assert(number_of_edges(A)==3)


def degree_distribution_test():
    A = np.array([[1,0,0,0,0],
                  [1,1,1,1,1],
                  [0,1,1,1,0],
                  [0,0,0,1,1],
                  [0,0,1,0,1]])
    distr = degree_distribution(A)
    assert(distr[0]==1)
    assert(distr[1]==2)
    assert(distr[2]==1)
    assert(distr[4]==1)
