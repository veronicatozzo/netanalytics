import numpy as np
import pandas as pd

from netanalytics.degree import degree_distribution
from netanalytics.clusters import graph_clustering_coefficient
from netanalytics.graphlets import graphlet_degree_vectors
from netanalytics.graphlets import _graphlet_degree_distribution
from netanalytics.graphlets import graphlet_correlation_matrix
from netanalytics.utils import _from_adjacency_to_edges_list, \
                               _from_edges_list_to_adjacency


class Graph:
    """

    """
    def __init__(self, nodes_list, adjacency=None, edges_list=None, values = None
                directed=False, axis=0):
        self.adjacency = adj
        self.nodes = nodes_list
        self.edges = edges_list
        self.values = values
        self.directed = directed
        self.axis = axis

    def fit(self, graphlet_size=5):
        if not (graphlet_size == 4 or graphlet_size == 5):
            raise ValueError("The maximum graphlet size must be either 4 or 5.")

        if self.adjacency is None and self.edges_list is None:
            raise ValueError("You must provide one between adjacency matrix or"
                              "edges list")

        # TODO: check values length
        # TODO: check nodes length correspond to all the rest

        if self.adjacency is None:
            self.adjacency = _from_edges_list_to_adjacency(self.edges,
                                                           self.values)
        elif self.edges_list is None:
            self.edges_list, self.values = \
                            _from_adjacency_to_edges_list(self.adjacency)

        self.nodes_degree, self.degree_distribution = \
            degree_distribution(self.adjacency, axis=self.axis)
        self.clustering_coefficient = \
            graph_clustering_coefficient(self.adjacency,self.directed,
                                         self.axis)
        self.GDV = graphlet_degree_vectors(nodes_list, edges_list,
                                           graphlet_size=graphlet_size)
        self.GDD = _graphlet_degree_distribution(self.GDV)
        self.GCM = graphlet_correlation_matrix(self.GDV)
        self.is_fitted = True
        return self
