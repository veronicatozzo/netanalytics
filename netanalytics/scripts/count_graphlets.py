"""
	Code to identify the graphlet signatures of all nodes in a network
	using the ORCA method.

    Usage example:
		./count_graphlets.py <tabular_edge_list_file> <graphlet_size>

    It writes the output on a ndump2 file.

    Re-arrangement of original code from Omer Nebil Yaveroglu.
"""
import pandas as pd

from netanalytics.io import from_edge_list
from netanalytics.graphlet import get_graphlet_orbits_count

if __name__ == "__main__":
    filename = str(sys.argv[1])
    if len(sys.argv)-1 == 1:
        graphlet_size = 5
    else:
        graphlet_size = int(sys.argv[2])
        if not (graphlet_size == 4 or graphlet_size == 5):
            raise ValueError("The maximum graphlet size must be either 4 or 5.")
    print(filename)
    _, nodes_list, edges_list = from_edge_list(filename, sep=' ',
											   only_adjacency=False)
    df = get_graphlet_orbits_count(nodes_list, edges_list, graphlet_size)
	df.to_csv(filename.split(".")[0]+".csv")
