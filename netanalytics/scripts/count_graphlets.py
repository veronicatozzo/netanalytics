"""
	Code to identify the graphlet signatures of all nodes in a network
	using the ORCA method.

    Usage example:
		./count_graphlets.py <tabular_edge_list_file> <graphlet_size>

    It writes the output on a ndump2 file.

    Re-arrangement of original code from Omer Nebil Yaveroglu.
"""

from netanalytics.io import read_edgelist, to_ORCA

if __name__ == "__main__":
	filename = sys.argv[1]
    if len(sys.argv)-1 == 1:
         graphlet_size = 5
    else:
        graphlet_size = sys.argv[2]
        if not (graphlet_size == 4 or graphlet_size == 5):
            raise ValueError("The maximum graphlet size must be either 4 or 5.")

	_, nodes_list, edges_list = read_edgelist(filename)

	aux1 = filename.split('.')[0] + '_orca.txt'
	to_ORCA(nodes_list, edges_list, aux)

	aux2 = filename.split('.')[0] + '_orca.txt'
	cmd = './../../_orca/orca '+ str(graphlet_size)  + aux1 + ' ' + aux2
	os.system(cmd)

	originalNdump2 = netFileName.rsplit('.', 1)[0] + '.ndump2'
	with open(aux, 'rb') as _from:
	       with open(filename.rsplit('.')[0] + '.ndump2', 'wb') as _to:
               i = 0
               for line in _from:
                   _to.write(str(nodes_list[i])+' '+ line)
                   i +=1
	os.remove(aux1)
	os.remove(aux2)
