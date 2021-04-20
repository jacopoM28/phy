#!/usr/bin/env python

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import argparse

parser = argparse.ArgumentParser(description='Plot of domain architecture networks.')
#Domain architecture input file input file
parser.add_argument('--in_domains', required=True, help='The domain structure file.')

args = parser.parse_args()

INPUT_FILE = args.in_domains

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import argparse

#Reading of the input file (output of module 1) and creation of dataframe for subsequent analyses :
Domain_arr = pd.read_csv(INPUT_FILE,sep='\t',header=0)
list_tmp = list(Domain_arr['Domain'].str.split('-'))
list_1 = []
list_2 = []

for sub_list in list_tmp :
    if len(sub_list) == 1 :
        list_1.append(sub_list[0])
        list_2.append("NONE")
    else :
        i = 1
        while i < len(sub_list) :
            list_1.append(sub_list[i-1])
            list_2.append(sub_list[i])
            i += 1
data_tuples = list(zip(list_1,list_2))
Domain_Network_tmp = pd.DataFrame(data_tuples, columns=['A','B']) #Can be overwritten 
Domain_Network = Domain_Network_tmp.groupby(Domain_Network_tmp.columns.tolist()).size().reset_index().rename(columns={0:'records'})

#Creation of a MultiDiGraph (Directionality, parallele branches and self connections allowed) :
G=nx.from_pandas_edgelist(Domain_Network, 'A', 'B', edge_attr=['records'],create_using=nx.MultiDiGraph())
G.remove_node("NONE")

#2. Create of a self-loop table in the same order as nodes in the network. 
#If a row is a self loop add at the column number of self-loops, else add 0.
Domain_Network
Domain_Network['Self_loop'] = np.where(Domain_Network['A']==Domain_Network['B'], Domain_Network['records'], 0)
Nodes_SelfLoops = pd.DataFrame(Domain_Network.drop(columns=['records','B']))
Nodes_SelfLoops = pd.DataFrame(Nodes_SelfLoops.groupby('A', as_index=False).max())
Nodes_SelfLoops
for index, row in Domain_Network.iterrows() :
    if not row[1] == "NONE" :
        boolean_finding = Nodes_SelfLoops['A'].str.contains(row[1]).any()
        if boolean_finding == False :
            missing_rows = pd.DataFrame({"A": [row[1]],"Self_loop": [0],})
            Nodes_SelfLoops = pd.concat([Nodes_SelfLoops, missing_rows])
Nodes_SelfLoops= Nodes_SelfLoops.set_index('A')
Nodes_SelfLoops=Nodes_SelfLoops.reindex(G.nodes())

#Creation of a MultiDiGraph (Directionality, parallele branches and self connections allowed) :
G=nx.from_pandas_edgelist(Domain_Network, 'A', 'B', edge_attr=['records'],create_using=nx.MultiDiGraph())
G.remove_node("NONE")

#---------MAIN---------#

#1. Summary Statistics :

#1.1 Check for singletons nodes and remove them from the plot :
#Rationale (Likely, they are miss-annotations)

Singletons = list(nx.isolates(G))
if not Singletons:
    print("---> Your gene families has not singleton domains!")
else :
    print("Your gene family has :" , len(list(nx.isolates(G))), "singleton domain(s). They will be removed from subsequent analyses but are stored in ...")
    G.remove_nodes_from(list(nx.isolates(G)))

#Extract widths attribute for plot (color or width of branches) :

colors = list(nx.get_edge_attributes(G,'records').values())

#Important: These plot cannot visualize self connections, for this purpose use the generated .graphml format and cytoscape

nx.draw(G, with_labels=True, node_color=Nodes_SelfLoops['Self_loop'].astype(int),arrows=True, edge_color=colors, cmap = plt.cm.Blues,font_size = 5,pos = nx.fruchterman_reingold_layout(G),node_size= 1000, edgecolors = "Black",edge_cmap=plt.cm.viridis, connectionstyle="arc3,rad=0.1")

ax = plt.gca()
ax.margins(0.2)
ax.collections[0] 
plt.savefig('random_geometric_graph.pdf')
plt.show()