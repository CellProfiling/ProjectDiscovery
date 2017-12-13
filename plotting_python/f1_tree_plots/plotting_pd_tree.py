import textwrap
import copy
import matplotlib.pyplot as plt
import matplotlib.colors
from tree_utils import hierarchical_tree_pos, draw_hierarchical_tree
from tree_utils import DEFAULT_HIERARCHY
from tree_plot import DISP_HIERARCHY
import networkx as nx
import csv
import numpy as np

if __name__ == '__main__':
    G = nx.Graph()
    nodes = []
    node_sizes = []
    node_colors = []
    edges = []
    f1s = {}
    norm = matplotlib.colors.Normalize(0.3,1.0)
    cmap = plt.cm.ScalarMappable(norm, cmap=plt.cm.plasma)
    row_names = []
    #cmap.set_clim(0.5,0.6)
    np.set_printoptions(threshold=np.nan)
    #for f in ['merge_level_notuning0.csv', 'merge_level_notuning1.csv', 'merge_level_notuning2.csv']:
    for f in ['cell_merge_level_v4_0.csv', 'cell_merge_level_v4_1.csv', 'cell_merge_level_v4_2.csv']:
        for row in csv.DictReader(open(f)):
            try:
                DEFAULT_HIERARCHY[row['Name']]
            except:
                print(row['Name'], 'No code!')
            
            code = DEFAULT_HIERARCHY[row['Name']]
            if code == '' or int(float(row['count'])) == 0:
                continue

            print(row['Name'], code)
            row_names.append(row['Name'])
            G.add_edge(code[:-1], code)
            nodes.append(code)
            edges.append((code[:-1], code))
            #node_colors.append(plt.cm.viridis(float(row['f1_score'])))
            #node_colors.append(float(row['f1_score']))
            node_colors.append(cmap.to_rgba(float(row['f1_score'])))
            num = int(float(row['count']))
            node_sizes.append(max(300, np.sqrt(num) * 6)) # /3.0)  # min(1000, lengths[code]))

            f1s[code] = float(row['f1_score'])
    print(row_names)
    #nmin = min(node_colors)
    #nmax = max(node_colors)
    #node_colors = list(map(plt.cm.viridis, map(lambda x: ((0.6-0.5) * (x - nmin) / (nmax - nmin)) + 0.5, node_colors)))
    print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    pos = hierarchical_tree_pos(G, '')
    print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    print(pos.keys())
    distances = np.ones([len(nodes),len(nodes)])
    print(len(nodes)) 
    for i,(nodei,currname) in enumerate(zip(nodes,row_names)):
        print(currname)
        for j,nodej in enumerate(nodes):
            distances[i][j] = nx.shortest_path_length(G,source=nodei,target=nodej)
    key_string = ','.join(row_names)
    np.savetxt('distances.csv',distances,delimiter=',',header=key_string)   
#    with open('distances.csv', 'w') as csvfile:
#        spamwriter = csv.writer(csvfile, delimiter=',',
#                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
#        spamwriter.writerow(DEFAULT_HIERARCHY)
#        for label in DEFAULT_HIERARCHY:
#            spamwriter.writerow(label+','+)

    #distances = get_tree_distances(G,nodes=nodes)
    fig, ax = draw_hierarchical_tree(G, pos, nodes=nodes, node_size=node_sizes,
                                     node_color=node_colors, edges=edges, sm=cmap)



    labels = {k: v for (k, v) in DISP_HIERARCHY.items()}
    labelpos = copy.deepcopy(pos)
    f_tex = "{:.2f}"
    l1 = len('(fibrillar center)')
    l2 = len('bodies/speckles')
    for c in labelpos:
        if c == '':
            continue
        originalpos = labelpos[c]
        if c in f1s:
            if f1s[c] <= 0.45:
                color = 'white'
            else:
                color = 'black'
            ax.text(originalpos[0]-0.010, originalpos[1]-0.0025,
                    f_tex.format(f1s[c]), fontsize=8, color=color)
            if len(c) > 2:
                ax.text(originalpos[0], originalpos[1]-0.032,
                        '\n'.join(textwrap.wrap(labels[c], l1)), rotation=-40, fontsize=8,
                        fontdict={'family': 'Sans-serif'})
            else:  # if len(c) <= 2
                ax.text(originalpos[0]+0.0260, originalpos[1]-0.005,
                        '\n'.join(textwrap.wrap(labels[c], l2)), rotation=0, fontsize=8,
                        fontdict={'family': 'Sans-serif'})
    
    fig.savefig('devin_tree2.pdf', transparent=True)
    #fig.savefig('devin_tree_prop.pdf', transparent=True)
