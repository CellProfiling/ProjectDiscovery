import matplotlib.pyplot as plt
import networkx as nx
from collections import OrderedDict
from pprint import pprint


def hierarchical_tree_pos(G, root, width=1., vgap=0.2, starty=20, startx=0.5):
    """
    Finds hierarchical positions of a graph and returns positions as
    a hierarcical tree.
    Args:
        G: The NetworkX graph (hierachical tree) to get the drawing positions
        for.
        root: The label of the root node of the tree.
        width: The width (counted from zero) that can be used for the figure.
        vgap: The vertical gap between nodes.
        starty: The y position for the root node.
        startx: The x position for the root node.
    Returns:
        A position dictionary to be used for NetworkX graph drawing functions.
    """
    pos = {root: (startx, starty)}
    neighbors = G.neighbors(root)
    parents = {root}
    print('sorting neighbors')
    print(neighbors)
    neighbors = sorted(neighbors)
    print(neighbors)

    while neighbors:
        dwidth = width / len(neighbors)
        base_x = startx - width/2 - dwidth/2
        new_neighbors = []
        starty = starty - vgap

        for neighbor in sorted(neighbors):
            base_x += dwidth
            pos[neighbor] = (base_x, starty)
            g_neighbors = G.neighbors(neighbor)
            new_neighbors.extend(g_neighbors)

        new_neighbors = [n for n in sorted(new_neighbors) if n not in parents]
        parents = set(neighbors)
        neighbors = new_neighbors
    return pos


def draw_hierarchical_tree(G, pos, node_size=500, nodes=None,
                           node_color='g', edges=None,
                           edge_width=1.0, edge_color='k',
                           axes=[0.1, 0.1, 0.9, 0.9],
                           figsize=(20, 10),
                           sm=plt.cm.ScalarMappable(cmap=plt.cm.viridis)):
    """
    Draws a hierarchical tree using the tree defined in G and the
    corresponding pos variable.
    Args:
        G: The NetworkX graph to draw.
        pos: A position dictionary, such as the one returned from
             the hierarhical_tree_pos function.
        node_size: A list of node sizes. Each item in the list should have a
                   corresponding item in `nodes`.
                   If this is an integer, each node will be drawn with that
                   size.
        nodes: A list of nodes to be drawn from the graph.
               If None, every node will be drawn.
        node_color: A list of colors. Each element should have a corresponding
                    item in `nodes`.
                    If this is a single color, each node will be drawn with
                    that colour.
        edges: The edges to be drawn from the graph.
               If None, every edge in the graph will be drawn.
        edge_width:  A list of floats describing the widths of the edges.
                     Each element should have a corresponding item in `edges`.
                     If this is a single number, every edge will have that
                     width.
        edge_color: A list of colors describing the colors of the edges.
                    Each element should have a corresponding item in `edges`.
                    If this is a single color, each edge will have that color.
        axes: The starting axes of the plot.
        figsize: A tuple of the size, in inches, of the figure.
    Returns:
        A pyplot figure and axes.
        The figure will be pyplots current figure.
    """
    fig = plt.figure(figsize=figsize)
    ax = plt.Axes(fig, axes)
    ax.set_axis_off()
    fig.add_axes(ax)

    if not nodes:
        nodes = G.nodes()
    if not edges:
        edges = G.edges()

    print('edges ')
    print(edges)
    print('pos: ')
    #print(nodes)
    print(pos)
    print('sorting pos')
    sort_pos = OrderedDict(pos)
    print(sort_pos)
    nx.draw_networkx_nodes(G, sort_pos, nodelist=nodes, node_size=node_size,
                           node_color=node_color)
    nx.draw_networkx_edges(G, sort_pos, edgelist=edges, width=edge_width,
                           edge_color=edge_color, arrows=False)

    sm._A = []
    fig.colorbar(sm, shrink=0.9)
    return fig, ax


# Hierarchy is defined as a string containing numbers
# Highest number indicates type, Nucleus:0, Cytoplasm:1, or Secretory:2
# Second highest indicates category within type
# Third indicates class within category
DEFAULT_HIERARCHY = {
        'ROOT': '',

        'Nuclear': '1',

        'Nucleus_1': '10',
        'Nucleus': '101',
        'Nucleoplasm': '102',

        'Nuclear bodies/speckles': '11',
        'Nuclear bodies': '111',
        'Nuclear bodies (many)': '112',
        'Nuclear speckles': '113',

        'Nucleoli_1': '12',
        'Nucleoli': '121',
        'Nucleoli (Fibrillar center)': '123',
        'Nucleoli (fib center)': '123',
        'Nucleoli fibrillar center': '123',
        'Nucleoli rim': '124',
        'Nucleoli (rim)': '124',

        'Nuclear membrane_1': '13',
        'Nuclear membrane': '131',

        'Cytoplasmic': '2',

        'Cytoplasmic structures': '20',
        'Cytoplasm': '201',
        'Cytosol': '201',
        'Aggresome': '202',
        'Mitochondria': '203',

        'Cytoskeleton_1': '21',
        'Cytoskeleton (Intermediate filaments)': '211',
        'Intermediate filaments': '211',
        'Cytoskeleton (Microtubule end)': '212',
        'Microtubule ends': '212',
        'Microtubule end': '212',
        'Cytoskeleton (Microtubules)': '213',
        'Microtubules': '213',
        'Cytoskeleton (Actin filaments)': '214',
        'Actin filaments': '214',
        'Cytoskeleton (Cytokinetic bridge)': '215',
        'Cytokinetic bridge': '215',
        'Cytoskeleton': '21',

        'MTOC': '22',
        'Microtubule organizing center': '221',
        'Centrosome': '222',

        'Secretory': '23',
        'Endoplasmic reticulum': '231',
        'Golgi apparatus': '232',
        'Vesicles': '233',
        'Lipid droplets': '234',

        'Periphery': '24',
        'Cell Junctions': '240',
        'Focal Adhesions': '241',
        'Focal adhesion sites': '241',
        'Plasma membrane': '242',

        'Unspecific/Negative': '9',
        'Neg/Unspec': '90',
        'Negative': '901',
        'Unspecific': '902',
        'NULL': '903',

        'R&R': '',
        'Rods & Rings': '',
        'CCD': '',
        'Periphery_1': '',

        'CCV': '',
        'CCV_2': '',
        }
