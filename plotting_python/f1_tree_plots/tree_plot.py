"""
This plots an f1-score tree for the DNN model of the project discovery/machine learning paper.
The code depends on networkx as well as matplotlib and numpy.

To run, make sure that the files
tree_level0-{1,2,3,4}, tree_level1-{1,2,3,4}, tree_level2-{1,2,3,4} and hpa_gold_standard_tree.csv is in the same folder
as this program and then run
python3 tree_plot.py

The program will output a pdf called tree_output5.pdf that will contain the result tree.
"""
import textwrap
import logging
import tree_utils
import matplotlib.colors
import csv
import networkx as nx
import matplotlib.pyplot as plt
import pickle
import collections
import numpy as np
import copy
from collections import OrderedDict

CUTOFF_VALUES = collections.defaultdict(lambda: 0.4, {
    # H0
    'Aggresome': 0.060,
    'Cell Junctions': 0.098,
    'Centrosome': 0.094,
    'Cytoplasm': 0.437,
    'Cytoskeleton (Actin filaments)': 0.134,
    'Cytoskeleton (Cytokinetic bridge)': 0.004,
    'Cytoskeleton (Intermediate filaments)': 0.085,
    'Cytoskeleton (Microtubule end)': 0.070,
    'Cytoskeleton (Microtubules)': 0.165,
    'Endoplasmic reticulum': 0.185,
    'Focal Adhesions': 0.126,
    'Golgi apparatus': 0.181,
    'Microtubule organizing center': 0.048,
    'Mitochondria': 0.264,
    'Negative': 0.306,
    'Nuclear bodies': 0.141,
    'Nuclear membrane': 0.125,
    'Nuclear speckles': 0.161,
    'Nucleoli (Fibrillar center)': 0.168,
    'Nucleoli': 0.183,
    'Nucleoplasm': 0.392,
    'Nucleus': 0.238,
    'Plasma membrane': 0.278,
    'Unspecific': 0.443,
    'Vesicles': 0.271,
    'NULL': 1.5,


    # H1
    'Cytoplasmic structures': 0.427,
    'Cytoskeleton_1': 0.213,
    'MTOC': 0.093,
    'Neg/Unspec': 0.307,
    'Nuclear bodies/speckles': 0.172,
    'Nuclear membrane_1': 0.129,
    'Nucleoli_1': 0.201,
    'Nucleus_1': 0.500,
    'Periphery': 0.282,
    'Secretory': 0.386,

    # H2
    'Cytoplasmic': 0.516,
    'Nuclear': 0.478,
    'Cytoskeleton': 0.216,
    'Unspecific/Negative': 0.296,
    })

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
        'Nucleoli (rim)': '124',

        'Nuclear membrane_1': '13',
        'Nuclear membrane': '131',

        'Cytoplasmic': '2',

        'Cytoplasmic structures': '20',
        'Cytoplasm': '201',
        'Aggresome': '202',
        'Mitochondria': '203',

        'Cytoskeleton_1': '21',
        'Intermediate filaments': '211',
        'Cytoskeleton (Intermediate filaments)': '211',
        'Microtubule end': '212',
        'Cytoskeleton (Microtubule end)': '212',
        'Microtubules': '213',
        'Cytoskeleton (Microtubules)': '213',
        'Actin filaments': '214',
        'Cytoskeleton (Actin filaments)': '214',
        'Cytokinetic bridge': '215',
        'Cytoskeleton (Cytokinetic bridge)': '215',
        'Cytoskeleton': '',

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
        'Plasma membrane': '242',

        'Unspecific/Negative': '9',
        'Neg/Unspec': '90',
        'Negative': '901',
        'Unspecific': '902',
        'NULL': '903',

        'R&R': '',
        'CCD': '',
        'Periphery_1': '',

        'CCV': '',
        'CCV_2': '',
        ' ': '',
        }
DISP_HIERARCHY = {
        '': 'ROOT',

        '1': 'Nuclear',

        '10': 'Nucleus',
        '101': 'Nucleus',
        '102': 'Nucleoplasm',

        '11': 'Subnuclear structures',
        '111': 'Nuclear bodies',
        '112': 'Nuclear bodies (many)',
        '113': 'Nuclear speckles',

        '12': 'Nucleoli',
        '121': 'Nucleoli',
        '123': 'Nucleoli (Fibrillar center)',
        '124': 'Nucleoli (rim)',

        '13': 'Nuclear membrane',
        '131': 'Nuclear membrane',

        '2': 'Cytoplasmic',

        '20':  'Cytoplasm',
        '201': 'Cytosol',
        '202': 'Aggresome',
        '203': 'Mitochondria',

        '21': 'Cytoskeleton',
        '211': 'Intermediate filaments',
        '212': 'Microtubule end',
        '213': 'Microtubules',
        '214': 'Actin filaments',
        '215': 'Cytokinetic bridge',
        '': 'Cytoskeleton',

        '22': 'MTOC',
        '221': 'Microtubule organizing center',
        '222': 'Centrosome',

        '23': 'Secretory',
        '231': 'Endoplasmic reticulum',
        '232': 'Golgi apparatus',
        '233': 'Vesicles',
        '234': 'Lipid droplets',

        '24': 'Cell Periphery',
        '240': 'Cell Junctions',
        '241': 'Focal Adhesions',
        '242': 'Plasma membrane',

        '': 'Unspecific/Negative',
        '': 'Neg/Unspec',
        '': 'Negative',
        '': 'Unspecific',
        '': 'NULL',

        '': 'R&R',
        '': 'CCD',
        '': 'Periphery_1',

        '': 'CCV',
        '': 'CCV_2',

        }


skip_classes_in_output_l0 = {4, 11, 14}


def read_single_fov_counts_layer(preds, plates, classes, skip_classes, location_string):
    counts = collections.defaultdict(lambda: collections.defaultdict(int))
    for fov, fov_dict in preds.items():
        fov = '_'.join(fov.split('_')[:-1])
        fov_prediction = fov_dict['prediction']
        fov_prediction = [fov_prediction[oc] for oc in range(len(fov_prediction)) if oc not in skip_classes]
        fov_prediction = np.asarray(fov_prediction)

        fov_locations = set(plates[fov][location_string].split(','))
        fov_classes = set()
        for i in range(int(len(classes)/2)):
            c = classes[i]
            if c not in CUTOFF_VALUES:
                print(c, 'not in cutoff values. Using 0.4')
            if fov_prediction[i] >= CUTOFF_VALUES[c]:
                fov_classes.add(c)

        for c in fov_classes:
            code = DEFAULT_HIERARCHY[c]
            if c in fov_locations:
                counts[code]['tp'] += 1
            else:
                counts[code]['fp'] += 1

        for c in fov_locations:
            code = DEFAULT_HIERARCHY[c]
            if c not in fov_classes:
                counts[code]['fn'] += 1
    return counts


def get_if_info_classes(if_info, locations_header):
    classes = set()
    for p in if_info:
        unsplit = if_info[p][locations_header]
        split = unsplit.split(',')
        classes.update(split)
    if locations_header == 'locations_h1':
        classes.add('NULL')
    if locations_header == 'locations_h2':
        classes.add('Neg/Unspec')
    classes = sorted(classes)
    ordered_classes = collections.OrderedDict()
    for i, c in enumerate(classes):
        if isinstance(c, int):
            logging.critical('Class CANNOT be an int')
            raise ValueError('{} is int when it should be str'.format(c))

        ordered_classes[c] = i
        ordered_classes[i] = c
    return ordered_classes


def read_if_file(if_file, include={}, non_include={}, identifiers=['if_plate_id', 'position', 'sample'],
                 split_to_well=True):
    dictreader = csv.DictReader(open(if_file))
    if_info = dict()
    for line in dictreader:
        skipping = False

        for constraint in include:
            actual_line = line[constraint]
            actual_line = actual_line.split(',')

            any_in = [actual in include[constraint] for actual in actual_line]
            if not any(any_in):
                skipping = True
                break

        if skipping:
            continue

        for constraint in non_include:
            actual_line = line[constraint]
            value = actual_line.split(',')

            any_in = [actual in non_include[constraint] for actual in value]
            if any(any_in):
                skipping = True
                break

        if skipping:
            continue

        id_ = ''
        for i in identifiers:
            id_ += line[i] + '_'

        if split_to_well:
            id_ = '_'.join(id_.split('_')[:-2])
        else:
            id_ = '_'.join(id_.split('_')[:-1])
        if_info[id_] = line
    return if_info


if __name__ == '__main__':
    counts = collections.defaultdict(lambda: collections.defaultdict(int))

    print('H0')
    reader = read_if_file('hpa_gold_standard_tree.csv')
    classes = get_if_info_classes(reader, 'locations')
    for f in ['tree_level0-1', 'tree_level0-2', 'tree_level0-3', 'tree_level0-4']:
        h0_preds = pickle.load(open(f, 'rb'))
        h0_counts = read_single_fov_counts_layer(h0_preds, reader, classes, skip_classes_in_output_l0, 'locations')
        for c in h0_counts:
            counts[c]['tp'] += h0_counts[c]['tp']
            counts[c]['fp'] += h0_counts[c]['fp']
            counts[c]['fn'] += h0_counts[c]['fn']

    print('H1')
    reader = read_if_file('hpa_gold_standard_tree.csv')
    classes = get_if_info_classes(reader, 'locations_h1')
    for f in ['tree_level1-1', 'tree_level1-2', 'tree_level1-3', 'tree_level1-4']:
        h1_preds = pickle.load(open(f, 'rb'))
        h1_counts = read_single_fov_counts_layer(h1_preds, reader, classes, [], 'locations_h1')
        for c in h1_counts:
            if c == DEFAULT_HIERARCHY['Nucleus_1']:
                print(c, counts[c]['tp'])
            counts[c]['tp'] += h1_counts[c]['tp']
            counts[c]['fp'] += h1_counts[c]['fp']
            counts[c]['fn'] += h1_counts[c]['fn']

    print('H2')
    reader = read_if_file('hpa_gold_standard_tree.csv')
    classes = get_if_info_classes(reader, 'locations_h2')
    for f in ['tree_level2-1', 'tree_level2-2', 'tree_level2-3', 'tree_level2-4']:
        h2_preds = pickle.load(open(f, 'rb'))
        h2_counts = read_single_fov_counts_layer(h2_preds, reader, classes, [], 'locations_h2')
        for c in h2_counts:
            counts[c]['tp'] += h2_counts[c]['tp']
            counts[c]['fp'] += h2_counts[c]['fp']
            counts[c]['fn'] += h2_counts[c]['fn']

    precisions = {}
    recalls = {}
    lengths = {}
    f1s = {}
    print('making precision_recalls')
    for c in counts:
        tp = counts[c]['tp']
        fp = counts[c]['fp']
        fn = counts[c]['fn']
        precisions[c] = tp/((tp + fp) or 1)
        recalls[c] = tp/((tp + fn) or 1)
        lengths[c] = tp + fn
        f1s[c] = 2 * (precisions[c] * recalls[c])/(
                (precisions[c] + recalls[c]) or 1)

    G = nx.DiGraph()
    nodes = []
    node_sizes = []
    edges = []
    edge_colors = []
    node_colors = []

    labels = {k: v for (k, v) in DISP_HIERARCHY.items()}

    norm = matplotlib.colors.Normalize(0.3, 1.0)
    cmap = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.plasma)
    print('sorting')
    sort_f1s = OrderedDict(sorted(f1s.items()))
    sort_recalls = OrderedDict(sorted(recalls.items()))
    sort_precisions = OrderedDict(sorted(precisions.items()))
    sort_lengths = OrderedDict(sorted(lengths.items()))
    label_sizes = {}
    for code in sort_precisions:
        if code == '' or code == '9' or code == '90' or code == '901' or code == '902' or code == '903':
            continue
        G.add_edge(code[0:-1], code)
        edges.append((code[0:-1], code))
        edge_colors.append(plt.cm.viridis(sort_recalls[code]))

        node_sizes.append(max(300, np.sqrt(sort_lengths[code]) * 12.7))
        label_sizes[code] = node_sizes[-1]
        print(code)
        nodes.append(code)
        node_colors.append(cmap.to_rgba(sort_f1s[code]))

    pos = tree_utils.hierarchical_tree_pos(G, '')
    print(pos)
    sort_edges = sorted(edges)
    fig, ax = tree_utils.draw_hierarchical_tree(
            G, pos, nodes=nodes, node_size=node_sizes, node_color=node_colors,
            edges=sort_edges, edge_color='k', sm=cmap)
    print(node_sizes)

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
            else:
                ax.text(originalpos[0]+0.0250, originalpos[1]-0.005,
                        '\n'.join(textwrap.wrap(labels[c], l2)), rotation=0, fontsize=8,
                        fontdict={'family': 'Sans-serif'})

    fig.savefig('tree_output5.pdf', dpi=100, transparent=True)
