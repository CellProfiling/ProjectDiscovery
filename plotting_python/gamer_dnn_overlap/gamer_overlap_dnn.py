"""
NOTE:
Depends on the forked version of matplotlib_venn as available at:
    https://github.com/cwinsnes/matplotlib-venn

Info:
Plotting the overlap between project discovery gamers and the DNN.

Example run:
    python3 gamer_overlap_dnn.py ../../data/D1_hpa_v14_gold_standard.csv

    Requires the file IF_gamer.csv in the current directory as well as the files dnn_predictions-{1,2,3,4}.

Arguments:
    First argument has to be a file containing the columns
        if_plate_id,position,sample,locations
    with the corresponding information in the correct places.
    Any other data available in the file is ignored.

Outputs:
    A pdf in the current directory called gamer_over_dnn.py containing the gamer overlap venn diagrams.
"""
#!/usr/bin/python3
import sys
import csv
import matplotlib.patches
import pickle
import numpy as np
import collections
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import logging
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


class keydefdict(collections.defaultdict):
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        else:
            ret = self[key] = self.default_factory(key)
            return ret

# Create a common dictionary between the gamer classes and the dnn
TRANSLATIONS = keydefdict(lambda x: x, {
    # DNN
    'Aggresome': 'Aggresome',
    'Cell Junctions': 'Cell Junctions',
    'Centrosome': 'Centrosome',
    'Cytoplasm': 'Cytoplasm',
    'Cytoskeleton (Actin filaments)': 'Actin filaments',
    'Cytoskeleton (Cytokinetic bridge)': 'Cytokinetic bridge',
    'Cytoskeleton (Intermediate filaments)': 'Intermediate filaments',
    'Cytoskeleton (Microtubule end)': 'Microtubule end',
    'Cytoskeleton (Microtubules)': 'Microtubules',
    'Endoplasmic reticulum': 'Endoplasmic reticulum',
    'Focal Adhesions': 'Focal adhesions',
    'Golgi apparatus': 'Golgi apparatus',
    'Microtubule organizing center': 'Microtubule organizing center',
    'Mitochondria': 'Mitochondria',
    'NULL': '',
    'Nuclear bodies': 'Nuclear bodies',
    'Nuclear membrane': 'Nuclear membrane',
    'Nuclear speckles': 'Nuclear speckles',
    'Nucleoli': 'Nucleoli',
    'Nucleoli (Fibrillar center)': 'Nucleoli fibrillar center',
    'Nucleoplasm': 'Nucleoplasm',
    'Nucleus': 'Nucleus',
    'Plasma membrane': 'Plasma membrane',
    'Vesicles': 'Vesicles',


    # Gamer
    'Cytosol': 'Cytoplasm',
    'NULL': '',
    'Endosomes': 'Vesicles',
    'Lysosomes': 'Vesicles',
    'Peroxisomes': 'Vesicles',
    'Focal adhesion sites': 'Focal adhesions',
    'Focal Adhesions': 'Focal adhesions',
    'Microtubule ends': 'Microtubule end',
    'Nucleoli fibrillar center': 'Nucleoli fibrillar center',
    'Nucleoli (fib center)': 'Nucleoli fibrillar center',
    'Rods & Rings': 'R&R',
    }
                          )


def get_if_info_classes(if_info, locations_header='locations'):
    classes = set()
    for p in if_info:
        unsplit = if_info[p][locations_header]
        split = unsplit.split(',')
        classes.update(split)
    classes = sorted(classes)
    ordered_classes = collections.OrderedDict()
    for i, c in enumerate(classes):
        # in the unlikely event that a class is an int
        # (not a str representation of int)
        # throw an error because that shit is wrong yo
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


def get_cutoffs(path):
    cutoffs = {}
    for line in open(path, 'r'):
        split = line.split(':')
        cutoff = float(split[1])
        if (cutoff < 0.0001):
            cutoff = 0.4
        cutoffs[split[0]] = cutoff
    return cutoffs

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Must have an IF file argument')
        sys.exit(1)
    cutoffs = get_cutoffs('cutoffs')
    if_info = read_if_file(sys.argv[1], split_to_well=True)
    model_classes = get_if_info_classes(if_info)
    for i in range(int(len(model_classes)/2)):
        print(model_classes[i])

    gamer_to_dnn = collections.defaultdict(lambda: {'same': 0, 'non_same': 0})
    dnn_to_gamer = collections.defaultdict(lambda: {'same': 0, 'non_same': 0})
    hpa_overlap = collections.defaultdict(lambda: {'same': 0, 'non_same': 0})

    dnn_correct = collections.defaultdict(int)
    gamer_correct = collections.defaultdict(int)
    both_correct = collections.defaultdict(int)

    num_of_class = collections.defaultdict(int)

    gamer_if_info = read_if_file('IF_gamer.csv', split_to_well=True)
    gamer_model_classes = get_if_info_classes(gamer_if_info)

    num_iters = 0
    print(len(model_classes))
    # cutoffs should have been tuned to data-0. do not include that one in this analysis
    for f in ['dnn_predictions-1', 'dnn_predictions-2', 'dnn_predictions-3', 'dnn_predictions-4']:
        predictions = pickle.load(open(f, 'rb'))
        for fov, fov_dict in predictions.items():
            fov = '_'.join(fov.split('_')[:-1])
            num_iters += 1
            fov_prediction = fov_dict['prediction']
            fov_prediction = [fov_prediction[i] for i in range(len(fov_prediction)) if i not in [4, 11, 14]]
            fov_prediction = np.asarray(fov_prediction)

            fov_classes = set()

            for i in range(int(len(model_classes)/2)):
                c = model_classes[i]
                if fov_prediction[i] > cutoffs[c]:
                    fov_classes.add(c)
            fov_classes = {TRANSLATIONS[x] for x in fov_classes}
            try:
                gamer_predictions = gamer_if_info[fov]['gamer_consensus'].split(',')
                gamer_predictions = {TRANSLATIONS[x] for x in gamer_predictions}
                hpa_predictions = {TRANSLATIONS[x] for x in if_info[fov]['locations'].split(',')}

            except KeyError as e:
                print('Keyerror:', e)
                continue

            print('#'*30)
            print(hpa_predictions)
            print(fov_classes)
            print(gamer_predictions)
            print('#'*30)

            for dnn_classification in fov_classes:
                if dnn_classification in hpa_predictions:
                    if dnn_classification in gamer_predictions:
                        both_correct[dnn_classification] += 1
                    else:
                        dnn_correct[dnn_classification] += 1

                if dnn_classification in gamer_predictions:
                    dnn_to_gamer[dnn_classification]['same'] += 1
                else:
                    dnn_to_gamer[dnn_classification]['non_same'] += 1

                if dnn_classification in hpa_predictions:
                    if dnn_classification in gamer_predictions:
                        hpa_overlap[dnn_classification]['same'] += 1
                    else:
                        hpa_overlap[dnn_classification]['non_same'] += 1

            for gamer_classification in gamer_predictions:
                if gamer_classification in hpa_predictions:
                    if gamer_classification not in fov_classes:
                        gamer_correct[gamer_classification] += 1

                if gamer_classification in fov_classes:
                    gamer_to_dnn[gamer_classification]['same'] += 1
                else:
                    gamer_to_dnn[gamer_classification]['non_same'] += 1

                if gamer_classification in hpa_predictions:
                    if gamer_classification not in fov_classes:
                        hpa_overlap[gamer_classification]['non_same'] += 1

            for hpa_classification in hpa_predictions:
                num_of_class[hpa_classification] += 1

    print(num_iters)
    print()

    for c in num_of_class:
        print(num_of_class[c])
    for c in dnn_to_gamer:
        if c not in gamer_to_dnn:
            print(c)

    for i in range(int(len(model_classes)/2)):
        c = model_classes[i]
        if TRANSLATIONS[c] not in gamer_to_dnn:
            print(c)

    print('Class', ' '*15, '|DNN Overlap|Gamer overlap|HPA Overlap')
    for c in sorted(gamer_to_dnn.keys()):
        gamer_same = gamer_to_dnn[c]['same']
        gamer_non_same = gamer_to_dnn[c]['non_same']
        gamer_overlap = gamer_same/((gamer_same + gamer_non_same) or 1)

        dnn_same = dnn_to_gamer[c]['same']
        dnn_non_same = dnn_to_gamer[c]['non_same']
        dnn_overlap = dnn_same/((dnn_same + dnn_non_same) or 1)

        hpa_same = hpa_overlap[c]['same']
        hpa_non_same = hpa_overlap[c]['non_same']
        hpa_overcount = hpa_same/((hpa_same + hpa_non_same) or 1)

        print('{:22}|{:11}|{:13}|{:13}'.format(c[:22], str(dnn_overlap)[:11],
                                               str(gamer_overlap)[:13],
                                               str(hpa_overcount)[:13]))

    print('dnn_correct', dnn_correct)
    print()
    print('gamer_correct', gamer_correct)
    print()
    print('both_correct', both_correct)
    print()

    fig = plt.figure(figsize=(16, 9))
    skipped = 0
    dnn_color = '#03256C'
    gamer_color = '#228B22'
    middle_color = '#F4C095'
    for i, c in enumerate(sorted(gamer_to_dnn.keys())):
        if num_of_class[c] < 1 or c == 'Unspecific':
            skipped += 1
            continue
        ax = fig.add_subplot(5, 5, i+1-skipped)
        ax.set_title(c + ' (' + str(num_of_class[c]) + ')', size=14, y=0.88, fontname='Arial')
        v = venn2(subsets=(dnn_correct[c], gamer_correct[c], both_correct[c]),
                  set_labels=('', ''), set_colors=(dnn_color, gamer_color, middle_color),
                  alpha=0.5, normalize_to=0.2, subset_label_fontargs={'fontname': 'Arial', 'fontsize': 11})

    c1 = matplotlib.patches.Circle((0, 0), 0.1, color=dnn_color, alpha=0.7)
    c2 = matplotlib.patches.Circle((0, 1), 0.1, color=gamer_color, alpha=0.7)
    fig.legend((c1, c2), ('DNN', 'Gamer'), loc='lower center')
    fig.savefig('gamer_overlap.pdf', dpi=1000, transparent=True)
