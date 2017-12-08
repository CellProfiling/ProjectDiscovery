#!/usr/bin/python3
"""
Example how to run this file:
    python3 f1_per_class.py ../data/D1_hpa_v14_gold_standard.csv ../data/bj_dnn/ ../data/bj_dnn/

    In this repo both D1_hpa_v14_gold_standard.csv and the bj_dnn folder are available under data/.

Arguments:
    First argument has to be a file containing the columns
        if_plate_id,position,sample,locations
    with the corresponding information in the correct places.
    Any other data available in the file is ignored.

    Cutoffs are first tuned on the 0th output set in second argument.
    It is important that the cutoff file ends in '-0' and that there is only one such file in the folder.

    Results are calculated on the items that do not end in '-0' from the folder in the third argument.

Outputs:
    Results are printed to stdout.
    The first set of results are the overall f1, precision, recall, and hamming score (jaccard similarity)
    for the result sets. The order in which they are printed are given from the order the files are printed
    before the outputs are printed.

    f1.pdf is a pdf of the f1 over cutoff for each class, given the cutoff dataset.
"""

import operator
import argparse
import pickle
import glob
import os
import numpy as np
import collections
import csv
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt


MINIMUM_CUTOFF = 0.07


skip_classes_in_output = {4, 11, 14}


def hamming_score(reals, predicts):
    numerator = 0
    denominator = 0
    for (r, p) in zip(reals, predicts):
        if len(r) != len(p):
            raise ValueError('Array lengths do not agree')

        true = set(np.where(r)[0])
        pred = set(np.where(p)[0])

        intersection = true.intersection(pred)
        union = true.union(pred)
        numerator += len(intersection)
        denominator += len(union)
    return numerator/denominator


def precision_recall(reals, predicts, instance='total'):
    if instance == 'total':
        axis = None
    elif instance == 'class':
        axis = 0
    else:
        raise ValueError('Not a valid instance type')

    reals = np.asarray(reals)
    predicts = np.asarray(predicts)

    truepos = np.logical_and(reals, predicts)
    false = reals - predicts
    falseneg = false > 0
    falsepos = false < 0

    truepos = np.sum(truepos, axis=axis)
    falseneg = np.sum(falseneg, axis=axis)
    falsepos = np.sum(falsepos, axis=axis)

    with np.errstate(divide='ignore', invalid='ignore'):
        precision = truepos/(truepos + falsepos)
        recall = truepos/(truepos + falseneg)
        if not np.isscalar(precision):
            precision[~np.isfinite(precision)] = 0
            recall[~np.isfinite(recall)] = 0
    return precision, recall


def convert_to_binary_classes(string_classes, if_info_classes):
    string_classes = string_classes.split(',')
    bin_classes = [0] * int(len(if_info_classes)/2)
    for c in string_classes:
        bin_classes[if_info_classes[c]] = 1
    return bin_classes


def get_if_info_classes(if_info, locations_header='locations'):
    classes = set()
    for p in if_info:
        unsplit = if_info[p][locations_header]
        split = unsplit.split(',')
        classes.update(split)
    classes = sorted(classes)
    ordered_classes = collections.OrderedDict()
    for i, c in enumerate(classes):
        if isinstance(c, int):
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


def calculate_cutoffs(predictions, if_info, if_info_classes, plot=False, locations_column='locations',
                      cutoff_step=0.001):
    best_results = collections.defaultdict(lambda: 0.01)
    best_cutoffs = collections.defaultdict(lambda: MINIMUM_CUTOFF)

    for i in range(int(len(if_info_classes)/2)):
        c = if_info_classes[i]
        best_cutoffs[c]  # To make sure every class is in there

    if plot:
        pp = matplotlib.backends.backend_pdf.PdfPages('f1.pdf')
        class_cutoff_f1s = collections.defaultdict(dict)
        num_instances = None

    for cutoff in np.arange(0, 1.1, cutoff_step):
        reals = []
        preds = []
        num_instances = collections.defaultdict(int)
        for id_, fov in predictions.items():
            id_ = '_'.join(id_.split('_')[:-1])
            real_classes = if_info[id_][locations_column]
            bin_classes = convert_to_binary_classes(real_classes, if_info_classes)
            prediction = fov['prediction']
            prediction = [prediction[oc] for oc in range(len(prediction)) if oc not in skip_classes_in_output]
            prediction = prediction > cutoff

            reals.append(bin_classes)
            preds.append(prediction)
        (p, r) = precision_recall(reals, preds, instance='class')

        for i in range(int(len(if_info_classes)/2)):
            c = if_info_classes[i]
            p_c = p[i]
            r_c = r[i]
            # Divide by 1 if both are 0, because the entire expression is 0 anyway
            f1_c = 2 * (p_c * r_c) / ((p_c + r_c) or 1.0)
            if best_results[c] < f1_c and cutoff > 0:
                best_results[c] = f1_c
                best_cutoffs[c] = cutoff

            if plot:
                class_cutoff_f1s[c][cutoff] = f1_c

    num_instances = [0] * len(reals[0])
    for r in reals:
        for j in range(len(r)):
            num_instances[j] += r[j]
    # for i in range(len(num_instances)):
    #     if num_instances[i] == 0:
    #         best_cutoffs[if_info_classes[i]] = 0.1
    if plot:
        for c in class_cutoff_f1s:
            if not num_instances[if_info_classes[c]]:
                continue
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel('Cutoff')
            ax.set_ylabel('F1-score')
            plt.ylim([0, 1.0])
            plt.xlim([0, 1.1])
            plt.title("%s (%d)" % (c, num_instances[if_info_classes[c]]))

            cutoff_values, f1_values = zip(*sorted(class_cutoff_f1s[c].items(), key=operator.itemgetter(0)))
            ind = np.argmax(f1_values)
            plt.axvline(x=cutoff_values[ind], ymin=0, ymax=f1_values[ind], color='r')
            plt.text(cutoff_values[ind]-0.05, f1_values[ind]+0.05,
                     '({:.3f}, {:.3f})'.format(cutoff_values[ind], f1_values[ind]))
            plt.plot(cutoff_values, f1_values)
            plt.savefig(pp, format='pdf')
            plt.close()
        pp.close()

    return best_cutoffs, num_instances


def get_precision_recall(predictions, if_info, if_info_classes, cutoffs, locations_column='locations'):
    num_instances = [0] * int(len(if_info_classes)/2)
    cutoff_array = [0] * int(len(if_info_classes)/2)
    for c in cutoffs:
        cutoff_array[if_info_classes[c]] = cutoffs[c]

    reals = []
    preds = []
    for id_, fov in predictions.items():
        id_ = '_'.join(id_.split('_')[:-1])
        real_classes = if_info[id_][locations_column]
        bin_classes = convert_to_binary_classes(real_classes, if_info_classes)
        num_instances = np.add(num_instances, bin_classes)
        prediction = fov['prediction']
        prediction = [prediction[oc] for oc in range(len(prediction)) if oc not in skip_classes_in_output]

        # prediction = np.greater_equal(prediction, 0.4)
        prediction[np.argmax(prediction)] = 1
        prediction = np.greater_equal(prediction, cutoff_array)

        reals.append(bin_classes)
        preds.append(prediction)
    (p, r) = precision_recall(reals, preds)
    hamming = hamming_score(reals, preds)
    print(2*(p*r)/(p+r), p, r, hamming)
    (p, r) = precision_recall(reals, preds, instance='class')

    prec_rec = {}
    for i in range(int(len(if_info_classes)/2)):
        c = if_info_classes[i]
        p_c = p[i]
        r_c = r[i]
        # Divide by 1 if both are 0, because the entire expression is 0 anyway
        f1_c = 2 * (p_c * r_c) / ((p_c + r_c) or 1.0)
        prec_rec[c] = (p_c, r_c, f1_c)
    return prec_rec, num_instances


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('if_file')
    parser.add_argument('cutoff_folder', help='Uses the 0th test set from the folder')
    parser.add_argument('results_folder', help='Uses the 1-4th test sets from the folder')
    parser.add_argument('--location-column', default='locations')
    args = parser.parse_args()

    if_file = args.if_file
    results_folder = args.results_folder
    results_files = glob.glob(os.path.join(results_folder, '*'))
    results_files = [f for f in results_files if not f.endswith('-0')]

    cutoff_folder = args.cutoff_folder
    cutoff_folder = glob.glob(os.path.join(cutoff_folder, '*-0'))
    cutoff_file = cutoff_folder[0]

    location_column = args.location_column

    if_info = read_if_file(if_file)
    if_info_classes = get_if_info_classes(if_info, location_column)
    print(len(if_info_classes)/2)
    for i in range(int(len(if_info_classes)/2)):
        print(i, if_info_classes[i])

    predictions = pickle.load(open(cutoff_file, 'rb'))
    cutoffs, _ = calculate_cutoffs(predictions, if_info, if_info_classes,
                                   plot=True, locations_column=location_column)

    num_instances = [0] * int(len(if_info_classes)/2)
    precs = [0] * int(len(if_info_classes)/2)
    recs = [0] * len(precs)
    f1s = [0] * len(precs)
    classes = set()
    print('#' * 30)
    print(results_files)
    print('#' * 30)
    print('Overall:')
    print('F1, precision, recall, hamming')
    for f in results_files:
        predictions = pickle.load(open(f, 'rb'))
        p_r, n_i = get_precision_recall(predictions, if_info, if_info_classes, cutoffs, location_column)
        num_instances = np.add(num_instances, n_i)

        for c in p_r:
            p, r, f = p_r[c]
            precs[if_info_classes[c]] += p
            recs[if_info_classes[c]] += r
            f1s[if_info_classes[c]] += f
            classes.add(c)
    print('#' * 30)

    i = 0
    avg_p = 0
    avg_r = 0
    avg_f = 0
    print('class, f1, precision, recall, cutoff')
    for c in sorted(classes):
        p = precs[if_info_classes[c]]/len(results_files)
        r = recs[if_info_classes[c]]/len(results_files)
        f = f1s[if_info_classes[c]]/len(results_files)
        print(c, '{:.2f}, {:.2f}, {:.2f}, [{:.3f}][{}]'.format(f, p, r, cutoffs[c],
                                                               num_instances[if_info_classes[c]]))

        if num_instances[if_info_classes[c]]:
            avg_f += f
            avg_p += p
            avg_r += r
            i += 1

    print('Avg F1:', avg_f/i)
    print('Avg Pr:', avg_p/i)
    print('Avg Re:', avg_r/i)
