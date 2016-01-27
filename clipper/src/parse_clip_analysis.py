__author__ = 'gpratt'

import pandas as pd
import cPickle as pickle
import os
import numpy as np
import argparse


def output_region_density(fn):
    with open(fn) as f_handle:
        result = pickle.load(f_handle)
    name = os.path.basename(fn).split(".")[0]
    regions_dict = {}
    for region in result['distributions']:
        region_type_dict = {}
        for region_type in ['total', 'individual']:
            count, bins = np.histogram(result['distributions'][region][region_type], bins=100, normed=True)
            region_type_dict[region_type] = {key: value for key, value in zip(np.arange(0,1, .01), count)}
        regions_dict[region] = pd.DataFrame(region_type_dict).T
    pd.concat(regions_dict).to_csv(name + ".distributions.csv")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Takes a pickle file from the output of clip_analysis and outputs the region distributions as csv in histogram form")
    parser.add_argument("-i", "--input", help="pickle file from clip analysis to parse")
    args = parser.parse_args()

    output_region_density(args.input)
