#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='stack multiple npz matrix')
parser.add_argument('-i', dest='input', nargs='+', help='input npz files')
parser.add_argument('-o', type=str, dest="output", help='output prefix')

args = parser.parse_args()

import sys
import copy
import numpy as np
from time import perf_counter as pc
from scipy import sparse
from scipy.sparse import save_npz, load_npz
from scipy.sparse import vstack

def main():
	input_lists = args.input
	outf = args.output
	start_time = pc()
	npz_list = []
	for i in input_lists:
		npz = load_npz(i)
		npz_list.append(npz)
	merged_npz = vstack(npz_list)
	outnpz = ".".join([outf,"npz"])
	save_npz(outnpz, merged_npz)
	end_time = pc()
	print('Used (secs): ', end_time - start_time)

if __name__ == "__main__":
	main()
