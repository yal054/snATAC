#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='plot sankey diagram')
parser.add_argument('--data', type=str, dest="data", help='summary of links for sankey diagram')
parser.add_argument('--idx', type=str, dest="idx", help='index for each elemnets')
parser.add_argument('-o', '--outPrefix', type=str, dest="outPrefix", help='output prefix')

args = parser.parse_args()

import pandas as pd
import numpy as np
import colorlover as cl
import plotly as py

from time import perf_counter as pc

def run():
	""" Run standard NMF on rank """
	start_time = pc()
	""" init input files """
	dataf = args.data
	idxf = args.idx
	prefix = args.outPrefix
	inF = read_files(dataf, idxf)
	data = inF['data']
	idx = inF['idx']
	mergedf = merge_df(data, idx)
	colors = map_color(idx)
	print("draw sankey diagram to html")
	plot_sankey(mergedf, idx, colors, prefix)
	end_time = pc()
	print('Used (secs): ', end_time - start_time)


def read_files(dataF, idxF):
	data = np.genfromtxt(dataF, dtype=None, names=True)
	idx = np.genfromtxt(idxF, dtype=None, names=True)
	return {'data':data, 'idx':idx}

def merge_df(data, idx):
	mergedf = pd.merge(pd.DataFrame(data),pd.DataFrame(idx),left_on='from', right_on='name', how="left")
	mergedf = mergedf.merge(pd.DataFrame(idx), left_on='to', right_on='name', how="left")
	return mergedf

def map_color(idx):
	color = cl.scales['12']['qual']['Set3'];
	color = cl.to_rgb(cl.interp(color, 27))
	colors = list(map(lambda i: color[i], list(idx['color'].astype(int))))
	return colors

def plot_sankey(mergedf, idx, colors, prefix):
	data = dict(
		type='sankey',
		node = dict(
			pad = 15,
			thickness = 20,
			line = dict(color = "black", width = 0.5),
			label = list(idx["name"].astype(str)),
			color = colors
		),
		link = dict(
			source = list(mergedf['idx_x']),
			target = list(mergedf['idx_y']),
			value = list(mergedf['count'])
		))

	layout =  dict(
		title = "Sankey Diagram for cell clusters at each rank",
		font = dict(size = 10)
		)

	fig = dict(data=[data], layout=layout)
	py.offline.plot(fig, filename = '.'.join([prefix, 'sankey_ranks.html']), auto_open=False)

if __name__ == "__main__":
	"""Draw sankey diagram for cell clusters at each rank"""
	run()
