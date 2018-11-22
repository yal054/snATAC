#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='Run NMF using sklearn.')
parser.add_argument('--normH', type=str, dest="normH", help='coefficient matrix H')
parser.add_argument('--statH', type=str, dest="statH", help='input statH matrix')
parser.add_argument('-p', '--perplexity', type=int, dest="perplexity", default=0.3, help='an int for perplexity')
parser.add_argument('-o', '--outPrefix', type=str, dest="outPrefix", help='output prefix')

args = parser.parse_args()

from os.path import dirname, abspath, join
from warnings import warn

import numpy as np
import scipy as sp
from scipy import io

from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_samples, silhouette_score, pairwise_distances
import umap

from time import perf_counter as pc

import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.switch_backend('agg')
try:
    from matplotlib.pyplot import savefig, imshow, set_cmap
except ImportError as exc:
    warn("Matplotlib must be installed.")

def run():
	""" Run standard NMF on rank """
	start_time = pc()
	""" init input files """
	normHf = args.normH
	statHf = args.statH
	outPrefix = args.outPrefix
	perplexity = args.perplexity
	inF = read_files(normHf, statHf)
	normH = inF['normH']
	o_stat_H = inF['statH']
	rank = inF['rank']
	print("draw silhouette & tsne plot")
	o_silhouette = cal_silhouette(normH,o_stat_H)
	X_dist = cal_pairwise_pearson(normH)
	"""
	X_transformed = cal_tSNE(X_dist, perplexity)
	np.savetxt('.'.join([outPrefix, "tsne.xy"]), X_transformed, fmt= "%g", delimiter="\t")
	"""
	X_transformed = cal_umap(X_dist)
	np.savetxt('.'.join([outPrefix, "umap.xy"]), X_transformed, fmt= "%g", delimiter="\t")
	
	plot_silhouette_tsne(o_silhouette, X_transformed, o_stat_H, rank, outPrefix)
	end_time = pc()
	print('Used (secs): ', end_time - start_time)

def read_files(normH, statH):
	normH = np.loadtxt(normH)
	statH = np.genfromtxt(statH, dtype=None, names=True)
	rank = normH.shape[0]
	return {'normH':normH, 'statH':statH, 'rank':rank}
	
def cal_silhouette(normH,o_stat_H):
	X = normH.T
	n_clusters = normH.shape[0]
	cluster_labels = o_stat_H['class0'].astype(float).astype(int)
	"""
	# The silhouette_score gives the average value for all the samples.
	# This gives a perspective into the density and separation of the formed clusters
	"""
	silhouette_avg = silhouette_score(X, cluster_labels)
	print("For n_clusters =", n_clusters, "The average silhouette_score is :", silhouette_avg)
	"""# Compute the silhouette scores for each sample"""
	sample_silhouette_values = silhouette_samples(X, cluster_labels)
	return {'silhouette':silhouette_avg, 'silhouette_values':sample_silhouette_values}

def cal_nmds(X):
	embedding = MDS(n_components=2, dissimilarity='precomputed', n_jobs=2, verbose=2)
	X_transformed = embedding.fit_transform(X)
	return X_transformed

def cal_pairwise_dist(normH):
	X = normH.T
	X_pairwise_dist = pairwise_distances(X)
	return X_pairwise_dist

def cal_pairwise_pearson(normH):
	X = normH.T
	X_corrcoef = np.corrcoef(X)
	X_corrcoef_dist = np.sqrt(2*(1-X_corrcoef))
	return X_corrcoef_dist

def cal_tSNE(X_dist, p):
	X_transformed = TSNE(n_components=2, perplexity=p, random_state=1, verbose=2, metric="precomputed").fit_transform(X_dist)
	return X_transformed

def cal_umap(X):
	embedding = umap.UMAP(metric='correlation')
	#embedding = umap.UMAP(n_neighbors=p, min_dist=0.1, metric='correlation')
	X_transformed = embedding.fit_transform(X)
	return X_transformed

def plot_silhouette_tsne(o_silhouette, X_transformed, o_stat_H, rank, prefix):
	n_clusters = rank
	silhouette_avg = o_silhouette['silhouette']
	sample_silhouette_values = o_silhouette['silhouette_values']
	cluster_labels = o_stat_H['class0'].astype(float).astype(int)

	"""# Create a subplot with 1 row and 2 columns"""
	fig, (ax1, ax2) = plt.subplots(1, 2)
	fig.set_size_inches(18, 7)
	fig.set_dpi(300)
	"""
	# The 1st subplot is the silhouette plot
	# The silhouette coefficient can range from -1, 1 but in this example all
	# lie within [-0.1, 1]
	"""
	ax1.set_xlim([-0.1, 1])
	"""
	# The (n_clusters+1)*10 is for inserting blank space between silhouette
	# plots of individual clusters, to demarcate them clearly.
	"""
	ax1.set_ylim([0, len(X_transformed) + (n_clusters + 1) * 10])

	y_lower = 10
	for i in range(n_clusters):
		"""# Aggregate the silhouette scores for samples belonging to cluster i, and sort them """
		ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
		ith_cluster_silhouette_values.sort()
		size_cluster_i = ith_cluster_silhouette_values.shape[0]
		y_upper = y_lower + size_cluster_i
		color = cm.nipy_spectral(float(i) / n_clusters)
		ax1.fill_betweenx(np.arange(y_lower, y_upper),
			0, ith_cluster_silhouette_values,
			facecolor=color, edgecolor=color, alpha=0.7)
		"""# Label the silhouette plots with their cluster numbers at the middle"""
		ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i+1))
		"""# Compute the new y_lower for next plot"""
		y_lower = y_upper + 10  # 10 for the 0 samples

	ax1.set_title("The silhouette plot for the various clusters.")
	ax1.set_xlabel(" ".join(["The silhouette coefficient values (mean:",str(silhouette_avg),")"]))
	ax1.set_ylabel("Cluster label")

	"""# The vertical line for average silhouette score of all the values """
	ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
	ax1.set_yticks([])  # Clear the yaxis labels / ticks
	ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

	"""# 2nd Plot showing the actual clusters formed"""
	colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
	ax2.scatter(X_transformed[:, 0], X_transformed[:, 1], marker='.', s=30, lw=0, alpha=0.7, c=colors, edgecolor='k')

	""" lable non-classified cells"""

	ax2.set_title("The visualization of the clustered data.")
	ax2.set_xlabel("Feature space for the 1st tsne")
	ax2.set_ylabel("Feature space for the 2nd tsne")
	plt.suptitle(("Silhouette analysis for tSNE clustering on coefficient matrix H "
		"with n_clusters = %d" % n_clusters),fontsize=14, fontweight='bold')
	#fig.savefig('.'.join([prefix, "silhouette_tsne", "png"]))
	fig.savefig('.'.join([prefix, "silhouette_umap", "png"]))


if __name__ == "__main__":
	"""Draw silhouette and tsne plot for coefficient matrix"""
	run()
