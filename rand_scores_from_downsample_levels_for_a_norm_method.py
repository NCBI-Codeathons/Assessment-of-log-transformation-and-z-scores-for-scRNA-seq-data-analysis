#Load libraries.

import pandas as pd
import numpy as np

from sklearn.decomposition import PCA
from sklearn.metrics.cluster import adjusted_rand_score

import phenograph

import matplotlib.pyplot as plt

from pylab import *

#Write function.
#Accept a dictionary of normalized matrices where the keys are downsample levels (0.1 to 1).
#Would run this method once per normalization method.
#Returns a single list.

def adjusted_rand_score_vector(normalized_matrices):
	PCA_model = PCA(n_components=1000,svd_solver='randomized')
	PC_column_names = ['PC' + str(i) for i in list(range(1,1001))]
	components_normed_data_full = pd.DataFrame(data = PCA_model.fit_transform(normalized_matrices[1]),columns = PC_column_names)
	full_communities, full_graph, full_Q = phenograph.cluster(components_normed_data_full)
	adj_rand_scores = []
	for split in list(np.array(range(1,10))/10):
		components_normed_data_downsample = pd.DataFrame(data = PCA_model.fit_transform(normalized_matrices[split]),columns = PC_column_names)
		downsample_communities,downsample_graph,downsample_Q = phenograph.cluster(components_normed_data_downsample)
		adj_rand_scores.append(adjusted_rand_score(full_communities,downsample_communities))
	return adj_rand_scores
