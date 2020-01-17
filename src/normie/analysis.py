import os
import sys
import scipy
import random
import umap
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import sparse
from sklearn.decomposition import PCA

import matplotlib
from matplotlib import pyplot as plt


def correlation_in_gene_groups(raw_data, normalized_data, save_dir=None, nbins = 5):
	# convert pandas dataframe to numpy array
	raw_data, normalized_data = raw_data.values, normalized_data.values

	# raw and normalized data must have cells along the rows and genes along the columns
	mean_before = np.asarray(np.mean(raw_data, axis = 0)).squeeze()
	log_mean_before = np.log2(np.asarray(np.mean(raw_data, axis = 0)).squeeze())

	mean_after = np.asarray(np.mean(normalized_data, axis = 0)).squeeze()
	log_mean_after = np.log2(np.asarray(np.mean(normalized_data, axis = 0)).squeeze())

	var_before = np.asarray(np.var(raw_data, axis = 0)).squeeze()
	log_var_before = np.log2(np.asarray(np.var(raw_data, axis = 0)).squeeze())

	var_after = np.asarray(np.var(normalized_data, axis = 0)).squeeze()
	log_var_after = np.log2(np.asarray(np.var(normalized_data, axis = 0)).squeeze())

	digitized = np.digitize(log_mean_before, np.linspace(min(log_mean_before), 
														 max(log_mean_before), nbins))

	log_mean_corr = []
	var_corr = []
	log_var_corr = []
	for j in np.unique(digitized):
		id_genes = np.where(digitized == j)[0]
		if len(id_genes) > 1:
			mean_corr.append(np.corrcoef(mean_before[id_genes], mean_after[id_genes])[0, 1])
			log_mean_corr.append(np.corrcoef(log_mean_before[id_genes], log_mean_after[id_genes])[0, 1])
			var_corr.append(np.corrcoef(var_before[id_genes], var_after[id_genes])[0, 1])
			log_var_corr.append(np.corrcoef(log_var_before[id_genes], log_var_after[id_genes])[0, 1])

	# plotting
	fig = plt.figure(figsize = (8*2, 6*2))
	ax = fig.add_subplot(2, 2, 1)
	ax.plot(range(len(mean_corr)), mean_corr)
	ax.set_xlabel('gene group')
	ax.set_ylabel('correlation of averages')

	ax = fig.add_subplot(2, 2, 2)
	ax.plot(range(len(log_mean_corr)), log_mean_corr)
	ax.set_xlabel('gene group')
	ax.set_ylabel('correlation of log averages')

	ax = fig.add_subplot(2, 2, 3)
	ax.plot(range(len(var_corr)), var_corr)
	ax.set_xlabel('gene group')
	ax.set_ylabel('correlation of variances')

	ax = fig.add_subplot(2, 2, 4)
	ax.plot(range(len(log_var_corr)), log_var_corr)
	ax.set_xlabel('gene group')
	ax.set_ylabel('correlation of log variances')

	if save_dir is not None:
		plt.savefig(save_dir)

	return mean_corr, log_mean_corr, var_corr, log_var_corr


def compare_statistics(raw_data, normalized_data, save_dir=None):
	# convert pandas dataframe to numpy array
	raw_data, normalized_data = raw_data.values, normalized_data.values

	# raw and normalized data must have cells along the rows and genes along the columns
	mean_before = np.asarray(np.mean(raw_data, axis = 0)).squeeze()
	mean_after = np.asarray(np.mean(normalized_data, axis = 0)).squeeze()

	mean_corr = np.corrcoef(mean_before, mean_after)[0, 1]
	log_mean_corr = np.corrcoef(np.log2(mean_before), np.log2(mean_after))[0, 1]

	var_before = np.asarray(np.var(raw_data, axis = 0)).squeeze()
	var_after = np.asarray(np.var(normalized_data, axis = 0)).squeeze()

	var_corr = np.corrcoef(var_before, var_after)[0, 1]
	log_var_corr = np.corrcoef(np.log2(var_before), np.log2(var_after))[0, 1]

	# plotting
	fig = plt.figure(figsize = (8*2, 6*2))
	ax = fig.add_subplot(2, 2, 1)
	ax.scatter(mean_before, mean_after)
	ax.set_xlabel('mean before normalization')
	ax.set_ylabel('mean after normalization')
	ax.set_title('Correlation = ' + str(mean_corr))

	ax = fig.add_subplot(2, 2, 2)
	ax.scatter(np.log2(mean_before), np.log2(mean_after))
	ax.set_xlabel('log2 mean before normalization')
	ax.set_ylabel('log2 mean after normalization')
	ax.set_title('Correlation = ' + str(log_mean_corr))

	ax = fig.add_subplot(2, 2, 3)
	ax.scatter(var_before, var_after)
	ax.set_xlabel('variance before normalization')
	ax.set_ylabel('variance after normalization')
	ax.set_title('Correlation = ' + str(var_corr))

	ax = fig.add_subplot(2, 2, 4)
	ax.scatter(np.log2(var_before), np.log2(var_after))
	ax.set_xlabel('log2 variance before normalization')
	ax.set_ylabel('log2 variance after normalization')
	ax.set_title('Correlation = ' + str(log_var_corr))

	if save_dir is not None:
		plt.savefig(save_dir)

	return mean_corr, log_mean_corr, var_corr, log_var_corr


def corr_lib_size(raw_data, normalized_data):
	# convert pandas dataframe to numpy array
	raw_data, normalized_data = raw_data.values, normalized_data.values

	# raw and normalized data must have cells along the rows and genes along the columns
	lib_size = np.asarray(np.sum(raw_data, axis = 1)).squeeze()
	corr_lib_size = [np.corrcoef(normalized_data[:, j], lib_size)[0, 1] for j in range(normalized_data.shape[1])]
	corr_log2_lib_size = [np.corrcoef(normalized_data[:, j], np.log2(lib_size))[0, 1] for j in range(normalized_data.shape[1])]

	return corr_lib_size, corr_log2_lib_size


def downsample_raw_data(raw_data, data_split, overlap_factor = 0.0, random_state_number = None):
	# convert pandas dataframe to numpy array
	umis = raw_data.values

	random_state = np.random.RandomState() if random_state_number is None \
										   else np.random.RandomState(random_state_number)

	umis_X_disjoint = random_state.binomial(umis, data_split - overlap_factor)
	umis_Y_disjoint = random_state.binomial(umis - umis_X_disjoint, (1 - data_split) / (1 - data_split + overlap_factor))
	overlap_factor = umis - umis_X_disjoint - umis_Y_disjoint
	umis_X = umis_X_disjoint + overlap_factor
	umis_Y = umis_Y_disjoint + overlap_factor

	return sparse.csr_matrix(umis_X), sparse.csr_matrix(umis_Y)


def co_embed(normalized_full_data, normalized_downsampled_data, save_dir=None):
	# convert pandas dataframe to numpy array
	normalized_full_data, normalized_downsampled_data = normalized_full_data.values, normalized_downsampled_data.values

	# combine data
	combined_data = np.concatenate((normalized_full_data, normalized_downsampled_data))

	# sample id
	sample_id = [0]*normalized_full_data.shape[0] + [1]*normalized_downsampled_data.shape[0]
	pca = PCA(n_components=200, svd_solver='randomized')
	principalComponents = pca.fit_transform(combined_data)
	umap_result = UMAP(n_neighbors=30, metric='euclidean', 
							random_state=7).fit_transform(principalComponents)
	color_pallete = ['g', 'r']
	file_names = ['full data', 'downsampled data']

	# plotting
	fig = plt.figure(figsize = (8*3, 6*1))
	ax = plt.add_subplot(1, 3, 1)
	for j in np.unique(sample_id):
		id_s = sample_id == j
		ax.scatter(umap_result[id_s, 0], umap_result[id_s, 1], 
					s = 2, 
					c = color_pallete[j], 
					label = file_names[j]);
	lgnd = plt.legend(loc='left', bbox_to_anchor= (1.01, 1.01), fontsize = 12, markerscale = 3);
	ax.set_xlabel('UMAP-1');
	ax.set_ylabel('UMAP-2');
	ax.set_xticks([], [])
	ax.set_yticks([], [])
	for item in np.unique(sample_id):
		ax = fig.add_subplot(1, 3, item + 2)
		id_s = sample_id == j
		ax.scatter(umap_result[:, 0], umap_result[:, 1], s = 0.1, c = 'k')
		ax.scatter(umap_result[id_s, 0], umap_result[id_s, 1], s = 5, c = color_pallete[item]);
		ax.set_xticks([], [])
		ax.set_yticks([], [])

	if save_dir is not None:
		plt.savefig(save_dir)

	return umap_result


def vst_hvg(count_data, n_top_genes = 50, span = 0.3):
	x, y = np.mean(count_data, axis = 0), np.var(count_data, axis = 0)
	id_non_const = (y != 0)

	# fit a loess regression
	w = sm.nonparametric.lowess(np.log10(y[id_non_const]), np.log10(x[id_non_const]), 
								frac=span, return_sorted = False)

	expected_variance = np.zeros(len(y))
	expected_variance[id_non_const] = 10**w
	for j in range(len(expected_variance)):
		if expected_variance[j] > np.sqrt(count_data.shape[0]):
			expected_variance[j] = np.sqrt(count_data.shape[0])

	standardized_count = np.zeros(shape = count_data.shape)
	for j in range(standardized_count.shape[1]):
		if id_non_const[j]:
			res =  (count_data[:, j] - x[j])/expected_variance[j]
			standardized_count[:, j] = res
		else:
			standardized_count[:, j] = count_data[:, j]

	vst_mean = x
	vst_var = np.var(standardized_count, axis = 0)
	gene_id = np.argsort(-1*vst_var)
	top_gene_id = gene_id[0:n_top_genes]

	return vst_mean, vst_var, top_gene_id


def hvg_before_after(raw_data, normalized_full_data, normalized_downsampled_data, save_dir=None):
	# convert pandas dataframe to numpy array
	raw_data = raw_data.values
	normalized_full_data = normalized_full_data.values
	normalized_downsampled_data = normalized_downsampled_data.values

	_, _, top_gene_id = vst_hvg(raw_data)
	mean_full = np.log2(np.asarray(np.mean(normalized_full_data, axis = 0)).squeeze())
	var_full = np.log2(np.asarray(np.var(normalized_full_data, axis = 0)).squeeze())
	mean_ds = np.log2(np.asarray(np.mean(normalized_downsampled_data, axis = 0)).squeeze())
	var_ds = np.log2(np.asarray(np.var(normalized_downsampled_data, axis = 0)).squeeze())

	# plotting
	fig = plt.figure(figsize = (8*2, 6*1))
	ax = fig.add_subplot(1, 2, 1)
	ax.scatter(mean_full, var_full, s = 2)
	ax.scatter(mean_full[top_gene_id], var_full[top_gene_id], s = 2, c = 'r')
	ax.set_title('normalized full data')
	ax.set_xlabel('log mean (full data)')
	ax.set_ylabel('log variance (full data)')

	ax = fig.add_subplot(1, 2, 2)
	ax.scatter(mean_ds, var_ds, s = 2)
	ax.scatter(mean_ds[top_gene_id], var_ds[top_gene_id], s = 2, c = 'r')
	ax.set_title('normalized downsampled data')
	ax.set_xlabel('log mean (full data)')
	ax.set_ylabel('log variance (full data)')

	if save_dir is not None:
		plt.savefig(save_dir)

