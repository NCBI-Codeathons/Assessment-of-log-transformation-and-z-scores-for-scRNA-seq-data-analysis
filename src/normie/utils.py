import os
import sys
import scipy
import glob
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from scipy.stats import zscore
from glmpca.glmpca import glmpca


def load_data(data_dir):
	matrix_dir = data_dir
	mat = scipy.io.mmread(glob.glob(data_dir+'/*.mtx*')[0]).tocsr().T

	genes_path = glob.glob(data_dir+'/*genes*')[0]
	gene_names = pd.read_csv(genes_path, index_col=0, header=None, sep='\t').iloc[:, 0].tolist()
	
	barcodes_path = glob.glob(data_dir+'/*barcodes*')[0]
	valid_bc = pd.read_csv(barcodes_path, header=None, sep='\t').iloc[:, 0].tolist()

	# remove features not detected in all observations
	data_df = pd.DataFrame(mat.todense(), index=valid_bc, columns=gene_names)
	data_df = data_df.sort_index()

	print(data_df.shape)

	return data_df


def clean_data(data_df):
	binarized = (data_df > 0)

	# filter out 15-cell genes
	gene_cell_counts = binarized.sum(axis=0)
	use_genes = gene_cell_counts.index[gene_cell_counts > 15]
	data_df = data_df.loc[:, use_genes]

	# filter out 15-gene cells
	cell_gene_counts = binarized.sum(axis=1)
	use_cells = cell_gene_counts.index[cell_gene_counts > 15]
	data_df = data_df.loc[use_cells, :]

	# remove zero-columns
	data_df = data_df.loc[:, (data_df > 0).any(axis=0)]

	return data_df


def load_indexes(data_dir):
	genes_path = glob.glob(directory+'/*genes*')[0]
	gene_names = pd.read_csv(genes_path, index_col=0, header=None, sep='\t').iloc[:, 0].tolist()
	
	barcodes_path = glob.glob(directory+'/*barcodes*')[0]
	valid_bc = pd.read_csv(barcodes_path, header=None, sep='\t').iloc[:, 0].tolist()

	return valid_bc, gene_names


def run_median(counts):
	ls = counts.sum(axis = 1)
	norm_counts = counts.div(ls, axis=0).mul(np.median(ls), axis=0)
	
	return norm_counts


def run_median_log(counts, pc=0.1):
	ls = counts.sum(axis = 1)
	norm_counts = counts.div(ls, axis=0).mul(np.median(ls), axis=0)
	log_norm_counts = np.log2(norm_counts + pc) - np.log2(pc)
	
	return log_norm_counts


def run_median_log_z(counts, pc=0.1):
	ls = counts.sum(axis = 1)
	norm_counts = counts.div(ls, axis=0).mul(np.median(ls), axis=0)
	log_norm_counts = np.log2(norm_counts + pc) - np.log2(pc)
	z_log_norm_counts = zscore(log_norm_counts)
	
	return z_log_norm_counts


def run_median_log_lr(counts, pc=0.1):
	ls = counts.sum(axis = 1)
	norm_counts = counts.div(ls, axis=0).mul(np.median(ls), axis=0)
	log_norm_counts = np.log2(norm_counts + pc) - np.log2(pc)
	
	# fit linear regression to 
	lr = LinearRegression()
	log_ls = np.log2(pc) - np.log2(pc)
	lr.fit(log_ls, log_norm_counts)
	lr_log_norm_counts = log_norm_counts - lr.predict(log_ls)
	
	return lr_log_norm_counts


def run_glmpca(counts, n_latent=10, likelihood='poi'):
	glm_res = glmpca(counts, n_latent, fam=likelihood)
	glcpca_norm_counts = np.dot(glm_res['factors'], glm_res['loadings'].T)
	
	return glcpca_norm_counts


def run_scvi(data):
	return


def run_normalization(data_df, method, kwargs=None, r_dir=None):
	supported_methods = ['median', 'median_log', 'median_log_z', 'median_log_lr', 
						 'glmpca', 'scvi', 'scran', 'sctransform', 'linnorm']

	assert method.lower() in supported_methods, \
		'method is not one of the supported: {}'.format(', '.join(supported_methods))

	method_map = {'median': run_median,
				  'median_log': run_median_log,
				  'median_log_z': run_median_log_z,
				  'median_log_lr': run_median_log_lr,
				  'glmpca': run_glmpca,
				  'scvi': run_scvi,
				  'scran': 'scran_normalization.R',
				  'sctransform': 'sctransform_normalization.R', 
				  'linnorm': 'linnorm_normalization.R'}

	if method in supported_methods[:-3]:
		# python method, load data and do work directly
		counts = data_df
		if kwargs is not None:
			normed_data = method_map[method](counts, **kwargs)
		else:
			normed_data = method_map[method](counts)
	else:
		# R method, feed data dir to Rscript
		assert r_dir is not None, \
			'Rscripts required for running R method'

		r_dir = os.path.expanduser(r_dir)
		script_path = r_dir + method_map[method]

		# write the dataframe as a mtx for R to load
		input_file_path = r_dir+'tmp.mtx'
		scipy.io.mmwrite(input_file_path, data_df.values.T)

		os.system('RScript {} {} {}'.format(script_path, input_file_path, r_dir))

		barcodes, gene_names = data_df.index, data_df.columns

		# load the R output
		if method == 'scran':
			normed_data = io.mmread(r_dir + "_scran_normalized.mtx").tocsr().T
			normed_data = pd.DataFrame(normed_data.todense(), index=barcodes, columns=gene_names)
		elif method == 'sctransform':
			normed_data = pd.read_feather(r_dir + "_sctransform_normalized").tocsr().T
			normed_data = pd.DataFrame(normed_data.todense(), index=barcodes, columns=gene_names)
		elif method == 'linnorm':
			normed_data1 = io.mmread(r_dir + "_linnorm_with_dot_norm.mtx").tocsr().T
			normed_data2 = io.mmread(r_dir + "_linnorm.mtx").tocsr().T
			normed_data = {'with_dot_norm': pd.DataFrame(normed_data1.todense(), index=barcodes, columns=gene_names), 
						   'without_dot_norm': pd.DataFrame(normed_data2.todense(), index=barcodes, columns=gene_names)}

	return normed_data



