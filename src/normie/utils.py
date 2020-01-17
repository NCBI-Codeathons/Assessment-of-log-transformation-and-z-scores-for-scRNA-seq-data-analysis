import os
import sys
import numpy as np
import pandas as pd
import scipy


def load_data(data_dir):
	return


def run_normalization(df, method):
	return


def run_median(counts):
    ls = counts.sum(axis = 1)
    norm_counts = counts.div(ls, axis=0).mul(np.median(ls), axis=0)
    
    return norm_counts


def run_median_log(counts, pc=0.1):
    ls = counts.sum(axis = 1)
    norm_counts = counts.div(ls, axis=0).mul(np.median(ls), axis=0)
    log_norm_counts = np.log2(norm_counts + pc) - np.log2(pc)
    
    return log_norm_counts


from scipy.stats import zscore
def run_median_log_z(counts, pc=0.1):
    ls = counts.sum(axis = 1)
    norm_counts = counts.div(ls, axis=0).mul(np.median(ls), axis=0)
    log_norm_counts = np.log2(norm_counts + pc) - np.log2(pc)
    z_log_norm_counts = zscore(log_norm_counts)
    
    return z_log_norm_counts


from sklearn.linear_model import LinearRegression
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


from glmpca.glmpca import glmpca
def run_glmpca(counts, n_latent=10, likelihood='poi'):
    glm_res = glmpca(counts, n_latent, fam=likelihood)
    glcpca_norm_counts = np.dot(glm_res['factors'], glm_res['loadings'].T)
    
    return glcpca_norm_counts


def run_scvi(data):
    return



