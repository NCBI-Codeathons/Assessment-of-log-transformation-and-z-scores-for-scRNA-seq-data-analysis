#Load libraries.

import matplotlib.pyplot as plt
import numpy as np
from pylab import *

#Accept as input a dictionary where each key is a normalization method.
#Value for the key is a list of adjusted rand scores comparing clusters with downsample levels 0.1 to 0.9 vs. the full data.

def adj_rand_index_plot(scores_vs_norm_method_dict):
	norm_methods = list(scores_vs_norm_method_dict.keys())

	splits_for_plot_forward_order = list(np.array(range(1,11))/10)
	splits_for_plot_reverse_order = list(np.array(range(1,11))/10)
	splits_for_plot_reverse_order.reverse()

	adj_rand_scores_method1 = scores_vs_norm_method_dict[norm_methods[0]]
	adj_rand_scores_method1 = adj_rand_scores_method1 + [1]
	adj_rand_scores_method1.reverse()

	ax = subplot(111)

	plt.plot(splits_for_plot_forward_order,adj_rand_scores_method1)
	plt.scatter(splits_for_plot_forward_order,adj_rand_scores_method1)

	plt.axis([0.1, 1, 0, 1])
	ax.set_xticklabels(splits_for_plot_reverse_order)
	plt.xlabel('Downsampling proportion')
	plt.ylabel('Adjusted rand score')

	for method in norm_methods:
		adj_rand_scores_this_method = scores_vs_norm_method_dict[method]
		adj_rand_scores_this_method = adj_rand_scores_this_method + [1]
		adj_rand_scores_this_method.reverse()

		plt.plot(splits_for_plot_forward_order,adj_rand_scores_this_method,label=method)
		plt.scatter(splits_for_plot_forward_order,adj_rand_scores_this_method)

	ax.legend()

	plt.show()

#Demo:

#adj_rand_scores_sctransform = [0.22684950098489198, 0.28455666308851585, 0.4112603869094376, 0.5483396194259293, 0.7209636382886637, 0.7314746314196946, 0.7318730262529415, 0.8373387272980577, 0.8704451666754582]
#adj_rand_scores_method2 = [i - 0.05 for i in adj_rand_scores_sctransform]
#adj_rand_index_plot({'sctransform':adj_rand_scores_sctransform,'dummy_method':adj_rand_scores_method2})
