import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scvi.dataset import Dataset10X
from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
import torch


# load data
data_dir, out_file = sys.argv[1], sys.argv[2]
gene_dataset = Dataset10X(save_path=os.path.expandusr(data_dir), measurement_names_column=1)

# set variables
save_path = '.'
n_epochs = 400
lr = 1e-3
use_batches = False
use_cuda = True

# define VAE
vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches)
trainer = UnsupervisedTrainer(
    vae,
    gene_dataset,
    train_size=0.75,
    use_cuda=use_cuda,
    frequency=5,
)

# train
trainer.train(n_epochs=n_epochs, lr=lr)
torch.save(trainer.model.state_dict(), '%s/vae.pkl' % save_path)

# extract normalized
full = trainer.create_posterior(trainer.model, local_10X_dataset, indices=np.arange(len(local_10X_dataset)))
normalized_values = full.sequential().get_sample_scale()
normalized_values = pd.DataFrame(normalized_values, index=data_df.index, columns=data_df.columns)

# save data
normalized_values.to_hdf(out_file, key='scvi_norm_df')


