import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
import random
from BatchCorrector import BatchCorrectorFactory
from PseudotimeInference import PseudotimeInferenceFactory
from PseudobulkClustering import PseudobulkClustererFactory
from NetworkinferenceTool import NetworkinferenceToolFactory

counts_file = 'hlca_core.h5ad'
bc_method = 'harmony'
filter_hvg = False
paga = False
subsets = []
out_afx = 'hlca_core_metacells_all'

data_dir = os.path.join('/nfs', 'data3', 'akett_data')
adata_all = sc.read_h5ad(os.path.join(data_dir, counts_file))
sc.settings.verbosity = 0

# use raw counts
adata_all.X = adata_all.raw.X.astype(np.float32)

out_dir = os.path.join(data_dir, 'raw_data_no_excl_ann2_13_new_nofilter')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# iterate through datasets to dodge batch effects in metacell computation
for dataset in adata_all.obs['dataset'].unique():
    with open('excluded_cells.txt', 'a') as f:
        f.write(f'\n{dataset}--------------------------\n')

    adata = adata_all[adata_all.obs['dataset'] == dataset]

    dataset = adata.obs['dataset'].unique()[0]
    assert len(adata.obs['dataset'].unique()) == 1
    adata.obs['dataset'] = dataset

    print(f'\n----------\nStarting processing of {dataset}.')
    if os.path.exists(os.path.join(out_dir, f'hlca_core_metacells_{dataset}.h5ad')):
        adata = sc.read_h5ad(os.path.join(out_dir, f'hlca_core_metacells_{dataset}.h5ad'))
    else:
        pb = PseudobulkClustererFactory().get_pseudobulk_clusterer(adata, method='metacells', target_metacell_size=13, out_dir=out_dir, verbose=False)
        adata = pb.cluster_pseudobulks()  

    # filter out uninformative metacells as per entropy comparison performed after metacell clustering
    adata = adata[adata.obs['informative_ann_level_2']]
    adata.obs['dataset'] = dataset

    subsets.append(adata)

# concat metacells of different datasets into one adata and store
metacells_all = ad.concat(subsets)
with open('excluded_cells.txt', 'a') as f:
    f.write(f'Number of stem metacells: {len(metacells_all[metacells_all.obs["ann_level_2"] == "Hematopoietic stem cells"])}\n')
    f.write(f'Total number of cells: {len(metacells_all.obs)}')
with open('excluded_cells.txt', 'a') as f:
    f.write(f'Total number of cells after removing uninformative metacells regarding ann_level_2: {len(metacells_all.obs)}')
metacells_all.write_h5ad(os.path.join(out_dir, f'hlca_core_metacells_all.h5ad'))

# transform gene fractions into pseudo counts by multiplying with arbitrary UMI count peer metacell
metacells_all.X = metacells_all.X * 100000 #np.dot(adata.X, np.diag(adata.obs['total_umis'].values))

# compute batch corrected PCA embedding
bc = BatchCorrectorFactory().get_batch_corrector(metacells_all, method='harmony', num_pcs=50, data_dir=data_dir, out_dir=out_dir, out_afx=out_afx, verbose=True)
bc.correct_batch()
metacells_all.write_h5ad(filename=os.path.join(out_dir, 'hlca_core_metacells_all_batch_corrected.h5ad'))

# compute pseudo time from batch corrected PCA embedding
pti = PseudotimeInferenceFactory().get_pseudotime_algorithm(metacells_all, method='dpt', num_pcs=50, paga=paga, data_dir=data_dir, out_dir=out_dir, out_afx=out_afx, verbose=True)
pti.infer_pseudotime()
metacells_all.write_h5ad(filename=os.path.join(out_dir, 'hlca_core_metacells_all_dpt.h5ad'))

tf_path = os.path.join(data_dir, 'human_known_tfs.txt')
network_cells = metacells_all[(metacells_all.obs['ann_level_2'] == 'Hematopoietic stem cells') |
                              (metacells_all.obs['ann_level_2'] == 'Mesothelium') |
                              (metacells_all.obs['ann_level_2'] == 'Lymphatic EC')]

nit = NetworkinferenceToolFactory().get_networkinference_tool(network_cells, method='grnboost2', tf_path=tf_path, out_dir=out_dir, verbose=True)
network = nit.infer_network()
network.to_csv(os.path.join(out_dir, f'network.csv'), sep='\t', index=False, header=False)





