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

counts_file = 'hlca_core.h5ad'
bc_method = 'harmony'
filter_hvg = False
paga = False

flavor = 'seurat_v3'
data_dir = os.path.join('/nfs', 'data3', 'akett_data')
adata_all = sc.read_h5ad(os.path.join(data_dir, counts_file))
sc.settings.verbosity = 0

if not os.path.exists(os.path.join(data_dir, 'results')):
    os.mkdir(os.path.join(data_dir, 'results'))

for sub in ['Alveolar epithelium']:#set(adata_all.obs['ann_level_2'].values[0]):

    print(f'\n----------\nStarting processing of {sub} subset.')
    adata = adata_all[(adata_all.obs['ann_level_2'] == sub) | (adata_all.obs['ann_level_2'] == 'Hematopoietic stem cells')]

    out_dir = os.path.join('/nfs', 'data3', 'akett_data', 'results', f'subset_{sub}_bc_{bc_method}_hvg_{filter_hvg}_paga_{paga}')
    out_afx = f'{sub}_{bc_method}_paga_{paga}'

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if os.path.exists(os.path.join(out_dir, out_afx+'_batch_corrected')):
        print('Found batch corrected data.')
        adata = sc.read_h5ad(os.path.join(out_dir, out_afx+'_batch_corrected.h5ad'))
    else:
        bc = BatchCorrectorFactory().get_batch_corrector(adata, method='harmony', num_pcs=50, data_dir=data_dir, out_dir=out_dir, out_afx=out_afx, verbose=True)
        bc.correct_batch()
        adata.write_h5ad(filename=os.path.join(out_dir, out_afx+'_batch_corrected.h5ad'))

    pb = PseudobulkClustererFactory().get_pseudobulk_clusterer(adata, method='metacells')
    pb.cluster_pseudobulks()

    #pti = PseudotimeInferenceFactory().get_pseudotime_algorithm(adata, method='dpt', num_pcs=50, paga=paga, data_dir=data_dir, out_dir=out_dir, out_afx=out_afx, verbose=True)
    #pti.infer_pseudotime()


