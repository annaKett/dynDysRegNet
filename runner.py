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

counts_file = 'hlca_core.h5ad'
bc_method = 'harmony'
filter_hvg = False
flavor = 'seurat_v3'
data_dir = os.path.join('/nfs', 'data3', 'akett_data')
adata_all = sc.read_h5ad(os.path.join(data_dir, counts_file))
sc.settings.verbosity = 0

for sub in ['Alveolar epithelium']:#set(adata_all.obs['ann_level_2'].values[0]):

    print(f'\n----------\nStarting processing of {sub} subset.')
    subset = sub
    adata = adata_all[(adata_all.obs['ann_level_2'] == subset) | (adata_all.obs['ann_level_2'] == 'Hematopoietic stem cells')]
    st = adata[adata.obs['ann_level_2'] == 'Hematopoietic stem cells']
    print(f'Number of stem cells: {st}')

    ct = set(adata.obs['cell_type'].values)
    print(f'Present cell types (cell_type): {ct}')

    out_dir = os.path.join('/nfs', 'data3', 'akett_data', f'subset_{subset}_bc_{bc_method}_hvg_{filter_hvg}')
    out_afx = f'subset_{subset}_bc_{bc_method}_hvg_{filter_hvg}'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    print(f'out_dir={out_dir}, counts_file={counts_file}, bc_method={bc_method}, subset={subset}, filter_hvg={filter_hvg}')

    print(f'Start batch correction with {bc_method}, generate 50 PCs.\n')
    bc = BatchCorrectorFactory().get_batch_corrector(adata, method='harmony', num_pcs=50, data_dir=data_dir, out_dir=out_dir, out_afx=out_afx, verbose=True)
    bc.correct_batch()
    #adata.write_h5ad(filename=os.path.join(out_dir, out_afx+'_batch_corrected'))

    print(f'Start pseudotime inference with DPT.\n')
    pti = PseudotimeInferenceFactory().get_pseudotime_algorithm(adata, method='dpt', num_pcs=50, data_dir=data_dir, out_dir=out_dir, out_afx=out_afx, verbose=True)
    pti.infer_pseudotime()

    print(f'Start pseudo bulk clustering with TODO')




    ## ------
    #if os.path.exists(os.path.join(out_dir, out_afx+'_batch_corrected')):
        #print(f'Corrected adata exists at {os.path.join(out_dir, out_afx)}_batch_corrected.h5ad')
        #adata = sc.read_h5ad(os.path.join(out_dir, out_afx+'_batch_corrected'))
    #else:
        #if filter_hvg:
            # like in harmony paper https://www.nature.com/articles/s41592-019-0619-0 filter top 1000 hvg. BUT: harmony paper uses ranking by coefficient of variation, here, we use normaized variance
            #if flavor == 'seurat_v3':
                #sc.pp.highly_variable_genes(adata.raw, n_top_genes=1000, span=0.3, n_bins=20, flavor=flavor, batch_key=None)
            #else:
                #sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25, flavor=flavor, batch_key=None) # TODO
            #adata = adata[:,adata.var.highly_variable]
            #sc.pp.scale(adata, max_value=10) # scale to unit variance and zero mean?
            #print(f'Left only highly varaible genes: flavor={flavor}\n')