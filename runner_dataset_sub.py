import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import os
from BatchCorrector import BatchCorrectorFactory
from PseudotimeInference import PseudotimeInferenceFactory
from PseudobulkClustering import PseudobulkClustererFactory
from NetworkinferenceTool import NetworkinferenceToolFactory

if __name__ == '__main__':
    sc.settings.verbosity = 0
    counts_file = 'hlca_core.h5ad'
    bc_method = 'harmony'
    filter_hvg = False
    paga = False
    subsets = []
    out_afx = 'hlca_core_metacells_all'
    data_dir = os.path.join('/nfs', 'data3', 'akett_data')
    
    out_dir = os.path.join('results')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_dir = os.path.join(data_dir, 'raw_data_scanvi')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    adata_all = sc.read_h5ad(os.path.join(data_dir, counts_file))

    # use raw counts
    adata_all.X = adata_all.raw.X.astype(np.float32)
    
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
        with open('excluded_cells.txt', 'a') as f:
            f.write(f'Will exclude {len(adata[~adata.obs["informative_ann_level_2"]])} cells in entropy check for dataset {dataset}.\n')
        adata = adata[adata.obs['informative_ann_level_2']]
        adata.obs['dataset'] = dataset

        subsets.append(adata)

    # concat metacells of different datasets into one adata and store
    metacells_all = ad.concat(subsets)
    with open('excluded_cells.txt', 'a') as f:
        f.write(f'Number of stem metacells: {len(metacells_all[metacells_all.obs["ann_level_2"] == "Hematopoietic stem cells"])}\n')
    with open('excluded_cells.txt', 'a') as f:
        f.write(f'Total number of cells after removing uninformative metacells regarding ann_level_2: {len(metacells_all.obs)}')    
    with open('excluded_cells.txt', 'a') as f:
        f.write(f'Number of metacell stem cells: {len(metacells_all.obs[metacells_all.obs["ann_level_2"] == "Hematopoietic stem cells"])}\n')
    
    metacells_all.obs_names = metacells_all.obs_names + metacells_all.obs['dataset'].astype(str)
    metacells_all.write_h5ad(os.path.join(out_dir, f'hlca_core_metacells_all.h5ad'))
    
    tmp = metacells_all

    # compute batch corrected PCA embedding
    if os.path.exists(os.path.join(out_dir, f'hlca_core_metacells_all_batch_corrected.h5ad')):
        metacells_all = sc.read_h5ad(os.path.join(out_dir, f'hlca_core_metacells_all_batch_corrected.h5ad'))
    else:
        bc = BatchCorrectorFactory().get_batch_corrector(
            metacells_all, 
            method='scanvi', 
            num_pcs=50, 
            data_dir=data_dir, 
            out_dir=out_dir, 
            out_afx=out_afx, 
            verbose=True
        ) 
        bc.correct_batch()
        metacells_all.write_h5ad(filename=os.path.join(out_dir, 'hlca_core_metacells_all_batch_corrected.h5ad'))
    
    # compute pseudo time from batch corrected PCA embedding
    if os.path.exists(os.path.join(out_dir, f'hlca_core_metacells_all_dpt.h5ad')):
        metacells_all = sc.read_h5ad(os.path.join(out_dir, f'hlca_core_metacells_all_dpt.h5ad'))
    else:
        pti = PseudotimeInferenceFactory().get_pseudotime_algorithm(
            metacells_all, 
            method='dpt', 
            num_pcs=None, 
            paga=paga, 
            data_dir=data_dir, 
            out_dir=out_dir, 
            out_afx=out_afx, 
            verbose=True, 
            rep='X_scANVI'
        )
        pti.infer_pseudotime()
        metacells_all.write_h5ad(filename=os.path.join(out_dir, 'hlca_core_metacells_all_dpt.h5ad'))

    tf_path = os.path.join(data_dir, 'human_known_tfs.txt')
    
    # continue with metacells and all genes
    metacells_all = tmp
    print(len(metacells_all.var_names))
    
    import scipy
    if scipy.sparse.issparse(metacells_all.X):
        dense_X = metacells_all.X.toarray()
    else:
        dense_X = metacells_all.X
    expr = pd.DataFrame(dense_X, index=metacells_all.obs_names, columns=metacells_all.var_names)
    
    network_cells = metacells_all[(metacells_all.obs['ann_level_2'] == 'Hematopoietic stem cells') |(metacells_all.obs['ann_level_2'] == 'Airway epithelium')]                                

    nit = NetworkinferenceToolFactory().get_networkinference_tool(network_cells, method='grnboost2', tf_path=tf_path, out_dir=out_dir, verbose=True)
    network = nit.infer_network()
    network = network.sort_values(by='importance', ascending=False)
    network.to_csv(os.path.join(out_dir, f'network_stem_airway_epithelium.csv'), sep='\t', index=False)
    adjacency_matrix = network.pivot(index='TF', columns='target', values='importance').fillna(0)
    adjacency_matrix.to_csv('network_stem_airway_epithelium_adjacency.csv')
    
        





