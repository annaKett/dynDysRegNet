import pandas as pd
import anndata as ad
import os
from NetworkinferenceTool import NetworkinferenceToolFactory

data_dir = os.path.join('/nfs', 'data3', 'akett_data')
out_dir = os.path.join(data_dir, 'raw_data_scanvi')
tf_path = os.path.join(data_dir, 'human_known_tfs.txt')

# important: from the output including the DPT values (output from the runner_dataset_sub.py), 
# you need to check which ann_level_2 categories have the lowest mean DPT. For template network inference, you use as few cells as
# as possible, belonging to a ann_level_2 type with low mean DPT. E.g., we used metacells of the types Hematopoietic stem cells and 
# Airway epithelium.

network_metacell_types = ['Hematopoietic stem cells', 'Airway epithelium']

# continue with metacells including all genes
metacells_all = ad.read_h5ad(os.path.join(out_dir, f'hlca_core_metacells_all.h5ad'))
print(len(metacells_all.var_names))

import scipy
if scipy.sparse.issparse(metacells_all.X):
    dense_X = metacells_all.X.toarray()
else:
    dense_X = metacells_all.X
expr = pd.DataFrame(dense_X, index=metacells_all.obs_names, columns=metacells_all.var_names)

# only use the selected metacell types for template network inference
network_cells = metacells_all[metacells_all.obs['ann_level_2'].isin(network_metacell_types)]                                

nit = NetworkinferenceToolFactory().get_networkinference_tool(network_cells, method='grnboost2', tf_path=tf_path, out_dir=out_dir, verbose=True)
network = nit.infer_network()
network = network.sort_values(by='importance', ascending=False)
network.to_csv(os.path.join(out_dir, f'network_stem_airway_epithelium.csv'), sep='\t', index=False)
adjacency_matrix = network.pivot(index='TF', columns='target', values='importance').fillna(0)
adjacency_matrix.to_csv('network_stem_airway_epithelium_adjacency.csv')

    



