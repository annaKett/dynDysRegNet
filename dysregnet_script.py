import dysregnet
import pandas as pd
import os
import anndata as ad
import scipy
import numpy as np

conCol = 'condition'
CatCov = []#['dataset']
ConCov = []#['dpt_pseudotime']

network_path = os.path.join('/nfs', 'data3', 'akett_data', 'raw_data_scanvi')
apx= '_'
output_file_path = os.path.join(network_path, f'dysregnet_counts_{apx}.csv')

expr = ad.read_h5ad(os.path.join(network_path, 'hlca_core_metacells_all.h5ad'))
print(f'LEN expr: {len(expr)}')
pheno = ad.read_h5ad(os.path.join(network_path, 'hlca_core_metacells_all_dpt.h5ad')).obs.copy()
print(f'LEN pheno: {len(pheno)}')
assert not pheno['dpt_pseudotime'].isna().any()
assert not pheno['dataset'].isna().any()

if scipy.sparse.issparse(expr.X):
    dense_X = expr.X.toarray()
else:
    dense_X = expr.X

expr = pd.DataFrame(dense_X, index=expr.obs_names, columns=expr.var_names)

assert expr.index.is_unique
assert pheno.index.is_unique
print(f'HEYO length expr genes: {len(expr.columns)}')

pheno['condition'] = 0
pheno.loc[~pheno['ann_level_2'].isin(['Airway epithelium', 'Hematopoietic stem cells']), 'condition'] = 1

pheno['sample_id'] = pheno.index
expr['sample_id'] = expr.index

control_expr = expr[expr['sample_id'].isin(pheno[pheno['condition'] == 0]['sample_id'])]
case_expr = expr[expr['sample_id'].isin(pheno[pheno['condition'] == 1]['sample_id'])]

# Calculate mean and std for control samples
mean_control = control_expr.drop(columns=['sample_id']).mean()
std_control = control_expr.drop(columns=['sample_id']).std()

control_expr_standardized = (control_expr.drop(columns=['sample_id']).sub(mean_control)).div(std_control)
case_expr_standardized = (case_expr.drop(columns=['sample_id']).sub(mean_control)).div(std_control)

standardized_expr = pd.concat([control_expr_standardized, case_expr_standardized])
standardized_expr['sample_id'] = standardized_expr.index

standardized_expr= standardized_expr.set_index('sample_id').reset_index()

for con in ConCov:
    mean_control_pheno = pheno.loc[pheno['condition'] == 0, con].mean()
    std_control_pheno = pheno.loc[pheno['condition'] == 0, con].std()

    pheno.loc[pheno['condition'] == 0, con] = (pheno.loc[pheno['condition'] == 0, con] - mean_control_pheno) / std_control_pheno
    pheno.loc[pheno['condition'] == 1, con] = (pheno.loc[pheno['condition'] == 1, con] - mean_control_pheno) / std_control_pheno

standardized_pheno = pheno
standardized_pheno['sample_id'] = standardized_pheno.index
standardized_pheno= standardized_pheno.set_index('sample_id').reset_index()

template_network = pd.read_csv(os.path.join(network_path, 'network_stem_airway_epithelium.csv'), header=0, index_col=None, sep='\t')
template_network = template_network.sort_values(by='importance', ascending=False).head(100000)
template_network = template_network[['TF', 'target']]

out = dysregnet.run(expression_data=standardized_expr,
                    meta=standardized_pheno,
                    GRN=template_network,
                    conCol=conCol,
                    CatCov=CatCov,
                    ConCov=ConCov,
                    direction_condition=False)

df = out.get_results()
model_stats = out.get_model_stats()

df.to_csv(os.path.join(network_path,f'model_results_{apx}.csv'), index=True, header=True)
model_stats.to_csv(os.path.join(network_path,f'model_stats_{apx}.csv'), index=True, header=True)

#df = pd.DataFrame(columns=['pos_higher_than_expected', 'pos_lower_than_expected', 'neg_higher_than_expected', 'neg_lower_than_expected', 'num_repressed', 'num_activated'], index=df.index)
df['pos_higher_than_expected'] = 0
df['pos_lower_than_expected'] = 0
df['neg_higher_than_expected'] = 0
df['neg_lower_than_expected'] = 0
df['num_repressed'] = 0
df['num_activated'] = 0
    
for coef in ['coef_TF']:
    positive_edges = model_stats[model_stats[coef] > 0].index.tolist()
    negative_edges = model_stats[model_stats[coef] < 0].index.tolist()

    pos_higher_than_expected_vals = df[positive_edges].apply(lambda row: sum(row[edge] > 0 for edge in positive_edges), axis=1)
    pos_lower_than_expected_vals = df[positive_edges].apply(lambda row: sum(row[edge] < 0 for edge in positive_edges), axis=1)
    neg_higher_than_expected_vals = df[negative_edges].apply(lambda row: sum(row[edge] > 0 for edge in negative_edges), axis=1)
    neg_lower_than_expected_vals = df[negative_edges].apply(lambda row: sum(row[edge] < 0 for edge in negative_edges), axis=1)

    num_repressed_vals = neg_lower_than_expected_vals + pos_lower_than_expected_vals
    num_activated_vals = pos_higher_than_expected_vals + neg_higher_than_expected_vals

    df['pos_higher_than_expected'] += pos_higher_than_expected_vals
    df['pos_lower_than_expected'] += pos_lower_than_expected_vals
    df['neg_higher_than_expected'] += neg_higher_than_expected_vals
    df['neg_lower_than_expected'] += neg_lower_than_expected_vals
    df['num_repressed'] += num_repressed_vals
    df['num_activated'] += num_activated_vals

df[['pos_higher_than_expected', 'pos_lower_than_expected', 'neg_higher_than_expected', 'neg_lower_than_expected', 'num_repressed', 'num_activated']].to_csv(os.path.join(network_path, f'num_act_repr_{apx}.csv'))
"""    
for coef in CatCov:
    # Filter edges based on their coefficients in the model stats
    positive_edges = model_stats[model_stats[coef] > 0].index.tolist()
    negative_edges = model_stats[model_stats[coef] < 0].index.tolist()

    # Calculate values for all cases without immediately adding them to the DataFrame
    pos_higher_than_expected_vals = df[positive_edges].apply(lambda row: sum(row[edge] > 0 for edge in positive_edges), axis=1)
    pos_lower_than_expected_vals = df[positive_edges].apply(lambda row: sum(row[edge] < 0 for edge in positive_edges), axis=1)
    neg_higher_than_expected_vals = df[negative_edges].apply(lambda row: sum(row[edge] > 0 for edge in negative_edges), axis=1)
    neg_lower_than_expected_vals = df[negative_edges].apply(lambda row: sum(row[edge] < 0 for edge in negative_edges), axis=1)

    # Calculate 'num_repressed' and 'num_activated' based on the values computed above
    num_repressed_vals = neg_lower_than_expected_vals + pos_lower_than_expected_vals
    num_activated_vals = pos_higher_than_expected_vals + neg_higher_than_expected_vals

    # Now, assign the computed values to the DataFrame to ensure these fields are not included in any subsequent computations
    df['pos_higher_than_expected'] += pos_higher_than_expected_vals
    df['pos_lower_than_expected'] += pos_lower_than_expected_vals
    df['neg_higher_than_expected'] += neg_higher_than_expected_vals
    df['neg_lower_than_expected'] += neg_lower_than_expected_vals
    df['num_repressed'] += num_repressed_vals
    df['num_activated'] += num_activated_vals

  """  
