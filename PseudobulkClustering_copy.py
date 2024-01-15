import anndata as ad
import scanpy as sc
import utils
from abc import ABC, abstractmethod
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import metacells as mc
from scipy.stats import entropy

class PseudobulkClustererFactory:
    """Factory class to find and return the required pseudo bulk clustering algorithm. TODO Thid could be changed into an enum later.
    """

    @staticmethod
    def get_pseudobulk_clusterer(adata, method, **kwargs):
        """Return instance of a concrete PseudobulkClusterer, depending on the handed method parameter.

        Parameters
        ----------
        method : str
            String identifier used to select the pseudo bulk clustering method. TODO change this into a selector at a later point.
        **kwargs : object
            Object containing additional algorithms.
        """
        if method == 'metacells':
            return MetacellsClusterer(adata, **kwargs)
        else:
            raise ValueError(method)

class PseudobulkClusterer(ABC):

    @abstractmethod
    def cluster_pseudobulks(self):
        pass

class MetacellsClusterer(PseudobulkClusterer):

    def __init__(self, adata, **kwargs):
        self.adata = adata
        self.target_metacell_size = kwargs['target_metacell_size'] if kwargs['target_metacell_size'] else 100
        self.verbose = kwargs['verbose'] if kwargs['verbose'] else True
        self.out_dir = kwargs['out_dir'] if kwargs['out_dir'] else os.getcwd()

    def cluster_pseudobulks(self):
        #revert_log = self.adata.copy()
        #revert_log.X = np.expm1(revert_log.X)
        #cells = revert_log
        cells = self.adata

        dataset = self.adata.obs['dataset'].unique()[0] # TODO change this later, put saving, annotation of adata outside
        assert len(self.adata.obs['dataset'].unique()) == 1

        mc.ut.set_name(cells, "hca_bm.one-pass.preliminary.cells")

        mc.pipeline.exclude.exclude_genes(cells, random_seed=12345)
        #mc.pipeline.exclude.exclude_cells(cells, properly_sampled_min_cell_total=, properly_sampled_max_cell_total=, properly_sampled_max_excluded_genes_fraction=)
        #cells.obs['excluded_cell'] = False
        #mc.pipeline.exclude.extract_clean_data(cells)

        # actually, now they remove lateral genes, meaning genes that are cell-cycle genes, to avoid that the clustering decisions
        # made by the algorithm are influenced by cell states too much. We have to check if we have to perform this step, too.
        LATERAL_GENE_NAMES = []
        LATERAL_GENE_PATTERNS = []
        # This will mark as "lateral_gene" any genes that match the above, if they exist in the clean dataset.
        mc.pl.mark_lateral_genes(
            cells,
            lateral_gene_names=LATERAL_GENE_NAMES,
            lateral_gene_patterns=LATERAL_GENE_PATTERNS,
        )

        lateral_gene_mask = mc.ut.get_v_numpy(cells, "lateral_gene")
        lateral_gene_names = set(cells.var_names[lateral_gene_mask])

        # remove noisy genes
        NOISY_GENE_NAMES = []
        # This will mark as "noisy_gene" any genes that match the above, if they exist in the clean dataset.
        mc.pl.mark_noisy_genes(cells, noisy_gene_names=NOISY_GENE_NAMES)

        max_parallel_piles = mc.pl.guess_max_parallel_piles(cells)
        mc.pl.set_max_parallel_piles(max_parallel_piles)

        mc.pl.divide_and_conquer_pipeline(cells, random_seed=123456, target_metacell_size=self.target_metacell_size)

        cells.write_h5ad(os.path.join(self.out_dir, f'hlca_core_cells_{dataset}.h5ad'))

        metacells = mc.pl.collect_metacells(cells, name="hca_bm.one-pass.preliminary.metacells", random_seed=123456)
        print(f"Preliminary: {metacells.n_obs} metacells, {metacells.n_vars} genes")
                
        metacells.obs['dataset'] = dataset

        print(f'Removing {len(cells[cells.obs["metacell"] == -1])} unclustered cells:')
        cells = cells[cells.obs['metacell'] != -1]

        # map metacell index in cells anndata to actual metacell identifier
        self.map_metacell_index(cells, metacells)
        
        for index in ['ann_level_2', 'ann_level_3', 'ann_level_4', 'ann_level_5']:
            
            # compute overall entropy
            e = entropy(cells.obs[index].value_counts(normalize=True).values)
            result_df, shannon_entropy = self.compute_entropy_per_metacell(cells, index)
            maj = self.get_majority_vote(result_df)
            self.map_metacell_annotation(maj, metacells, index)

        metacells.write_h5ad(os.path.join(self.out_dir, f'hlca_core_metacells_{dataset}.h5ad'))

        self.adata = metacells

    def compute_entropy_per_metacell(self, cells, index):
        # pivot table such that we have ann_level_2 on rows, metacells on columns, and counts in cells
        result_df = cells.obs.pivot_table(index=index, columns='metacell', aggfunc='size', fill_value=0)
        result_df = result_df.div(result_df.sum(axis=0), axis=1)
        result_df.reset_index(inplace=True)
        result_df.index = result_df[index]
        result_df = result_df.drop(index, axis=1)
        
        # compute entropy and relative entropy compared to the overall entropy per metacell
        entropy_per_metacell = result_df.apply(lambda col: entropy(col, base=2), axis=0)
        shannon_entropy = entropy_per_metacell
        return result_df, shannon_entropy

    def get_majority_vote(self, result_df):
        return result_df.idxmax().to_dict()

    def map_metacell_index(self, cells, metacells):
        cells.obs['metacell'] = cells.obs['metacell'].map(lambda i: metacells.obs_names[i])

    def map_metacell_annotation(self, maj, metacells, index):
        metacells.obs[index] = metacells.obs.index.map(lambda i: maj[i])
