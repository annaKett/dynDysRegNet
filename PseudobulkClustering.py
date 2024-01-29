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
    """Factory class to find and return the required batch correction algorithm.
    TODO This could be changed into an enum later.
    """

    @staticmethod
    def get_pseudobulk_clusterer(adata, method, **kwargs):
        """Return instance of a concrete PseudobulkClusterer, depending on the handed method parameter.

        Parameters
        ----------
        adata : AnnData
            Annotated data object containing single-cell RNA-seq data.
        method : str
            String identifier used to select the pseudo bulk clustering method. 
            Currently supports 'metacells'.
        **kwargs : object
            Object containing additional algorithms.

        Returns
        -------
        PseudobulkClusterer
            Instance of the selected PseudobulkClusterer.
        """
        if method == 'metacells':
            return MetacellsClusterer(adata, **kwargs)
        else:
            raise ValueError(method)

class PseudobulkClusterer(ABC):

    @abstractmethod
    def cluster_pseudobulks(self):
        """Abstract method to be implemented by subclasses for pseudobulk clustering."""
        pass

class MetacellsClusterer(PseudobulkClusterer):

    def __init__(self, adata, **kwargs):
        """Initialize MetacellsClusterer.

        Parameters
        ----------
        adata : AnnData
            Annotated data object containing single-cell RNA-seq data.
        **kwargs : dict
            Additional keyword arguments.
            - target_metacell_size: int, optional
              Target size for metacells. Default is 100.
            - verbose: bool, optional
              Verbosity flag. Default is True.
            - out_dir: str, optional
              Output directory. Default is the current working directory.
        """
        self.adata = adata
        self.target_metacell_size = kwargs['target_metacell_size'] if kwargs['target_metacell_size'] else 100
        self.verbose = kwargs['verbose'] if kwargs['verbose'] else True
        self.out_dir = kwargs['out_dir'] if kwargs['out_dir'] else os.getcwd()

    def cluster_pseudobulks(self):
        """Cluster pseudobulks using the Metacells algorithm.

        The process includes log reversion, quality control, metacells clustering,
        annotation assignment, and filtering based on entropy calculations.
        Resulting metacells data is saved to H5AD files.

        Returns
        -------
        AnnData
            Annotated data object containing metacells clustering results.
        """
        cells = self.adata.copy()

        mc.ut.set_name(cells, "hca_bm.one-pass.preliminary.cells")

        cells = self.set_lateral_genes(cells) # this doesn't do anything, but the according masks have to be set for the algorithm to work

        with open('excluded_cells.txt', 'a') as f:
            f.write(f'Starting clustering with number of stem cells: {len(cells.obs[cells.obs["ann_level_2"] == "Hematopoietic stem cells"])}\n')

        self.set_parallelity(cells)
        mc.pl.divide_and_conquer_pipeline(cells, target_metacell_size=13, random_seed=123456)
        
        # as metacells authors recommend:
        max_metacell = cells.obs['metacell'].max()
        cells.obs.loc[cells.obs['ann_level_2'] == 'Hematopoietic stem cells', 'metacell'] = max_metacell+1
        
        with open('excluded_cells.txt', 'a') as f:
            f.write(f'Masking {len(cells[cells.obs["metacell"] == -1])} unclustered cells with is_clustered=False in .obs:\n')
        
        metacells = mc.pl.collect_metacells(cells, name="hca_bm.iteration-1.metacells", random_seed=123456)

        with open('excluded_cells.txt', 'a') as f:
            f.write(f"{metacells.n_obs} metacells, {metacells.n_vars} genes\n")
        
        # remove the 'unclustered'-metacell and map metacell index in cells anndata to actual metacell identifier
        with open('excluded_cells.txt', 'a') as f:
            f.write(f'Removing {len(cells[cells.obs["metacell"] == -1])} unclustered cells.\n')
        cells = cells[cells.obs['metacell'] != -1]
        cells = self.map_metacell_index(cells, metacells)

        result_df, shannon_entropy = self.compute_entropy_per_metacell(cells, 'ann_level_2')
        maj = self.get_majority_vote(result_df)
        metacells = self.map_metacell_annotation(maj, metacells, 'ann_level_2')

        
        for index_before, index in zip([None, 'ann_level_1', 'ann_level_2', 'ann_level_3', 'ann_level_4'], ['ann_level_1', 'ann_level_2', 'ann_level_3', 'ann_level_4', 'ann_level_5']):
            # remember original annotation
            cells.obs['original_' + index] = cells.obs[index].copy()

            e = np.inf
            if index_before:
                # inherit annotation from coarser annotation level if annotation is None on current annotation level
                cells = self.inherit_level_annotation_if_none(cells, index_before, index)
                # compute overall entropy from original data (don't filter out unclustered cells)
                e = entropy(cells.obs[index].value_counts(normalize=False).values)
            else:
                e = entropy(cells.obs[index].value_counts(normalize=False).values)

            with open('excluded_cells.txt', 'a') as f:
                f.write(f'Shannon entropy of {index} variable across all metacells: {e}\n')
            
            # compute entropy per metacell
            result_df, shannon_entropy = self.compute_entropy_per_metacell(cells, index)
            maj = self.get_majority_vote(result_df)
            metacells = self.map_metacell_annotation(maj, metacells, index)

            # mask metacells that are less informative on the current level than the unclustered cells
            is_informative = shannon_entropy <= e
            metacells.obs['informative_' + index] = None
            metacells.obs = metacells.obs.astype({'informative_' + index: bool})
            metacells.obs.loc[is_informative, 'informative_' + index] = True
            metacells.obs.loc[~is_informative, 'informative_' + index] = False
        
        # save both cells and metacells anndata
        dataset = cells.obs['dataset'].unique()[0]
        assert len(cells.obs['dataset'].unique()) == 1

        metacells.write_h5ad(os.path.join(self.out_dir, f'hlca_core_metacells_{dataset}.h5ad'))
        cells.write_h5ad(os.path.join(self.out_dir, f'hlca_core_cells_{dataset}.h5ad'))

        # return metacells
        self.adata = metacells
        return self.adata

    def set_lateral_genes(self, cells):
        """Mark lateral genes and remove noisy genes from the given cell data. Does nothing, only called to set empty masks on cells data.

        Parameters
        ----------
        cells : AnnData
            Annotated data object containing single-cell RNA-seq data.
        """
        # they now would remove lateral genes, meaning genes that are cell-cycle genes, to avoid that the clustering decisions
        # made by the algorithm are influenced by cell states too much. We have to check if we have to perform this step, too.
        LATERAL_GENE_NAMES = []
        LATERAL_GENE_PATTERNS = []
        # This will mark as "lateral_gene" any genes that match the above, if they exist in the clean dataset.
        mc.pl.mark_lateral_genes(
            cells,
            lateral_gene_names=LATERAL_GENE_NAMES,
            lateral_gene_patterns=LATERAL_GENE_PATTERNS,
        )

        NOISY_GENE_NAMES = []
        # This will mark as "noisy_gene" any genes that match the above, if they exist in the clean dataset.
        mc.pl.mark_noisy_genes(cells, noisy_gene_names=NOISY_GENE_NAMES)
        return cells

    def set_parallelity(self, cells):
        """Set the maximum number of parallel piles for the metacells algorithm.

        Parameters
        ----------
        cells : AnnData
            Annotated data object containing single-cell RNA-seq data.
        """
        max_parallel_piles = mc.pl.guess_max_parallel_piles(cells)
        mc.pl.set_max_parallel_piles(max_parallel_piles)

    def compute_entropy_per_metacell(self, cells, index):
        """Compute entropy and relative entropy per metacell for a specified annotation index.

        Parameters
        ----------
        cells : AnnData
            Annotated data object containing single-cell RNA-seq data.
        index : str
            Annotation index for which entropy is computed.

        Returns
        -------
        Tuple
            Tuple containing a DataFrame with entropy calculations and a Series
            representing the Shannon entropy.
        """
        # pivot table such that we have ann_level_2 on rows, metacells on columns, and counts in cells
        result_df = cells.obs.pivot_table(index=index, columns='metacell', aggfunc='size', fill_value=0)
        result_df.reset_index(inplace=True)
        result_df.index = result_df[index]
        result_df = result_df.drop(index, axis=1)
        
        # compute entropy and relative entropy compared to the overall entropy per metacell
        shannon_entropy = result_df.apply(lambda col: entropy(col, base=2), axis=0)
        return result_df, shannon_entropy

    def inherit_level_annotation_if_none(self, cells, level_before, cur_level):
        """Inherit level annotation from a coarser level if annotation is None on the current level.

        Parameters
        ----------
        cells : AnnData
            Annotated data object containing single-cell RNA-seq data.
        level_before : str
            Annotation level before the current level.
        cur_level : str
            Current annotation level.
        """
        cells.obs.loc[cells.obs[cur_level] == 'None', cur_level] = cells.obs[level_before].astype(cells.obs[cur_level].dtype)
        return cells

    def get_majority_vote(self, result_df):
        """Compute majority vote for each metacell based on the given result DataFrame.

        Parameters
        ----------
        result_df : DataFrame
            DataFrame containing metacell data.

        Returns
        -------
        dict
            Dictionary mapping metacell indices to majority votes.
        """
        return result_df.idxmax().to_dict()

    def map_metacell_index(self, cells, metacells):
        """Map metacell indices in cells AnnData to actual metacell identifiers.

        Parameters
        ----------
        cells : AnnData
            Annotated data object containing single-cell RNA-seq data.
        metacells : AnnData
            Annotated data object containing metacell clustering results.
        """
        max_ = cells.obs['metacell'].max()
        cells.obs['metacell'] = cells.obs['metacell'].map(lambda i: metacells.obs_names[i])
        return cells


    def map_metacell_annotation(self, maj, metacells, index):
        """Map metacell annotations based on majority votes.

        Parameters
        ----------
        maj : dict
            Dictionary containing majority votes for metacells.
        metacells : AnnData
            Annotated data object containing metacell clustering results.
        index : str
            Annotation index for mapping.
        """
        metacells.obs[index] = metacells.obs.index.map(lambda i: maj[i])
        return metacells
