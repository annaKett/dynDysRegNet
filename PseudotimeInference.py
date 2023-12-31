import anndata as ad
import scanpy as sc
import utils
from abc import ABC, abstractmethod
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

class PseudotimeInferenceFactory:
    """Factory class to find and return the required pseudo time inference algorithm. TODO Thid could be changed into an enum later.
    """

    @staticmethod
    def get_pseudotime_algorithm(adata, method, **kwargs):
        """Return instance of a concrete BatchCorrector, depending on the handed method parameter.

        Parameters
        ----------
        method : str
            String identifier used to select the pseudotime inference algorithm. TODO change this into a selector at a later point.
        **kwargs : object
            Object containing additional algorithms.
        """
        if method == 'dpt':
            return DPT(adata, **kwargs)
        elif method == 'paga':
            return Paga(adata, **kwargs)
        else:
            raise ValueError(method)


class PseudotimeInference(ABC):
    """Abstract pseudotime inference method."""

    @abstractmethod
    def infer_pseudotime(self):
        """Abstract batch correction method."""
        pass


class DPT(PseudotimeInference):
    """Concrete class implementing pseudotime inference with DPT. TODO insert source."""

    def __init__(self, adata, **kwargs):
        self.adata = adata
        self.num_pcs = kwargs['num_pcs'] if kwargs['num_pcs'] else 50
        self.data_dir = kwargs['data_dir'] if kwargs['data_dir'] else os.getcwd()
        self.verbose = kwargs['verbose'] if kwargs['verbose'] else True
        self.out_dir = kwargs['out_dir'] if kwargs['out_dir'] else os.getcwd()
        self.out_afx = kwargs['out_afx'] if kwargs['out_afx'] else 'unspecified'
        self.paga = kwargs['paga'] if kwargs['paga'] else False

    def infer_pseudotime(self):
        """Infer diffusion pseudotime (DPT)."""
        #  first, use paga to infer a graph of cell paths
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=self.num_pcs, use_rep='X_pca')
        sc.tl.diffmap(self.adata) # The width of the connectivity kernel is implicitly determined by the number of neighbors used to compute the single-cell graph in neighbors()
        sc.pp.neighbors(self.adata, n_neighbors=10, use_rep='X_diffmap') # compute neigbors based on diffmap representation
        
        if self.paga:
            sc.tl.draw_graph(self.adata) # uses connectivities from neighbors
        
            fig, ax = plt.subplots()
            ax.set_title('1: Force directed graph based on neighbors computed from diffmap colored by cell_type')
            sc.pl.draw_graph(self.adata, color='cell_type', legend_loc='on data', show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_1_fa_diffmap_neigh_graph.pdf'), bbox_inches='tight')

            sc.tl.leiden(self.adata) # uses connectivities from neighbors
            sc.tl.paga(self.adata, groups='leiden') # uses leiden clustering

            fig, ax = plt.subplots()
            ax.set_title('2.1&2.2: Force directed graph colored by cell_type and Leiden clusters')
            sc.pl.draw_graph(self.adata, color=['cell_type', 'leiden'], legend_loc='on data', show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_21_fa_diffmap_neigh_graph.pdf'), bbox_inches='tight')

            fig, ax = plt.subplots()
            ax.set_title('2: Paga graph colored by Leiden clusters')
            sc.pl.paga(self.adata, color=['leiden'], show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_2_paga.pdf'), bbox_inches='tight')

            sc.tl.draw_graph(self.adata, init_pos='paga')

            fig, ax = plt.subplots()
            ax.set_title('3: Force directed graph colored by Leiden clusters and cell_type')
            sc.pl.draw_graph(self.adata, color=['leiden', 'cell_type'], legend_loc='on data', show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_3_paga_draw_graph.pdf'), bbox_inches='tight')

        # then, use dpt to link paga paths with a pseudotime
        self.adata.uns['iroot'] = np.flatnonzero(self.adata.obs['ann_level_2']  == 'Hematopoietic stem cells')[0]
        sc.tl.dpt(self.adata) # uses neighbors and diffmap

        if self.paga:
            fig, ax = plt.subplots()
            ax.set_title('4: Force directed graph colored by pseudo time')
            sc.pl.draw_graph(self.adata, color=['dpt_pseudotime'], legend_loc='on data', show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_4_dpt.pdf'), bbox_inches='tight')

        fig, ax = plt.subplots()
        ax.set_title('5: PCA embedding colored by pseudo time')
        sc.pl.scatter(
            self.adata,
            basis="pca",
            color=["dpt_pseudotime"],
            color_map="gnuplot2",
            title='pca_dpt_pseudotime',
            ax=ax)
        fig.savefig(os.path.join(self.out_dir, self.out_afx + '_5_dpt.pdf'), bbox_inches='tight')

        fig, ax = plt.subplots()
        ax.set_title('6: Diffmap colored by pseudo time')
        sc.pl.scatter(
            self.adata,
            basis="diffmap",
            color=["dpt_pseudotime"],
            color_map="gnuplot2",
            components=[2, 3],
            title='diffmap_dpt_pseudotime',
            ax=ax)
        fig.savefig(os.path.join(self.out_dir, self.out_afx + '_6_dpt.pdf'), bbox_inches='tight')


        fig, ax = plt.subplots()  
        ax.set_title('7: Diffmap colored by cell_type') 
        sc.pl.scatter(
            self.adata,
            basis="diffmap",
            color=["cell_type"],
            color_map="gnuplot2",
            components=[2, 3],
            title='cell_type_diffmap_dpt_pseudotime',
            ax=ax)
        fig.savefig(os.path.join(self.out_dir, self.out_afx + '_7_dpt.pdf'), bbox_inches='tight')

        #sc.pl.dpt_groups_pseudotime(self.adata)
        #sc.pl.dpt_timeseries(self.adata)

        self.adata.write_h5ad(filename=os.path.join(self.out_dir, self.out_afx+'_dpt.h5ad'))
