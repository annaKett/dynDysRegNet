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
        else:
            raise ValueError(method)


class PseudotimeInference(ABC):
    """Abstract pseudotime inference method."""

    @abstractmethod
    def infer_pseudotime(self):
        """Abstract batch correction method."""
        pass

    @staticmethod
    def plot_scatter(adata, save, title, basis, colors, show=False):
        fig, ax = plt.subplots()
        ax.set_title(title)
        sc.pl.scatter(
            adata,
            basis=basis,
            color=colors,
            color_map="gnuplot2",
            ax=ax)
        if save:
            fig.savefig(save, bbox_inches='tight')
        if show:
            plt.show()


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
        self.rep = kwargs['rep'] if kwargs['rep'] else 'X'

    def infer_pseudotime(self):
        """Infer diffusion pseudotime (DPT)."""
        #  first, get neighbor graph from specified embedding
        sc.pp.neighbors(self.adata, n_neighbors=15, n_pcs=self.num_pcs, use_rep=self.rep)
        
        # The width of the connectivity kernel is implicitly determined by the number of neighbors used to compute the single-cell graph in neighbors()
        sc.tl.diffmap(self.adata)

        if self.paga:
            self._apply_paga()

        stem = self.adata[self.adata.obs['ann_level_2']  == 'Hematopoietic stem cells']
        largest_stem = stem.obs['grouped'].idxmax()

        largest_stem = self.adata.obs.index.get_loc(largest_stem)
        self.adata.uns['iroot'] = largest_stem
        print(f'Root cell index for dpt inference: {self.adata.uns["iroot"]}')
        
        # use neighbors and diffmap to compute dpt pseudotime
        sc.tl.dpt(self.adata)


    def _apply_paga(self):
            sc.tl.draw_graph(self.adata) # uses connectivities from neighbors
        
            fig, ax = plt.subplots()
            ax.set_title('1: Force directed graph based on neighbors computed from diffmap colored by cell_type')
            sc.pl.draw_graph(self.adata, color='cell_type', legend_loc='on data', show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_1_fa_diffmap_neigh_graph.pdf'), bbox_inches='tight')

            sc.tl.leiden(self.adata) # uses connectivities from neighbors
            sc.tl.paga(self.adata, groups='leiden') # uses leiden clustering

            fig, ax = plt.subplots()
            ax.set_title('2.1: Force directed graph colored by cell_type clusters')
            sc.pl.draw_graph(self.adata, color=['cell_type'], legend_loc='on data', show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_21_fa_diffmap_neigh_graph.pdf'), bbox_inches='tight')

            fig, ax = plt.subplots()
            ax.set_title('2.2: Force directed graph colored by Leiden clusters')
            sc.pl.draw_graph(self.adata, color=['leiden'], legend_loc='on data', show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_22_fa_diffmap_neigh_graph.pdf'), bbox_inches='tight')

            fig, ax = plt.subplots()
            ax.set_title('2: Paga graph colored by Leiden clusters')
            sc.pl.paga(self.adata, color=['leiden'], show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_2_paga.pdf'), bbox_inches='tight')

            sc.tl.draw_graph(self.adata, init_pos='paga')

            fig, ax = plt.subplots()
            ax.set_title('3: Force directed graph colored by cell_type clusters')
            sc.pl.draw_graph(self.adata, color=['cell_type'], legend_loc='on data', show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_31_paga_draw_graph.pdf'), bbox_inches='tight')

            fig, ax = plt.subplots()
            ax.set_title('3: Force directed graph colored by Leiden clusters')
            sc.pl.draw_graph(self.adata, color=['leiden'], legend_loc='on data', show=False, ax=ax)
            fig.savefig(os.path.join(self.out_dir, self.out_afx + '_32_paga_draw_graph.pdf'), bbox_inches='tight')
