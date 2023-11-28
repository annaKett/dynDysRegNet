import anndata as ad
import scanpy as sc
import utils
from abc import ABC, abstractmethod
import pandas as pd
import matplotlib.pyplot as plt
import os

class BatchCorrectorFactory:
    """Factory class to find and return the required batch correction algorithm. TODO Thid could be changed into an enum later.
    """

    @staticmethod
    def get_batch_corrector(adata, method, **kwargs):
        """Return instance of a concrete BatchCorrector, depending on the handed method parameter.

        Parameters
        ----------
        method : str
            String identifier used to select batch correction method. TODO change this into a selector at a later point.
        **kwargs : object
            Object containing additional algorithms.
        """
        if method == 'harmony':
            return HarmonyBatchCorrector(adata, **kwargs)
        else:
            raise ValueError(method)


class BatchCorrector(ABC):
    """Abstract batch corrector."""

    @abstractmethod
    def correct_batch(self):
        """Abstract batch correction method."""
        pass

    @staticmethod
    def plot_pcs(adata, save=False, title='PC1, PC2', show=False):
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.suptitle(title)
        g = sc.pl.pca(adata, color='cell_type', show=show, ax=ax1, title='colored by cell_type')
        f = sc.pl.pca(adata, color='dataset', show=show, ax=ax2, title='colored by dataset')
        fig.savefig(save, bbox_inches='tight')


class HarmonyBatchCorrector(BatchCorrector):
    """Concrete class implementing batch correction with harmony. TODO insert source."""

    def __init__(self, adata, **kwargs):
        # check args and set default if argument is not available
        self.adata = adata
        self.num_pcs = kwargs['num_pcs'] if kwargs['num_pcs'] else 50
        self.data_dir = kwargs['data_dir'] if kwargs['data_dir'] else os.getcwd()
        self.verbose = kwargs['verbose'] if kwargs['verbose'] else True
        self.out_dir = kwargs['out_dir'] if kwargs['out_dir'] else os.getcwd()
        self.out_afx = kwargs['out_afx'] if kwargs['out_afx'] else 'unspecified'

    def correct_batch(self): #https://support.parsebiosciences.com/hc/en-us/articles/7704577188500-How-to-analyze-a-1-million-cell-data-set-using-Scanpy-and-Harmony
        """Apply harmony batch correction."""
        sc._settings.ScanpyConfig.n_jobs = 5
        utils.u_print(f'Computing PCA with {self.num_pcs} PCs as basis for batch effect correction.', verbose=self.verbose)
        sc.tl.pca(self.adata, n_comps=50)

        # plot before batch correction
        super().plot_pcs(self.adata, save=os.path.join(self.out_dir, self.out_afx+'_before_bc.pdf'), title='PC1 and PC2 embedding before batch correction')

        # apply batch correction; afterwards, use (batch corrected) X_pca_harmony as default pca embedding
        sc.external.pp.harmony_integrate(self.adata, key='dataset', basis='X_pca')
        utils.u_print('Storing batch corrected embedding in adata.obsm[\'X_pca\'].', verbose=self.verbose)
        utils.u_print(f'Storing adata with batch corrected PCA embedding in {self.data_dir}/{self.out_dir}/X_pca_from_harmony.h5ad.', verbose=self.verbose)
        self.adata.obsm['X_pca'] = self.adata.obsm['X_pca_harmony']

        # plot after batch correction
        super().plot_pcs(self.adata, save=os.path.join(self.out_dir, self.out_afx+'_after_bc.pdf'), title='PC1 and PC2 embedding after batch correction')


