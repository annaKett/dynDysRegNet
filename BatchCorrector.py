import anndata as ad
import scanpy as sc
import utils
from abc import ABC, abstractmethod
import pandas as pd
import matplotlib.pyplot as plt
import os
import scvi
import torch
from scvi.model.utils import mde

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
        elif method == 'scanvi':
            return SCANVIBatchCorrector(adata, **kwargs)
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

        # apply batch correction; afterwards, use (batch corrected) X_pca_harmony as default pca embedding
        sc.external.pp.harmony_integrate(self.adata, key='dataset', basis='X_pca')
        utils.u_print('Storing batch corrected embedding in adata.obsm[\'X_pca\'].', verbose=self.verbose)
        utils.u_print(f'Storing adata with batch corrected PCA embedding in {self.data_dir}/{self.out_dir}/X_pca_from_harmony.h5ad.', verbose=self.verbose)
        self.adata.obsm['X_pca'] = self.adata.obsm['X_pca_harmony']


class SCANVIBatchCorrector(BatchCorrector):
    """Concrete class implementing batch correction with scANVI. TODO insert source."""

    def __init__(self, adata, **kwargs):
        # check args and set default if argument is not available
        self.adata = adata
        self.num_pcs = kwargs['num_pcs'] if kwargs['num_pcs'] else 50
        self.data_dir = kwargs['data_dir'] if kwargs['data_dir'] else os.getcwd()
        self.verbose = kwargs['verbose'] if kwargs['verbose'] else True
        self.out_dir = kwargs['out_dir'] if kwargs['out_dir'] else os.getcwd()
        self.out_afx = kwargs['out_afx'] if kwargs['out_afx'] else 'unspecified'

    def correct_batch(self): #https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/harmonization.html
        # pretrained model: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/query_hlca_knn.html
        # training a new model here
        """Apply scANVI batch correction."""
        scvi.settings.seed = 12345
        sc.pp.neighbors(self.adata)
        sc.tl.umap(self.adata)
        sc.pl.umap(
            self.adata, 
            title='UMAP embedding colored by dataset before applying batch correction', 
            color=['dataset', 'ann_level_2'],
        )
        sc.set_figure_params(figsize=(4, 4))
        torch.set_float32_matmul_precision("high")
        sc.pp.highly_variable_genes(
            self.adata,
            flavor="seurat_v3",
            n_top_genes=2000,
            batch_key="dataset",
            subset=True,
        )
        with open('excluded_cells.txt', 'a') as f:
            f.write('Highly variable genes scanvi:\n')
            f.write(len(self.adata.var[self.adata.var['highly_variable']]))
        if os.path.exists(os.path.join(self.out_dir, 'raw_hlca_core.pt')):
            scanvi_model = scvi.load(os.path.join(self.out_dir, 'raw_hlca_core_scanvi_model.pt'), self.adata)
        else:
            scvi.model.SCVI.setup_anndata(self.adata, batch_key="dataset")
            model = scvi.model.SCVI(self.adata, n_layers=2, n_latent=30, gene_likelihood="nb")
            model.train()
            self.adata.obsm["X_scVI"] = model.get_latent_representation()
            self.adata.obsm["X_mde"] = mde(self.adata.obsm["X_scVI"])
            scanvi_model = scvi.model.SCANVI.from_scvi_model(
                model,
                adata=self.adata,
                labels_key="ann_level_2",
                unlabeled_category="None",
            )
            scanvi_model.train(max_epochs=20, n_samples_per_label=100)
            #scanvi_model.save_state_dict(os.path.join(self.out_dir, 'raw_hlca_core.pt'))
            #torch.save(scanvi_model.state_dict(), os.path.join(self.out_dir, 'raw_hlca_core_scanvi_model.pt'))
            scanvi_model.save(os.path.join(self.out_dir, 'raw_hlca_core_scanvi_model.pt'))

        SCANVI_LATENT_KEY = "X_scANVI"
        self.adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(self.adata)
        SCANVI_MDE_KEY = "X_scANVI_MDE"
        self.adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(self.adata.obsm[SCANVI_LATENT_KEY])
        sc.pl.embedding(
            self.adata, 
            basis=SCANVI_MDE_KEY, 
            color=["ann_level_2"], 
            ncols=1, 
            frameon=False
        )
  
        

