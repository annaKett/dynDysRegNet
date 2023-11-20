import anndata as ad

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

    def __init__(self, **kwargs):
        self.adata = adata
        self.num_pcs = kwargs['num_pcs'] if kwargs['num_pcs'] else 50
        self.data_dir = kwargs['data_dir'] if kwargs['data_dir'] else os.getcwd()
        self.verbose = kwargs['verbose'] if kwargs['verbose'] else True
        self.out_dir = kwargs['out_dir'] if kwargs['out_dir'] else os.getcwd()
        self.out_afx = kwargs['out_afx'] if kwargs['out_afx'] else 'unspecified'

    def infer_pseudotime(self): # TODO check again: does this touch the uncorrected data .X somehow?
        """Infer diffusion pseudotime (DPT)."""
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=self.num_pcs, use_rep='X_pca')
        sc.tl.diffmap(self.adata) # uses connectivities from neighbors
        sc.pp.neighbors(self.adata, n_neighbors=10, use_rep='X_diffmap') # recompute neighbors
        sc.tl.draw_graph(self.adata) # uses connectivities from neighbors
        f = sc.pl.draw_graph(self.adata, color='cell_type', legend_loc='on data')

        sc.tl.leiden(self.adata) # uses connectivities from neighbors
        sc.tl.paga(self.adata, groups='leiden') # uses leiden clustering
        sc.pl.paga(adata, color=['leiden']) # ""

        adata.uns['iroot'] = np.flatnonzero(adata.obs['ann_level_2']  == 'Hematopoietic stem cells')[0]
        sc.tl.dpt(self.adata) # uses neighbors and diffmap
        adata.write_h5ad(filename=os.path.join(out_dir, out_afx+'_dpt'))

class Paga(PseudotimeInference):
    """Concrete class implementing pseudotime inference with Paga. TODO insert source."""

    def __init__(self, **kwargs):
        pass

    def infer_pseudotime(self):
        """Infer Paga pseudotime."""
        pass