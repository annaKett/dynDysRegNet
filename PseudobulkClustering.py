import os
import anndata as ad 
import abc

class PseudobulkClustering(abc):

    @abstractmethod
    def __init__(method='none'):
        self.method = method

    @abstractmethod
    def cluster_pseudobulks(adata):
        pass