import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from abc import ABC, abstractmethod
import numpy as np
import os
import anndata as ad

class NetworkinferenceToolFactory:
    """Factory class to find and return the required network inference algorithm.
    """

    @staticmethod
    def get_networkinference_tool(adata, method, **kwargs):
        """Return instance of a concrete NetworkinferenceTool, depending on the handed method parameter.

        Parameters
        ----------
        method : str
            String identifier used to select the network inference method.
        **kwargs : object
            Object containing additional algorithms.
        """
        if method == 'grnboost2':
            return GRNBoost2InferenceTool(adata, **kwargs)
        else:
            raise ValueError(method)

class NetworkinferenceTool(ABC):

    @abstractmethod
    def infer_network(self):
        pass

class GRNBoost2InferenceTool(NetworkinferenceTool):
    
    def __init__(self, adata, **kwargs):
        self.adata = adata
        self.verbose = kwargs['verbose'] if kwargs['verbose'] else True
        self.out_dir = kwargs['out_dir'] if kwargs['out_dir'] else os.getcwd()
        self.tf_path = kwargs['tf_path']

    def infer_network(self):
        std = self.adata.X.toarray().std(axis=0)
        mean = self.adata.X.toarray().mean(axis=0)
        with open('excluded_cells.txt', 'a') as f:
            f.write(f'Mean: {mean}\n')
            f.write(f'Std: {std}\n')
        arr = (self.adata.X-mean)/std
        ex_matrix = pd.DataFrame(arr, columns=self.adata.var_names, index=self.adata.obs_names)
        tf_names = load_tf_names(self.tf_path)
        network = grnboost2(expression_data=ex_matrix, tf_names=tf_names)
        return network

