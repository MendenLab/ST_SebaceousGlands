import scanpy as sc
import numpy as np


def normalize_by_sizefactor(adata):
    """Normalising the count matrix using sizefactors
    The basic preprocessing includes assuming all size factors are equal
    (library size normalization to counts per million - CPM) and log-transforming the count data

    Parameters
    ----------
    adata : annData

    Returns
    -------
    adata : annData
        The count data has been normalized and log-transformed with an offset of 1.
        The offset of 1 ensures that zero counts map to zeros. Keep this data in '.raw' part of the AnnData object
        as it will be used to visualize gene expression and perform statistical tests such as computing marker genes
        for clusters

    """

    # Keep the count data in a counts layer
    adata.layers["counts"] = adata.X.copy()

    # Normalize adata
    adata.X /= adata.obs['size_factors'].values[:, None]

    # log-transforms and updates adata
    # and log or Box-Cox transformation (lambda=0)
    # because count distribution follows a power-law and afterwards a normal distribution -> easier to apply stat-tests
    sc.pp.log1p(adata)

    # modify resulting matrix
    adata.X = np.asarray(adata.X)

    # Store the full data set in 'raw' as log-normalised data for statistical testing
    adata.raw = adata

    return adata
