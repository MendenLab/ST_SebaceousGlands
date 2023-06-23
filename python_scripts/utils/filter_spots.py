import scanpy as sc


def qc_filter(adata, min_counts_gene, max_counts_gene, mt_threshold, min_genes, min_shared_counts,
              max_gene_count_threshold=False, min_counts_spots=None, max_counts_spots=None):
    """Filter count matrix according to identified and given QC thresholds --> sort out outliers

    Parameters
    ----------
    adata : annData
    min_counts_gene : int
    max_counts_gene : int
    mt_threshold : float
    min_genes : int
    min_shared_counts : int
    max_gene_count_threshold : bool
        determine if in spatial transcriptomics max gene count threshold is useful
    min_counts_spots : int
    max_counts_spots : int

    Returns
    -------
    adata : annData
        filtered adata object

    """
    # Filter spots:
    # Minimum number of counts required for a spot to pass filtering.
    sc.pp.filter_cells(adata, min_counts=min_counts_spots)
    print('Number of spots after min count filter: {:d}'.format(adata.n_obs))

    # Maximum number of counts required for a spot to pass filtering.
    sc.pp.filter_cells(adata, max_counts=max_counts_spots)
    print('Number of spots after max count max_counts_spots: {:d}'.format(adata.n_obs))

    # Filter out spots which have a low number of genes expressed
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print('Number of spots after gene filter: {:d}'.format(adata.n_obs))

    # --> we have bulk resolution
    if max_gene_count_threshold:
        # Filters out genes by maximum number of counts required for a gene to pass filtering
        # -> removes potential doublets
        sc.pp.filter_cells(adata, max_counts=max_counts_gene)

    # Filter MT-fraction:
    print('Total number of spots: {:d}'.format(adata.n_obs))
    # Threshold for MT-fraction is 20% to 25%
    # ==> filter out spots with an overall high proportion of dying or highly stressed cells
    adata = adata[adata.obs['mt_frac'] < mt_threshold]
    print('Number of spots after MT filter: {:d}'.format(adata.n_obs))

    # Filter genes:
    print('Total number of genes: {:d}'.format(adata.n_vars))

    # Minimum number of cells expressed required for a gene to pass filtering
    sc.pp.filter_genes(adata, min_cells=min_shared_counts)
    print('Number of genes after spots filter: {:d}'.format(adata.n_vars))

    # Min. 20 UMI-counts - filters out 0 count genes by minimum number of counts required for a gene to pass filtering.
    sc.pp.filter_genes(adata, min_counts=min_counts_gene)

    return adata
