from python_scripts.utils import remove_observables
import numpy as np
import pandas as pd
import os


def side_notes(adata, single_adata):
    """
    Side notes to know of adata object but can also be looked up in the summary created by 10x Genomics Spaceranger:
    1. No. spots under tissue
    2. Total No. genes
    3. Median genes per spot
    4. Total No. UMI counts
    5. Median No. UMI counts per Spot

    :param adata: [annData]
    :param single_adata: [bool] if False multiple adatas stored in one adata
    :return:
    """
    if single_adata:
        # All samples
        # get barcodes of each sample and therefore the number of spots under tissue for each sample
        model = adata.obs[['sample'] + []]
        batch_info = model.groupby('sample').groups.values()
        n_batches = np.array([len(v) for v in batch_info])
        print("Sorted No. spots under tissue for each sample: ", n_batches)
    else:
        unique_samples = np.unique(adata.obs['sample'])
        for sample in unique_samples:
            adata_sup = adata[adata.obs['sample'] == sample]
            print("\nSide notes of {} ".format(adata_sup.obs['sample'].values[1]))
            number_spots_under_tissue = len(
                np.where(adata_sup.obs['sample'].values == str(adata_sup.obs['sample'].values[1]))[0])

            print("No. spots under tissue: ", number_spots_under_tissue)
            # count number of expressed genes (count one gene over all spots)
            counts_gene = adata_sup[:number_spots_under_tissue].X.sum(0)
            counts_gene_sorted = np.sort(counts_gene)
            print("Total No. genes detected: ", np.count_nonzero(counts_gene_sorted))

            # Calculate median genes per spot
            copy_sample_1 = adata_sup[:number_spots_under_tissue].X.copy()
            mask = copy_sample_1 > 0
            zero_array = np.zeros_like(copy_sample_1)
            # count numbers of True == numbers of gene overall spots
            zero_array[mask] = 1
            median_genes_per_spot = np.median(zero_array.sum(1))
            median_umi_counts_per_spot = np.median(copy_sample_1.sum(1))
            print("Median genes per spot: ", median_genes_per_spot)
            print("Total No. UMI counts: ", sum(copy_sample_1.sum(1)))
            print("Median No. UMI counts per Spot: ", median_umi_counts_per_spot)

            # Second option: load from hdf5 files
            # If adata is read as h5 file otherwise remove todense() ..
            # copy_sample_1 = adata_obj[:number_spots_under_tissue].X.todense().copy()
            # mask = copy_sample_1 > 0
            # zero_array = np.zeros_like(copy_sample_1)
            # # count numbers of True == numbers of gene overall spots
            # zero_array[mask] = 1
            # median_genes_per_spot = np.median(zero_array.sum(1), axis=0)
            # median_umi_counts_per_spot = np.median(copy_sample_1.sum(1), axis=0)
            # print("Median genes per spot: ", median_genes_per_spot)
            # print("Total No. of UMI Counts: ", sum(copy_sample_1.sum(1)))
            # print("Median UMI Counts per Spot: ", median_umi_counts_per_spot)
    print("\n")


def minmax_spots_assigned_label(adata, output_folder):
    """Print labels which have the most and least assigned spots and save number of spots for each label in a data frame

    Parameters
    ----------
    adata : annData
    output_folder : str

    Returns
    -------

    """
    # remove cell cycle annotations from tissue_cell_labels list
    adata = remove_observables.remove_obs(
        adata=adata, colname=["G1_score", "G2M_score", "S_score", "M_score", 'ANNOTATOR'])

    # Use annotations from pathologist instead of clusters
    obs_keys = list(adata.obs_keys())

    # get all manual annotations by extracting all keys with upper case characters
    annotations = [char for char in obs_keys if any(c.isupper() for c in char)]

    spots_per_label = dict()
    for label in annotations[:-15]:
        spots_per_label[label] = len(adata.obs[label][adata.obs[label] == 1])

    # find label with highest and lowest number of spots
    max_spots_label = max(spots_per_label.keys(), key=(lambda k: spots_per_label[k]))
    print("The most assigned spots has the category ", max_spots_label)
    min_spots_label = min(spots_per_label.keys(), key=(lambda k: spots_per_label[k]))
    print("The least assigned spots has the category ", min_spots_label)

    # save num spots per label in df
    df = pd.DataFrame.from_dict(spots_per_label, orient='index', columns=['Number_spots'])

    df.to_csv(os.path.join(output_folder, "Number_of_spots_per_annotation.csv"))
    df.to_excel(os.path.join(output_folder, "Number_of_spots_per_annotation.xlsx"))
