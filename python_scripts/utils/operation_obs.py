import numpy as np
import pandas as pd
import math
import copy
from operator import itemgetter


def remove_obs(adata, colname):
    # check if all columns are keys in adata.obs
    obs = list(set(adata.obs_keys()) & set(colname))
    adata.obs = adata.obs.drop(columns=obs)
    return adata


def orderofmagnitude(number):
    return math.floor(math.log(number, 10))


def map_sample_batch_list(adata, num_samples_patient):
    """

    :param adata:
    :param num_samples_patient: [list] [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2]
    :return:
    """

    metadata_batch = pd.DataFrame(columns=["sample", "batch"])
    metadata = dict()
    list_project_samples = dict()
    for ind, project in enumerate(np.unique(adata.obs['project'])):
        list_project_samples[project] = np.unique(adata.obs['sample'][adata.obs['project'] == project])
        sample_names = [int(i.split('_')[1]) for i in list_project_samples[project]]

        temp_num_oom = 10 ** orderofmagnitude(np.amax(sample_names))
        max_samplenumber = np.amax(sample_names) - temp_num_oom

        rest = max_samplenumber % num_samples_patient[ind]
        if rest != 0:
            max_samplenumber = max_samplenumber + rest

        # determine no batches
        num_batches = 0
        list_batches = num_samples_patient.copy()
        number_samples = copy.copy(max_samplenumber)
        while True:
            list_batches.remove(num_samples_patient[num_batches])
            if ((number_samples - num_samples_patient[num_batches]) == 0) | \
                    (num_batches == (len(num_samples_patient) - 1)):
                num_batches += 1
                break
            number_samples = number_samples - num_samples_patient[num_batches]
            num_batches += 1

        batches = []
        for i_c in range(num_batches):
            batches.extend([i_c + 1 + ind] * num_samples_patient[i_c])

        array_sampleids = np.arange(1 + temp_num_oom, temp_num_oom + max_samplenumber + 1)

        for i_c in range(len(array_sampleids)):
            metadata_batch.loc[i_c] = ["_".join([project, str(array_sampleids[i_c])]), batches[i_c]]

        metadata[project] = pd.DataFrame(columns=["sample", "batch"])
        for i_c in range(len(sample_names)):
            sample_batch = metadata_batch[(metadata_batch ==
                                           list_project_samples[project][i_c]).any(1)].stack().unique()
            metadata[project].loc[i_c] = [sample_batch[0], sample_batch[1]]
        num_samples_patient = list_batches

    return metadata


def get_tissue_annot(adata):
    # Use annotations from pathologist instead of clusters
    obs_keys = list(adata.obs_keys())
    # get all manual annotations by extracting all keys with upper case characters
    annotations = [char for char in obs_keys if any(c.isupper() for c in char)]

    if "ANNOTATOR" in annotations:
        annotations.remove("ANNOTATOR")

    # split annotations by disease and tissue / cell types
    elements = ["LESONAL", "NON LESIONAL"]
    try:
        index = []
        for el in elements:
            index.append(annotations.index(el))
        target_index = np.amax(index) + 1
        disease_index = np.amin(index)
    except ValueError:
        target_index = None
        print("REMINDER: Insert 'LESONAL' and 'NON LESIONAL' in your excel file")

    tissue_cell_labels = annotations[target_index:]
    disease_labels = np.array(annotations[:disease_index])
    lesion_labels = np.array(itemgetter(*index)(annotations))

    # remove cell cycle annotations from tissue_cell_labels list
    for cc in ["G1_score", "G2M_score", "S_score", "M_score"]:
        try:
            tissue_cell_labels.remove(cc)
        except ValueError:
            continue

    return tissue_cell_labels, disease_labels, lesion_labels


def add_task_observable(adata, task_name, mask_conditions):
    """

    Parameters
    ----------
    adata
    task_name
    mask_conditions

    Returns
    -------

    """
    # assign to observable cond 1 cond 2 ..
    default_label = np.array(['default'] * adata.shape[0], dtype='<U24')
    adata.obs[task_name] = default_label

    for ind, mask_cond in enumerate(np.asarray(mask_conditions)):
        adata.obs[task_name][mask_cond] = str(ind)
    adata.obs[task_name] = adata.obs[task_name].astype('category')

    return adata


def get_dge_analysis_tasks(path_xlsx_file):
    """Get DGE analysis tasks provided in .xlsx file
    2: Or operation
    1: category included in condition and AND operation
    0: category not considered in condition
    -1: category excluded in condition and DEMUX or AND with NOT operation
    -2: XOR
    Example:
                    upper EPIDERMIS middle EPIDERMIS basal EPIDERMIS Hairfollicle
    condition 1:        1               0               0               -1
    condition 1:        0               1               0               0

    How to read:
    condition 1: Get all spots with label upper EPIDERMIS and exclude those with label upper EPIDERMIS and Hairfollicle
    condition 2: Get all spots with label middle EPIDERMIS

    Parameters
    ----------
    path_xlsx_file : str
        path to .xlsx file

    Returns
    -------
    sheet_to_df_map : dict
        dictionary with DGE analysis tasks

    """
    # get all sheet names in .xlsx file
    xls = pd.ExcelFile(path_xlsx_file)  # xls.sheet_names

    # loop through sheet_names
    sheet_to_df_map = dict()
    for sheet_name in xls.sheet_names:
        sheet_to_df_map[sheet_name] = xls.parse(sheet_name)

    return sheet_to_df_map
