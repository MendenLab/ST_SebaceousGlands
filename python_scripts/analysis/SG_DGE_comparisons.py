import python_scripts.utils.operation_obs as opobs
from datetime import date
import os
from fnmatch import fnmatch
import glob
import sys
from operator import itemgetter
import numpy as np
import pandas as pd
import scanpy as sc

"""
['<command_to_run>', '<path_to_script>', 'arg1' , 'arg2', 'arg3', 'arg4']
import subprocess
import rpy2.robjects as robjects
Define command and arguments
command = 'Rscript'
Call rscript:      
r_source = robjects.r['source']         
r_source(os.path.join(rscript_path, "DESeq2_DGE_analysis.R"))       
OR
args = (csv_filename, "test", False, True, True)
cmd = [command, os.path.join(rscript_path, "DESeq2_DGE_analysis.R")] + args  
x = subprocess.check_output(cmd, universal_newlines=True)

Call function in rscript:   r.r.source(os.path.join(rscript_path, "DESeq2_DGE_analysis.R"))
                            output = r.r["main"](csv_filename, "test", False, True, True)
"""

minor_conditions = ["PSO", "AE", "LICHEN", "PRP", "LESONAL", "NON LESIONAL"]


def get_spots_per_condition(adata, tasks, cond_filename, save_folder, file_shortcuts, save_condition_masks=True):
    """Condition handling:
    - disease * lesion label * tissue labels are connected via AND operator
    - diseases are connected via OR operator
    - lesion labels are connected via OR operator
    - tissue labels are connected via OR, AND, DEMUX or XOR operator

    Parameters
    ----------
    adata : annData
    tasks : dict
    cond_filename : str
    save_folder : str
    file_shortcuts : str
        path to shortcuts for current manual annotations (not automized yet)
    save_condition_masks : bool
        if mask for conditions for pairwise DGE analysis shall be saved or not (default is True)

    Returns
    -------

    """
    mask_conditions = []

    # read in short_cuts
    sc_categories = pd.read_csv(os.path.join(os.getcwd(), file_shortcuts), sep=";")

    # task keys as list
    tasks_keys = list(tasks.keys())

    # get patient -- Attention: currently manually specified for towards our dataset
    df_patients = opobs.map_sample_batch_list(
        adata=adata, num_samples_patient=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 4])
    # get batches assigned to each patient (currently we have 2-4 samples per patient)
    df_batches = opobs.map_sample_batch_list(
        adata=adata, num_samples_patient=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4])

    for ind_tasks, name_task in enumerate(tasks_keys):
        if "default" not in name_task:
            # task file
            df_tmp = tasks[name_task].set_index("Unnamed: 0")
            df_condition = pd.DataFrame({"geneNames": adata.var.index})
            # meta Data
            df_metadata = pd.DataFrame(columns=["sample_id", "barcode", "condition", "label", "patient", "disease",
                                                "biopsy_type", "num_spots", "sizefactor", 'batch'])
            for cond in df_tmp.index:
                # Save .h5 object and reset condition mask
                if ("1" in cond) & (len(mask_conditions) > 0):
                    # task from before
                    task_previous = "{}_{}".format(cond_filename, tasks_keys[ind_tasks - 1])
                    adata = opobs.add_task_observable(
                        adata=adata, task_name=task_previous, mask_conditions=mask_conditions)

                    # reset mask for new task
                    mask_conditions = []

                # Load categories (columns in task sheet)
                categories = df_tmp.loc[cond]
                if sum(categories.values) > 0:
                    # get spots of observables with 1 and -1 in condition
                    working_conditions = categories[categories != 0]
                    condition_name = ""
                    # Incorporate abbreviations
                    for ind, c in enumerate(working_conditions):
                        # get abbreviation of condition_name
                        sc_cond_name = sc_categories["Abbreviation"][sc_categories['Explanation'] ==
                                                                     working_conditions.index[ind]].values[0]
                        if c == -1:
                            condition_name = "_".join([str(condition_name), sc_cond_name, "no_"]).replace(" ", "_")
                        else:
                            condition_name = "_".join([condition_name, sc_cond_name]).replace(" ", "_")
                    # Group conditions into disease, non-lesonial/lesional, tissue parts

                    # 1. Get all manual annotations by extracting all keys with upper case characters
                    annotations = [char for char in categories.index if any(c.isupper() for c in char)]

                    elements = ["LESIONAL", "LESONAL", "NON LESIONAL"]
                    index = []
                    for el in elements:
                        if el in annotations:
                            index.append(annotations.index(el))
                    try:
                        target_index = np.amax(index) + 1
                        disease_index = np.amin(index)
                    except ValueError:
                        print("REMINDER: Insert 'LESONAL' and 'NON LESIONAL' in your excel file")
                        sys.exit()

                    bool_list = list(map(bool, categories.values))
                    # 2. Save and group disease, lesion and tissue labels and create bool arrays for labels
                    tissue_cell_labels = annotations[target_index:]
                    sign_assignment_tissue = categories.values[target_index:]
                    bool_tissue_labels = bool_list[target_index:]
                    disease_labels = np.array(annotations[:disease_index])
                    bool_disease_labels = np.array(bool_list[:disease_index])
                    lesion_labels = np.array(itemgetter(*index)(annotations))
                    bool_lesion_labels = np.array(itemgetter(*index)(bool_list))

                    # 3.1 read out spots belonging to disease
                    disease_spots = disease_labels[bool_disease_labels]
                    # 3.2 read out spots being labeled as lesion, non lesion or both
                    lesion_spots = lesion_labels[bool_lesion_labels]
                    # 3.3 read out tissue part
                    tissue_spots = np.array(tissue_cell_labels)[bool_tissue_labels]
                    sign_assignment_tissue = np.array(sign_assignment_tissue)[bool_tissue_labels]
                    # reverse sorting to get the order -2, -1, 0, 1, 2 - todo check
                    index_signs = np.argsort(sign_assignment_tissue)[::-1]
                    tissue_spots = tissue_spots[index_signs]
                    # reverse sorting to get the order 2, 1, 0, -1, -2
                    sign_assignment_tissue[::-1].sort()

                    # 4. Get bool arrays for all groups
                    # If comparison are done with specific tissue spots or disease lesion vs non-lesional
                    if len(tissue_spots) > 0:
                        bool_tissue_spots = adata.obs[tissue_spots].values.astype(bool)
                        # get first column from bool array
                        bool_array_tissue = bool_tissue_spots[:, 0]
                        try:
                            for i_c in range((bool_tissue_spots.shape[1] - 1)):
                                if sign_assignment_tissue[i_c + 1] == 2:
                                    # OR
                                    bool_array_tissue = np.logical_or(bool_array_tissue, bool_tissue_spots[:, i_c + 1])
                                elif sign_assignment_tissue[i_c + 1] == 1:
                                    # AND with OR
                                    bool_array_tissue = np.logical_and(bool_tissue_spots[:, i_c + 1], bool_array_tissue)
                                elif sign_assignment_tissue[i_c + 1] == -1:
                                    # AND with NOT
                                    bool_array_tissue = np.logical_and(
                                        bool_array_tissue, np.logical_not(bool_tissue_spots[:, i_c + 1]))
                                elif sign_assignment_tissue[i_c + 1] == -2:
                                    # XOR
                                    bool_array_tissue = np.logical_xor(bool_array_tissue, bool_tissue_spots[:, i_c + 1])
                        except ValueError:
                            print("Please use one of the following denotations in your excel file")
                            print("OR = 2\n AND = 1\n DEMUX[Y0] (AND with NOT) = -1\n XOR = -2")

                        # concat biopsy type with OR
                        bool_lesion_spots = adata.obs[lesion_spots].values.astype(bool)
                        bool_array_lesion = np.any(bool_lesion_spots, axis=1)

                        # concat diseases with OR
                        bool_disease_spots = adata.obs[disease_spots].values.astype(bool)
                        bool_array_disease = np.any(bool_disease_spots, axis=1)

                        # multiply (= AND operator) groups of diseases, lesion/non-lesion and tissue parts
                        bool_array_spots = bool_array_disease * bool_array_lesion * bool_array_tissue
                    else:
                        # If only disease and lesion/non lesional spots shall be considered
                        # concat biopsy type with OR
                        bool_lesion_spots = adata.obs[lesion_spots].values.astype(bool)
                        bool_array_lesion = np.any(bool_lesion_spots, axis=1)

                        # concat diseases with OR
                        bool_disease_spots = adata.obs[disease_spots].values.astype(bool)
                        bool_array_disease = np.any(bool_disease_spots, axis=1)

                        # multiply groups of diseases and lesion/non-lesion
                        bool_array_spots = bool_array_disease * bool_array_lesion

                    # spots in samples having all required conditions fulfilled
                    sample_names = np.unique(adata.obs['sample'][bool_array_spots])

                    # store mask in object -> read out sub adata object later
                    if save_condition_masks:
                        mask_conditions.append(bool_array_spots)

                    # per sample: get count from spots and add counts to dataframe
                    for ind, sample in enumerate(sample_names):
                        index_sample = np.where(adata.obs['sample'] == sample)  # index array
                        adata_sample = adata[adata.obs['sample'] == sample].copy()  # copy sample adata

                        # get mask of spots involved in condition
                        sample_bool_array_spots = bool_array_spots[index_sample]

                        if 'counts' in adata_sample.layers.keys():
                            sample_count_matrix = adata_sample.layers['counts']
                        else:
                            # Read in normalised counts
                            sample_count_matrix = adata_sample.X

                        # Important information for DGE Analysis
                        no_spots = sample_count_matrix[sample_bool_array_spots].shape[0]
                        # save to metaData
                        # "sample_id", "condition", "label", "batch", "disease", "num_spots"
                        project = sample.split("_")[0]

                        """ ATTENTION: This works only if we have unique tissue biopsies per capture area !! """
                        disease_columns = adata_sample.obs[disease_labels] == 1
                        disease_columns = np.asarray(disease_columns.idxmax(1).values)
                        type_lesional = adata_sample.obs[lesion_labels] == 1
                        type_lesional = np.asarray(type_lesional.idxmax(1).values)
                        index_batch_sample = df_batches[project].index[df_batches[project]['sample'] == sample][0]
                        index_patient_sample = df_patients[project].index[df_patients[project]['sample'] == sample][0]

                        # read out spots and label them in metaData with condition name
                        sample_barcodes = adata_sample.obs.index
                        matrix = sample_count_matrix[sample_bool_array_spots]
                        for index, spot in enumerate(matrix):
                            df_condition[".".join(
                                [sample, sample_barcodes[index], condition_name, cond]).replace(" ", "_")] = spot

                            # save to metaData
                            # Size Factors were calculated on whole data set
                            df_temp = pd.DataFrame(
                                {"sample_id": [sample],
                                 "barcode": [sample_barcodes[index]],
                                 "condition": [cond.split("_")[-1]],
                                 "label": [adata_sample.obs['tissue_type'].values[sample_bool_array_spots][index]],
                                 "patient": [int(df_patients[project]['batch'][index_patient_sample])],
                                 "disease": np.unique(disease_columns),
                                 "biopsy_type": np.unique(type_lesional),
                                 "num_spots": [no_spots],
                                 "sizefactor": [
                                     adata_sample.obs['size_factors'].values[sample_bool_array_spots][index]],
                                 "batch": [int(df_batches[project]['batch'][index_batch_sample])]})  # Slide

                            df_metadata = df_metadata.append(df_temp, ignore_index=True)
                elif sum(categories.values) < 0:
                    print("Seems like you only want to exclude spots.")
                    print("=> Please check your Excel files <=")
                    sys.exit(0)
                else:
                    # Do this if condition contains only 0's:
                    # spots having all required conditions fulfilled
                    sample_names = np.unique(adata.obs['sample'])

                    for ind, sample in enumerate(sample_names):
                        if 'counts' in adata.layers.keys():
                            sample_count_matrix = adata.layers['counts'][np.where(adata.obs['sample'] == sample)]
                        else:
                            # Read in normalised counts
                            sample_count_matrix = adata.X[np.where(adata.obs['sample'] == sample)]
                        df_condition["".join([sample, cond])] = np.zeros(sample_count_matrix.sum(axis=0).shape)

            # df_metadata.to_csv(os.path.join(save_folder, "{}_{}.csv".format("metaData", name_task)))
            # df_condition.to_csv(os.path.join(save_folder, "{}.csv".format(name_task)))

    # Save task to .h5 object for the case only one task is in the .xlsx sheet
    if len(tasks_keys) == 1:
        task_previous = "{}_{}".format(cond_filename, tasks_keys[0])
        adata = opobs.add_task_observable(
            adata=adata, task_name=task_previous, mask_conditions=mask_conditions)

    return adata


def main(adata, save_folder,  name_analysis_folder):
    """DGE, GO-term, Pathway Enrichment and Gene Co-expression Analysis with manual Annotations from Alex and Kilian
    We assume that we have biological (technical) replicates.
    --> lesioned or healthy tissue biopsies from a patient are considered as biological replicates
    This means that we pretend to have multiple sequencing runs of the same library.

    Parameters
    ----------
    adata : annData
    save_folder :str
    name_analysis_folder : str

    Returns
    -------

    """
    # get abbreviations of categories
    file_shortcuts = os.path.join("..", "..", "input", "Abbreviations_categories.csv")

    # path to DGE analysis excel file
    wd_tasks = os.path.join("..", "..", "input", 'DGE_tasks', name_analysis_folder)

    # Remove unwanted observables from adata
    adata = opobs.remove_obs(adata, colname=['ANNOTATOR', "G1_score", "G2M_score", "S_score", "M_score"])

    # Remove all spots (~304) without a tissue layer annotation
    adata = adata[adata.obs['tissue_type'] != 'Unknown'].copy()

    # 1. Analysis tasks
    # 1.1 load tasks file
    all_xlsx_files = [file for path, subdir, files in os.walk(wd_tasks)
                      for file in glob.glob(os.path.join(path, '*.xlsx'))]

    pattern = ['Abkuuerzungen.xlsx', 'Abkürzungen.xlsx']
    # 1.2 get tasks and metaData
    for xlsx_file in all_xlsx_files:
        if (not fnmatch(xlsx_file, pattern[0])) & (not fnmatch(xlsx_file, pattern[1])):
            # create folder for tasks inlcuded in an .xlsx file
            cond_filename = xlsx_file.split(os.sep)[-1].split('.')[0]
            save_folder_task = os.path.join(save_folder, cond_filename)
            os.makedirs(save_folder_task, exist_ok=True)

            # get .csv files for Analysis DGE analysis comparisons
            tasks = opobs.get_dge_analysis_tasks(path_xlsx_file=xlsx_file)
            adata = get_spots_per_condition(adata=adata, tasks=tasks, cond_filename=cond_filename,
                                            save_folder=save_folder_task, file_shortcuts=file_shortcuts)

    sc.write(os.path.join(save_folder, "".join(["annData", name_analysis_folder, ".h5"])), adata)

    # path to R scripts
    # rscript_path = os.path.join(Path(__file__).parent.parent, "R_scripts")
    # # .csv file names generated in R scripts
    # csv_filename = 'joined_adata_R_dge.csv'
    # dge_csv_filename = "w_replicates_healthy_DGE_all_genes.csv"

    # 2. Analysis tasks in R
    print("-------- Start Analysis in R --------")
    # for conditions in tasks.keys():
    #     # 2.1 DGE Analysis
    #     # args = (sample, save_folder, replicate_type, logfc_factor, fdr_value, p_value, test_method)
    #     r.r.source(os.path.join(rscript_path, "glmGamPoi_DGE_analysis.R"))
    #     r.r["main"](save_folder, 'biological', 1.2, 0.05, 0.05, 'BH')
    #
    #     # 2.2 Gene Co-expression Analysis
    #     # args = (csv_filename, sample, name_tissue, save_folder, replicates, healthy, p_value)
    #     r.r.source(os.path.join(rscript_path, "CEMiTool_Co-expression_analysis.R"))
    #     # TODO provide also all genes found in data set
    #     r.r["main"](csv_filename, conditions, save_folder, True, True, 0.05)
    #
    #     # 2.3 GO-term enrichment Analysis with TopGo
    #     # args = (csv_filename, sample, name_comparison, DGE_analysis, save_folder, go_term, replicates, healthy,
    #     # test_methods, algorithms)
    #     r.r.source(os.path.join(rscript_path, "TopGO_GO_term_analysis.R"))
    #     r.r["main"](dge_csv_filename, conditions, "edgeR", "test", True, True, True, ["t", "ks"],
    #                 ["weight01", "elim"])

    print("-------- Finished Analysis in R --------")
    # --------------------- End DGE, Gene Co-expression, Go-term and Pathway enrichment Analysis --------------------- #


if __name__ == '__main__':
    today = date.today()

    # parameter
    spatial_cluster_label = 'tissue_type'
    # 20210728_Tasksforfigures/20210728_Tasksforfigures.xlsx
    adata_filename = "2021-07-29_Visium_Data_QC_BC_clustered.h5"
    analysis_folder = '20210616_Tasksforfigures'  # '20210507_PilosebaceousUnitApproaches_PS'

    # create saving folder in current project path
    savepath = os.path.join('..', '..', "output", "dge_tasks__pilosebaceous_unit", str(today))
    os.makedirs(savepath, exist_ok=True)

    # path to anndata object
    spatial_adata = sc.read(os.path.join('/Users/christina.hillig/Documents/Projects/adata_storage', adata_filename))

    main(adata=spatial_adata, save_folder=savepath, name_analysis_folder=analysis_folder)
