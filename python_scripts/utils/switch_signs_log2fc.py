# 1. Read in DGE results file
# 2. Read in MetaData
# 3. Set condition 1 to reference and condition 2 to control
# 4. Save as .csv file

from datetime import date
import pandas as pd
import numpy as np
import os
import glob
import re


def check_condition(df):
    unique_labels = list(np.unique(df['label']))
    sorted_alphabetically_labels = sorted(unique_labels)

    # reference == positive log2fc and control == negative log2fc
    df['mask'] = np.where(
        ((df['label'] == sorted_alphabetically_labels[0]) & (df['condition'] == 'reference')) |
        ((df['label'] == sorted_alphabetically_labels[1]) & (df['condition'] == 'control')),
        np.ones(df['label'].values.shape, dtype=bool),
        np.zeros(df['label'].values.shape, dtype=bool))

    return df


def main():
    input_dir = os.path.join(
        basepath, 'input', 'switch_signs', 'ST_cdr_annotation_condition',
        '5_DGE_Approaches_Glands_modifiedPS_20210507')

    # get DGE results files in directory
    all_csv_files = glob.glob(input_dir + "/**/*DGE_all_genes.csv", recursive=True)
    # pattern = re.compile("(DGE|metaData)")

    # get metaData  files
    all_metadata_files = glob.glob(input_dir + "/**/*metaData.csv", recursive=True)

    # read out dge results file together with metaData file
    for ind_files in range(0, len(all_csv_files)):
        # Get comparison name
        name_comparison = "/".join(all_csv_files[ind_files].split(os.sep)[-3:-1])

        print(all_csv_files[ind_files])

        # Read out only those genes which are specific for a comparison
        df_dgeres = pd.read_csv(all_csv_files[ind_files], error_bad_lines=False)
        # Remove column Unnamed: 0
        if 'Unnamed: 0' in df_dgeres.keys():
            df_dgeres = df_dgeres.drop(['Unnamed: 0'], axis=1)

        # find file which belongs to dge
        matching = [s for s in all_metadata_files if name_comparison in s]
        df_metadata = pd.read_csv(matching[0], error_bad_lines=False)
        if 'Unnamed: 0' in df_metadata.keys():
            df_metadata = df_metadata.drop(['Unnamed: 0'], axis=1)

        # check whether signs need to be switched
        sub_df_metadata = df_metadata[['label', 'condition']]
        sub_df_metadata = check_condition(df=sub_df_metadata)

        # if mask False -> switch signs of log2fc and change reference and control labels in condition
        if ~np.all(sub_df_metadata['mask'].values):
            df_dgeres['log2fc'] = -df_dgeres['log2fc']
            # get index of control and replace with reference
            mask_ctrl = sub_df_metadata['condition'] == 'control'
            df_metadata['condition'][mask_ctrl] = 'reference'
            # get index of reference and replace with control
            mask_reference = sub_df_metadata['condition'] == 'reference'
            df_metadata['condition'][mask_reference] = 'control'

        # check again :)
        sub_df_metadata = df_metadata[['label', 'condition']]
        sub_df_metadata = check_condition(df=sub_df_metadata)
        if ~np.all(sub_df_metadata['mask'].values):
            print('Something went wrong ...')
        else:
            print("Successfully changed :)")

        # set paths
        new_dge_path = re.sub('{}plots_DGE_analysis'.format(os.sep), "",
                              all_csv_files[ind_files].replace(os.path.join("..", "..", 'input'), savepath))
        os.makedirs(os.path.dirname(os.path.abspath(new_dge_path)), exist_ok=True)
        new_metadata_path = re.sub('{}plots_DGE_analysis'.format(os.sep), "",
                                   matching[0].replace(os.path.join("..", "..", 'input'), savepath))
        os.makedirs(os.path.dirname(os.path.abspath(new_metadata_path)), exist_ok=True)
        # save as .csv and .xlsx
        df_dgeres.to_csv(new_dge_path)
        df_dgeres.to_excel(new_dge_path.replace(".csv", ".xlsx"))

        df_metadata.to_csv(new_metadata_path)
        df_metadata.to_excel(new_metadata_path.replace(".csv", ".xlsx"))


if __name__ == '__main__':
    today = date.today()
    basepath = os.path.join("..", "..")

    # create saving folder in current project path
    savepath = os.path.join(basepath, "output", "switch_signs_log2fc", str(today))
    os.makedirs(savepath, exist_ok=True)

    main()
