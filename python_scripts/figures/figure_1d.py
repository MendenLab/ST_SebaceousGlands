import scanpy as sc
import numpy as np
import pandas as pd
import math

import matplotlib.pyplot as plt

import os
from datetime import date


fontsize_xylabel = 14
fontsize_xyticks = 12
fontsize_legend_title = 12
fontsize_text = 10
fontsize_legend = 10
fontsize_axis_title = 8


def plot_aeh(coord_sample, histology_results, patterns, num_patterns, specimen, biopsy_type,
             save_folder, invert, disease, figsize):

    fig = plt.figure(figsize=figsize)
    # Plot AEH patterns
    for i in range(1, num_patterns + 1):
        ax = fig.add_subplot(3, math.ceil(num_patterns / 3), i)
        if specimen == 'SN-V11J13-122_B_20':
            img = ax.scatter(coord_sample['y'], coord_sample['x'], c=patterns[i], s=2)
        else:
            img = ax.scatter(coord_sample['x'], coord_sample['y'], c=patterns[i], s=2)
        ax.invert_xaxis()
        if invert:
            ax.invert_yaxis()
        ax.set_title('Pattern {} - {} genes'.format(i, histology_results.query('pattern == @i').shape[0]),
                     fontsize=fontsize_axis_title)
        cb = plt.colorbar(img, ax=ax)
        cb.ax.tick_params(labelsize=6)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'AEH_Pattern_{}_{}_{}.pdf'.format(specimen, disease, biopsy_type)))
    plt.close()


def main(path_adata, save_folder, path_to_df):
    n_pattern = 9
    specimens = []
    invert = {"2-V19S23-004-V4_2": False}
    # (width, height)
    figure_sizes = {"2-V19S23-004-V4_2": (5, 5)}
    filename = "2-V19S23-004-V4_2_AD_LESIONAL"

    if os.path.isfile(path_to_df):
        df = pd.read_excel(path_to_df, sheet_name='coordinates_pattern')
        specimen = df['specimen'].unique()[0]
        biopsy_type = df['biopsy_type'].unique()[0]
        disease = df['disease'].unique()[0]

        coord_sample = df.loc[:, df.columns.isin(['x', 'y'])]
        patterns = df.loc[:, ~df.columns.isin(['x', 'y', 'disease', 'biopsy_type', 'specimen'])]

        histology_results = pd.read_excel(path_to_df, sheet_name='AEH', index_col=0)

        plot_aeh(coord_sample=coord_sample, histology_results=histology_results,
                 patterns=patterns, num_patterns=n_pattern,
                 specimen=specimen, biopsy_type=biopsy_type, save_folder=save_folder, invert=invert[specimen],
                 disease=disease, figsize=figure_sizes[specimen])

    else:
        h5files = [name for name in os.listdir(path_adata) if os.path.isdir(os.path.join(path_adata, name))]
        try:
            h5files.remove('Pathway_enrichment_analysis')
        except ValueError:
            print('No Pathway_enrichment_analysis folder in path')

        adata = sc.read(os.path.join(path_adata, 'st_QC_normed_BC_project_PsoAD.h5'))

        print('File: {}'.format(filename))
        specimen = "_".join(filename.split(sep='_')[:-2])
        adata_sample = adata[adata.obs['specimen'] == specimen].copy()
        specimens.append(specimen)
        biopsy_type = adata_sample.obs['biopsy_type'].cat.categories[0]
        disease = adata_sample.obs['DISEASE'].cat.categories[0]

        coord_sample = adata_sample.obsm['spatial']
        coord_sample = pd.DataFrame.from_dict({'x': coord_sample[:, 0], 'y': coord_sample[:, 1]})
        coord_sample.index = adata_sample.obs.index

        histology_results = pd.read_csv(os.path.join(path_adata, filename, 'AEH__{}.csv'.format(specimen)), index_col=0)

        patterns = pd.DataFrame(index=adata_sample.obs.index, columns=np.arange(1, n_pattern + 1))
        for pattern in range(1, n_pattern + 1):
            patterns.loc[:, pattern] = adata_sample.obs.loc[:, 'Pattern_intensity_{}'.format(pattern)]

        plot_aeh(coord_sample=coord_sample, histology_results=histology_results, patterns=patterns, num_patterns=n_pattern,
                 specimen=specimen, biopsy_type=biopsy_type, save_folder=save_folder, invert=invert[specimen],
                 disease=disease, figsize=figure_sizes[specimen])

        df = patterns.copy()
        df['x'] = coord_sample['x']
        df['y'] = coord_sample['y']
        num_patterns = []
        for i in range(1, n_pattern + 1):
            num_patterns.append(histology_results.query('pattern == @i').shape[0])
        # num_patterns.extend([np.nan, np.nan])
        # df.loc['number of genes'] = num_patterns
        df['disease'] = disease
        df['biopsy_type'] = biopsy_type
        df['specimen'] = specimen

        # writer = pd.ExcelWriter(os.path.join(save_folder, "Plots_AEH_Pattern.xlsx"), engine='xlsxwriter')
        writer = pd.ExcelWriter(path_to_df, engine="xlsxwriter")
        # Save figure parameters to excel file
        # df.to_excel(writer, sheet_name="Plot_{}_{}_{}".format(
        #     specimen, disease, "".join(next(zip(*biopsy_type.split(' '))))), index=False)
        df.to_excel(writer, index=False, sheet_name='coordinates_pattern')
        histology_results.to_excel(writer, sheet_name='AEH', index=False)

        writer.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "figure_1d", str(today))
    os.makedirs(savepath, exist_ok=True)

    adata_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE/2023-09-18_paper_figures'

    path_df = os.path.join("/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer",
                           "output", "figure_1d", "Figure_1d.xlsx")
    main(path_adata=adata_path, save_folder=savepath, path_to_df=path_df)
