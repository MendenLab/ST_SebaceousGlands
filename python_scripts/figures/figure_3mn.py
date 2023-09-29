import scanpy as sc
import pandas as pd
from operator import itemgetter

import seaborn as sns
import matplotlib.pyplot as plt


import os
from datetime import date


def plot_image(img_rotated, df, disease, specimen, biopsy_type, save_folder, invert, colors):
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(data=df, x='spatial_x', y='spatial_y', hue='spot_type', palette=colors, ax=ax, s=8, zorder=2)
    ax.imshow(img_rotated, zorder=1)

    ax.invert_xaxis()
    if invert:
        ax.invert_yaxis()

    sns.despine(fig=fig, ax=ax, top=False, right=False, left=False, bottom=False)
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    plt.tight_layout()

    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'HE_image_{}_{}_{}.pdf'.format(specimen, disease, biopsy_type)))
    plt.close()


def main(path_adata, save_folder):
    h5files = ['2-V19S23-004-V4_2_AD_LESIONAL', '10-V19T12-025-V4_10_Pso_LESIONAL']

    adata = sc.read(os.path.join(path_adata, 'st_QC_normed_BC_project_PsoAD.h5'))

    # remove JUNCTION
    mask = adata.obs["spot_type"] == "JUNCTION"
    mask_upper_epidermis = adata.obs["upper EPIDERMIS"] == 1
    mask_middle_epidermis = adata.obs["middle EPIDERMIS"] == 1
    mask_basal_epidermis = adata.obs["basal EPIDERMIS"] == 1
    mask_dermis = adata.obs["DERMIS"] == 1
    adata.obs.loc[mask & mask_upper_epidermis, "spot_type"] = "upper EPIDERMIS"
    adata.obs.loc[mask & mask_middle_epidermis, "spot_type"] = "middle EPIDERMIS"
    adata.obs.loc[mask & mask_basal_epidermis, "spot_type"] = "basal EPIDERMIS"
    adata.obs.loc[mask & mask_dermis, "spot_type"] = "DERMIS"
    adata.obs["spot_type"] = adata.obs["spot_type"].cat.remove_unused_categories()

    colors_spottype = dict(zip(adata.obs['spot_type'].cat.categories.to_list(), [
        "#1f77b4", "#ff7f0e", "#279e68", '#e377c2', '#8c564b', '#aa40fc', '#b5bd61', '#17becf', '#aec7e8']))

    specimens = []
    writer = pd.ExcelWriter(os.path.join(save_folder, "Plots_HE_images_spottypes.xlsx"), engine='xlsxwriter')

    invert_x = {'2-V19S23-004-V4_2': True, '10-V19T12-025-V4_10': True}
    invert_y = {'2-V19S23-004-V4_2': True, '10-V19T12-025-V4_10': True}

    # Read out specimen
    for filename in h5files:
        print('File: {}'.format(filename))

        specimen = "_".join(filename.split(sep='_')[:-2])
        adata_sample = adata[adata.obs['specimen'] == specimen].copy()
        specimens.append(specimen)
        biopsy_type = adata_sample.obs['biopsy_type'].cat.categories[0]
        disease = adata_sample.obs['DISEASE'].cat.categories[0]
        sample = adata_sample.obs['sample'].cat.categories[0]

        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.spatial(adata=adata_sample, color='spot_type', library_id=sample, ax=ax, title='', show=False,
                      palette=itemgetter(*adata_sample.obs['spot_type'].cat.categories.to_list())(colors_spottype))
        if invert_x[specimen]:
            ax.invert_xaxis()
        if invert_y[specimen]:
            ax.invert_yaxis()
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, 'HE_image_{}_{}_{}.pdf'.format(specimen, disease, biopsy_type)))
        plt.close()

        # convert spot types to colors
        colors = pd.DataFrame(adata_sample.obs['spot_type'].to_list(), columns=['spot_type'])
        colors = colors.replace({"spot_type": colors_spottype})

        df = pd.DataFrame.from_dict({'spatial_x': adata_sample.obsm[
                'spatial'][:, 1] * adata_sample.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'],
                                     'spatial_y': adata_sample.obsm[
                'spatial'][:, 0] * adata_sample.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'],
                                     'color': colors.values.T[0],
                                     'spot_type': adata_sample.obs['spot_type'].to_list()})

        df['spot_type'] = df['spot_type'].astype('category')
        df['spot_type'] = df['spot_type'].cat.reorder_categories(list(adata_sample.obs['spot_type'].cat.categories))

        # Save figure parameters to Excel file
        df.to_excel(writer, sheet_name="Plot_{}_{}_{}".format(
            specimen, disease, "".join(next(zip(*biopsy_type.split(' '))))), index=False)

    writer.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "figure_3mn__spatialDE_HE_images", str(today))
    os.makedirs(savepath, exist_ok=True)

    adata_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE/2023-04-12_paper_figures'

    main(path_adata=adata_path, save_folder=savepath)
