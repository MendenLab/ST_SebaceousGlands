import scanpy as sc
import pandas as pd

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
    adata = sc.read(os.path.join(path_adata, 'st_QC_normed_BC_project_PsoAD.h5'))

    writer = pd.ExcelWriter(os.path.join(save_folder, "Plots_HE_images_spottypes.xlsx"), engine='xlsxwriter')

    invert_x = {'11-V19T12-012-V2_11': True}
    invert_y = {'11-V19T12-012-V2_11': False}

    filename = '11-V19T12-012-V2_11_AD_NON LESIONAL'
    print('File: {}'.format(filename))

    specimen = "_".join(filename.split(sep='_')[:-2])
    adata_sample = adata[adata.obs['specimen'] == specimen].copy()
    biopsy_type = adata_sample.obs['biopsy_type'].cat.categories[0]
    disease = adata_sample.obs['DISEASE'].cat.categories[0]
    sample = adata_sample.obs['sample'].cat.categories[0]

    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.spatial(adata=adata_sample, color='spot_type', library_id=sample, ax=ax, title='', show=False)
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
    color_mapping = dict(
        zip(adata_sample.obs['spot_type'].cat.categories, list(adata_sample.uns['spot_type_colors'])))
    colors = pd.DataFrame(list(adata_sample.obs['spot_type'].values), columns=['spot_type'])
    colors = colors.replace({"spot_type": color_mapping})

    df = pd.DataFrame.from_dict({
        'spatial_x': adata_sample.obsm['spatial'][:, 1] * adata_sample.uns['spatial'][sample]['scalefactors'][
            'tissue_hires_scalef'],
        'spatial_y': adata_sample.obsm['spatial'][:, 0] * adata_sample.uns['spatial'][sample]['scalefactors'][
            'tissue_hires_scalef'],
        'color': colors.values.T[0],
        'spot_type': list(adata_sample.obs['spot_type'].values)})

    df['spot_type'] = df['spot_type'].astype('category')
    df['spot_type'] = df['spot_type'].cat.reorder_categories(list(adata_sample.obs['spot_type'].cat.categories))

    # Save figure parameters to excel file
    df.to_excel(writer, sheet_name="Plot_{}_{}_{}".format(
        specimen, disease, "".join(next(zip(*biopsy_type.split(' '))))), index=False)

    writer.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "figure_2b__spatialDE_HE_images", str(today))
    os.makedirs(savepath, exist_ok=True)

    adata_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE/2023-04-12_paper_figures'

    main(path_adata=adata_path, save_folder=savepath)
