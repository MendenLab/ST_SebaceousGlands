import scanpy as sc
import pandas as pd
import numpy as np

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


def export_legend(legend, save_folder, ncol, key, filename="legend", expand=[-5, -5, 5, 5]):
    # https://stackoverflow.com/questions/4534480/get-legend-as-a-separate-picture-in-matplotlib
    fig = legend.figure
    sns.despine(left=True, bottom=True, fig=fig)
    fig.canvas.draw()
    # bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    bbox = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    # fig.tight_layout()
    fig.savefig(os.path.join(save_folder, "{}_{}_{}{}".format(filename, ncol, key, '.pdf')), dpi='figure',
                bbox_inches=bbox)
    plt.close(fig=fig)


def main(path_adata, save_folder):
    h5files = ['10-V19T12-025-V3_10_Pso_LESIONAL', '11-V19T12-012-V1_11_AD_NON LESIONAL', '2-V19S23-004-V3_2_AD_LESIONAL']

    adata = sc.read(os.path.join(path_adata, 'st_QC_normed_BC_project_PsoAD.h5'))
    adata.uns['spot_type_colors'] = np.asarray(
        ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#e377c2', '#8c564b',
         '#aa40fc', '#b5bd61', '#17becf', '#aec7e8'])

    specimens = []
    writer = pd.ExcelWriter(os.path.join(save_folder, "Plots_HE_images_spottypes.xlsx"), engine='xlsxwriter')

    invert_x = {'10-V19T12-025-V3_10': True,  '11-V19T12-012-V1_11': True, '2-V19S23-004-V3_2': True}
    invert_y = {'10-V19T12-025-V3_10': True,  '11-V19T12-012-V1_11': False, '2-V19S23-004-V3_2': True}

    # Read out specimen
    for filename in h5files:
        print('File: {}'.format(filename))

        specimen = "_".join(filename.split(sep='_')[:-2])
        adata_sample = adata[adata.obs['specimen'] == specimen].copy()
        specimens.append(specimen)
        biopsy_type = adata_sample.obs['biopsy_type'].cat.categories[0]
        disease = adata_sample.obs['DISEASE'].cat.categories[0]
        sample = adata_sample.obs['sample'].cat.categories[0]

        fig, ax = plt.subplots(figsize=(10, 10))
        sc.pl.spatial(adata=adata_sample, color='spot_type', library_id=sample, ax=ax, title='', show=False,
                      legend_loc='best')
        if invert_x[specimen]:
            ax.invert_xaxis()
        if invert_y[specimen]:
            ax.invert_yaxis()
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        # # Put a legend below current axis
        # leg = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), title='',
        #                 fancybox=True, shadow=True, ncol=10, prop={'size': 16}, frameon=False, title_fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, 'HE_image_{}_{}_{}.pdf'.format(specimen, disease, biopsy_type)))
        plt.close()

        # Plot legend separately - Needs to be run in debugging mode
        labels = list(adata.obs['spot_type'].cat.categories)
        f = lambda m, c: plt.plot([], [], marker='o', color=c, ls="none")[0]
        handles = [f("s", adata.uns['spot_type_colors'][i]) for i in range(len(adata.uns['spot_type_colors']))]

        legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=False, ncol=10, prop={'size': 22})
        export_legend(legend, save_folder=save_folder, ncol=10, key='spot_type')

        # convert spot types to colors
        color_mapping = dict(zip(adata_sample.obs['spot_type'].cat.categories,
                                 list(adata_sample.uns['spot_type_colors'])))
        colors = pd.DataFrame(list(adata_sample.obs['spot_type'].values), columns=['spot_type'])
        colors = colors.replace({"spot_type": color_mapping})

        df = pd.DataFrame.from_dict({'spatial_x': adata_sample.obsm[
                'spatial'][:, 1] * adata_sample.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'],
                                     'spatial_y': adata_sample.obsm[
                'spatial'][:, 0] * adata_sample.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'],
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
        "Fig_S2abc", str(today))
    os.makedirs(savepath, exist_ok=True)

    adata_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE/2023-04-12_paper_figures'

    main(path_adata=adata_path, save_folder=savepath)
