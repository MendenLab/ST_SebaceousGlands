import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import os
import numpy as np
import pandas as pd
from operator import itemgetter
from datetime import date


def main(adata, save_folder, path_to_df=None):
    if os.path.isfile(path_to_df):
        df = pd.read_excel(path_to_df, index_col=0)

        palette = dict(zip(df['spot_type'], df['color']))

        fig, ax = plt.subplots(figsize=(8, 8))
        sns.scatterplot(data=df, x='UMAP1', y='UMAP2', hue='spot_type', s=18,
                        ax=ax, linewidth=0, edgecolor='k', palette=itemgetter(*list(df['spot_type'].unique()))(palette))
        sns.despine(fig=fig, ax=ax)
        ax.set_ylabel("UMAP2", fontsize=18)
        ax.set_xlabel("UMAP1", fontsize=18)
        # Hide X and Y axes tick marks
        ax.set_xticks([])
        ax.set_yticks([])

        # Put a legend below current axis
        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.05),
            title="",
            fancybox=True,
            shadow=True,
            ncol=3,
            prop={"size": 12},
            frameon=False,
            title_fontsize=14,
        )
        plt.tight_layout()
        
        plt.savefig(
            os.path.join(save_folder, "UMAP_spot_type_LNL_AD_Pso.pdf"), bbox_inches="tight"
        )
        plt.close()
    else:
        # Remove Lesion LP samples
        adata = adata[~((adata.obs['DISEASE'] == 'LP') & (adata.obs['biopsy_type'] == 'LESIONAL'))].copy()

        # drop all samples from newest cohort
        new_samples = list(adata.obs['sample'].cat.categories[adata.obs['sample'].cat.categories.str.contains('P21093')])
        adata = adata[~adata.obs['sample'].isin(new_samples)].copy()

        adata.uns['spot_type_colors'] = np.asarray(
            ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#e377c2', '#8c564b',
             '#aa40fc', '#b5bd61', '#17becf', '#aec7e8'])

        # separate adata into adata dermis + epidermis and tissue structures
        mask = adata.obs['spot_type'].isin(
            ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'JUNCTION', 'DERMIS'])

        # Remove Epidermis and Dermis
        adata = adata[~mask].copy()

        color = dict(zip(adata.obs['spot_type'].cat.categories, adata.uns['spot_type_colors']))
        df = pd.DataFrame.from_dict({'UMAP1': adata.obsm['X_umap'][:, 0], 'UMAP2': adata.obsm['X_umap'][:, 1],
                                     'spot_type': adata.obs['spot_type'],
                                     'color': itemgetter(*list(adata.obs['spot_type']))(color)})
        df.to_excel(path_to_df)

        fig, ax = plt.subplots(figsize=(8, 8))
        sc.pl.umap(adata=adata, color='spot_type', ax=ax, show=False, title='')
        sns.despine(fig=fig, ax=ax)
        ax.set_ylabel('UMAP2', fontsize=18)
        ax.set_xlabel('UMAP1', fontsize=18)

        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), title='',
                  fancybox=True, shadow=True, ncol=3, prop={'size': 12}, frameon=False, title_fontsize=14)
        plt.tight_layout()

        plt.savefig(os.path.join(
            save_folder, 'UMAP_spot_type_LNL_AD_Pso.pdf'), bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "figure_1c", str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')

    path_df = os.path.join("/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer",
                           "output", "figure_1c", "Figure_1c.xlsx")

    main(adata=spatial_adata, save_folder=savepath, path_to_df=path_df)
