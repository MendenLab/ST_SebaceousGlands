import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import os
import pandas as pd
from operator import itemgetter
from datetime import date


def main(adata, save_folder, path_to_df):
    if os.path.isfile(path_to_df):
        df = pd.read_excel(path_to_df, index_col=0)
        palette = dict(zip(df["biopsy_type"], df["color"]))
        
        fig, ax = plt.subplots(figsize=(8, 8))
        sns.scatterplot(
            data=df,
            x="UMAP1",
            y="UMAP2",
            hue="biopsy_type",
            s=18,
            ax=ax,
            linewidth=0,
            edgecolor="k",
            palette=itemgetter(*list(df["biopsy_type"].unique()))(palette),
        )
        sns.despine(fig=fig, ax=ax)
        ax.set_ylabel("UMAP2", fontsize=18)
        ax.set_xlabel("UMAP1", fontsize=18)
        # Hide X and Y axes tick marks
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Legend below plot
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), title='',
                  fancybox=True, shadow=True, ncol=2, prop={'size': 15}, frameon=False, title_fontsize=14)
        plt.tight_layout()
        
        plt.savefig(
            os.path.join(save_folder, "UMAP_biopsy_type_AD_Pso_legend_bottom.pdf"), bbox_inches="tight"
        )
        plt.close()
    else:
        # drop all samples from newest cohort
        new_samples = list(adata.obs['sample'].cat.categories[adata.obs['sample'].cat.categories.str.contains('P21093')])
        adata = adata[~adata.obs['sample'].isin(new_samples)].copy()

        # Remove Epidermis and Dermis
        adata = adata[~adata.obs['spot_type'].isin(
            ['DERMIS', 'upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'JUNCTION'])].copy()

        color = dict(zip(adata.obs["biopsy_type"].cat.categories, ['lightcoral', 'teal']))
        df = pd.DataFrame.from_dict(
            {
                "UMAP1": adata.obsm["X_umap"][:, 0],
                "UMAP2": adata.obsm["X_umap"][:, 1],
                "biopsy_type": adata.obs["biopsy_type"],
                "color": itemgetter(*list(adata.obs["biopsy_type"]))(color),
            }
        )
        df.to_excel(path_to_df)

        fig, ax = plt.subplots(figsize=(6, 6))
        sc.pl.umap(adata=adata, color='biopsy_type', ax=ax, show=False, title='', palette=['lightcoral', 'teal'])
        sns.despine(fig=fig, ax=ax)
        ax.set_ylabel('UMAP2', fontsize=18)
        ax.set_xlabel('UMAP1', fontsize=18)

        # Legend below plot
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), title='',
                  fancybox=True, shadow=True, ncol=2, prop={'size': 15}, frameon=False, title_fontsize=14)
        plt.tight_layout()

        plt.savefig(os.path.join(
            save_folder, 'UMAP_biopsy_type_AD_Pso_legend_bottom.pdf'), bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata=adata, color='biopsy_type', ax=ax, show=False, title='', palette=['lightcoral', 'teal'])
        sns.despine(fig=fig, ax=ax)
        ax.set_ylabel('UMAP2', fontsize=18)
        ax.set_xlabel('UMAP1', fontsize=18)

        # Legend right, next to plot
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='',
                  fancybox=True, shadow=True, ncol=1, prop={'size': 15}, frameon=False, title_fontsize=14)
        plt.tight_layout()

        plt.savefig(os.path.join(
            save_folder, 'UMAP_biopsy_type_AD_Pso_legend_right.pdf'), bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "figure_3a", str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Remove LP
    mask = (spatial_adata.obs['DISEASE'] == 'LP') & (spatial_adata.obs['biopsy_type'] == 'LESIONAL')
    spatial_adata = spatial_adata[~mask].copy()

    path_df = os.path.join("/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer",
                           "output", "figure_3a", "Figure_3a.xlsx")

    main(adata=spatial_adata, save_folder=savepath, path_to_df=path_df)
