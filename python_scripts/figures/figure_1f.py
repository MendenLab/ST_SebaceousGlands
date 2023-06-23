import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import os
from datetime import date


def main(adata, save_folder):
    # Remove Epidermis and Dermis
    adata = adata[~adata.obs['spot_type'].isin(
        ['DERMIS', 'upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'JUNCTION'])].copy()

    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata=adata, color='spot_type', ax=ax, show=False, title='')
    sns.despine(fig=fig, ax=ax)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)

    # Legend right, next to plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='',
              fancybox=True, shadow=True, ncol=1, prop={'size': 15}, frameon=False, title_fontsize=14)
    plt.tight_layout()

    plt.savefig(os.path.join(
        save_folder, 'UMAP_spot_type_LP_AD_Pso.pdf'), bbox_inches='tight')
    plt.close()

    adata = adata[adata.obs['DISEASE'] != 'LP'].copy()

    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata=adata, color='spot_type', ax=ax, show=False, title='')
    sns.despine(fig=fig, ax=ax)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)

    # Legend right, next to plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='',
              fancybox=True, shadow=True, ncol=1, prop={'size': 15}, frameon=False, title_fontsize=14)
    plt.tight_layout()

    plt.savefig(os.path.join(
        save_folder, 'UMAP_spot_type_AD_Pso.pdf'), bbox_inches='tight')
    plt.close()

    adata = adata[adata.obs['biopsy_type'] != 'NON LESIONAL'].copy()

    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata=adata, color='spot_type', ax=ax, show=False, title='')
    sns.despine(fig=fig, ax=ax)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)

    # Legend right, next to plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='',
              fancybox=True, shadow=True, ncol=1, prop={'size': 15}, frameon=False, title_fontsize=14)
    plt.tight_layout()

    plt.savefig(os.path.join(
        save_folder, 'UMAP_spot_type_NL_AD_Pso.pdf'), bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "figure_1f", str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Remove LP
    # spatial_adata = spatial_adata[spatial_adata.obs['DISEASE'] != 'LP'].copy()

    main(adata=spatial_adata, save_folder=savepath)
