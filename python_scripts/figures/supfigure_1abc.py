import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import os
import numpy as np
from datetime import date


def main(adata, save_folder):
    adata = adata[adata.obs['DISEASE'] != 'LP'].copy()
    adata.uns['spot_type_colors'] = np.asarray(
        ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#e377c2', '#8c564b',
         '#aa40fc', '#b5bd61', '#17becf', '#aec7e8'])

    # a)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata=adata, color='DISEASE', ax=ax, show=False, title='', palette=['darkred', 'darkblue'], s=15)
    sns.despine(fig=fig, ax=ax)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), title='',
              fancybox=True, shadow=True, ncol=3, prop={'size': 15}, frameon=False, title_fontsize=14)
    plt.tight_layout()

    plt.savefig(os.path.join(
        save_folder, 'UMAP_DISEASE_LNL_AD_Pso.pdf'), bbox_inches='tight')
    plt.close()

    # b)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata=adata, color='biopsy_type', ax=ax, show=False, title='', palette=['darkorange', 'darkgreen'], s=15)
    sns.despine(fig=fig, ax=ax)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), title='',
              fancybox=True, shadow=True, ncol=3, prop={'size': 15}, frameon=False, title_fontsize=14)
    plt.tight_layout()

    plt.savefig(os.path.join(
        save_folder, 'UMAP_biopsy_type_LNL_AD_Pso.pdf'), bbox_inches='tight')
    plt.close()

    # separate adata into adata dermis + epidermis and tissue structures
    mask = adata.obs['spot_type'].isin(
        ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'JUNCTION', 'DERMIS'])

    # c)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata=adata[mask], color='spot_type', ax=ax, show=False, title='', s=5)
    sc.pl.umap(adata=adata[~mask], color='spot_type', ax=ax, show=False, title='', s=25)
    sns.despine(fig=fig, ax=ax)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), title='',
              fancybox=True, shadow=True, ncol=3, prop={'size': 15}, frameon=False, title_fontsize=14)
    plt.tight_layout()

    plt.savefig(os.path.join(
        save_folder, 'UMAP_spot_type_EpidermisDermis_LNL_AD_Pso.pdf'), bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "Fig_S1abc", str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')

    main(adata=spatial_adata, save_folder=savepath)
