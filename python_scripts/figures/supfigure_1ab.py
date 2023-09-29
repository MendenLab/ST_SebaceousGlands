import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import os
from datetime import date
from operator import itemgetter


def main(adata, save_folder):
    adata = adata[~((adata.obs['DISEASE'] == 'LP') & (adata.obs['biopsy_type'] == 'LESIONAL'))].copy()

    # drop all samples from the newest cohort
    new_samples = list(adata.obs['sample'].cat.categories[adata.obs['sample'].cat.categories.str.contains('P21093')])
    adata = adata[~adata.obs['sample'].isin(new_samples)].copy()

    # remove JUNCTION
    mask = adata.obs["spot_type"] == "JUNCTION"
    mask_middle_epidermis = adata.obs["middle EPIDERMIS"] == 1
    mask_basal_epidermis = adata.obs["basal EPIDERMIS"] == 1
    mask_upper_epidermis = adata.obs["upper EPIDERMIS"] == 1
    mask_dermis = adata.obs["DERMIS"] == 1
    adata.obs.loc[mask & mask_middle_epidermis, "spot_type"] = "middle EPIDERMIS"
    adata.obs.loc[mask & mask_basal_epidermis, "spot_type"] = "basal EPIDERMIS"
    adata.obs.loc[mask & mask_upper_epidermis, "spot_type"] = "upper EPIDERMIS"
    adata.obs.loc[mask & mask_dermis, "spot_type"] = "DERMIS"
    adata.obs["spot_type"] = adata.obs["spot_type"].cat.remove_unused_categories()

    adata.obs["spot_type"] = adata.obs['spot_type'].cat.reorder_categories(
        ["upper EPIDERMIS", "middle EPIDERMIS", "basal EPIDERMIS", "DERMIS", 'MUSCLE', "VESSEL", "HAIR FOLLICLE",
         "SEBACEOUS GLAND", "SWEAT GLAND"])

    colors_spottype = dict(zip(adata.obs['spot_type'].cat.categories.to_list(), [
        "#1f77b4", "#ff7f0e", "#279e68", '#e377c2', '#8c564b', '#aa40fc', '#b5bd61', '#17becf', '#aec7e8']))

    adata.obs['NonLesion_disease'] = adata.obs['DISEASE'].copy()
    adata.obs['NonLesion_disease'] = adata.obs['NonLesion_disease'].astype(str)
    adata.obs.loc[adata.obs['biopsy_type'] == 'NON LESIONAL', "NonLesion_disease"] = 'NON LESIONAL'
    adata.obs['NonLesion_disease'] = adata.obs['NonLesion_disease'].astype('category')
    # set order of categories
    adata.obs["NonLesion_disease"] = adata.obs["NonLesion_disease"].cat.reorder_categories(
        ['NON LESIONAL', 'AD', 'Pso'])

    # a)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata=adata, color='NonLesion_disease', ax=ax, show=False, title='',
               palette=['grey', 'darkred', 'darkblue'], s=20)
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

    # separate adata into adata dermis + epidermis and tissue structures
    mask = adata.obs['spot_type'].isin(
        ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'DERMIS'])

    # b)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata=adata[mask], color='spot_type', ax=ax, show=False, title='', s=5,
               palette=itemgetter(*adata.obs.loc[
                   mask, 'spot_type'].cat.remove_unused_categories().cat.categories.to_list())(colors_spottype))
    sc.pl.umap(adata=adata[~mask], color='spot_type', ax=ax, show=False, title='', s=25,
               palette=itemgetter(*adata.obs.loc[
                   ~mask, 'spot_type'].cat.remove_unused_categories().cat.categories.to_list())(colors_spottype))
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
        "Fig_S1ab", str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')

    main(adata=spatial_adata, save_folder=savepath)
