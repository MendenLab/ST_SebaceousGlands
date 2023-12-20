import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import os
import pandas as pd
from datetime import date
from operator import itemgetter


def plot_scatter(df, hue, palette, save_folder):
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="UMAP1",
        y="UMAP2",
        hue=hue,
        s=18,
        ax=ax,
        linewidth=0,
        edgecolor="k",
        palette=itemgetter(*list(df[hue].unique()))(palette),
    )
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
        os.path.join(save_folder, "UMAP_{}_LNL_AD_Pso.pdf".format(hue)), bbox_inches="tight"
    )
    plt.close()


def main(adata, save_folder, path_to_df=None):
    if os.path.isfile(path_to_df):
        df = pd.read_excel(path_to_df, index_col=0)

        for hue_tmp in ['NonLesion_disease', 'spot_type']:
            palette = dict(zip(df["hue_tmp"], df["color_{}".format(hue_tmp)]))
            plot_scatter(df=df, hue=hue_tmp, palette=palette, save_folder=save_folder)

    else:
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
            'limegreen', 'mediumseagreen', 'darkgreen', '#e377c2', '#8c564b', '#aa40fc', '#b5bd61', '#17becf', '#aec7e8']))

        adata.obs['NonLesion_disease'] = adata.obs['DISEASE'].copy()
        adata.obs['NonLesion_disease'] = adata.obs['NonLesion_disease'].astype(str)
        adata.obs.loc[adata.obs['biopsy_type'] == 'NON LESIONAL', "NonLesion_disease"] = 'NON LESIONAL'
        adata.obs['NonLesion_disease'] = adata.obs['NonLesion_disease'].astype('category')
        # set order of categories
        adata.obs["NonLesion_disease"] = adata.obs["NonLesion_disease"].cat.reorder_categories(
            ['NON LESIONAL', 'AD', 'Pso'])

        color_NL_disease = dict(zip(adata.obs['NonLesion_disease'].cat.categories, ['teal', 'darkred', 'darkblue']))
        df = pd.DataFrame.from_dict({
            'UMAP1': adata.obsm['X_umap'][:, 0], 'UMAP2': adata.obsm['X_umap'][:, 1],
            'spot_type': adata.obs['spot_type'], 'NonLesion_disease': adata.obs["NonLesion_disease"],
            'color_spot_type': itemgetter(*list(adata.obs['spot_type']))(colors_spottype),
            'color_NonLesion_disease': itemgetter(*list(adata.obs['NonLesion_disease']))(color_NL_disease)})
        df.to_excel(path_to_df)

        # a)
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata=adata, color='NonLesion_disease', ax=ax, show=False, title='',
                   palette=['teal', 'darkred', 'darkblue'], s=20)
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

    path_df = os.path.join("/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer",
                           "output", "Fig_S1ab", "SFigure_1ab.xlsx")

    main(adata=spatial_adata, save_folder=savepath, path_to_df=path_df)
