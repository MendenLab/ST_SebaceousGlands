import scanpy as sc
import pandas as pd
import numpy as np
from operator import itemgetter

import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator

import os
from datetime import date


def plot_boxplot(df_melt, palette, box_pairs, plotting_parameters, formatted_pvalues, save_folder):
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    sns.boxplot(
        data=df_melt,
        x="gene_name",
        y="value",
        hue="condition",
        fliersize=0,
        ax=ax,
        palette=palette,
    )
    sns.stripplot(
        data=df_melt,
        x="gene_name",
        y="value",
        hue="condition",
        jitter=True,
        dodge=True,
        marker=".",
        color="black",
        ax=ax
    )
    # Add annotations
    annotator = Annotator(ax, box_pairs, **plotting_parameters)
    annotator.set_custom_annotations(formatted_pvalues)
    # annotator.get_configuration()
    annotator.pvalue_format.fontsize = 12
    annotator.annotate()
    sns.despine(fig=fig, ax=ax)
    ax.set_ylabel("normed counts", fontsize=18)
    ax.set_xlabel("", fontsize=18)
    ax.tick_params(labelsize=14)
    # Put a legend below current axis
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[0:2], labels[0:2],
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
    leg.get_texts()[1].set_text("rest of skin")
    leg.get_texts()[0].set_text("SG")
    plt.tight_layout()

    plt.savefig(
        os.path.join(save_folder, "Boxplots_normed_AD_Pso.pdf"), bbox_inches="tight"
    )
    plt.close()


def main(adata, save_folder, path_to_df):
    genes = ["S100A8", "KRT79", "SAA1", "FABP7", "HSD3B1"]
    box_pairs = [((p, 0), (p, 1)) for p in genes]

    if os.path.isfile(path_to_df):
        df_melt = pd.read_excel(path_to_df, sheet_name='boxplot')
        df_formatted_pvalues = pd.read_excel(path_to_df, sheet_name="formatted_pvalues")
        formatted_pvalues = df_formatted_pvalues['formatted_pvalues'].to_list()

        palette = dict(zip(df_melt['condition'], df_melt['color']))
        palette = list(itemgetter(*list(df_melt['condition'].unique()))(palette))
        plotting_parameters = {'data': df_melt, 'x': 'gene_name', 'y': 'value', 'hue': 'condition',
                               'palette': palette}

        plot_boxplot(df_melt, palette, box_pairs, plotting_parameters, formatted_pvalues, save_folder)

    else:
        # Comparison AD_L_SEB_vs_RestofSkin
        adata.uns['spot_type_colors'] = np.asarray(
            ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#e377c2', '#8c564b',
             '#aa40fc', '#b5bd61', '#17becf', '#aec7e8'])
        palette = [adata.uns['spot_type_colors'][
                adata.obs['spot_type'].cat.categories.str.contains('SEBACEOUS GLAND')][0], 'grey']

        df_deg = pd.read_excel(os.path.join(
            "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
            'switch_signs_log2fc/2021-09-04/switch_signs/ST_cdr_patient_annotation_condition',
            '13_DGE_Approaches_Signature_AE_modifiedPS_20210507/L, SEB G__condition 1_vs_condition 2',
            'ST_condition 1_vs_condition 2_glmGamPoi_DGE_all_genes.xlsx'), index_col=0)
        df_deg = df_deg.loc[df_deg['gene_symbol'].isin(genes), :][['gene_symbol', 'padj', 'log2fc']]
        df_deg.index = df_deg['gene_symbol']
        # reorder rows
        df_deg = df_deg.loc[genes, :]

        df = adata.to_df()
        df = df[genes].copy()
        df = np.log2(df)
        df['condition'] = adata.obs['SEBACEOUS GLAND']

        df_melt = pd.melt(df, id_vars='condition')
        df_melt['condition'] = df_melt['condition'].astype('category')
        df_melt['condition'] = df_melt['condition'].cat.reorder_categories([1, 0])

        plotting_parameters = {'data': df_melt, 'x': 'gene_name', 'y': 'value', 'hue': 'condition',
                               'palette': palette}
        # Transform each p-value to "p=" in scientific notation
        formatted_pvalues = ["".join([r"log$_2$fc=", "{:.2f}".format(df_deg.iloc[value]['log2fc']),
                                     "\npadj={:.2e}".format(df_deg.iloc[value]['padj'])])
                             if df_deg.iloc[value]['padj'] < 0.05 else "ns" for value in range(df_deg.shape[0])]

        plot_boxplot(df_melt, palette, box_pairs, plotting_parameters, formatted_pvalues, save_folder)

        color = dict(zip(df_melt['condition'].cat.categories, palette))
        df_melt['color'] = itemgetter(*list(df_melt['condition']))(color)
        df_formatted_pvalues = pd.DataFrame(formatted_pvalues, columns=['formatted_pvalues'])

        writer = pd.ExcelWriter(path_to_df, engine="xlsxwriter")
        df_melt.to_excel(writer, sheet_name='boxplot')
        df_formatted_pvalues.to_excel(writer, sheet_name='formatted_pvalues')
        writer.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "figure_1e", str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Remove LP & Pso
    spatial_adata = spatial_adata[spatial_adata.obs['DISEASE'] != 'LP'].copy()

    path_df = os.path.join("/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer",
                           "output", "figure_1e", "Figure_1e.xlsx")
    main(adata=spatial_adata, save_folder=savepath, path_to_df=path_df)
