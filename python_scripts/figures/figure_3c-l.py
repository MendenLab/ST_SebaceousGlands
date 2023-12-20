from python_scripts.utils import gene_lists

import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import os
from datetime import date


fontsize_xylabel = 14
fontsize_xyticks = 12
fontsize_legend_title = 12
fontsize_text = 10
fontsize_legend = 10
fontsize_axis_title = 8


def plot_svgs(
        adata, auto_corr_genes, disease, sample, biopsy_type, specimen, save_folder,
        invert_x, invert_y, key, figure_sizes):

    for gene in auto_corr_genes:
        fig, ax = plt.subplots(figsize=figure_sizes)
        sc.pl.spatial(
            adata, color=gene, size=1, cmap='Reds', library_id=sample, show=False, ax=ax,
            layer='counts')
        ax.set_title(gene, fontsize=16)
        if invert_x:
            ax.invert_xaxis()
        if invert_y:
            ax.invert_yaxis()

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        # plt.tight_layout()
        plt.savefig(os.path.join(save_folder, 'spatialDE_{}_SVG_{}_{}_{}_{}.pdf'.format(
            key, gene, specimen, disease, biopsy_type)))
        plt.close('all')


def plot_svgs_from_df(
        df, auto_corr_genes, disease, biopsy_type, specimen, save_folder,
        invert_x, invert_y, key, figure_sizes):

    for gene in auto_corr_genes:
        fig, ax = plt.subplots(figsize=figure_sizes)
        sns.scatterplot(
            data=df,
            x="spatial_x",
            y="spatial_y",
            hue=gene,
            s=60,
            ax=ax,
            linewidth=0,
            edgecolor="k", cmap='Reds'
        )
        ax.set_title(gene, fontsize=16)
        if invert_x:
            ax.invert_xaxis()
        if invert_y:
            ax.invert_yaxis()

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        # plt.tight_layout()
        plt.savefig(os.path.join(save_folder, 'spatialDE_{}_SVG_{}_{}_{}_{}.pdf'.format(
            key, gene, specimen, disease, biopsy_type)))
        plt.close('all')


def main(path_adata, input_dir, save_folder, path_to_df):
    h5files = [name for name in os.listdir(path_adata) if os.path.isdir(os.path.join(path_adata, name))]
    h5files.remove('Pathway_enrichment_analysis')

    invert_x = {'2-V19S23-004-V4_2': True, 'SN-V11J13-122_B_20': False, '11-V19T12-012-V2_11': True,
                '10-V19T12-025-V3_10': True, '10-V19T12-025-V4_10': True, '11-V19T12-012-V1_11': True,
                '2-V19S23-004-V3_2': True}

    invert_y = {'2-V19S23-004-V4_2': True, 'SN-V11J13-122_B_20': True, '11-V19T12-012-V2_11': False,
                '10-V19T12-025-V3_10': True, '10-V19T12-025-V4_10': True, '11-V19T12-012-V1_11': False,
                '2-V19S23-004-V3_2': True}
    # (width, height)
    figure_sizes = {'2-V19S23-004-V4_2': (6, 6), 'SN-V11J13-122_B_20': (6, 6), '11-V19T12-012-V2_11': (6, 5),
                    '10-V19T12-025-V3_10': (9, 7), '10-V19T12-025-V4_10': (8, 8), '11-V19T12-012-V1_11': (5, 4),
                    '2-V19S23-004-V3_2': (6, 6)}
    # Get common SVGs between Pso and AD
    df = pd.read_excel(os.path.join(input_dir, "SVGs_shared_AD_PSO.xlsx"))
    com_genes = df["Shared"]
    # Get unique SVGs of AD and Pso
    df_unique = pd.read_excel(os.path.join(input_dir, "SVGs_unique_AD_PSO.xlsx"))

    if os.path.isfile(path_to_df):
        for filename in h5files:
            print("File: {}".format(filename))
            specimen = "_".join(filename.split(sep="_")[:-2])

            df_counts = pd.read_excel(path_to_df, index_col=0, sheet_name=specimen)
            specimen = df_counts['specimen'].unique()[0]
            biopsy_type = df_counts['biopsy_type'].unique()[0]
            disease = df_counts['disease'].unique()[0]

            if biopsy_type != "NON LESIONAL":
                save_folder_tmp = os.path.join(save_folder, specimen)
                os.makedirs(save_folder_tmp, exist_ok=True)

                df_svg_sig = pd.read_excel(
                    os.path.join(path_adata, filename, '{}_{}_{}.xlsx'.format(
                        specimen, disease, biopsy_type)), index_col=0)
                df_svg_sig = df_svg_sig[df_svg_sig['g'].isin(df_unique[disease])]
                # Plot most significant SVGs among the SVGs which belong to patterns enriched in SGs
                df_svg_sig_sorted = df_svg_sig.sort_values(['qval'], ascending=True)

                var_names = df.columns[~df.columns.isin(
                    ['spatial_x', 'spatial_y', 'disease', 'biopsy_type', 'specimen'])]
                # Most significant SVGs
                svg_sig_sorted = list(np.intersect1d(list(df_svg_sig_sorted.loc[:, 'g'][:10]), var_names))
                plot_svgs_from_df(
                    df=df_counts, auto_corr_genes=svg_sig_sorted,
                    disease=disease, biopsy_type=biopsy_type, specimen=specimen,
                    save_folder=save_folder_tmp,
                    invert_x=invert_x[specimen], invert_y=invert_y[specimen], key='Most_significant_SVGs',
                    figure_sizes=figure_sizes[specimen])

                # Golden standard disease specific SVGs
                gs_genes = gene_lists.goldenstandard_sebcaeous_glands_genes()
                gs_genes = list(np.intersect1d(list(gs_genes[disease]), var_names))
                plot_svgs_from_df(
                    df=df_counts, auto_corr_genes=gs_genes,
                    disease=disease, biopsy_type=biopsy_type, specimen=specimen, save_folder=save_folder_tmp,
                    invert_x=invert_x[specimen], invert_y=invert_y[specimen],
                    key='Golden_Standard_Disease_specific_SVGs', figure_sizes=figure_sizes[specimen])

                # Shared SVGs
                com_genes = list(np.intersect1d(list(com_genes), var_names))
                plot_svgs_from_df(
                    df=df_counts, auto_corr_genes=com_genes,
                    disease=disease, biopsy_type=biopsy_type, specimen=specimen, save_folder=save_folder_tmp,
                    invert_x=invert_x[specimen], invert_y=invert_y[specimen],
                    key='Shared_SVGs', figure_sizes=figure_sizes[specimen])
    else:
        adata = sc.read(os.path.join(path_adata, 'st_QC_normed_BC_project_PsoAD.h5'))
        writer = pd.ExcelWriter(os.path.join(path_to_df), engine='xlsxwriter')

        for filename in h5files:
            print('File: {}'.format(filename))

            specimen = "_".join(filename.split(sep='_')[:-2])
            adata_sample = adata[adata.obs['specimen'] == specimen].copy()
            biopsy_type = adata_sample.obs['biopsy_type'].cat.categories[0]
            disease = adata_sample.obs['DISEASE'].cat.categories[0]
            sample = adata_sample.obs['sample'].cat.categories[0]

            if biopsy_type != 'NON LESIONAL':
                save_folder_tmp = os.path.join(save_folder, specimen)
                os.makedirs(save_folder_tmp, exist_ok=True)

                df_counts = pd.DataFrame.from_dict({'spatial_x': adata_sample.obsm[
                    'spatial'][:, 0] * adata_sample.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'],
                                         'spatial_y': adata_sample.obsm[
                    'spatial'][:, 1] * adata_sample.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef']})
                df_counts["disease"] = disease
                df_counts["biopsy_type"] = biopsy_type
                df_counts["specimen"] = specimen

                df_svg_sig = pd.read_excel(
                    os.path.join(path_adata, filename, '{}_{}_{}.xlsx'.format(
                        specimen, disease, biopsy_type)), index_col=0)
                df_svg_sig = df_svg_sig[df_svg_sig['g'].isin(df_unique[disease])]

                # Plot most significant SVGs among the SVGs which belong to patterns enriched in SGs
                df_svg_sig_sorted = df_svg_sig.sort_values(['qval'], ascending=True)

                # Check if genes are expressed on specimen
                svg_sig_sorted = list(np.intersect1d(list(df_svg_sig_sorted.loc[:, 'g'][:10]), adata.var_names))
                plot_svgs(adata=adata_sample, auto_corr_genes=svg_sig_sorted,
                          disease=disease, sample=sample, biopsy_type=biopsy_type, specimen=specimen,
                          save_folder=save_folder_tmp,
                          invert_x=invert_x[specimen], invert_y=invert_y[specimen], key='Most_significant_SVGs',
                          figure_sizes=figure_sizes[specimen])

                # Golden standard disease specific SVGs
                gs_genes = gene_lists.goldenstandard_sebcaeous_glands_genes()
                gs_genes = list(np.intersect1d(list(gs_genes[disease]), adata.var_names))
                plot_svgs(adata=adata_sample, auto_corr_genes=gs_genes,
                          disease=disease, sample=sample, biopsy_type=biopsy_type, specimen=specimen,
                          save_folder=save_folder_tmp,
                          invert_x=invert_x[specimen], invert_y=invert_y[specimen],
                          key='Golden_Standard_Disease_specific_SVGs', figure_sizes=figure_sizes[specimen])

                df_svg_sig_sorted.to_excel(os.path.join(save_folder_tmp, 'Disease_specific_SVGs_sorted.xlsx'))

                # Shared SVGs
                com_genes = list(np.intersect1d(list(com_genes), adata.var_names))
                plot_svgs(adata=adata_sample, auto_corr_genes=com_genes,
                          disease=disease, sample=sample, biopsy_type=biopsy_type, specimen=specimen,
                          save_folder=save_folder_tmp,
                          invert_x=invert_x[specimen], invert_y=invert_y[specimen],
                          key='Shared_SVGs', figure_sizes=figure_sizes[specimen])

                # Add raw gene counts to dataframe
                for genes in [svg_sig_sorted, gs_genes, com_genes]:
                    df_gene_counts_tmp = adata_sample.to_df(layer='counts')[genes]
                    df_gene_counts_tmp = df_gene_counts_tmp.reset_index(drop=True)
                    df_counts = pd.concat([df_counts, df_gene_counts_tmp], axis=1)

                # Save figure parameters to Excel file
                df_counts.to_excel(writer, sheet_name=specimen, index=False)
        writer.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "Figure_3cl__spatialDE_SVGs", str(today))
    os.makedirs(savepath, exist_ok=True)

    input_path = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "Analysis_3cl__spatialDE_SVGs")

    adata_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE/2023-09-18_paper_figures'

    path_df = os.path.join("/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer",
                           "output", "Figure_3cl__spatialDE_SVGs", "Figure_3cl.xlsx")

    main(path_adata=adata_path, save_folder=savepath, input_dir=input_path, path_to_df=path_df)
