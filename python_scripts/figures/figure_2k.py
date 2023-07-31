import scanpy as sc
import numpy as np
import pandas as pd
import xlsxwriter

import seaborn as sns
import matplotlib.pyplot as plt
import statannot
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator


import os
from datetime import date


fontsize_xylabel = 22
fontsize_xyticks = 20
fontsize_legend_title = 20
fontsize_text = 20
fontsize_legend = 18


def plot_boxplot(adata, df_melt, num_patterns, specimen, disease, biopsy_type, save_folder, obs='SEBACEOUS GLAND'):

    box_pairs = [((p, 0), (p, 1)) for p in range(num_patterns)]

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=df_melt, x='variable', y='value', hue=obs, ax=ax,
        palette=[adata.uns['spot_type_colors'][
            adata.obs['spot_type'].cat.categories.str.contains(obs)][0], 'grey'])
    ax, test_result = statannot.add_stat_annotation(
        ax, data=df_melt, x='variable', y='value', hue=obs, box_pairs=box_pairs,
        test='Mann-Whitney-gt', loc='inside', verbose=2, perform_stat_test=True, fontsize=fontsize_text)
    sns.despine(fig=fig, ax=ax)
    ax.set_xlabel('Pattern', fontsize=fontsize_xylabel)
    ax.set_ylabel('Relative expression levels', fontsize=fontsize_xylabel)
    ax.tick_params(axis='both', labelsize=fontsize_xyticks)
    # Put a legend to the right of the current axis
    leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fontsize_legend, frameon=False,
                    fancybox=False, shadow=False, ncol=1, title="".join(next(zip(*obs.split()))),
                    title_fontsize=fontsize_legend_title)
    leg.get_texts()[1].set_text('No')
    leg.get_texts()[0].set_text('Yes')
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Boxplot_Pattern_{}__{}_{}_{}.pdf'.format(
        num_patterns, specimen, disease, "".join(next(zip(*biopsy_type.split(' ')))), obs)))
    plt.close('all')

    plotting_parameters = {'data': df_melt, 'x': 'variable', 'y': 'value', 'hue': obs,
                           'palette': [adata.uns['spot_type_colors'][
                                           adata.obs['spot_type'].cat.categories.str.contains(obs)][0], 'grey']}
    # Transform each p-value to "p=" in scientific notation
    test_result_pvals = [result.__dict__['pval'] for result in test_result]
    formatted_pvalues = ["{:.2e}".format(pvalue) if pvalue < 0.05 else "ns" for pvalue in test_result_pvals]

    # adjust order of statsannot pvals to match normal order
    test_result_dicts = [result.__dict__ for result in test_result]
    df = pd.DataFrame(test_result_dicts)
    statsannot_order = np.asarray(list(df['box2'])).T[0]
    ind_neworder = np.argsort(statsannot_order)
    formatted_pvalues = np.asarray(formatted_pvalues)[ind_neworder]

    # Plot with formatted p-vales
    fig, ax = plt.subplots(figsize=(14, 6))
    sns.boxplot(**plotting_parameters)
    # Add annotations
    annotator = Annotator(ax, box_pairs, **plotting_parameters)
    annotator.set_custom_annotations(formatted_pvalues)
    # annotator.get_configuration()
    annotator.pvalue_format.fontsize = fontsize_text
    annotator.annotate()
    sns.despine(fig=fig, ax=ax)
    ax.set_xlabel('Pattern', fontsize=fontsize_xylabel)
    ax.set_ylabel('Relative expression levels', fontsize=fontsize_xylabel)
    ax.tick_params(axis='both', labelsize=fontsize_xyticks)
    # Put a legend to the right of the current axis
    leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fontsize_legend, frameon=False,
                    fancybox=False, shadow=False, ncol=1, title="".join(next(zip(*obs.split()))),
                    title_fontsize=fontsize_legend_title)
    leg.get_texts()[1].set_text('No')
    leg.get_texts()[0].set_text('Yes')
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'ReformatedBoxplot_Pattern_{}__{}_{}_{}.pdf'.format(
        num_patterns, specimen, disease, "".join(next(zip(*biopsy_type.split(' ')))), obs)))
    plt.close('all')

    return test_result


def calculate_enrichment_pattern_sebaceousglands(adata_sample, patterns, obs='SEBACEOUS GLAND'):
    df_patterns_annot = patterns.copy()
    df_patterns_annot[obs] = adata_sample.obs[obs]

    pvals = []
    # df_patterns_annot[obs] == 1 -> 1 if spot_type is in spot 0 if its not
    mask_sebaceousglands = df_patterns_annot[obs] == 1
    for pattern in patterns.columns:
        stats, pval = mannwhitneyu(
            x=df_patterns_annot.loc[mask_sebaceousglands, pattern],
            y=df_patterns_annot.loc[~mask_sebaceousglands, pattern],
            use_continuity=True, alternative='greater')
        pvals.append(pval)

    df_melt = pd.melt(df_patterns_annot, id_vars=[obs])

    return pvals, df_melt


def main(path_adata, save_folder):

    adata = sc.read(os.path.join(path_adata, 'st_QC_normed_BC_project_PsoAD.h5'))

    specimens = []
    writer = pd.ExcelWriter(os.path.join(save_folder, "Boxplots_Pattern.xlsx"), engine='xlsxwriter')

    # Read out specimen
    filename = '11-V19T12-012-V2_11_AD_NON LESIONAL'

    print('File: {}'.format(filename))

    specimen = "_".join(filename.split(sep='_')[:-2])
    adata_sample = adata[adata.obs['specimen'] == specimen].copy()
    specimens.append(specimen)
    biopsy_type = adata_sample.obs['biopsy_type'].cat.categories[0]
    disease = adata_sample.obs['DISEASE'].cat.categories[0]

    patterns = pd.DataFrame(index=adata_sample.obs.index, columns=np.arange(0, 9))
    for pattern in range(0, 9):
        patterns.loc[:, pattern] = adata_sample.obs.loc[:, 'Pattern_intensity_{}'.format(pattern)]

    pvals, df_melt = calculate_enrichment_pattern_sebaceousglands(
        adata_sample=adata_sample, patterns=patterns, obs='SEBACEOUS GLAND')

    df_melt['SEBACEOUS GLAND'] = df_melt['SEBACEOUS GLAND'].astype('category')
    df_melt['SEBACEOUS GLAND'] = df_melt['SEBACEOUS GLAND'].cat.reorder_categories([1, 0])
    df_melt['color'] = [
        adata.uns['spot_type_colors'][
            adata.obs['spot_type'].cat.categories.str.contains(
                'SEBACEOUS GLAND')][0] if val == 1 else "grey" for val in df_melt['SEBACEOUS GLAND']]
    test_result = plot_boxplot(adata, df_melt, num_patterns=9, specimen=specimen, disease=disease,
                               biopsy_type=biopsy_type, save_folder=save_folder, obs='SEBACEOUS GLAND')

    test_result = [result.__dict__ for result in test_result]
    df = pd.DataFrame(test_result)

    df['DeltaMean'] = 0
    df['mean_values_box1'] = 0
    df['mean_values_box2'] = 0
    for ind, (box1, box2) in enumerate(zip(df['box1'], df['box2'])):
        mask_1 = (df_melt['variable'] == box1[0]) & (df_melt['SEBACEOUS GLAND'] == box1[1])
        mask_2 = (df_melt['variable'] == box2[0]) & (df_melt['SEBACEOUS GLAND'] == box2[1])
        # Read out values
        box1_values = df_melt.loc[mask_1, 'value'].values
        box2_values = df_melt.loc[mask_2, 'value'].values
        # log2fc = log2(x/y) Fold change: (Y - X)/X
        # fold_change = (box2_values.mean() - box1_values.mean()) / box1_values.mean()
        mean_box1 = box1_values.mean()
        mean_box2 = box2_values.mean()
        # log2fc = np.log2(abs(box1_values.mean() / box2_values.mean()))
        delta_mean = (box1_values.mean() - box2_values.mean())
        df.loc[ind, 'DeltaMean'] = delta_mean
        df.loc[ind, 'mean_values_box1'] = mean_box1
        df.loc[ind, 'mean_values_box2'] = mean_box2

    # Save figure parameters to excel file
    df_melt.to_excel(writer, sheet_name="Plot_{}_{}_{}".format(
        specimen, disease, "".join(next(zip(*biopsy_type.split(' '))))), index=False)
    df.to_excel(writer, sheet_name="Stats_{}_{}_{}".format(
        specimen, disease, "".join(next(zip(*biopsy_type.split(' '))))), index=False)

    writer.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "figure_2k_spatialDE_SG_enrichment", str(today))
    os.makedirs(savepath, exist_ok=True)

    adata_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE/2023-04-12_paper_figures'

    main(path_adata=adata_path, save_folder=savepath)
