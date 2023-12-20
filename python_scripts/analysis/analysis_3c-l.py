import scanpy as sc
import numpy as np
import pandas as pd
from collections import Counter

import seaborn as sns
import matplotlib.pyplot as plt
import statannot

import os


def prepare_dataframe(adata_sample, patterns, obs='SEBACEOUS GLAND'):
    df_patterns_annot = patterns.copy()
    df_patterns_annot[obs] = adata_sample.obs[obs]

    df_melt = pd.melt(df_patterns_annot, id_vars=[obs])

    return df_melt


def get_enrichment(adata, df_melt, num_patterns, obs='SEBACEOUS GLAND'):
    box_pairs = [((p, 0), (p, 1)) for p in range(num_patterns)]

    fig, ax = plt.subplots()
    sns.boxplot(
        data=df_melt, x='variable', y='value', hue=obs, ax=ax,
        palette=['grey', adata.uns['spot_type_colors'][
            adata.obs['spot_type'].cat.categories.str.contains(obs)][0]])
    ax, test_result = statannot.add_stat_annotation(
        ax, data=df_melt, x='variable', y='value', hue=obs, box_pairs=box_pairs,
        test='Mann-Whitney-ls', loc='inside', verbose=2, perform_stat_test=True, comparisons_correction='bonferroni')

    plt.close(fig=fig)

    # adjust order of statsannot pvals to match normal order
    test_result_pval_corrected = [result.__dict__['pval'] for result in test_result]
    test_result_dicts = [result.__dict__ for result in test_result]
    df = pd.DataFrame(test_result_dicts)
    statsannot_order = np.asarray(list(df['box2'])).T[0]
    ind_neworder = np.argsort(statsannot_order)
    formatted_pval_corrected = np.asarray(test_result_pval_corrected)[ind_neworder]

    return formatted_pval_corrected


def get_commonly_regulated_genes(a, b):
    print("Get commonly regulated genes between 2 groups")
    intersection_genes = set(a) & set(b)

    return intersection_genes


def get_uniquely_regulated_genes(a, rem_b):
    print("Get uniquely regulated genes in a group")
    unique_reg_genes = set(a) - set(rem_b)

    return unique_reg_genes


def main(path_adata, save_folder):
    h5files = [name for name in os.listdir(path_adata) if os.path.isdir(os.path.join(path_adata, name))]
    h5files.remove('Pathway_enrichment_analysis')

    adata = sc.read(os.path.join(path_adata, 'st_QC_normed_BC_project_PsoAD.h5'))

    svg = {'AD': [], 'Pso': []}

    # Read out specimen
    for filename in h5files:
        print('File: {}'.format(filename))
        # filename = 'SN-V11J13-122_B_20_AD_LESIONAL'

        specimen = "_".join(filename.split(sep='_')[:-2])
        adata_sample = adata[adata.obs['specimen'] == specimen].copy()
        biopsy_type = adata_sample.obs['biopsy_type'].cat.categories[0]
        disease = adata_sample.obs['DISEASE'].cat.categories[0]

        patterns = pd.DataFrame(index=adata_sample.obs.index, columns=np.arange(0, 9))
        for pattern in range(0, 9):
            patterns.loc[:, pattern] = adata_sample.obs.loc[:, 'Pattern_intensity_{}'.format(pattern)]

        df_melt = prepare_dataframe(adata_sample=adata_sample, patterns=patterns, obs='SEBACEOUS GLAND')
        pvals = get_enrichment(adata, df_melt, num_patterns=len(patterns.columns), obs='SEBACEOUS GLAND')
        pvals = np.asarray(pvals)
        mask_ind = np.argwhere(pvals < 0.05).T[0]
        enriched_patterns = list(patterns.columns[mask_ind])

        # save to excel sheet
        df_melt.to_excel(os.path.join(save_folder, 'Pattern_information_boxplot.xlsx'.format(specimen)))

        # Get pattern specific SVGs
        histology_results = pd.read_csv(os.path.join(path_adata, filename, 'AEH__{}.csv'.format(specimen)), index_col=0)
        histology_results.loc[(histology_results['SEBACEOUS GLAND'] == 1) & ~(
            histology_results['pattern'].isin(enriched_patterns)), 'SEBACEOUS GLAND'] = 0
        histology_results.to_csv(os.path.join(path_adata, filename, 'AEH__{}_corrected.csv'.format(specimen)))
        histology_results.to_excel(os.path.join(path_adata, filename, 'AEH__{}_corrected.xlsx'.format(specimen)))

        # We need SVGs per pattern which are enriched in SGs
        svgs = list(histology_results.loc[(histology_results['SEBACEOUS GLAND'] == 1) & (
            histology_results['pattern'].isin(enriched_patterns)), 'g'].values)

        if biopsy_type != 'NON LESIONAL':
            svg[disease].extend(svgs)

    # Keep only genes which can be found in all samples of same disease
    for diag in svg.keys():
        counts_genes = Counter(svg[diag])
        max_occurence = np.amax(list(counts_genes.values()))
        svg[diag] = np.asarray(list(counts_genes.keys()))[list(counts_genes.values()) == max_occurence]

    # Get common SVGs between Pso and AD
    com_genes = get_commonly_regulated_genes(a=svg['Pso'], b=svg['AD'])
    # Get unique SVGs
    unique_disease = {'AD': list(get_uniquely_regulated_genes(a=svg['AD'], rem_b=com_genes)),
                      'Pso': list(get_uniquely_regulated_genes(a=svg['Pso'], rem_b=com_genes))}
    df_unique = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in unique_disease.items()]))

    svgs_partioned = {'Shared': list(com_genes),
                      'AD': list(get_uniquely_regulated_genes(a=svg['AD'], rem_b=com_genes)),
                      'Pso': list(get_uniquely_regulated_genes(a=svg['Pso'], rem_b=com_genes))}
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in svgs_partioned.items()]))

    # Save results
    df_unique.to_excel(os.path.join(save_folder, 'SVGs_unique_AD_PSO.xlsx'))
    df.to_excel(os.path.join(save_folder, 'SVGs_shared_AD_PSO.xlsx'))


if __name__ == '__main__':
    # create saving folder in current project path
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "Analysis_3cl__spatialDE_SVGs")
    os.makedirs(savepath, exist_ok=True)

    adata_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE/2023-09-18_paper_figures'

    main(path_adata=adata_path, save_folder=savepath)
