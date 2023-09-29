import NaiveDE
import SpatialDE
import scanpy as sc
from itertools import combinations
from scipy.spatial import distance
from sklearn.cluster import KMeans

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import statannot
from scipy.stats import mannwhitneyu
import math

import os
from datetime import date

# Tutorial: https://github.com/Teichlab/SpatialDE


def plot_geneexpression(coord: pd.DataFrame, norm_expr: pd.DataFrame, genes: list,
                        specimen: str, biopsy_type: str, save_folder: str):
    fig, ax = plt.subplots(1, len(genes), figsize=(6, 4))
    for i, g in enumerate(genes):
        ax[i].scatter(coord['x'], coord['y'], c=norm_expr[g], s=1)
        ax[i].set_title(g)
        ax[i].axis('equal')
        ax[i].invert_xaxis()
    plt.tight_layout()

    if len(genes) > 1:
        fig.savefig(os.path.join(save_folder, 'GeneExpression__{}_{}_{}.pdf'.format(
            specimen, biopsy_type, "_".join(genes))))
    else:
        fig.savefig(os.path.join(save_folder, 'GeneExpression__{}_{}_{}.pdf'.format(
            specimen, biopsy_type, genes)))
    plt.close(fig=fig)


def plot_spatial_variance(df_results: pd.DataFrame, specimen: str, biopsy_type: str, save_folder: str):
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.scatterplot(data=df_results, x="FSV", y="qval", hue="Gene_class", s=5,
                    hue_order=['Not significant', 'PER', 'SE', 'linear'],
                    palette=['black', 'darkorange', 'darkblue', 'forestgreen'], ax=ax, linewidth=0, alpha=1)
    ax.axhline(0.05, xmin=0, xmax=np.max(df_results['FSV']) + 1, c='black', lw=1, ls='--')
    ax.set_xlabel('Fraction spatial variance')
    ax.set_ylabel('Adj. P-value')
    ax.set_yscale('log')
    ax.invert_yaxis()
    sns.despine(fig=fig, ax=ax)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False,
              fancybox=False, shadow=False, ncol=1, title='Gene')
    fig.tight_layout()
    fig.savefig(os.path.join(save_folder, 'SpatialVariance__{}_{}.pdf'.format(specimen, biopsy_type)))
    plt.close()


def plot_mean_vs_variance(df: pd.DataFrame, specimen: str, biopsy_type: str, save_folder: str):
    plt.loglog()
    plt.scatter(df.mean(), df.var())
    xy = [1e-2, 1e3]
    plt.plot(xy, xy, c='C1')
    plt.plot(xy, np.square(xy), c='C2')
    plt.xlabel('Mean')
    plt.ylabel('Variance')
    plt.savefig(os.path.join(save_folder, 'Mean_vs_Variance__{}_{}'.format(specimen, biopsy_type)))
    plt.close()


def plot_lengthscales(results, mean_lengthscale, specimen, biopsy_type, disease, save_folder):

    df_plot = results[['qval', 'l', 'Gene_class']].query('qval < 0.05')
    df_plot = df_plot.sort_values('l')
    df_plot = df_plot.reset_index()
    df_plot['x'] = df_plot.index
    df_plot['Gene_class'] = df_plot['Gene_class'].cat.remove_unused_categories()
    df_plot['size'] = 1
    df_plot['size'] = df_plot.loc[df_plot['Gene_class'] == 'linear', 'l'].value_counts()

    # Per gene pattern get mean length scale
    mean_lengthscale_linear = df_plot.loc[df_plot['Gene_class'] == 'linear', 'l'].mean()
    mean_lengthscale_per = df_plot.loc[df_plot['Gene_class'] == 'PER', 'l'].mean()
    mean_lengthscale_se = df_plot.loc[df_plot['Gene_class'] == 'SE', 'l'].mean()

    # TODO add number of counts as size
    fig, ax = plt.subplots()
    sns.histplot(
        data=df_plot, x='l', color='grey', ax=ax, log_scale=(False, True))
    ax.axvline(mean_lengthscale, c='lightgrey', ls='--')
    sns.histplot(
        data=df_plot, x='l', hue='Gene_class', color=['blue', 'orange', 'green'],
        ax=ax, multiple="stack", log_scale=(False, True))
    ax.axvline(mean_lengthscale_linear, c='darkgreen', ls='--')
    ax.axvline(mean_lengthscale_se, c='darkorange', ls='--')
    ax.axvline(mean_lengthscale_per, c='darkblue', ls='--')
    ax.set_ylabel('counts')
    ax.set_xlabel('lengthscale')
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Lengthscale__{}_{}_{}.pdf'.format(specimen, disease, biopsy_type)))
    plt.close()

    fig, ax = plt.subplots()
    sns.boxplot(data=df_plot, y='l', x='Gene_class', palette=['blue', 'orange', 'green'], ax=ax)
    ax.set_ylabel('lengthscale')
    ax.set_xlabel('Gene pattern class')
    ax.axhline(mean_lengthscale_linear, c='darkgreen', ls='--')
    ax.axhline(mean_lengthscale_se, c='darkorange', ls='--')
    ax.axhline(mean_lengthscale_per, c='darkblue', ls='--')
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Lengthscale_boxplot__{}_{}_{}.pdf'.format(
        specimen, disease, biopsy_type)))
    plt.close()


def plot_boxplot(adata, df_melt, num_patterns, specimen, biopsy_type, save_folder, obs='SEBACEOUS GLAND'):

    box_pairs = [((p, 0), (p, 1)) for p in range(1, num_patterns + 1)]

    fig, ax = plt.subplots()
    sns.boxplot(
        data=df_melt, x='variable', y='value', hue=obs, ax=ax,
        palette=['grey', adata.uns['spot_type_colors'][
            adata.obs['spot_type'].cat.categories.str.contains(obs)][0]])
    ax, test_result = statannot.add_stat_annotation(
        ax, data=df_melt, x='variable', y='value', hue=obs, box_pairs=box_pairs,
        test='Mann-Whitney-ls', loc='inside', verbose=2, perform_stat_test=True)
    sns.despine(fig=fig, ax=ax)
    ax.set_xlabel('Pattern')
    ax.set_ylabel('Expression levels')
    # Put a legend to the right of the current axis
    leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                    fancybox=False, shadow=False, ncol=1, title=obs)
    leg.get_texts()[0].set_text('No')
    leg.get_texts()[1].set_text('Yes')
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Boxplot_Pattern_{}__{}_{}_{}.pdf'.format(
        num_patterns, specimen, biopsy_type, obs)))
    plt.close()

    # adjust order of statsannot pvals to match normal order
    test_result_pvals = [result.__dict__['pval'] for result in test_result]
    test_result_dicts = [result.__dict__ for result in test_result]
    df = pd.DataFrame(test_result_dicts)
    statsannot_order = np.asarray(list(df['box2'])).T[0]
    ind_neworder = np.argsort(statsannot_order)
    formatted_pvalues = np.asarray(test_result_pvals)[ind_neworder]

    return formatted_pvalues


def plot_aeh(coord_sample, histology_results, patterns, num_patterns, specimen, biopsy_type, save_folder, key=""):

    fig = plt.figure(figsize=(7, 3 + math.ceil(num_patterns / 3)))
    # Plot AEH patterns
    for i in range(num_patterns):
        ax = fig.add_subplot(3, math.ceil(num_patterns / 3), i + 1)
        img = ax.scatter(coord_sample['x'], coord_sample['y'], c=patterns[i + 1], s=2)
        ax.invert_yaxis()
        ax.set_title('Pattern {} - {} genes'.format(i + 1, histology_results.query('pattern == @i').shape[0]),
                     fontsize=8)
        cb = plt.colorbar(img, ax=ax)
        cb.ax.tick_params(labelsize=6)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'AEH_Pattern_{}__{}_{}_{}.pdf'.format(
        num_patterns, specimen, biopsy_type, key)))
    plt.close()


def plot_aeh_clustered(coord_sample, histology_results, patterns, num_patterns, specimen,
                       biopsy_type, save_folder, key=""):

    for i in range(num_patterns):
        km = KMeans(n_clusters=4, random_state=0)
        labels = km.fit_predict(np.asarray(patterns[i + 1]).reshape(-1, 1))

        fig = plt.figure(figsize=(5, 2))
        # Plot AEH pattern intensity
        ax = fig.add_subplot(1, 2, 1)
        img = ax.scatter(coord_sample['x'], coord_sample['y'], c=patterns[i + 1], s=5)
        ax.invert_yaxis()
        ax.set_title('Pattern {} - {} genes'.format(i + 1, histology_results.query('pattern == @i').shape[0]),
                     fontsize=6)
        cb = plt.colorbar(img, ax=ax)
        cb.ax.tick_params(labelsize=6)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        fig.tight_layout()

        # Plot KMeans Clusters
        ax = fig.add_subplot(1, 2, 2)
        unique_labels = np.unique(labels)
        for n in unique_labels:
            ax.scatter(coord_sample['x'][labels == n],
                       coord_sample['y'][labels == n], label=n, s=5)
        ax.invert_yaxis()
        ax.set_title('Clustering of genes', fontsize=6)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False, ncol=1,
                  fontsize=6)

        fig.tight_layout()
        plt.savefig(os.path.join(save_folder, 'AEH_Cluster_Pattern_{}__{}_{}_{}.pdf'.format(
            i + 1, specimen, biopsy_type, key)))
        plt.close()


def calculate_enrichment_pattern_sebaceousglands(adata_sample, patterns, obs='SEBACEOUS GLAND'):
    df_patterns_annot = patterns.copy()
    df_patterns_annot[obs] = adata_sample.obs[obs]

    pvals = []
    # df_patterns_annot[obs] == 1 -> 1 if spot_type is in spot 0 if it's not
    mask_sebaceousglands = df_patterns_annot[obs] == 1
    for pattern in patterns.columns:
        stats, pval = mannwhitneyu(
            x=df_patterns_annot.loc[mask_sebaceousglands, pattern],
            y=df_patterns_annot.loc[~mask_sebaceousglands, pattern],
            use_continuity=True, alternative='greater')
        pvals.append(pval)

    df_melt = pd.melt(df_patterns_annot, id_vars=[obs])

    return pvals, df_melt


def normalize_data(data):
    min_data = np.min(data)
    max_data = np.max(data)
    normed_data = (data - min_data) / (max_data - min_data)
    return normed_data


def similarity_pattern(patterns, n_pattern, highest_similarity):

    if n_pattern > 1:
        similarity_list = []
        comb_patterns = list(combinations(range(1, n_pattern + 1), r=2))
        for tuple_pattern in comb_patterns:
            # Normalise pattern intensities between 0 and 1
            normed_values_0 = normalize_data(data=np.asarray(patterns[tuple_pattern[0]].values))
            normed_values_1 = normalize_data(data=np.asarray(patterns[tuple_pattern[1]].values))

            # Calculate similarity between patterns
            dst = distance.euclidean(np.asarray(normed_values_0),
                                     np.asarray(normed_values_1))
            similarity_tmp = np.divide(1, 1 + dst)
            print("Similarity between {} and {} for #pattern {}: ".format(
                tuple_pattern[0], tuple_pattern[1], n_pattern), similarity_tmp)
            similarity_list.append(similarity_tmp)

        argmax_sim = np.nanargmax(similarity_list)
        print('Most similar patterns: {}, MSE: {}, #Pattern {}'.format(
            comb_patterns[argmax_sim], similarity_list[argmax_sim], n_pattern))
        max_similarity = similarity_list[argmax_sim]
        highest_similarity.append(max_similarity)
    else:
        highest_similarity.append(0)

    return highest_similarity


def get_aeh(adata, adata_sample, coord_sample, res, ms_results, n_pattern, mean_lengthscale, specimen,
            highest_similarity, histology_results_list):

    histology_results_tmp, patterns = SpatialDE.aeh.spatial_patterns(
        X=coord_sample, exp_mat=res[ms_results['g']], DE_mll_results=ms_results,
        C=n_pattern, l=mean_lengthscale, verbosity=1)
    # g: Gene, pattern: Gene assigned to pattern,
    # membership: Maximal posterior pattern assignment probability
    print(histology_results_tmp.head())

    # add + 1 to pattern to start counting at 1
    histology_results_tmp['pattern'] = histology_results_tmp['pattern'] + 1
    patterns.columns = patterns.columns + 1

    # Boxplot of Pattern Intensities with stats test
    pvals, df_melt = calculate_enrichment_pattern_sebaceousglands(
        adata_sample=adata_sample, patterns=patterns, obs='SEBACEOUS GLAND')

    # Add Patterns to adata.obs
    for ind_pval, pval in enumerate(pvals):
        adata.obs.loc[(adata.obs['specimen'] == specimen),
                      'Pattern_intensity_{}'.format(ind_pval + 1)] = patterns[ind_pval + 1]

    # It is usually interesting to see what the coexpressed genes determining a
    # histological pattern are:
    for i in histology_results_tmp.sort_values('pattern').pattern.unique():
        print('Pattern {}'.format(i))
        print(', '.join(
            histology_results_tmp.query('pattern == @i').sort_values('membership')['g'].tolist()))
        print()

    # Calculate MSE to make informed, objective decision on increasing number of patterns
    highest_similarity = similarity_pattern(
        patterns=patterns, n_pattern=n_pattern, highest_similarity=highest_similarity)

    histology_results_list.append(histology_results_tmp)

    return adata, histology_results_list, patterns, highest_similarity, df_melt, pvals


def main(path_adata, save_folder):
    h5files = list(filter(lambda fname: '.h5' in fname, os.listdir(path_adata)))
    h5files.remove('st_QC_normed_BC_project_PsoAD.h5')

    adata = sc.read(os.path.join(adata_path, 'st_QC_normed_BC_project_PsoAD.h5'))
    adata = adata[adata.obs['DISEASE'] != 'LP'].copy()

    n_pattern = 9
    for pattern_val in range(0, n_pattern):
        adata.obs['Pattern_intensity_{}'.format(pattern_val)] = 0

    # Store SVG and AEH results in dataframe
    specimens = []

    # Read out specimen
    for filename in h5files:
        print('File: {}'.format(filename))
        adata_sample = sc.read(os.path.join(path_adata, filename))
        # specimen = '2-V19S23-004-V3_2'  # P15509_1003
        specimen = "_".join(filename.split(sep='_')[:-1])
        specimens.append(specimen)
        biopsy_type = adata_sample.obs['biopsy_type'].cat.categories[0]
        disease = adata_sample.obs['DISEASE'].cat.categories[0]

        save_folder_tmp = os.path.join(save_folder, specimen)
        os.makedirs(save_folder_tmp, exist_ok=True)

        # Plot H&E image with spot types
        fig, ax = plt.subplots()
        sc.pl.spatial(adata=adata_sample, color='spot_type',
                      library_id=adata_sample.obs['sample'].cat.categories[0],
                      show=False, ax=ax)
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder_tmp, 'spot_type__{}_{}.pdf'.format(specimen, biopsy_type)))
        plt.close()

        # read out dataframe
        df = adata_sample.to_df(layer='counts')
        adata_sample.obs['total_counts'] = df.sum(axis=1)
        df.index = adata_sample.obs.index  # Align count matrix with metadata table
        coord_sample = adata_sample.obsm['spatial']
        coord_sample = pd.DataFrame.from_dict({'x': coord_sample[:, 0], 'y': coord_sample[:, 1]})
        sample_info = coord_sample.copy()
        sample_info['total_counts'] = adata_sample.obs['total_counts'].values

        # Plot data Mean vs Variance
        plot_mean_vs_variance(df=df, specimen=specimen, biopsy_type=biopsy_type, save_folder=save_folder_tmp)

        # 2. SpatialDE assumes normally distributed noise
        # -> approximately transform the data to normal distributed noise
        # Correct for library size or sequencing depth
        # total_counts must be a columns in .obs
        dfm = NaiveDE.stabilize(df.T).T
        res = NaiveDE.regress_out(sample_info=adata_sample.obs, expression_matrix=dfm.T,
                                  covariate_formula='np.log(total_counts)', design_formula='1', rcond=-1).T

        # =================================================================== #
        #                                 spatialDE                           #
        # =================================================================== #
        # Note: Takes about 15 minutes to run
        results = SpatialDE.run(X=coord_sample, exp_tab=res, kernel_space=None)
        # The most important columns in results are:
        # g - The name of the gene
        # pval - The P-value for spatial differential expression
        # qval - Significance after correcting for multiple testing
        # l - A parameter indicating the distance scale a gene changes expression over

        # Top spatialDE genes
        plot_geneexpression(coord=coord_sample, norm_expr=res, genes=results['g'].values[0:3],
                            save_folder=save_folder_tmp, specimen=specimen, biopsy_type=biopsy_type)
        # Not spatialDE genes
        plot_geneexpression(coord=coord_sample, norm_expr=res,
                            genes=results.sort_values('qval').tail(10)['g'].values[0:3],
                            save_folder=save_folder_tmp, specimen=specimen, biopsy_type=biopsy_type)

        if np.all(results['qval'] != 1):
            sign_results = results.query('qval < 0.05')
            m_sig = results.qval < 0.05
            # Group SVGs in periodic and linear by classifying DE genes to interpretable function classes
            ms_results = SpatialDE.model_search(
                X=coord_sample, exp_tab=res, DE_mll_results=sign_results, kernel_space=None)
            # 1. reorder ms_results
            sorter = results.loc[m_sig, 'g'].values
            df1 = ms_results.set_index('g')
            df1 = df1.reindex(sorter)
            df1 = df1.reset_index()
            df1.to_excel(os.path.join(save_folder_tmp, '{}_{}_{}.xlsx'.format(specimen, disease, biopsy_type)))
            # 2. add model to results
            results['Gene_class'] = 'Not significant'
            results.loc[m_sig, 'Gene_class'] = df1['model'].values
            results['Gene_class'] = results['Gene_class'].astype('category')

            # In regular differential expression analysis, we usually investigate the relation between
            # significance and effect size by 'so called' volcano plots. We don't have the concept of fold change
            # in our case, but we can investigate the fraction of variance explained by spatial variation.
            plot_spatial_variance(
                df_results=results, specimen=specimen, save_folder=save_folder_tmp, biopsy_type=biopsy_type)

            # =================================================================== #
            #                   Automatic expression histology (AEH)              #
            # =================================================================== #
            # AEH requires two parameters:
            # 1. Number of patterns C
            # 2. Characteristic lengthscale l for histological patterns
            # Guidance for picking the optimal lengthscale look at significant spatialDE genes
            mean_lengthscale = sign_results['l'].mean()

            # Plot lengthscales
            plot_lengthscales(results=results, mean_lengthscale=mean_lengthscale, specimen=specimen,
                              biopsy_type=biopsy_type, disease=disease, save_folder=save_folder_tmp)

            # Investigate linear, periodic and exponential decrease of pattern intensity
            for pattern_type in ['linear', 'PER', 'SE']:
                ms_results_type = ms_results[ms_results['model'] == pattern_type].copy()
                # patterns: The posterior mean underlying expression for genes in given spatial patterns.
                histology_results_type, patterns_type = SpatialDE.aeh.spatial_patterns(
                    X=coord_sample, exp_mat=res[ms_results_type['g']], DE_mll_results=ms_results_type,
                    C=1, l=mean_lengthscale, verbosity=1)

                # add + 1 to pattern to start counting at 1
                histology_results_type["pattern"] = histology_results_type["pattern"] + 1
                patterns_type.columns = patterns_type.columns + 1

                ms_results_type.to_excel(os.path.join(save_folder_tmp, '{}_{}_{}_{}.xlsx'.format(
                    specimen, disease, biopsy_type, pattern_type)))

                # Plot AEH
                plot_aeh(coord_sample=coord_sample, histology_results=histology_results_type,
                         patterns=patterns_type, num_patterns=1, specimen=specimen, biopsy_type=biopsy_type,
                         save_folder=save_folder_tmp, key=pattern_type)
                plot_aeh_clustered(
                    coord_sample=coord_sample, histology_results=histology_results_type, patterns=patterns_type,
                    num_patterns=1, specimen=specimen, biopsy_type=biopsy_type, save_folder=save_folder_tmp,
                    key=pattern_type)

            # Init parameters
            histology_results_list = []
            highest_similarity = []
            pvals_list = []

            # Run AEH
            adata, histology_results_list, patterns, highest_similarity, df_melt, pvals = get_aeh(
                adata, adata_sample, coord_sample, res, ms_results, n_pattern, mean_lengthscale,
                specimen, highest_similarity, histology_results_list)

            # Boxplot of Pattern Intensities with stats test
            padj_vals = plot_boxplot(
                adata=adata, df_melt=df_melt, num_patterns=n_pattern, specimen=specimen,
                biopsy_type=biopsy_type, save_folder=save_folder_tmp, obs='SEBACEOUS GLAND')
            pvals_list.append(padj_vals)

            # Plot AEH
            plot_aeh(coord_sample=coord_sample, histology_results=histology_results_list[0],
                     patterns=patterns, num_patterns=n_pattern, specimen=specimen, biopsy_type=biopsy_type,
                     save_folder=save_folder_tmp)
            plot_aeh_clustered(
                coord_sample=coord_sample, histology_results=histology_results_list[0],
                patterns=patterns, num_patterns=n_pattern, specimen=specimen, biopsy_type=biopsy_type,
                save_folder=save_folder_tmp)

            similarity_step = highest_similarity[0]

            print("Highest similarity of {} for {} pattern".format(similarity_step, n_pattern))

            # Store histology_results
            histology_results = histology_results_list[0]
            del histology_results_list
            # Append histology_results
            histology_results['specimen'] = specimen

            # Add information of whether pattern is enriched in Sebaceous Gland
            histology_results['SEBACEOUS GLAND'] = 0
            for ind_pval, pval in enumerate(pvals_list[0]):
                if pval < 0.05:
                    histology_results.loc[histology_results['pattern'] == ind_pval + 1, 'SEBACEOUS GLAND'] = 1

            filename_aeh = os.path.join(save_folder_tmp, 'AEH__{}.csv'.format(specimen))
            histology_results.to_csv(filename_aeh)

            # Save results of SVG
            results['specimen'] = specimen
        else:
            results['specimen'] = specimen

        # Save dataframe
        filename_svg = os.path.join(save_folder_tmp, 'SVG__{}.csv'.format(specimen))
        results.to_csv(filename_svg)

        print('=============================\n')

    # Save adata
    sc.write(os.path.join(save_folder, 'st_QC_normed_BC_project_PsoAD.h5'), adata=adata)


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "spatialDE", str(today))
    os.makedirs(savepath, exist_ok=True)

    adata_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE'

    main(path_adata=adata_path, save_folder=savepath)
