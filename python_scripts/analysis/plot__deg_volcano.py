#!/usr/bin/env python
"""Create Volcano and Violin plots using the DGE Analysis result
    File name: Fig3B__ST_Volcanoplot_cytokine.py
    Author: Christina Hillig
    Date created: November/03/2020
    Date last modified: April/30/2021
    Python Version: 3.7
"""

from python_scripts.utils import add_observables, anndata_information, gene_lists

import os
import glob
from fnmatch import fnmatch
from datetime import date

import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from adjustText import adjust_text
import plotly.graph_objects as go

# source: https://www.genomics-online.com/resources/16/5049/housekeeping-genes/
_house_keeping_genes = ['ACTB', 'GAPDH', 'PGK1', 'PPIA', 'RPLP0', 'ARBP', 'B2M', 'YWHAZ', 'SDHA', 'TFRC', 'GUSB',
                        'HMBS', 'HPRT1', 'TBP']

# plot params
dotsize = 8
figure_size = (8, 8)
xy_fontsize = 16
xy_ticks = 12
title_fontsize = 18
legend_fontsize = 14
text_fontsize = 12
file_format = '.pdf'
img_key = 'hires'


def load_adata(type_dataset, cluster_label):
    """Load annData object

    Parameters
    ----------
    type_dataset : str
        data set type (ST or scRNA-seq)
    cluster_label : str
        name of observable
    Returns
    -------
    adata : annData

    """
    adata_path = os.path.join("..", 'adata_storage/2021-05-18')
    if type_dataset == 'ST':
        adata = sc.read(os.path.join(adata_path, 'annData_PS_tasks.h5'))
        # remove all spots without a tissue label
        adata = adata[adata.obs[cluster_label] != 'Unknown']
    elif type_dataset == 'SC':
        adata = sc.read(os.path.join(adata_path, '2020-12-04_SC_Data_QC_clustered.h5'))
    else:
        raise ValueError("Please provide either 'SC' or 'ST' as input for type_dataset")

    return adata


def get_expression_values(adata, gene):
    """Add gene count, label and cluster annotation as column to annData

    Parameters
    ----------
    adata : annData
    gene : str or 'list' ['str']

    Returns
    -------

    """

    # Get observable for gene
    adata, obs_name = add_observables.convert_variable_to_observable(
        adata=adata, gene_names=gene, task='cell_gene', label='celltype', condition=None)

    return adata, obs_name


def exclude_cytokine_dp(adata, condition_resps_dict):
    """Exclusively label double positive cells as double positive

    Parameters
    ----------
    adata : annData
    condition_resps_dict : dict

    Returns
    -------
    adata : annData
    obs_name : str

    """
    # check if all genes are in adata
    for skin_layer in condition_resps_dict.keys():
        genes = list(set(adata.var.index) & set(condition_resps_dict[skin_layer]))
        condition_resps_dict[skin_layer] = genes

    # add label - todo
    adata, obs_name = add_observables.convert_variable_to_observable(
        adata=adata, gene_names=condition_resps_dict, task='annotate_cells',
        label='gsg_epidermis', condition=[np.all, np.all, np.all])

    return adata, obs_name


def create_obs_cytopos_cytoneg(adata, cyto, gene, observable):
    """Add observable to annData object

    Parameters
    ----------
    adata : AnnData
    cyto : str
    gene : str
    observable : str

    Returns
    -------

    """
    if gene in adata.var.index:
        counts = np.copy(adata.layers['counts'])[:, np.where(adata.var.index == gene)[0]][:, 0]
        m_gene = counts > 0
        m_gene_cytopos = np.logical_and(m_gene, adata.obs[observable].values == cyto)
        obs_name = "_".join([gene, 'group'])
        adata.obs[obs_name] = counts
        adata.obs[obs_name][m_gene_cytopos] = counts[m_gene_cytopos]
        adata.obs[obs_name][~m_gene_cytopos] = counts[~m_gene_cytopos]
        adata.obs[obs_name] = adata.obs[obs_name]

    return adata


def get_labels(df, annotations, log2fc='log2fc', pval='pval', log2fc_cut=1.0, threshold=0.05, xlim=None, ylim=None):
    """Get label positions as log2fc and p-value and colors, blue for Responder genes and purple for Dirver genes

    Parameters
    ----------
    df : pandas.Dataframe
    annotations : dict or list
    log2fc : str
    pval : str
    log2fc_cut : float
    threshold : float
    xlim : list, optional
    ylim : list, optional

    Returns
    -------

    """
    label_info = pd.DataFrame()
    if isinstance(annotations, dict):
        driver_group = annotations['Driver_genes']
        responder_group = annotations['Responder_genes']
        merged_genes = driver_group.copy()
        merged_genes.extend(responder_group)
        ind_genes = np.where(df['gene_symbol'][np.newaxis, :] == np.array(merged_genes)[:, np.newaxis])[1]
        label_df = df.iloc[ind_genes]

        # get position and color for genes of interest
        label_info['x'], label_info['y'] = label_df[log2fc], -np.log10(label_df[pval])
        label_info['gene_symbol'] = label_df['gene_symbol'].values
        label_info['color'] = ['mediumblue'] * len(label_info.index)
        ind_drivers = np.where(label_info['gene_symbol'][np.newaxis, :] == np.array(driver_group)[:, np.newaxis])[1]
        label_info['color'].values[ind_drivers] = ['purple'] * len(ind_drivers)
    else:
        # Using loop + keys()
        res_list = []
        for sublist in list(annotations):
            res_list.extend(sublist)

        ind_genes = np.where(df['gene_symbol'][np.newaxis, :] == np.array(res_list)[:, np.newaxis])[1]
        label_df = df.iloc[ind_genes]
        if (xlim is not None) & (ylim is not None):
            m_text = (abs(label_df[log2fc]) < max(xlim)) & (label_df[pval] < max(ylim))
            # get position and color for genes of interest
            label_info['x'], label_info['y'] = label_df[log2fc][m_text], -np.log10(label_df[pval][m_text])
            label_info['gene_symbol'] = label_df['gene_symbol'][m_text]
            label_info['color'] = ['k'] * len(label_info.index)
        else:
            m_text = (abs(label_df[log2fc]) > log2fc_cut) & (label_df[pval] < threshold)
            # get position and color for genes of interest
            label_info['x'], label_info['y'] = label_df[log2fc][m_text], -np.log10(label_df[pval][m_text])
            label_info['gene_symbol'] = label_df['gene_symbol'][m_text]
            label_info['color'] = ['k'] * len(label_info.index)

    return label_info


def set_axislimits(cytokine):
    """Set axis limits

    Parameters
    ----------
    cytokine : str
        Name of cytokine

    Returns
    -------

    """
    if cytokine == 'IFNG':
        xlim = [-15, 35]
        ylim = [0, 50]
    elif cytokine == 'IL13':
        xlim = [-15, 30]
        ylim = [0, 12.5]
    elif cytokine == 'IL17A':
        xlim = [-15, 30]
        ylim = [0, 42]
    else:
        xlim = None
        ylim = None
    return xlim, ylim


def point_plot(df, col_log2fc, col_pval, mask, label, colour, axs, order=1, alpha=1.0):
    """Plot scatter plot of genes meeting the conditions by a mask

    Parameters
    ----------
    df : pandas.Dataframe
    col_log2fc : str
    col_pval : str
    mask : numpy.array
    label : str
    colour : str
    axs : matplotlib.pyplot
    order : int
    alpha : float

    Returns
    -------

    """
    # dotsize was 10 before
    xs, ys = df[col_log2fc][mask], -np.log10(df[col_pval][mask])
    axs.scatter(xs, ys, c=colour, edgecolor='face', label=label, alpha=alpha, s=dotsize, zorder=order, linewidth=0.0)


def volcano_plot(df, df_keys, cytokine, label_genes, title, save_folder, adjust=None, log2fc_cut=1.0, threshold=0.05):
    """Plot DGE Analysis result in Volcano plot

    Parameters
    ----------
    df : pandas.Dataframe
    df_keys : 'list' ['str']
        containing dataframe keys from inbetweener, up- and down-regulated genes
    cytokine : str
    label_genes : 'dict' or 'list'
    title : str
    save_folder : str
    adjust : bool
    log2fc_cut : float
    threshold : float

    Returns
    -------

    """

    log2fc = df_keys[0]
    pval = df_keys[1]

    if adjust is None:
        adjust = True

    # Get axis limits
    xlim, ylim = set_axislimits(cytokine=cytokine)

    # Get positions and colors of Driver and Responder genes
    label_info = get_labels(df=df, annotations=label_genes, log2fc=log2fc, pval=pval,
                            log2fc_cut=log2fc_cut, threshold=threshold, xlim=xlim, ylim=ylim)

    # create masks for all scenarios
    m_sig_log2fc, m_sig_log2fc_pval, m_sig_pval, m_inbetween, index_hkg = get_updowninbetween_masks(
        df, pval, log2fc, pval_cut=threshold, log2fc_cut=log2fc_cut)

    # Create Plot
    fig, ax = plt.subplots(nrows=1, ncols=1, facecolor='w', edgecolor='k', figsize=figure_size)
    ax.set_axisbelow(True)
    # Plot points
    point_plot(df=df, col_log2fc=log2fc, col_pval=pval, mask=m_sig_log2fc_pval,
               colour='darkred', label=r"p-value $\leq$ 0.05 and |log$_2$FC| $\geq$ 1", axs=ax, order=3)
    point_plot(df=df, col_log2fc=log2fc, col_pval=pval, mask=m_sig_log2fc,
               colour='darkblue', label=r"p-value $>$ 0.05 and |log$_2$FC| $>$ 1", axs=ax, order=2, alpha=0.5)
    point_plot(df=df, col_log2fc=log2fc, col_pval=pval, mask=m_sig_pval,
               colour='darkorange', label=r"p-value $<$ 0.05 and |log$_2$FC| $<$ 1", axs=ax, order=2, alpha=0.5)
    point_plot(df=df, col_log2fc=log2fc, col_pval=pval, mask=m_inbetween,
               colour='lightgrey', label=r"p-value $>$ 0.05 and |log$_2$FC| $<$ 1", axs=ax, order=1, alpha=0.5)

    # draw legends
    ax.legend(loc='best', fontsize=legend_fontsize)

    # Add lines to plot
    ax.axhline(y=-np.log10(threshold), linestyle='--', color='k', linewidth=1.5, zorder=4)
    ax.axvline(x=log2fc_cut, linestyle='--', color='k', linewidth=1.5, zorder=4)
    ax.axvline(x=-log2fc_cut, linestyle='--', color='k', linewidth=1.5, zorder=4)

    # Axes properties
    ax.set_xlabel(r'log$_2$(FC)', fontsize=xy_fontsize)
    ax.set_ylabel(r'-log$_{10}$(p-value)', fontsize=xy_fontsize)
    # # sub region of the original image
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(labelsize=xy_ticks)
    # remove upper and right border
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["right"].set_visible(False)

    # # set y-axis to the right
    # ax.yaxis.set_ticks_position("right")

    # Add annotations: Driver and Responder genes
    texts = []
    for ind_label, (x, y) in enumerate(zip(label_info['x'], label_info['y'])):
        texts.append(ax.text(x, y, label_info['gene_symbol'].values[ind_label], size=text_fontsize,
                             color=label_info['color'].values[ind_label], zorder=5))
    if adjust:
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1, zorder=5), ax=ax, lim=500,
                    expand_text=(2, 2))

    # Add FDR cut-off to Volcano plot
    if xlim is None:
        ax.text(np.amax(df[log2fc]) - 3, -np.log10(threshold) + 0.2, '5% FDR',
                size=text_fontsize, color='k', zorder=5)
    else:
        ax.text(xlim[1] - 5, -np.log10(threshold) + 0.2, '5% FDR', size=text_fontsize, color='k', zorder=5)

    # remove legend box
    ax.legend(frameon=False)

    fig.savefig(os.path.join(save_folder, "".join([title, file_format])))
    plt.close()


def plotly_interactive_volcano(df, df_keys, save_folder, key, x_lab, y_lab, log2fc_cut=1.2, pval_cut=0.05):
    """Plot interactive Volcano plot using plotly

    Parameters
    ----------
    df : pandas.Dataframe
        DGE Analysis results
    df_keys : 'list' ['str']
        containing data frame keys from inbetween, up- and down-regulated genes
    save_folder : str
        path to save folder
    key : str
        control vs reference
    x_lab : str
        x-axis label
    y_lab : str
        y-axis label
    log2fc_cut : float
    pval_cut : float

    Returns
    -------

    """
    pval = df_keys[1]
    log2fc = df_keys[0]

    # create masks for all scenarios
    m_sig_log2fc, m_sig_log2fc_pval, m_sig_pval, m_inbetween, index_hkg = get_updowninbetween_masks(
        df, pval, log2fc, log2fc_cut=log2fc_cut, pval_cut=pval_cut)

    # plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df[df_keys[0]][m_inbetween], y=-np.log10(df[df_keys[1]][m_inbetween]),
                             mode='markers', marker=dict(size=dotsize, color='grey'), opacity=0.5,
                             text=df[df_keys[2]][m_inbetween],  # hover text goes here
                             name='inbetween'))

    fig.add_trace(go.Scatter(x=df[df_keys[0]][m_sig_log2fc_pval], y=-np.log10(df[df_keys[1]][m_sig_log2fc_pval]),
                             mode='markers', marker=dict(size=dotsize, color='darkred'),
                             text=df[df_keys[2]][m_sig_log2fc_pval], name='Significant'))
    fig.add_trace(go.Scatter(x=df[df_keys[0]][m_sig_pval], y=-np.log10(df[df_keys[1]][m_sig_pval]),
                             mode='markers', marker=dict(size=dotsize, color='darkorange'), opacity=0.5,
                             text=df[df_keys[2]][m_sig_pval], name=r'p-value < 0.05 and |log$_2$FC| < 1'))
    fig.add_trace(go.Scatter(x=df[df_keys[0]][m_sig_log2fc], y=-np.log10(df[df_keys[1]][m_sig_log2fc]),
                             mode='markers', marker=dict(size=dotsize, color='darkblue'), opacity=0.5,
                             text=df[df_keys[2]][m_sig_log2fc], name=r"p-value > 0.05 and |log$_2$FC| > 1"))
    fig.add_trace(go.Scatter(x=df[df_keys[0]].iloc[index_hkg], y=-np.log10(df[df_keys[1]].iloc[index_hkg]),
                             marker=dict(size=dotsize, color='darkgreen'), mode='markers',
                             text=df[df_keys[2]].iloc[index_hkg], name='Housekeeping genes'))

    fig.update_layout(title="DEx " + key, xaxis_title=x_lab, yaxis_title=y_lab,
                      font=dict(family="Courier New, monospace", size=18, color="#7f7f7f"))

    # fig.write_image(os.path.join(save_folder, key + "_non_interactive.png"))
    fig.write_html(os.path.join(save_folder, key + "_interactive.html"))


def get_updowninbetween_masks(df, pval, log2fc, log2fc_cut=1.0, pval_cut=0.05):
    """Get masks for Volcano plot marking the significant and inbetween genes

    Parameters
    ----------
    df : pandas.Dataframe
    pval : str
    log2fc : str
    log2fc_cut : float
    pval_cut : float

    Returns
    -------

    """

    # create masks for scenarios:
    # I): p-value > pval_cut and |log2FC| > log2fc_cut (grey)
    m_sig_log2fc = (df[pval] > pval_cut) & (abs(df[log2fc]) > log2fc_cut)
    # II): p-value <= pval_cut and |log2FC| >= log2fc_cut (darkred)
    m_sig_log2fc_pval = (df[pval] <= pval_cut) & (abs(df[log2fc]) >= log2fc_cut)
    # III): p-value < pval_cut and |log2FC| < log2fc_cut (orange)
    m_sig_pval = (df[pval] < pval_cut) & (abs(df[log2fc]) < log2fc_cut)
    # IV): p-value > pval_cut and |log2FC| < log2fc_cut (blue)
    m_inbetween = (df[pval] > pval_cut) & (abs(df[log2fc]) < log2fc_cut)

    # Check if you read out the same number of genes
    length_maskgenes = \
        len(df[pval][m_sig_log2fc_pval]) + len(df[pval][m_sig_log2fc]) + \
        len(df[pval][m_sig_pval]) + len(df[pval][m_inbetween])
    print(len(df[pval]) == length_maskgenes)

    # Index House keeping genes
    index_hkg = np.where(df['gene_symbol'][:, np.newaxis] == np.array(_house_keeping_genes)[np.newaxis, :])[0]

    return m_sig_log2fc, m_sig_log2fc_pval, m_sig_pval, m_inbetween, index_hkg


def _write_dataframe(adata, df, output_folder, method):
    """Read and save counts of DGE Analysis groups in a data frame and to a .xlsx / .csv file

    Parameters
    ----------
    adata : annData
    df : pandas.Dataframe
    output_folder : str
    method : str

    Returns
    -------

    """
    # Replace . with - in adata because some genes are written with - in adata but with . in df
    # res = list(filter(lambda x: '-' in x, adata.var_names))
    only_dge_genes = set(df['gene_symbol']) - set(adata.var_names)
    for ind_symbol, gene_name in enumerate(df['gene_symbol']):
        if "." in gene_name and gene_name in only_dge_genes:
            df['gene_symbol'][ind_symbol] = gene_name.replace('.', '-')

    # add counts of each group (condition 1 and condition 2) to df
    # varindex_genes = np.where(adata.var.index[np.newaxis, :] == np.array(df['gene_symbol'].values)[:, np.newaxis])[1]
    # cyto_adata = adata.copy()[:, varindex_genes]
    # m_cyto_counts = adata.obs[observable].values == condition
    # df['condition_1__counts'] = cyto_adata.layers['counts'][~m_cyto_counts].sum(axis=0)
    # df['condition_2__counts'] = cyto_adata.layers['counts'][m_cyto_counts].sum(axis=0)

    # add number of spots a gene is found in per group
    # df['condition_1__number_spots'] = np.count_nonzero(np.copy(cyto_adata.layers['counts'])[~m_cyto_counts], axis=0)
    # df['condition_2__number_spots'] = np.count_nonzero(np.copy(cyto_adata.layers['counts'])[m_cyto_counts], axis=0)

    # Save as .xlsx and .csv file
    df.to_excel(os.path.join(
        output_folder, "_".join([method, "DGE_analysis_results.xlsx"])), sheet_name=method, index=False)
    df.to_csv(os.path.join(output_folder, "_".join([method, "DGE_analysis_results.csv"])), index=False)

    return df


def main(dataset_type, save_folder, df_keys, dge_results_folder):
    """

    Parameters
    ----------
    dataset_type : str
    save_folder : str
    df_keys : list
    dge_results_folder : str

    Returns
    -------

    """
    # Determine name of cluster observable
    if dataset == 'SC':
        cluster_label = 'cluster_labels'
    else:
        cluster_label = 'tissue_type'

    print("# ------ Load data ------ #")
    skin_layers, skin_genes_dict = gene_lists.goldenstandard_genes_epidermis()

    # # 1. Load adata
    adata = load_adata(type_dataset=dataset_type, cluster_label=cluster_label)
    # # 2. Assign labels to spots using golden standard epidermis genes
    for skin_layer in skin_layers:
        adata = add_observables.add_columns_genes(adata, genes=skin_genes_dict[skin_layer], label=skin_layer)

    # read out number of spots for each manual annotation
    anndata_information.minmax_spots_assigned_label(adata=adata, output_folder=save_folder)

    # get files in directory
    all_csv_files = [file
                     for path, subdir, files in os.walk(dge_results_folder)
                     for file in glob.glob(os.path.join(path, '*.csv'))]

    pattern = "".join(['*', '*all_genes.csv'])
    for file in all_csv_files:
        if fnmatch(file, pattern):
            # Read out only those genes which are specific for a comparison
            allgenes_df = pd.read_csv(file, error_bad_lines=False)
            # Remove column Unnamed: 0
            allgenes_df = allgenes_df.drop(['Unnamed: 0'], axis=1)

            # Check if column names and row names are unique
            print("Unique Genes:", allgenes_df['gene_symbol'].is_unique)
            print("Unique Columns:", allgenes_df.columns.is_unique)
            # remove duplicated rows
            allgenes_df = allgenes_df.loc[~allgenes_df['gene_symbol'].duplicated(), :]

            # Name of used design function
            design = file.split(os.sep)[-4]
            # Name of used DGE Analysis method
            method = file.split(os.sep)[-1].split("_")[-4]

            # Create output folder
            output_folder = os.path.join(save_folder, design, file.split(os.sep)[-3], file.split(os.sep)[-2])
            os.makedirs(output_folder, exist_ok=True)

            allgenes_df = _write_dataframe(adata=adata, df=allgenes_df, output_folder=output_folder, method=method)

            print("# ------ Volcano plot ------ #")
            # 3. Volcano plot interactive plot
            plotly_interactive_volcano(df=allgenes_df, df_keys=df_keys, save_folder=output_folder,
                                       key="_".join([method, "Volcano_plot_zoom"]),
                                       x_lab=r'log$_2$(FC)', y_lab=r'-log$_{10}$(pvalue)',
                                       log2fc_cut=1, pval_cut=0.05)

            volcano_plot(df=allgenes_df, df_keys=df_keys, adjust=True, cytokine="",
                         label_genes=skin_genes_dict.values(),
                         title="_".join([method, "Volcano_plot_zoom"]),
                         save_folder=output_folder, log2fc_cut=1.0, threshold=0.05)


if __name__ == '__main__':
    today = date.today()
    dataset = 'ST'
    columns = ['log2fc', 'pval', 'gene_symbol']

    input_folder = os.path.join("..", 'input', 'DGE_analysis_output', '2021-06-01')

    # create output path
    output_path = os.path.join("..", "output", "plots_volcano", str(today))
    os.makedirs(output_path, exist_ok=True)

    main(dataset_type=dataset, save_folder=output_path, df_keys=columns, dge_results_folder=input_folder)
