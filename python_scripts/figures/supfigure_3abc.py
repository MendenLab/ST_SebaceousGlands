import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np
import random
import pandas as pd
import os
from datetime import date
import pingouin as pg
from adjustText import adjust_text


def plot_scatter(df, save_folder, pattern_name_1, pattern_name_2):
    df, df_tmp_shared, df_tmp_unique_1, df_tmp_unique_2 = get_infos(
        df=df, pattern_name_1=pattern_name_1, pattern_name_2=pattern_name_2)

    corr_val = pg.corr(x=df['log_p.adjust_x'], y=df['log_p.adjust_y'])

    texts = []
    for ind in range(df_tmp_unique_2[:4].shape[0]):
        texts.append(plt.text(df_tmp_unique_2.iloc[ind, :]['log_p.adjust_x'],
                              df_tmp_unique_2.iloc[ind, :]['log_p.adjust_y'],
                              df_tmp_unique_2.iloc[ind, :]['Description_y_adjusted']))

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(data=df, x='log_p.adjust_x', y='log_p.adjust_y', ax=ax, c=['black'], edgecolor='black', s=8, zorder=10)
    # Shared points
    sns.scatterplot(data=df[df['Description_x'].isin(list(df_tmp_shared.loc[:, 'Description_x'].values))],
                    x='log_p.adjust_x', y='log_p.adjust_y', ax=ax, c=['red'], edgecolor='red', s=20, zorder=10)
    # Unique pattern 1
    sns.scatterplot(data=df[df['Description_x'].isin(list(df_tmp_unique_1.loc[:, 'Description_x'].values))],
                    x='log_p.adjust_x', y='log_p.adjust_y', ax=ax, c=['red'], edgecolor='red', s=20, zorder=10)
    # Unique pattern 2
    sns.scatterplot(data=df[df['Description_y'].isin(list(df_tmp_unique_2.loc[:, 'Description_y'].values))],
                    x='log_p.adjust_x', y='log_p.adjust_y', ax=ax, c=['red'], edgecolor='red', s=20, zorder=10)
    ax.axvline(ymin=0, ymax=30, x=-np.log10(0.05), ls='dashed', c='grey', lw=1)
    ax.axhline(xmin=0, xmax=30, y=-np.log10(0.05), ls='dashed', c='grey', lw=1)
    ax.text(df['log_p.adjust_x'].max() / 4, df['log_p.adjust_y'].max(),
            'r: {:.2e}; pval: {:.2e}'.format(corr_val['r'][0], corr_val['p-val'][0]), fontsize=14)
    # Highlight shared pathways
    for i in range(df_tmp_shared.shape[0]):
        ax.annotate(df_tmp_shared.iloc[i, :]['Description_x_adjusted'],
                    xy=(df_tmp_shared.iloc[i, :]['log_p.adjust_x'], df_tmp_shared.iloc[i, :]['log_p.adjust_y']),
                    xytext=(df_tmp_shared.iloc[i, :]['log_p.adjust_x'] + random.uniform(0.2, 0.8),
                            df_tmp_shared.iloc[i, :]['log_p.adjust_y'] + random.uniform(0.2, 0.4)),
                    c='red', fontsize=12, zorder=5)
    for i in range(df_tmp_unique_1[:4].shape[0]):
        ax.annotate(df_tmp_unique_1.iloc[i, :]['Description_x_adjusted'],
                    xy=(df_tmp_unique_1.iloc[i, :]['log_p.adjust_x'], df_tmp_unique_1.iloc[i, :]['log_p.adjust_y']),
                    xytext=(df_tmp_unique_1.iloc[i, :]['log_p.adjust_x'] + 10 - 7*i,
                            df_tmp_unique_1.iloc[i, :]['log_p.adjust_y'] - 4),
                    fontsize=12, arrowprops=dict(arrowstyle='-|>', lw=1.5, ls='-'), zorder=1)
    for i in range(df_tmp_unique_2[:4].shape[0]):
        ax.annotate(df_tmp_unique_2.iloc[i, :]['Description_y_adjusted'],
                    xy=(df_tmp_unique_2.iloc[i, :]['log_p.adjust_x'], df_tmp_unique_2.iloc[i, :]['log_p.adjust_y']),
                    xytext=(df_tmp_unique_2.iloc[i, :]['log_p.adjust_x'],
                            df_tmp_unique_2.iloc[i, :]['log_p.adjust_y'] + 8 - 2*i),
                    fontsize=12, arrowprops=dict(arrowstyle='-|>', lw=1.5, ls='-'), zorder=1)
        # adjust_text(texts, only_move={'points': 'y', 'texts': 'y'}, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
    # Highlight unique pathways
    ax.set_xlabel(''.join([r'-log$_{10}$(p.adj) ', 'of Pattern {} pathways'.format(pattern_name_1)]), fontsize=18)
    ax.set_ylabel(''.join([r'-log$_{10}$(p.adj) ', 'of Pattern {} pathways'.format(pattern_name_2)]), fontsize=18)
    ax.tick_params(axis='both', labelsize=16)
    sns.despine()
    if pattern_name_2 == '7':
        ax.set_xlim([-4, df['log_p.adjust_x'].max() + 1])
        ax.set_ylim([-1, df['log_p.adjust_y'].max() + 1])
    elif pattern_name_1 == '7':
        ax.set_xlim([-1, df['log_p.adjust_x'].max() + 1])
        ax.set_ylim([-4, df['log_p.adjust_y'].max() + 1])
    else:
        ax.set_xlim([-4, df['log_p.adjust_x'].max() + 1])
        ax.set_ylim([-4, df['log_p.adjust_y'].max() + 1])
    # ax.set_xlim([-5, df['log_p.adjust_x'].max() + 1])
    # ax.set_ylim([-5, df['log_p.adjust_y'].max() + 1])
    plt.savefig(os.path.join(
        save_folder, 'Scatterplot__Patter_{}_vs_{}.pdf'.format(pattern_name_1, pattern_name_2)), bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(data=df, x='log_p.adjust_x', y='log_p.adjust_y', ax=ax, c=['black'], edgecolor='black', s=15,
                    zorder=10)
    # Shared points
    sns.scatterplot(data=df[df['Description_x'].isin(list(df_tmp_shared.loc[:, 'Description_x'].values))],
                    x='log_p.adjust_x', y='log_p.adjust_y', ax=ax, c=['red'], edgecolor='red', s=30, zorder=10)
    # Unique pattern 1
    sns.scatterplot(data=df[df['Description_x'].isin(list(df_tmp_unique_1.loc[:, 'Description_x'].values))],
                    x='log_p.adjust_x', y='log_p.adjust_y', ax=ax, c=['red'], edgecolor='red', s=30, zorder=10)
    # Unique pattern 2
    sns.scatterplot(data=df[df['Description_y'].isin(list(df_tmp_unique_2.loc[:, 'Description_y'].values))],
                    x='log_p.adjust_x', y='log_p.adjust_y', ax=ax, c=['red'], edgecolor='red', s=30, zorder=10)
    ax.axvline(ymin=0, ymax=30, x=-np.log10(0.05), ls='dashed', c='grey', lw=1)
    ax.axhline(xmin=0, xmax=30, y=-np.log10(0.05), ls='dashed', c='grey', lw=1)
    ax.text(df['log_p.adjust_x'].max() / 4, np.ceil(df['log_p.adjust_y'].max()),
            'r: {:.2e}; pval: {:.2e}'.format(corr_val['r'][0], corr_val['p-val'][0]), fontsize=14)
    ax.set_xlabel(''.join([r'-log$_{10}$(p.adj) ', 'of Pattern {} pathways'.format(pattern_name_1)]), fontsize=18)
    ax.set_ylabel(''.join([r'-log$_{10}$(p.adj) ', 'of Pattern {} pathways'.format(pattern_name_2)]), fontsize=18)
    ax.tick_params(axis='both', labelsize=16)
    sns.despine()
    if pattern_name_2 == '7':
        ax.set_xlim([-4, df['log_p.adjust_x'].max() + 1])
        ax.set_ylim([-1, df['log_p.adjust_y'].max() + 1])
    elif pattern_name_1 == '7':
        ax.set_xlim([-1, df['log_p.adjust_x'].max() + 1])
        ax.set_ylim([-4, df['log_p.adjust_y'].max() + 1])
    else:
        ax.set_xlim([-4, df['log_p.adjust_x'].max() + 1])
        ax.set_ylim([-4, df['log_p.adjust_y'].max() + 1])
    plt.savefig(os.path.join(
        save_folder, 'Scatterplot_unannotated__Patter_{}_vs_{}.pdf'.format(
            pattern_name_1, pattern_name_2)), bbox_inches='tight')
    plt.close()

    print('Unique pathways of Pattern {}: {}\n'.format(pattern_name_1, df_tmp_unique_1.shape[0]))
    print(df_tmp_unique_1['Description_x'])
    print('\n\nUnique pathways of Pattern {}: {}\n'.format(pattern_name_2, df_tmp_unique_2.shape[0]))
    print(df_tmp_unique_2['Description_y'])
    print('\n=======================\n')

    # Save information to excel file
    df.to_excel(os.path.join(save_folder, 'Pathways_{}_vs_{}.xlsx'.format(pattern_name_1, pattern_name_2)))


def get_infos(df, pattern_name_1, pattern_name_2):
    # Masks for text
    mask = (df.loc[:, 'log_p.adjust_x'] > -np.log10(0.05)) & (df.loc[:, 'log_p.adjust_y'] > -np.log10(0.05))
    df_tmp_shared = df.loc[mask, :]

    # Insert after every third word a newline
    df_tmp_shared = insert_linebreak(df=df_tmp_shared, colname='Description_x')

    mask_unique_12_pattern_1 = (df.loc[:, 'log_p.adjust_x'] > -np.log10(0.05)) & (
            df.loc[:, 'log_p.adjust_y'] == 0)
    df_tmp_unique_1 = df.loc[mask_unique_12_pattern_1, :]
    # Insert after every third word a newline
    df_tmp_unique_1 = insert_linebreak(df=df_tmp_unique_1, colname='Description_x')

    mask_unique_12_pattern_2 = (df.loc[:, 'log_p.adjust_x'] == 0) & (
            df.loc[:, 'log_p.adjust_y'] > -np.log10(0.05))
    df_tmp_unique_2 = df.loc[mask_unique_12_pattern_2, :]
    # Insert after every third word a newline
    df_tmp_unique_2 = insert_linebreak(df=df_tmp_unique_2, colname='Description_y')

    df['Pathway_type'] = 'unspecified'
    df.loc[mask, 'Pathway_type'] = 'shared_significant'
    df.loc[mask_unique_12_pattern_1, 'Pathway_type'] = 'unique_pattern_{}'.format(pattern_name_1)
    df.loc[mask_unique_12_pattern_2, 'Pathway_type'] = 'unique_pattern_{}'.format(pattern_name_2)

    return df, df_tmp_shared, df_tmp_unique_1, df_tmp_unique_2


def insert_linebreak(df, colname='Description_x'):
    colname_new = "{}_adjusted".format(colname)
    df[colname_new] = df[colname]
    # Insert after every third word a newline
    for i in df.index:
        df[colname_new][i] = ''.join([df[colname][i].split()[n] + ' ' if n % 3 != 2 else
                                      df[colname][i].split()[n] + '\n' for n in
                                      range(len(df[colname][i].split()))])

    return df


def main(input_folder, save_folder):

    filename = 'Patterns_11-V19T12-012-V2_Reactome_ORA_.xlsx'

    df_patter_1 = pd.read_excel(os.path.join(input_folder, filename), sheet_name='1')
    df_patter_1['log_p.adjust'] = -np.log10(df_patter_1['p.adjust'])
    df_patter_7 = pd.read_excel(os.path.join(input_folder, filename), sheet_name='7')
    df_patter_7['log_p.adjust'] = -np.log10(df_patter_7['p.adjust'])
    df_patter_9 = pd.read_excel(os.path.join(input_folder, filename), sheet_name='9')
    df_patter_9['log_p.adjust'] = -np.log10(df_patter_9['p.adjust'])

    # Scatterplot of padj val of pathways of Pattern 1 vs padj val of pathways of Pattern 9
    # Note: Description not unique but ID is
    df_19 = df_patter_1[['ID', 'Description', 'log_p.adjust']].merge(
        df_patter_9[['ID', 'Description', 'log_p.adjust']], on=['ID'], how='outer')
    df_19['log_p.adjust_x'] = df_19['log_p.adjust_x'].fillna(0)
    df_19['log_p.adjust_y'] = df_19['log_p.adjust_y'].fillna(0)

    plot_scatter(df=df_19, save_folder=save_folder, pattern_name_1='1', pattern_name_2='9')

    df_17 = df_patter_1[['ID', 'Description', 'log_p.adjust']].merge(
        df_patter_7[['ID', 'Description', 'log_p.adjust']], on=['ID'], how='outer')
    df_17['log_p.adjust_x'] = df_17['log_p.adjust_x'].fillna(0)
    df_17['log_p.adjust_y'] = df_17['log_p.adjust_y'].fillna(0)
    plot_scatter(df=df_17, save_folder=save_folder, pattern_name_1='1', pattern_name_2='7')

    df_79 = df_patter_7[['ID', 'Description', 'log_p.adjust']].merge(
        df_patter_9[['ID', 'Description', 'log_p.adjust']], on=['ID'], how='outer')
    df_79['log_p.adjust_x'] = df_79['log_p.adjust_x'].fillna(0)
    df_79['log_p.adjust_y'] = df_79['log_p.adjust_y'].fillna(0)
    plot_scatter(df=df_79, save_folder=save_folder, pattern_name_1='7', pattern_name_2='9')


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    input_dir = '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output'
    savepath = os.path.join(input_dir, "figure_S3abc", str(today))
    os.makedirs(savepath, exist_ok=True)

    input_path = os.path.join(input_dir, 'spatialDE/2023-09-18_paper_figures/Pathway_enrichment_analysis',
                              'Figure_2lmn/patient/ORA_REACTOME/Patterns/2023-09-18/11-V19T12-012-V2')

    main(input_folder=input_path, save_folder=savepath)
