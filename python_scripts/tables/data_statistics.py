import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import os
from datetime import date


def save_all_data_stats(df, save_folder):
    df.groupby(['DISEASE', 'patient', 'spot_type', 'biopsy_type']).size().to_excel(
        os.path.join(save_folder, 'Overview_spot_annotations.xlsx'))
    df_gb = df.groupby(['DISEASE', 'spot_type', 'biopsy_type']).size().unstack(level=2)

    df_gb.plot(kind='bar')
    plt.xlabel('disease', fontsize=18)
    plt.xticks(fontsize=16, rotation=0)
    plt.ylabel('counts', fontsize=18)
    plt.yticks(fontsize=16, rotation=0)
    plt.title('Spot types', fontsize=20)
    sns.despine()
    plt.legend(frameon=False)
    plt.tight_layout()

    plt.savefig(os.path.join(save_folder, 'spottypes__counts.pdf'))
    plt.close('all')


def save_sg_stats(df, save_folder):
    df_sg = df[df['spot_type'] == 'SEBACEOUS GLAND']
    df_sg.groupby(['DISEASE', 'patient', 'biopsy_type']).size().to_excel(
        os.path.join(save_folder, 'Overview_SG.xlsx'))

    df_gb = df_sg.groupby(['DISEASE', 'biopsy_type']).size().unstack(level=1)
    df_gb.plot(kind='bar')
    plt.xlabel('disease', fontsize=18)
    plt.xticks(fontsize=16, rotation=0)
    plt.ylabel('counts', fontsize=18)
    plt.yticks(fontsize=16, rotation=0)
    plt.title('Sebaceous glands', fontsize=20)
    sns.despine()
    plt.legend(frameon=False)
    plt.tight_layout()

    plt.savefig(os.path.join(save_folder, 'SG__counts.pdf'))
    plt.close('all')


def main(adata, save_folder):
    # keep used sample in analysis
    # adata = adata[adata.obs['sample'].isin(
    #     ['P16357_1035', 'P16357_1036', 'P16357_1003', 'P16357_1004', 'P16357_1037', 'P16357_1038'])].copy()  # 212

    df = adata.obs[['DISEASE', 'spot_type', 'biopsy_type', 'patient']]

    save_all_data_stats(df, save_folder)
    # SEBACEOUS GLAND
    save_sg_stats(df, save_folder)


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "data_statistics", str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Remove LP
    mask = (spatial_adata.obs['DISEASE'] == 'LP') & (spatial_adata.obs['biopsy_type'] == 'LESIONAL')
    spatial_adata = spatial_adata[~mask].copy()

    # drop all samples from newest cohort
    new_samples = list(spatial_adata.obs['sample'].cat.categories[
                           spatial_adata.obs['sample'].cat.categories.str.contains('P21093')])
    adata = spatial_adata[~spatial_adata.obs['sample'].isin(new_samples)].copy()

    main(adata=spatial_adata, save_folder=savepath)
