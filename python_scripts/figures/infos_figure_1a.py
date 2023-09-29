import scanpy as sc
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

import os
from datetime import date


def main(adata, save_folder):
    # drop all samples from newest cohort
    new_samples = list(adata.obs['sample'].cat.categories[adata.obs['sample'].cat.categories.str.contains('P21093')])
    adata_spots = adata[~adata.obs['sample'].isin(new_samples)].copy()
    # Number of spots without new and P16357_1039: 14866
    print('Total number of transcriptomes: ', adata_spots.shape[0])
    print(pd.crosstab(adata_spots.obs['DISEASE'], adata_spots.obs['biopsy_type']))

    # Pso: 6 NL + 6 L
    # AD: 8 NL + 10 L
    # LP: 10 NL

    # ======================
    # Read out number of patient which have SG spots
    adata = adata[adata.obs['SEBACEOUS GLAND'] == 1].copy()
    # drop unuseable sample
    adata = adata[adata.obs['sample'] != 'P16357_1039'].copy()  # 230

    # Number of Spots
    # No. spots AD: 191; So. spots Pso: 39
    pd.crosstab(adata.obs['DISEASE'], adata.obs['SEBACEOUS GLAND'])
    # No. spots AD L: 121, NL: 70; No. spots Pso L 39
    pd.crosstab(adata.obs['DISEASE'], adata.obs['biopsy_type'])

    # Number of samples per disease
    # 3 samples - 'P16357_1035', 'P16357_1036', 'P21093_21L008959'
    adata[adata.obs['DISEASE'] == 'Pso'].obs['sample'].cat.categories
    # 5 samples - 'P16357_1003', 'P16357_1004', 'P16357_1037', 'P16357_1038', 'P21093_21L008971'
    adata[adata.obs['DISEASE'] == 'AD'].obs['sample'].cat.categories

    # Drop newest samples -  removes only 18 spots ..
    adata = adata[~adata.obs['sample'].isin(['P21093_21L008959', 'P21093_21L008971'])].copy()  # 212
    # Number of Spots
    # No. spots AD: 177; So. spots Pso: 35
    pd.crosstab(adata.obs['DISEASE'], adata.obs['SEBACEOUS GLAND'])
    # No. spots AD L: 107, NL: 70; No. spots Pso L 35
    pd.crosstab(adata.obs['DISEASE'], adata.obs['biopsy_type'])

    # Number of samples per disease
    # 2 samples - 'P16357_1035', 'P16357_1036'
    adata[adata.obs['DISEASE'] == 'Pso'].obs['sample'].cat.categories
    # 4 samples - 'P16357_1003', 'P16357_1004', 'P16357_1037', 'P16357_1038'
    adata[adata.obs['DISEASE'] == 'AD'].obs['sample'].cat.categories


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output",
        "figure_1e", str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Remove LP
    mask = (spatial_adata.obs['DISEASE'] == 'LP') & (spatial_adata.obs['biopsy_type'] == 'LESIONAL')
    spatial_adata = spatial_adata[~mask].copy()

    main(adata=spatial_adata, save_folder=savepath)
