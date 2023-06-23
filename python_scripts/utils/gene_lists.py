from collections import OrderedDict


def goldenstandard_genes_epidermis():
    goldenstandard_genes = OrderedDict()
    # Source LCE: https://www.genenames.org/data/genegroup/#!/group/627
    goldenstandard_genes["uE"] = ["FLG", "HRN", "LCE1A", "LCE1B", "LCE1C", "LCE1D", "LCE1E", "LCE1F",
                                  "LCE2A", "LCE2B", "LCE2C", "LCE2D",
                                  "LCE3A", "LCE3B", "LCE3C", "LCE3D", "LCE3E",
                                  "LCE4A", "LCE5A", "LCE6A"]
    goldenstandard_genes["mE"] = ["FLG", "HRN", "LCE1A", "LCE1B", "LCE1C", "LCE1D", "LCE1E", "LCE1F",
                                  "LCE2A", "LCE2B", "LCE2C", "LCE2D",
                                  "LCE3A", "LCE3B", "LCE3C", "LCE3D", "LCE3E",
                                  "LCE4A", "LCE5A", "LCE6A"]
    # Source S100: https://www.genenames.org/data/genegroup/#!/group/459
    goldenstandard_genes["bE"] = ["KRT10", "DEFB4", "S100A1", "S100A2", "S100A3", "S100A4", "S100A5",
                                  "S100A6", "S100A7", "S100A7A", "S100A7L2", "S100A7P1", "S100A7P2", "S100A8", "S100A9",
                                  "S100A10", "S100A11", "S100A12", "S100A13", "S100A14", "S100A15A", "S100A16",
                                  "S100B", "S100G", "S100P", "S100Z"]

    return list(goldenstandard_genes.keys()), goldenstandard_genes


def goldenstandard_sebcaeous_glands_genes():
    goldenstandard_genes = OrderedDict()
    goldenstandard_genes['AD'] = ['EBP', 'HSD3B1', 'COX8A', 'COX5B']
    # goldenstandard_genes['PsO'] = ['any S100.., any KRT..., any FABP.., any SERPIN.., any IL17.., any CARD..',
    #                                'IL36', 'IL22', 'NOS2', 'CXCL8', 'DEFB4', 'CARD18', 'SERPINF1']
    goldenstandard_genes['Pso'] = ['ACOT4', 'S1PR3', 'SERPINF1']

    # common_genes = ['KRT17', 'KRT79', 'SAA1', 'FABP7', 'SERPINE2', 'ACAT2', 'HSD11B1', 'RARRES1']

    return goldenstandard_genes
