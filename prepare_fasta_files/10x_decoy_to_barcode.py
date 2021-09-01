import pandas as pd
import numpy as np

decoys = pd.read_csv("./B10x/decoy_summary.csv", index_col=0)
sequence_info = pd.read_csv("./B10x/B10x_all_sequences_with_10xID_barcode_binding_UMI.csv", index_col = 0)
binding_UMIs = pd.read_csv("./B10x/peptide_UMI_matrix.csv", index_col = 0)

decoy_bc = pd.DataFrame()

print(decoys)

for i in decoys.index.values:
    tcr_id = decoys.loc[i, "id_TCR"]
    peptide_id = decoys.loc[i, "id_pMHC"]
    decoy_id = decoys.loc[i, "decoy_name"]

    bc = sequence_info[sequence_info["10x_ID"] == tcr_id].iloc[0]["barcode"]
    peptide = sequence_info[sequence_info["10x_ID"] == peptide_id].iloc[0]["peptide"]
    umi = binding_UMIs.loc[bc, peptide]
    if umi > binding_UMIs.loc[bc, "neg_thresh"] and umi > 10:
        binding = 1
    else:
        binding = 0
    decoy_bc.loc[decoy_id, "barcode"] = bc
    decoy_bc.loc[decoy_id, "peptide"] = peptide
    decoy_bc.loc[decoy_id, "binder"] = binding
    decoy_bc.loc[decoy_id, "binding_UMIs"] = umi

decoy_bc["10x_ID"] = decoy_bc.index.values
decoy_bc.to_csv("./B10x/B10x_decoy_sequences_with_10xID_barcode_binding_UMI.csv")
