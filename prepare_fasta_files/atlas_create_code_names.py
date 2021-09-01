import pandas as pd
import numpy as np

atlas_info = pd.read_excel("./tcr_atlas/ATLAS.xlsx", index_col=0)
atlas_info["structure_id"] = ""
print(atlas_info)

for i in atlas_info.index.values:
    pdb = atlas_info.loc[i, "template_PDB"]
    if not isinstance(pdb, str):
        atlas_info.loc[i, "structure_id"] = atlas_info.loc[i, "true_PDB"].lower()
        continue
    else:
        mhc_mut = ".".join(atlas_info.loc[i, "MHC_mut"].split(" | "))
        mhc_mut_chain = ".".join(str(atlas_info.loc[i, "MHC_mut_chain"]).split(" | "))
        tcr_mut = ".".join(atlas_info.loc[i, "TCR_mut"].split(" | "))
        tcr_mut_chain = ".".join(str(atlas_info.loc[i, "TCR_mut_chain"]).split(" | "))
        pep_mut = ".".join(atlas_info.loc[i, "PEP_mut"].split(" | "))
        atlas_info.loc[i, "structure_id"] = pdb + "_" + mhc_mut + "_" + mhc_mut_chain + "_" + tcr_mut + "_" + tcr_mut_chain + "_" + pep_mut

print(atlas_info)
atlas_info.to_csv("./tcr_atlas/atlas_with_names.csv")
