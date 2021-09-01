import pandas as pd
import numpy as np

binarized_info = pd.read_csv("./B10x/vdj_v1_hs_aggregated_donor1_binarized_matrix.csv", index_col=0)
print(binarized_info)

# col_subset = [x for x in binarized_info.columns.values if "binder" in x]
# col_names = [x.split("_")[1] for x in col_subset]
# binary = binarized_info[col_subset]
# binary.columns = col_names
# binary = binary.astype(int)
# print(binary)
#
# binary.to_csv("./B10x/peptide_binding_matrix.csv")

col_subset = binarized_info.columns.values[17:67]
neg_ctrls = [x for x in col_subset if "NC" in x]
negatives = binarized_info[neg_ctrls]
print(negatives)
neg_thresh = {}
for bc in negatives.index.values:
    neg_thresh[bc] = 5*max(negatives.loc[bc])
print(neg_thresh)
print(col_subset)
# col_names = [x.split("_")[1] for x in col_subset]
binary = binarized_info[col_subset]
# binary.columns = col_names
binary = binary.astype(int)
print(binary)
thresh = pd.DataFrame.from_dict(neg_thresh, orient="index", columns=["neg_thresh"])
print(thresh)
new_binary = pd.merge(binary, thresh, left_index=True, right_index=True)
print(new_binary)

new_binary.to_csv("./B10x/peptide_UMI_matrix.csv")
