import pandas as pd
import numpy as np

binding_matrix = pd.read_csv("./new_10X/tcrs_with_binding_info.csv", index_col = 0)
sequences = pd.read_csv("./new_10X/10X_all_sequences_from_fasta.csv", index_col = 0)

print(sequences)

for id in sequences.index.values:
    alpha = sequences.loc[id, "alpha"]
    beta = sequences.loc[id, "beta"]
    peptide = sequences.loc[id, "peptide"]

    binding_info = binding_matrix.loc[(binding_matrix.alpha == alpha) & (binding_matrix.beta == beta)]
    if np.shape(binding_info)[0]>1:
        print("multiple hits for this sequence")
        break
    neg_thresh = binding_info.iloc[0]["neg_thresh"]
    pep_col = binding_info.columns[binding_info.columns.to_series().str.contains(peptide)][0]
    umi_count = binding_info.iloc[0][pep_col]
    # print(umi_count, neg_thresh)
    sequences.loc[id, "umi_count"] = umi_count
    sequences.loc[id, "binding"] = int((umi_count > neg_thresh) & (umi_count > 10))

classification_results = pd.read_csv("./post-bug/B10x/10X_classification_results_atchley-dist-atr.csv", index_col = 0)
print(classification_results)

classification_results.index = [x.split("_")[0] for x in classification_results.index.values]

merged = pd.merge(sequences, classification_results, left_index=True, right_index=True)
print(merged)

merged = merged.drop("real_answer", axis = 1)
# merged.to_csv("./new_10X/10X_results_crosschecked_with_UMIs_trained_on_all_templ.csv")

sequences.to_csv("./new_10X/10X_all_sequences_from_fasta_with_binding.csv")
