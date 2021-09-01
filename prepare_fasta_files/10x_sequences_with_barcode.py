import pandas as pd
import numpy as np
from Bio import SeqIO, SeqUtils

data10x = pd.read_csv("./B10x/10x_refined.csv", index_col = 0)
aa_seqs = pd.read_csv("./B10x/VJ_aa_seq_human.csv", index_col=0)
mhc_aa_seqs = pd.read_csv("./B10x/mhc_aa_seq_human.csv", index_col=0)
binding_info = pd.read_csv("./B10x/peptide_UMI_matrix.csv", index_col=0)

complete_10x = pd.DataFrame(columns = ["alpha", "beta", "peptide", "mhc", "barcode", "binder"])

for x in range(0, np.shape(data10x)[0]):
    print(x)
    data = data10x.iloc[x]
    barcode = data["barcode"]
    V_A_gene = data.V_A + "*01"
    J_A_gene = data.J_A + "*01"
    V_B_gene = data.V_B + "*01"
    J_B_gene = data.J_B + "*01"
    mhc_gene = data["mhc"]
    binding_UMI = binding_info.loc[barcode, data["pep_seq"]]
    ctrl_UMI = binding_info.loc[barcode, "neg_thresh"]
    binder = int((binding_UMI > ctrl_UMI) & (binding_UMI > 10))
    if V_A_gene in aa_seqs.index.values and J_A_gene in aa_seqs.index.values and V_B_gene in aa_seqs.index.values and J_B_gene in aa_seqs.index.values and  mhc_gene in mhc_aa_seqs.index.values:
        V_A = aa_seqs.loc[V_A_gene, "aa_seq"]
        J_A = aa_seqs.loc[J_A_gene, "aa_seq"]
        CDR3_A = data.CDR3A
        i = len(V_A) - 2
        while i > 0:
            CDR_begin = CDR3_A[0:2]
            if V_A[i:i+2] == CDR_begin:
                V_A = V_A[:i]
                break
            else:
                i -= 1
        j = 0
        while j < len(J_A):
            CDR_end = CDR3_A[len(CDR3_A)-2:]
            if J_A[j:j+2] == CDR_end:
                J_A = J_A[j+2:]
                break
            else:
                j += 1
        alpha = V_A+CDR3_A+J_A

        V_B = aa_seqs.loc[V_B_gene, "aa_seq"]
        J_B = aa_seqs.loc[J_B_gene, "aa_seq"]
        CDR3_B = data.CDR3B
        i = len(V_B) - 2
        while i > 0:
            CDR_begin = CDR3_B[0:2]
            if V_B[i:i+2] == CDR_begin:
                V_B = V_B[:i]
                break
            else:
                i -= 1
        j = 0
        while j < len(J_B):
            CDR_end = CDR3_B[len(CDR3_B)-2:]
            if J_B[j:j+2] == CDR_end:
                J_B = J_B[j+2:]
                break
            else:
                j += 1
        beta = V_B+CDR3_B+J_B

        mhc = mhc_aa_seqs.loc[mhc_gene, "aa_seq"]

        peptide = data["pep_seq"]

        chain_names = ["alpha", "beta", "peptide", "mhc"]
        for j, chain in enumerate([alpha, beta, peptide, mhc]):
            chain_name = chain_names[j]
            complete_10x.loc[x, chain_name] = chain

        complete_10x.loc[x, "barcode"] = barcode
        complete_10x.loc[x, "binder"] = binder

complete_10x.to_csv("./B10x/B10x_sequences_with_barcodes_UMI.csv")
