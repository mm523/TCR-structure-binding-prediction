import pandas as pd
import numpy as np
from Bio import SeqIO, SeqUtils

Dash_info = pd.read_csv("./Dash/Dash2017-fromVDJDb.tsv", sep = "\t", index_col = 0)
print(Dash_info)

beta_chains = Dash_info[Dash_info.Gene == "TRB"]
beta_chains.columns = [x + "_b" for x in beta_chains.columns]
alpha_chains = Dash_info[Dash_info.Gene == "TRA"]
alpha_chains.columns = [x + "_a" for x in alpha_chains.columns]

Dash_info2 = pd.merge(alpha_chains, beta_chains, on = "complex.id")
print(Dash_info2)

Dash_info2 = Dash_info2.drop_duplicates().reset_index(drop = True)
Dash_info2.to_csv("./Dash/Dash2017_worked.csv")

VJ_human = pd.read_csv("./B10x/VJ_aa_seq_human.csv", index_col=0)
VJ_mouse = pd.read_csv("./Dash/VJ_aa_seq_mouse.csv", index_col=0)
mhc_human = pd.read_csv("./B10x/mhc_aa_seq_human.csv", index_col=0)
mhc_mouse = pd.read_csv("./Dash/mhc_aa_seq_mouse.csv", index_col=0)
PDB_seqs = pd.read_csv("./post_bug/PDB_set/PDB_all_sequences_from_fasta.csv")
B10x_seqs = pd.read_csv("./post_bug/B10x/10X_all_sequences_from_fasta.csv")
existing_seqs = pd.concat([PDB_seqs, B10x_seqs])
complete_Dash = pd.DataFrame(columns = ["alpha", "beta", "peptide", "mhc"])

# for i in Dash_info2.index.values:
#     id = "Dash-" + "%04d" %i
#     print(i, id)
#     data = Dash_info2.loc[i]
#     print(data)
#     V_A_gene = data.V_a
#     J_A_gene = data.J_a
#     V_B_gene = data.V_b
#     J_B_gene = data.J_b
#     mhc_gene = data["MHC A_a"].replace("HLA-", "").replace("*", "").replace(":", "")
#     if data.Species_a == "HomoSapiens":
#         print("homo")
#         VJ_aa_seqs = VJ_human
#         mhc_aa_seqs = mhc_human
#     else:
#         VJ_aa_seqs = VJ_mouse
#         mhc_aa_seqs = mhc_mouse
#     print(V_A_gene, V_B_gene, J_A_gene, J_B_gene, mhc_gene)
#     if V_A_gene in VJ_aa_seqs.index.values and J_A_gene in VJ_aa_seqs.index.values and V_B_gene in VJ_aa_seqs.index.values and J_B_gene in VJ_aa_seqs.index.values and  mhc_gene in mhc_aa_seqs.index.values:
#         V_A = VJ_aa_seqs.loc[V_A_gene, "aa_seq"]
#         J_A = VJ_aa_seqs.loc[J_A_gene, "aa_seq"]
#         CDR3_A = data.CDR3_a
#         i = len(V_A) - 2
#         while i > 0:
#             CDR_begin = CDR3_A[0:2]
#             if V_A[i:i+2] == CDR_begin:
#                 V_A = V_A[:i]
#                 break
#             else:
#                 i -= 1
#         j = 0
#         while j < len(J_A):
#             CDR_end = CDR3_A[len(CDR3_A)-2:]
#             if J_A[j:j+2] == CDR_end:
#                 J_A = J_A[j+2:]
#                 break
#             else:
#                 j += 1
#         alpha = V_A+CDR3_A+J_A
#
#         V_B = VJ_aa_seqs.loc[V_B_gene, "aa_seq"]
#         J_B = VJ_aa_seqs.loc[J_B_gene, "aa_seq"]
#         CDR3_B = data.CDR3_b
#         i = len(V_B) - 2
#         while i > 0:
#             CDR_begin = CDR3_B[0:2]
#             if V_B[i:i+2] == CDR_begin:
#                 V_B = V_B[:i]
#                 break
#             else:
#                 i -= 1
#         j = 0
#         while j < len(J_B):
#             CDR_end = CDR3_B[len(CDR3_B)-2:]
#             if J_B[j:j+2] == CDR_end:
#                 J_B = J_B[j+2:]
#                 break
#             else:
#                 j += 1
#         beta = V_B+CDR3_B+J_B
#
#         mhc = mhc_aa_seqs.loc[mhc_gene, "aa_seq"]
#
#         peptide = data["Epitope_a"]
#
#         if peptide in list(set(existing_seqs.peptide)):
#             print("peptide found, checking for complex...")
#             to_check = existing_seqs.loc[existing_seqs.peptide == peptide]
#             for i in range(0, np.shape(to_check)[0]):
#                 alpha1 = to_check.iloc[i]["alpha"]
#                 beta1 = to_check.iloc[i]["beta"]
#                 if alpha == alpha1 and beta == beta1:
#                     print("this sequence is already in the dataset, discarding...")
#                     break
#
#         print("saving sequence")
#
#         f = open("./post_bug/Dash/March2021/fasta_files/" + id + ".fasta", "w+")
#
#         chain_names = ["alpha", "beta", "peptide", "mhc"]
#         for i, chain in enumerate([alpha, beta, peptide, mhc]):
#             chain_name = chain_names[i]
#             f.write(">" + chain_name + "\n")
#             f.write(chain + "\n")
#             complete_Dash.loc[id, chain_name] = chain
#
#         f.close()
#
# complete_Dash.to_csv("./post_bug/Dash/March2021/Dash_sequences.csv")

Dash_sequences = pd.read_csv("./post_bug/Dash/March2021/Dash_sequences.csv", index_col = 0)

peptides = Dash_sequences[['peptide', 'mhc']].drop_duplicates().reset_index(drop = True)
print(peptides)
peptide_mhc_dict = {}

for pep in peptides.peptide:
    mhc = peptides[peptides.peptide == pep].reset_index().loc[0, "mhc"]
    print(mhc)
    peptide_mhc_dict[pep] = mhc

decoy_sequences = pd.DataFrame(columns = ["alpha", "beta", "peptide", "mhc"])
j = 0

for a, b in zip(Dash_sequences.alpha, Dash_sequences.beta):
    tcr_binding = Dash_sequences[(Dash_sequences.alpha == a) & (Dash_sequences.beta == b)]
    binding_peptides = set(tcr_binding.peptide)
    for pep1 in peptide_mhc_dict.keys():
        if pep1 in binding_peptides:
            continue
        else:
            mhc1 = peptide_mhc_dict[pep1]
            decoy_name = "Dash-decoy" + "{0:0=4d}".format(j)
            f = open("./post_bug/Dash/March2021/fasta_files/" + decoy_name + ".fasta", "w+")
            chain_names = ["alpha", "beta", "peptide", "mhc"]
            for i, chain in enumerate([a, b, pep1, mhc1]):
                chain_name = chain_names[i]
                f.write(">" + chain_name + "\n")
                f.write(chain + "\n")
                decoy_sequences.loc[decoy_name, chain_name] = chain
            f.close()
            j+=1
decoy_sequences.to_csv("./post_bug/Dash/March2021/decoy_sequences.csv")
