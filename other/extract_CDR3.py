import pandas as pd
from datetime import datetime
from sys import argv
d = datetime.today().strftime('%Y-%m-%d')

s = argv[1]
f = "path_to_folder/"
f1 = "path_to_folder/" + s
anarci_table_path = f1 + "/anarci/"

_10x_sequences = pd.read_csv(f + "B10x/10x_sequences_with_binding.csv", index_col = 0)
dash_sequences = pd.read_csv(f + "Dash/Dash_all_sequences_with_binding.csv", index_col = 0)
expt_sequences = pd.read_csv(f + "expt/expt_constructs_sequences.csv", index_col = 0)
atlas_sequences = pd.read_csv(f + "atlas/atlas_all_sequences.csv", index_col = 0)
pdb_sequences =  pd.read_csv(f + "PDB_set/PDB_set_all_sequences.csv", index_col = 0)
# zhang_sequences = pd.read_csv(f + "zhang/zhang_all_sequences.csv", index_col = 0)
# os_sequences = pd.read_csv(f + "OS/OS_all_sequences_with_binding.csv", index_col = 0)
newVdj_sequences = pd.read_csv(f + "newVdj/newVdj_all_sequences.csv", index_col = 0)

sequences = {"B10x" : _10x_sequences, "Dash" : dash_sequences,
             "expt" : expt_sequences, "atlas":atlas_sequences, "PDB_set":pdb_sequences,
             "newVdj":newVdj_sequences} # "zhang":zhang_sequences, "OS":os_sequences,

list_of_structures = sequences[s]

def extract_CDR3beta(anarci_table_path, structure_id):
    anarci_table = pd.read_csv(anarci_table_path + structure_id + "_IMGT_numbering_beta.csv", index_col=0)
    CDR3 = []
    for i in anarci_table.index.values:
        j = int(str(i)[0:3])
        # print(i)
        if j >= 104 and j <= 118:
            if anarci_table.loc[i, " Res"] != "-":
                CDR3.append(anarci_table.loc[i, " Res"])
    CDR3_seq = "".join(CDR3)
    return CDR3_seq

def extract_CDR3alpha(anarci_table_path, structure_id):
    anarci_table = pd.read_csv(anarci_table_path + structure_id + "_IMGT_numbering_alpha.csv", index_col=0)
    CDR3 = []
    for i in anarci_table.index.values:
        j = int(str(i)[0:3])
        # print(i)
        if j >= 104 and j <= 118:
            if anarci_table.loc[i, " Res"] != "-":
                CDR3.append(anarci_table.loc[i, " Res"])
    CDR3_seq = "".join(CDR3)
    return CDR3_seq

CDR3_beta_results = pd.DataFrame(columns = ["CDR3_beta", "peptide"])
skipped = 0

for structure_id in list_of_structures.index.values:
    try:
        print(structure_id)
        CDR3beta = extract_CDR3beta(anarci_table_path, structure_id)
        CDR3alpha = extract_CDR3alpha(anarci_table_path, structure_id)
        CDR3_beta_results.loc[structure_id, "CDR3_alpha"] = CDR3alpha
        CDR3_beta_results.loc[structure_id, "CDR3_beta"] = CDR3beta
        CDR3_beta_results.loc[structure_id, "peptide"] = list_of_structures.loc[structure_id, "peptide"]
    except Exception as e:
        print("skipping")
        print(e)
        skipped += 1
        continue
    # print(CDR3_beta_results)
print("skipped: ", skipped)
CDR3_beta_results.to_csv(f1 + "/" + s + "_CDR3AB_with_key.csv")
