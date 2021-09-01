from Bio import SeqUtils
import pandas as pd
import numpy as np
from difflib import SequenceMatcher
from Bio.PDB import *
from itertools import product
from Bio.SeqUtils import seq1

res_numbers = list(range(27, 39)) + list(range(56, 66)) + list(range(105, 112)) + ["111A", "111B",
                "111C", "111D", "111E", "111F", "112E", "112F", "112D", "112C", "112B", "112A", "112"] + list(range(113, 118))
res_numbers = [str(x) for x in res_numbers]
p_nums = [str(x) for x in list(range(1,21))]
CDR1_start = 27
CDR1_end = 38
CDR2_start = 56
CDR2_end = 65
CDR3_start = 105
CDR3_end = 117

def extract_sequence(chain, structure):
    seq_res = []
    seq_num = []
    for res in structure[0][chain]:
        if SeqUtils.seq1(res.resname) != 'X':
            seq_res.append(SeqUtils.seq1(res.resname))
            seq_num.append(str(res.id[1]) + str(res.id[2]))
    sequence = "".join(seq_res) #sequence as string
    seq_num = [x.strip() for x in seq_num]
#    print seq_num, seq_res
    return sequence, seq_num, seq_res

def IMGT_vs_PDB_numbering(IMGT_dict, seq_res, seq_num):
    IMGT = list()
    IMGT_seq = []

    for element in IMGT_dict:
        print(element, element[1])
        IMGT_num = "%03d"%element[0][0] + element[0][1]
        if element[1] != "-":
            IMGT_seq.append(element[1])

    seq_IMGT = "".join(IMGT_seq)
    seq_PDB = "".join(seq_res)
    print("PDB: ", seq_PDB)
    print("IMGT: ", seq_IMGT)
    print(seq_IMGT in seq_PDB)
    print("PDB num: ", seq_num)
    # seq_match = SequenceMatcher()
    # seq_match.set_seq1(seq_IMGT)
    # seq_match.set_seq2(seq_PDB)
    seq_match = SequenceMatcher(autojunk=False, a = seq_IMGT, b = seq_PDB)
    match = seq_match.find_longest_match(0, len(seq_IMGT), 0, len(seq_PDB))
    # print(match)
    l,j,k = match
    assert seq_IMGT[l : l+k] == seq_PDB[j : j+k]
    # print(l, j, k)
    i = j #j is the beginning of the match of PDB residues
    # print(i, j)

    for element in IMGT_dict:
        IMGT_num = "%03d"%element[0][0] + element[0][1]
        res = element[1]
        if res != "-":
            if i < len(seq_num):
                pdb_num = seq_num[i]
                i += 1
        else:
            pdb_num = ""
        new_line = [IMGT_num, res, pdb_num]
        IMGT.append(new_line)
    IMGT = pd.DataFrame(IMGT, columns = ["num_IMGT", "res", "num_pdb"])
    print(IMGT)
    return IMGT

def standard_renumbering(sequence_P, alpha_IMGT, beta_IMGT, folder, filename, seq_num_P, pep_chain = "P", alpha_chain = "A", beta_chain = "B"):
    new_pdb = list()
    i = -1
    first_2 = [1,2] #min is 6, so I want some in the middle always
    last_2 = [19,20] #min is 6, so I want some in the middle always

    len_central = len(sequence_P) - 4 #how many residues excluding the start and end?
    initial = int(11.5 - (len_central - (len_central % 2)) / 2) #number from which to number middle chunk (8 is middle of 15)
    middle_chunk = list(range(initial, initial + len_central))
    numbering = first_2 + middle_chunk + last_2

    final_pdb = []
    alpha_IMGT['num_pdb'] = alpha_IMGT.num_pdb.astype(str)
    beta_IMGT['num_pdb'] = beta_IMGT.num_pdb.astype(str)

    with open(folder + filename) as f:
        pdb_file = f.readlines()
        pdb_file = [x.strip() for x in pdb_file]

    for line in pdb_file:
        if line[0:4] == "ATOM" and SeqUtils.seq1(line[17:20]) != "X":
            line = list(line)
            if line[21] == pep_chain:
                original_num = "".join(line[23:27]).strip()
                index = seq_num_P.index(original_num.strip())
                res_num = numbering[index]
                #print(res_num)
                line[23:26] = "%03d"%res_num
                new_line = "".join(line)
            elif line[21] == alpha_chain:
                original_num = "".join(line[23:27]).strip()
                if original_num in list(alpha_IMGT['num_pdb']):
    #                        print "yes"
                    new_num = alpha_IMGT.loc[alpha_IMGT["num_pdb"] == original_num]["num_IMGT"].iloc[0]
    #                        print original_num, new_num
                    if len(str(new_num)) <= 3:
                        line[23:26] = str(new_num)
                    else:
                        line[23:27] = str(new_num)
                    new_line = "".join(line)
                else:
                    new_line = ""
            elif line[21] == beta_chain:
                original_num = "".join(line[23:27]).strip()
                if original_num in list(beta_IMGT['num_pdb']):
                    new_num = beta_IMGT.loc[beta_IMGT["num_pdb"] == original_num]["num_IMGT"].iloc[0]
                    if len(str(new_num)) <= 3:
                        line[23:26] = str(new_num)
                    else:
                        line[23:27] = str(new_num)
                    new_line = "".join(line)
                else:
                    new_line = ""
            else:
                new_line = "".join(line)
            final_pdb.append(new_line)

    return final_pdb

def calculate_pairwise_distances(chain, chain_name, structure, pep_chain = "P"):
    res_numbers = list(range(27, 39)) + list(range(56, 66)) + list(range(105, 112)) + ["111A", "111B",
                "111C", "111D", "111E", "111F", "112E", "112F", "112D", "112C", "112B", "112A", "112"] + list(range(113, 118))
    res_numbers = [str(x) for x in res_numbers]
    p_nums = [str(x) for x in list(range(1,21))]
    all_atom_distances = list()
#    print(chain, chain_name)
    match_table = np.empty((20, len(res_numbers))) #35 is the sum of lengths of all CDRs, pep length*2 so I can plot alpha and beta together
    match_table[:] = np.nan
    match_table = pd.DataFrame(match_table, index = p_nums, columns = res_numbers)
    min_res1_res2 = list()
    for res1 in structure[0][chain]:
#         print(res1.id)
        if (res1.id[1] >= CDR1_start and res1.id[1] <= CDR1_end) or (res1.id[1] >= CDR2_start and res1.id[1] <= CDR2_end) or (res1.id[1] >= CDR3_start and res1.id[1] <= CDR3_end) :
            if res1.resname != 'HOH':
                res1id = "".join([str(res1.id[1]), str(res1.id[2])]).strip()#
#                 print(res1id, res1id in res_numbers)
                for res2 in structure[0][pep_chain]:
                    res2id = "".join([str(res2.id[1]), str(res2.id[2])]).strip()
                    distances = list()
                    if res2.resname != 'HOH':
                        res2_distances = list()
                        res1_res2_min = list()
                        for atom1 in res1:
                            for atom2 in res2:
                                distance = atom1 - atom2
                                distances.append(distance)
                                new_line = [str(res1id), str(res1.resname),
                                            str(atom1.name),str(res2id), str(res2.resname),
                                            str(atom2.name), distance]
                                all_atom_distances.append(new_line)

                            min_res1_res2.append([str(res1id), str(res1.resname), str(res2id), str(res2.resname), min(distances)])
                            match_table.at[res2id, res1id] = min(distances)
#                             print(match_table)
    all_atom_distances = pd.DataFrame(all_atom_distances, columns = ["res1_num", "res1", "atom1", "res2_num", "res2", "atom2", "distance"])
    min_res1_res2 = pd.DataFrame(min_res1_res2, columns = ["res1_num", "res1", "res2_num", "res2", "min_distance"])

    return all_atom_distances, match_table


def get_energy_table(folder1, structure_id, energy_type):
    energy_table_alpha = pd.read_csv(folder1 +  structure_id + "_residue_pair_energy_alpha.csv", index_col = 0)
    energy_per_pair_alpha = pd.DataFrame(energy_table_alpha[energy_type])
    energy_table_layer_alpha = pd.DataFrame(columns = res_numbers, index = p_nums)

    for i in range(0, np.shape(energy_per_pair_alpha)[0]): #make a table in which energies are plotted + track contributions of types of energy
        match = energy_per_pair_alpha.index[i].split("-")
        CDR = match[0]
        pep = match[1]
        energy = energy_per_pair_alpha.iloc[i][0]
        energy_table_layer_alpha.loc[pep, CDR] = energy

    energy_table_beta = pd.read_csv(folder1 +  structure_id + "_residue_pair_energy_beta.csv", index_col = 0)
    energy_per_pair_beta = pd.DataFrame(energy_table_beta[energy_type])
    energy_table_layer_beta = pd.DataFrame(columns = res_numbers, index = p_nums)

    for i in range(0, np.shape(energy_per_pair_beta)[0]): #make a table in which energies are plotted + track contributions of types of energy
        match = energy_per_pair_beta.index[i].split("-")
        CDR = match[0]
        pep = match[1]
        energy = energy_per_pair_beta.iloc[i][0]
        energy_table_layer_beta.loc[pep, CDR] = energy

#    print(np.shape(energy_table_layer_alpha))
#    print(energy_table_layer_alpha)
    alpha_linear = np.array(energy_table_layer_alpha).reshape(-1)
    beta_linear = np.array(energy_table_layer_beta).reshape(-1)

    col_names = []
    for row in p_nums:
        row = str(row)
        for col in res_numbers:
            col = str(col)
            col_names.append(str(row + "-" + col))

    alpha_col_names = list([str(x + "_alpha") for x in col_names])
    # print(len(alpha_col_names))
    # print(energy_table_layer_alpha.shape)
    # print(alpha_linear.shape)
    beta_col_names = list([str(x + "_beta") for x in col_names])
#    print(np.shape(alpha_linear), np.shape(alpha_col_names))
    alpha_linear = pd.DataFrame(alpha_linear.reshape(-1, len(alpha_linear)), columns = alpha_col_names, index = [structure_id])
    beta_linear = pd.DataFrame(beta_linear.reshape(-1, len(beta_linear)), columns = beta_col_names, index = [structure_id])
#    print(np.shape(alpha_linear))
    linearised_df = pd.concat([alpha_linear, beta_linear], axis = 1)

    return energy_table_layer_alpha, energy_table_layer_beta, linearised_df


def extract_CDR_and_peptide_sequences(structure_id, structure, chains):
    col_names1 = [str(x) + "_" + str(y) for x, y in product(["alpha", "beta"], res_numbers)]
    col_names2 = [str(x) + "_" + str(y) for x, y in product(["peptide"], list(range(1, 21)))]
    col_names = col_names1 + col_names2
    sequence_results = pd.DataFrame(columns = col_names)
    index = structure_id
    for chain_name in ["alpha", "beta"]:
        chain = chains[chain_name]
        for res in structure[0][chain]:
            if (res.id[1] >= CDR1_start and res.id[1] <= CDR1_end):
                sequence_results.loc[index, chain_name + "_" + "".join([str(res.id[1]), str(res.id[2])]).strip()] = seq1(res.resname)
            elif (res.id[1] >= CDR2_start and res.id[1] <= CDR2_end):
                sequence_results.loc[index, chain_name + "_" + "".join([str(res.id[1]), str(res.id[2])]).strip()] = seq1(res.resname)
            elif (res.id[1] >= CDR3_start and res.id[1] <= CDR3_end):
                sequence_results.loc[index, chain_name + "_" + "".join([str(res.id[1]), str(res.id[2])]).strip()] = seq1(res.resname)
#                print(sequence_results)
    for res in structure[0][chains["peptide"]]:
        sequence_results.loc[index, "peptide_" + "".join([str(res.id[1]), str(res.id[2])]).strip()] = seq1(res.resname)

    return sequence_results
