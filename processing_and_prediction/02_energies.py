## this script does not contain functions because PyRosetta gives a segFault if used...

from sys import argv
import os
import _pairwise_energies
import pandas as pd
from Bio.PDB import *
import pyrosetta
pyrosetta.init("-corrections::restore_talaris_behavior") ##this to revert to scorefxn talaris2014 (rather than ref15)
#from pyrosetta import * #for pose_from_pdb
#from pyrosetta.rosetta import * #for EMapVector
#from pyrosetta.rosetta.core.scoring import EMapVector
from pyrosetta.toolbox import cleanATOM
from pyrosetta import get_fa_scorefxn
from pyrosetta.rosetta.core.scoring import *
from pyrosetta import pose_from_pdb
#from pyrosetta.rosetta.core.scoring import EMapVector
"""
From rosetta forum support:

there used to be only one top-level namespace, now there's both `pyrosetta` and
`rosetta` - the latter being things which are from Rosetta's C++ library, and
the former being code which is PyRosetta specific.
(Mainly a number of convenience functions.)
In your case, things like fa_atr and EMapVector can be found in the rosetta.core.scoring
namespace. They may or may not be loaded with a `from rosetta import *` command -
you may need to import them (or the nested namespace) specifically.
"""

class TCRpep_select(Select): #function to create TCR + peptide PDB
    def accept_chain(self, chain, alpha_chain = "A", beta_chain = "B", pep_chain = "P"):
        if chain.full_id[2] == str(alpha_chain) or chain.full_id[2] == str(beta_chain) or chain.full_id[2] == str(pep_chain):
            return 1
        else:
            return 0

structure_id = argv[1].split("_TCR-pMHC.pdb")[0]
filename = "renumbered_pdb/" + structure_id + "_renumbered.pdb"
folder = "/home/regmili/Scratch/" + argv[2]
if not os.path.exists(folder+"Results_by_pdb/energies/"):
    os.mkdir(folder+"Results_by_pdb/energies/")
if not os.path.exists(folder+"TCR_peptide_only/"):
    os.mkdir(folder+"TCR_peptide_only/")

res_numbers = list(range(27, 39)) + list(range(56, 66)) + list(range(105, 112)) + ["111A", "111B",
                "111C", "111D", "111E", "111F", "112E", "112F", "112D", "112C", "112B", "112A", "112"] + list(range(113, 118))
res_numbers = [str(x) for x in res_numbers]
p_nums = [str(x) for x in list(range(1,21))]

Parser = PDBParser()
structure = Parser.get_structure(structure_id, folder + filename)
print("--------------------------------")
#print("initialising score function...")
scorefxn = get_fa_scorefxn()
#print("initialised")

energy_terms = [fa_atr, fa_elec, fa_rep, fa_sol]
column_names1 = ["fa_atr", "fa_elec", "fa_rep", "fa_sol"]

pep_chain = "P"
alpha_chain = "A"
beta_chain = "B"

#print()
#print("------------------------------------")
print()
print("calculating pairwise energies...")
##make PDB with only peptide and TCR to localise the PyRosetta search
#print("save structure with only peptide and TCR")
io = PDBIO()
io.set_structure(structure)
io.save(folder + "TCR_peptide_only/" + structure_id +"_TCR_peptide.pdb", TCRpep_select())
os.chdir(folder + "TCR_peptide_only/")
cleanATOM(structure_id +"_TCR_peptide.pdb")
#print("--------------------------------")
#print("load pose from pdb")
pdb_TCR_pep = pose_from_pdb(str(folder + "TCR_peptide_only/" + structure_id +"_TCR_peptide.clean.pdb"))
#print("loaded")
#print("--------------------------------")
os.chdir("path_to_folder/")

print("scoring pose...")
scorefxn.show(pdb_TCR_pep)

print("alpha calculations...")

energies = list()
res_couple = list()

# pdb_TCR_pep.pdb_info().icode(i) tells you whether there is an "A" after the res number
# pdb_TCR_pep.pdb_info().number(i) gives you the res number

for i in range(1, pdb_TCR_pep.total_residue() + 1):
    # print(pdb_TCR_pep.pdb_info().chain(i))
    # print(pdb_TCR_pep.pdb_info().number(i), pdb_TCR_pep.pdb_info().icode(i))
    if (pdb_TCR_pep.pdb_info().chain(i) == alpha_chain) and (str(pdb_TCR_pep.pdb_info().number(i)) in res_numbers):
        rsd1 = pdb_TCR_pep.residue(i)
        for j in range(1, pdb_TCR_pep.total_residue() + 1):
            print(j)
            if pdb_TCR_pep.pdb_info().chain(j) == pep_chain:
                emap = EMapVector()
                energy_vector = list()
                rsd2 = pdb_TCR_pep.residue(j)
                scorefxn.eval_ci_2b(rsd1, rsd2, pdb_TCR_pep, emap)
                a = str(pdb_TCR_pep.pdb_info().number(i))
                # print(a)
                # print(pdb_TCR_pep.pdb_info().number(i))
                a1 = str(pdb_TCR_pep.pdb_info().icode(i))
                b = str(pdb_TCR_pep.pdb_info().number(j))
                b1 = str(pdb_TCR_pep.pdb_info().icode(j))
                res_couple.append(("".join([a, a1]).strip(), "".join([b, b1]).strip()))
                # print(res_couple)
                for term in energy_terms:
                    energy_vector.append(emap[term])
                energies.append(energy_vector)
print("energies calculated, formatting...")
residue_energies_alpha = pd.DataFrame(columns = column_names1)

for i in range(0, len(res_couple)):
    pair = res_couple[i]
    print(pair)
    energy_pair = energies[i]
    residue_energies_alpha.loc["-".join(pair)] = energy_pair

residue_energies_alpha.to_csv(folder + "Results_by_pdb/energies/" + structure_id + "_residue_pair_energy_alpha.csv")

print("beta calculations...")

energies = list()
res_couple = list()

for i in range(1, pdb_TCR_pep.total_residue() + 1):
    if pdb_TCR_pep.pdb_info().chain(i) == beta_chain and str(pdb_TCR_pep.pdb_info().number(i)) in res_numbers:
        rsd1 = pdb_TCR_pep.residue(i)
        for j in range(1, pdb_TCR_pep.total_residue() + 1):
            if pdb_TCR_pep.pdb_info().chain(j) == pep_chain:
                emap = EMapVector()
                energy_vector = list()
                rsd2 = pdb_TCR_pep.residue(j)
                scorefxn.eval_ci_2b(rsd1, rsd2, pdb_TCR_pep, emap)
                a = str(pdb_TCR_pep.pdb_info().number(i))
                a1 = str(pdb_TCR_pep.pdb_info().icode(i))
                b = str(pdb_TCR_pep.pdb_info().number(j))
                b1 = str(pdb_TCR_pep.pdb_info().icode(j))
                res_couple.append(("".join([a, a1]).strip(), "".join([b, b1]).strip()))
                # print(res_couple)
                for term in energy_terms:
                    energy_vector.append(emap[term])
                energies.append(energy_vector)
print("energies calculated, formatting...")
residue_energies_beta = pd.DataFrame(columns = column_names1)

for i in range(0, len(res_couple)):
    pair = res_couple[i]
    energy_pair = energies[i]
    residue_energies_beta.loc["-".join(pair)] = energy_pair

residue_energies_beta.to_csv(folder + "Results_by_pdb/energies/" + structure_id + "_residue_pair_energy_beta.csv")

fa_atr, fa_elec, fa_rep, fa_sol = _pairwise_energies.energies_of_interest(argv[1], argv[2])

print("done")
#print(fa_atr, fa_elec)
print()
print("------------------------------------")
