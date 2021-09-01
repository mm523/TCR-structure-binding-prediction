from sys import argv
import _renumber
import _distances
import _extract_sequence_AF
import pandas as pd
from Bio.PDB import *

structure_id = argv[1].split("_TCR-pMHC.pdb")[0]
filename = "renumbered_pdb/" + structure_id + "_renumbered.pdb"
folder = "path_to_folder/" + argv[2]

print("------------------------------------")
print()
print("renumbering input structure...")
_renumber.renumbering(argv[1], argv[2])
print("done")
print()
print("------------------------------------")
print()
print("load input structure...")
Parser = PDBParser()
structure = Parser.get_structure(structure_id, folder + filename)
print()
print("calculating parwise distances...")
distances = _distances.calculate_distances(argv[1], argv[2], structure)
print("done")
print()
print("------------------------------------")
print()
print("extracting sequence information...")
atchley = _extract_sequence_AF.atchley_factors(argv[1], argv[2]) #atchley factors
print("success!")
#print(distances, atchley)
print()
