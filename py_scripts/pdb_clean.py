#function to parse pdb structures and write new files with only protein atoms

import pandas as pd
import MDAnalysis
import rdkit
from tqdm import tqdm
import os


def clean_pdbs(pdb_dir_path, index_file_path):

    pdbs = []
    with open(index_file_path, "r") as f:
        for line in f.readlines()[6:]:
            pdbs.append(line.split()[0])

    for pdb in tqdm(pdbs):
        missing = 0
        pdb_path = pdb_dir_path + pdb + f"/{pdb}_protein.pdb"
        if os.path.exists(pdb_path):
            u = MDAnalysis.Universe(pdb_path)
            protein = u.select_atoms("protein")
            protein.write(pdb_dir_path + pdb + f"/{pdb}_protein_cleaned")
        else:
            missing += 1
    print(missing)


clean_pdbs(
    "Data/pdbbind_2020_general/",
    "Data/pdbbind_2020_general/index/INDEX_general_PL_data.2020",
)
