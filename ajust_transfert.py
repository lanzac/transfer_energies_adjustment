#!/usr/bin/env python3
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'
import os

residueName3To1 = {
        "ALA" : "A",
        "ARG" : "R",
        "ASN" : "N",
        "ASP" : "D",
        "CYS" : "C",
        "GLU" : "E",
        "GLN" : "Q",
        "GLY" : "G",
        "HIS" : "H",
        "ILE" : "I",
        "LEU" : "L",
        "LYS" : "K",
        "MET" : "M",
        "PHE" : "F",
        "PRO" : "P",
        "SER" : "S",
        "THR" : "T",
        "TRP" : "W",
        "TYR" : "Y",
        "VAL" : "V"
    }

transfer_exp = {
        "ALA" : -1.3094,
        "ARG" : 4.2870,
        "ASN" : 2.5471,
        "ASP" : 3.2825,
        "CYS" : -6.4753,
        "GLN" : 0.95067,
        "GLU" : 2.7265,
        "GLY" : 0,
        "HIS" : -0.55605,
        "ILE" : -7.6413,
        "LEU" : -7.1928,
        "LYS" : 4.1973,
        "MET" : -5.1300,
        "PHE" : -7.5336,
        "PRO" : 0,
        "SER" : 0.19731,
        "THR" : -1.0762,
        "TRP" : -9.4709,
        "TYR" : -4.0179,
        "VAL" : -5.1300
}

# --------------------------------------------------------------------------------------------------
type_to_ducarmeTypes = {}
# Csp3 single-bonded carbons
type_to_ducarmeTypes["CT"] = ("Csp3", -0.105)

# Csp2 double-bonded and aromatic carbons
for csp2 in ["C", "CA", "CB", "CC", "CD", "CK", "CM", "CN", "CQ", "CR", "CV", "CW", "C*"]:
    type_to_ducarmeTypes[csp2] = ("Csp2", -0.105)

# H(=0) noncharged hydrogen (bound to C)
for Hnc in ["HC","H1","H2","H3","HA","H4","H5","HP","HZ"]:
    type_to_ducarmeTypes[Hnc] = ("Hnc", -0.0397)
# H(/0), charged hydrogen
for Hc in ["H", "HO", "HS", "HW"]:
    type_to_ducarmeTypes[Hc] = ("Hc", 0.0362)
# O, oxygen
for O in ["O", "O2", "OW", "OH", "OS"]:
    type_to_ducarmeTypes[O] = ("O", 0.0403)
# N, nitrogen
for N in ["N", "NA", "NB", "NC", "N2", "N3", "NT", "N*", "NY"]:
    type_to_ducarmeTypes[N] = ("N", 0.112)
# S, sulfur.
for S in ["S", "SH"]:
    type_to_ducarmeTypes[S] = ("S", -0.108)

for charged_null in ["CY", "CZ", "C0", "F", "Cl", "Br", "I", "IM", "IB", "MG", 
                     "P", "CU", "FE", "Li", "IP", "Na", "K", "Rb", "Cs", "Zn"]:
    type_to_ducarmeTypes[charged_null] = ("nc", 0.000)
# --------------------------------------------------------------------------------------------------


res_name3 = "ALA"


# Read nametotype file of all residues
nametotype_path = "/Users/andrelanrezac/Dev/BS_Project/biospring/data/forcefield/nametotype.txt"
nametotype_df = pd.read_csv(nametotype_path, sep="\t", names=["name", "type"], skip_blank_lines=True).dropna()

# Keep only name and type column
nametoype_all_df = nametotype_df[["name", "type"]]

# Add residue code 1 column
nametoype_all_df["res1"] = nametoype_all_df["name"].apply(lambda x: x[0][0])

# Group by res1
nametoype_all_groups = nametoype_all_df.groupby("res1")


for res_name3 in residueName3To1:

    res_name3 = "CYS"

    # Initialize INPUT VALUES
    # typetotransfer_path = "/Users/andrelanrezac/Dev/BS_Project/biospring/data/forcefield/typetotransfer.txt"
    # typetotransfert_df = pd.read_csv(typetotransfer_path, sep="\t", names=["type", "transfer"])
    ##typetotransfert_df = typetotransfert_df.set_index("type")

    # Read nametotype file of all residues
    # nametotype_path = "/Users/andrelanrezac/Dev/BS_Project/biospring/data/forcefield/nametotype.txt"
    # nametotype_df = pd.read_csv(nametotype_path, sep="\t", names=["name", "type"], skip_blank_lines=True).dropna()

    # # Keep only the current residue atom types
    # nametotype_df = nametotype_df[nametotype_df.name.str.startswith(residueName3To1[res_name3])]

    # Keep only name and type column
    #nametoype_all_df = nametotype_df[["name", "type"]]

    # 
    # nametotype_df = nametotype_df.drop("name", axis=1)
    # nametotype_df = nametotype_df.drop_duplicates(subset=["type"])
    # nametotype_df["tr_diff"] = 0.0
    # nametotype_df = nametotype_df.merge(typetotransfert_df, how='left')

    # ----------------------------------------------------------------------------------------------
    # INPUT
    # ----------------------------------------------------------------------------------------------
    residue_group = nametoype_all_groups.get_group(residueName3To1[res_name3])
    residue_types = residue_group
    
    # Filter by unique atom type, each atom type id associated by one transfer energy.
    # Different atom types can share the same transfer energy.
    # All transfert energy ajustments are made in accordance with the 7 atomic types of Ducarme1998.
    # Add references types (Ducarme1998) column -> filter by unique atomic type -> drop type name & res1 columns
    residue_types[["types_ref", "tr_ref"]] = residue_types.apply(lambda x: type_to_ducarmeTypes[x["type"]], axis=1, result_type='expand')
    residue_types = residue_types.drop_duplicates(subset=["tr_ref"])
    residue_types = residue_types.drop(["name", "type", "res1"], axis=1)
    # Add tr_diff column
    residue_types["tr_diff"] = 0.0
    # Reset index from 0
    residue_types = residue_types.reset_index(drop=True)
    print(residue_types)

    # Initialize residue_data
    residue_data = residue_group[["name"]].copy()

    print(residue_data)

    quit()

    amber_ff_path = "/Users/andrelanrezac/Dev/BS_Project/biospring/data/forcefield/amber.ff"
    headers = ["name", "charge", "radius", "epsilon", "mass", "tr_ref"]
    amber_ff_df = pd.read_csv(amber_ff_path, sep="\t", names=headers, skiprows=1)

    # tr_ref and av_sasa are constant values
    # tr_diff and eimp are variable values
    df = amber_ff_df[["name"]]

    df = df.merge(nametoype_all_df, how='left')
    df["av_sasa"] = 0.0
    df["eimp"] = 0.0

    # Filter to keep only the residue considered
    df = df[df.name.str.startswith(residueName3To1[res_name3])]


    # Set particle name as index
    df = df.set_index("name")

    # Exclude backbone atoms
    bb_atoms_names = ["CA","C","O","H","HA","N","H2","H3"]
    df = df.drop([residueName3To1[res_name3]+a for a in bb_atoms_names], errors='ignore')

    # Compute the average sasa value for each particle
    res_pdbs_dir_path = "/Users/andrelanrezac/Dev/ims_systems/impala/transfer_energies_validation/amino_acids_pdbs"
    pdb_count = 0
    for filename in os.listdir(res_pdbs_dir_path):
        if filename.startswith(res_name3):
            pdb_count += 1
            with open(os.path.join(res_pdbs_dir_path, filename), 'r') as f:
                for line in f:
                    atom_name = line[12:16].strip()
                    res_name3 = line[17:20]
                    sasa = line[60:66] # wrote in temp factor by freesasa
                    particle_name = residueName3To1[res_name3] + atom_name
                    # Add sasa of each particle
                    if (atom_name in bb_atoms_names): continue
                    df.loc[particle_name, "av_sasa"] += float(sasa)
                    
    df.av_sasa /= pdb_count


    # Adjuts transfert by modifying tr_diff
    success=False
    for i in range(1000):
        nametotype_df["tr_diff"] = np.random.uniform(-0.05, 0.05, size=nametotype_df.shape[0])
        nametotype_df["tr_new"] = nametotype_df["tr_diff"] + nametotype_df["transfer"]
        df = df.merge(nametotype_df[["type", "tr_new"]], how='left')
        # ----
        #print(nametotype_df) # <- tr_diff column is the input

        # Compute Eimp
        df.eimp = -df.av_sasa * df.tr_new * -0.5 + (-0.018 * df.av_sasa * -0.5)

        
        eimp_tot = sum(df.eimp)
        abs_diff_imp = abs(eimp_tot - transfer_exp[res_name3])
        
        if (abs_diff_imp < 0.01):
            print("*"*80)
            print(res_name3)
            print(nametotype_df)
            #print(df)
            print("Diff=",abs_diff_imp)
            print("Eimp=",eimp_tot)
            print("Eimp_exp=", transfer_exp[res_name3])
            success=True
            break
        df = df.drop("tr_new", axis=1)
    
    if (not success):
        print("="*80)
        print("adjusted transfert not found for res ", res_name3)
        print("="*80)

