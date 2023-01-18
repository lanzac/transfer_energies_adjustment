#!/usr/bin/env python3
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'
import os

res_pdbs_dir_path = "/Users/andrelanrezac/Dev/ims_systems/impala/transfer_energies_validation/amino_acids_pdbs"

residueName3To1 = {
        "ALA" : "A",
        "ARG" : "R",
        "ASN" : "N",
        "ASP" : "D",
        "CYS" : "C",
        "GLN" : "Q",
        "GLU" : "E",
        "GLY" : "G",
        "HIS" : "H",
        "ILE" : "I",
        "LEU" : "L",
        "LYS" : "K",
        "MET" : "M",
        "PHE" : "F",
        #"PRO" : "P",
        "SER" : "S",
        "THR" : "T",
        "TRP" : "W",
        "TYR" : "Y",
        "VAL" : "V"
    }

residueName1To3 = {
        "R" : "ARG",
        "N" : "ASN",
        "A" : "ALA",
        "D" : "ASP",
        "C" : "CYS",
        "E" : "GLU",
        "Q" : "GLN",
        "G" : "GLY",
        "H" : "HIS",
        "I" : "ILE",
        "L" : "LEU",
        "K" : "LYS",
        "M" : "MET",
        "F" : "PHE",
        #"P" : "PRO",
        "S" : "SER",
        "T" : "THR",
        "W" : "TRP",
        "Y" : "TYR",
        "V" : "VAL"
    }

transfer_exp = {
        "ALA" : -1.3094,
        "ARG" : 4.2870,
        "ASN" : 2.5471,
        "ASP" : 3.2825,
        "CYS" : -6.4753,
        "GLN" : 0.95067,
        "GLU" : 2.7265,
        "GLY" : 0.0001,
        "HIS" : -0.55605,
        "ILE" : -7.6413,
        "LEU" : -7.1928,
        "LYS" : 4.1973,
        "MET" : -5.1300,
        "PHE" : -7.5336,
        #"PRO" : 0,
        "SER" : 0.19731,
        "THR" : -1.0762,
        "TRP" : -9.4709,
        "TYR" : -4.0179,
        "VAL" : -5.1300
}

all_new_Csp3_tr = {
        "ALA" :  "failed",
        "ARG" :  "failed",
        "ASN" :  "failed",
        "ASP" :  "failed",
        "CYS" :  "failed",
        "GLN" :  "failed",
        "GLU" :  "failed",
        "GLY" :  "failed",
        "HIS" :  "failed",
        "ILE" :  "failed",
        "LEU" :  "failed",
        "LYS" :  "failed",
        "MET" :  "failed",
        "PHE" :  "failed",
        "PRO" :  "failed",
        "SER" :  "failed",
        "THR" :  "failed",
        "TRP" :  "failed",
        "TYR" :  "failed",
        "VAL" :  "failed"
}

bb_atoms_names = ["CA","C","O","H","HA","N","H2","H3"]

# --------------------------------------------------------------------------------------------------

def update_transfert_energies(res_type, diff_ener):

    # Update residue_types with new tranfert energies
    global residue_types
    residue_types["tr_new"] = new_Csp3 * residue_types["tr_ratio"]

    # Update residue_data tr_new from residue_types and keep atomname index
    global residue_data
    if ("tr_new" in residue_data.columns):
        residue_data = residue_data.drop("tr_new", axis=1)
    residue_data = residue_data.merge(residue_types[["types_ref", "tr_new"]], how='left').set_index(residue_data.index)

    # Update backbone_data & sidechain_data tr_new from residue_data
    global backbone_data
    backbone_data  = residue_data[is_bb]
    global sidechain_data
    sidechain_data = residue_data[~is_bb]

    return new_Csp3


def compute_average_sasa(res3, res_data):

    res_data.set_index("atomName", inplace=True)

    pdb_count = 0
    for filename in os.listdir(res_pdbs_dir_path):
        if filename.startswith(res3):
            pdb_count += 1
            with open(os.path.join(res_pdbs_dir_path, filename), 'r') as f:
                for line in f:
                    atom_name = line[12:16].strip()
                    res3 = line[17:20]
                    sasa = line[60:66] # wrote in temp factor by freesasa
                    # Add sasa of each particle
                    #if (atom_name in bb_atoms_names): continue
                    try:
                        res_data.loc[atom_name, "av_sasa"] += float(sasa)  
                    except:
                        # TO DOT : create a log error file
                        continue
                        #print("atom: " + atom_name + " not found in file " + filename + " in line \n:" + line)        
    res_data["av_sasa"] /= pdb_count


# --------------------------------------------------------------------------------------------------

ducarmeType_to_atomType = {
    "Csp3": ["CT"],
    "Csp2": ["C", "CA", "CB", "CC", "CD", "CK", "CM", "CN", "CQ", "CR", "CV", "CW", "C*"],
    "Hnc": ["HC","H1","H2","H3","HA","H4","H5","HP","HZ"],
    "Hc": ["H", "HO", "HS", "HW"],
    "O": ["O", "O2", "OW", "OH", "OS"],
    "N": ["N", "NA", "NB", "NC", "N2", "N3", "NT", "N*", "NY"],
    "S": ["S", "SH"],
    "nc": ["CY", "CZ", "C0", "F", "Cl", "Br", "I", "IM", "IB", "MG", "P", "CU", "FE", "Li", "IP", "Na", "K", "Rb", "Cs", "Zn"]
}

def get_ducarmeType(atomType):
    for ducarmeType, atomTypes in ducarmeType_to_atomType.items():
        if atomType in atomTypes:
            return ducarmeType
    return None

ducarmeType_to_transferEnergies = {"Csp3" : -0.105, 
                                    "Csp2" : -0.0134,
                                    "Hnc"  : -0.0397,
                                    "Hc"   : 0.0362,
                                    "O"    : 0.0403,
                                    "N"    : 0.112,
                                    "S"    : -0.108,
                                    #"nc"   : 0.000
                                    }




# --------------------------------------------------------------------------------------------------
# Final table
# Initialize the DataFrame with the given data
final_residues_ducarmeType = pd.DataFrame(columns=residueName3To1.keys(), index= list(ducarmeType_to_transferEnergies.keys()) + ["eimp_ref", "eimp_new"])

# Add all experimental Ducarme sidechains transfer energies values (Fig3 Ducarme1998) to eimp_ref row
final_residues_ducarmeType.loc["eimp_ref", transfer_exp.keys()] = list(transfer_exp.values())

eimp_new_all = dict.fromkeys(residueName3To1.keys(), None)


print(final_residues_ducarmeType)

# --------------------------------------------------------------------------------------------------
# Configure residues_ducarmeType
residues_ducarmeType = pd.DataFrame.from_dict(ducarmeType_to_transferEnergies, orient='index', columns=['tr_ref'])
residues_ducarmeType = residues_ducarmeType.rename_axis('ducarmeType')
# Set 'tr_diff' column to 0.0
residues_ducarmeType = residues_ducarmeType.assign(tr_diff=0.0)
# Set 'tr_new' column to values in 'tr_ref' column
residues_ducarmeType = residues_ducarmeType.assign(tr_new=residues_ducarmeType['tr_ref'])


# --------------------------------------------------------------------------------------------------
# Read nametotype file of all residues
residues_data_path = "/Users/andrelanrezac/Dev/BS_Project/biospring/data/forcefield/nametotype.txt"
residues_data_df = pd.read_csv(residues_data_path, sep='\t', header=None, names=['particleName', 'atomType'])
# Get residue code 1 from particleName column and add to a new res1 column
residues_data_df['res1'] = residues_data_df['particleName'].str.slice(0, 1)
# Get atom name from particleName column and add to a new atomName column
residues_data_df['atomName'] = residues_data_df['particleName'].str.slice(1)
# Drop particleName column
residues_data_df.drop(columns=["particleName"], inplace=True)
# Add average sasa column
residues_data_df = residues_data_df.assign(av_sasa=0.0)
# Add eimp column
residues_data_df = residues_data_df.assign(eimp=0.0)
# Get Ducarme type from each atomType and add it to a new column called ducarmeType
residues_data_df['ducarmeType'] = residues_data_df['atomType'].apply(get_ducarmeType)
residues_data_df.set_index('ducarmeType', inplace=True)

# Merge residues_ducarmeType to residues_data_df
residues_data_df = residues_data_df.join(residues_ducarmeType, how="right")

# Group data by res1
residues_data_df_grouped = residues_data_df.groupby("res1")


# --------------------------------------------------------------------------------------------------
# Compute all average sasa

for res3 in residueName3To1:
    print("-"*80)
    print("->", res3)
    res_data = residues_data_df_grouped.get_group(residueName3To1[res3]).reset_index()
    compute_average_sasa(res3, res_data)

    # Clean a bit
    res_data.reset_index(inplace=True)
    res_data.drop(columns=["res1"], inplace=True)
    res_data.set_index('ducarmeType', inplace=True)

    is_bb = res_data.atomName.isin(bb_atoms_names) 
    backbone_data  = res_data[is_bb]
    sidechain_data = res_data[~is_bb]



    # Get ducarmeTypes that transfert energie will not be ajusted (all types present (in bb and not in sidechain) OR (not in bb AND not in sidechain)).
    # Group the data by ducarmeType
    ducarmeType_grouped = is_bb.groupby('ducarmeType')
    # Filter the groups that have all True values
    fixed_ducarmeTypes_indices = ducarmeType_grouped.filter(lambda x: x.all()).index

    missing_type = list(set(ducarmeType_to_transferEnergies.keys()) - (set(res_data.index)))
    print("Missing types: {}".format(missing_type))

    fixed_ducarmeTypes_indices = fixed_ducarmeTypes_indices.union(pd.Index(missing_type))

    print("Fixed ducarmeTypes: {}".format(fixed_ducarmeTypes_indices.values))



    diff_ener = 0.0001
    # update_transfert_energies(residues_ducarmeType, diff_ener)


    
    

    # residues_ducarmeType = residues_ducarmeType.assign(gradient=0.0)
    # residues_ducarmeType = residues_ducarmeType.assign(direction=0.0)

    residues_ducarmeType = residues_ducarmeType.assign(norm_ref=(residues_ducarmeType.tr_ref-residues_ducarmeType.tr_ref.min())/
                                                                (residues_ducarmeType.tr_ref.max() - residues_ducarmeType.tr_ref.min()) -0.5)
    # [-1;1] : norm_ref=2*(residues_ducarmeType.tr_ref-residues_ducarmeType.tr_ref.min())/  (residues_ducarmeType.tr_ref.max() - residues_ducarmeType.tr_ref.min()) -1
    # [-0.5;0.5] norm_ref=(residues_ducarmeType.tr_ref-residues_ducarmeType.tr_ref.min())/  (residues_ducarmeType.tr_ref.max() - residues_ducarmeType.tr_ref.min()) -0.5

    # prev_gradient = None
    success = False
    sign = 1
    N = 100000
    ind = 0
    for i in range(N):

        # update sidechain.tr_new with residues_ducarmeType.tr_new
        # sidechain_data = sidechain_data.reset_index().merge(residues_ducarmeType[["tr_new"]], how='left').set_index('ducarmeType')
        # sidechain_data.tr_new = residues_ducarmeType.tr_new.map(residues_ducarmeType.index)
        sidechain_data.update(residues_ducarmeType)
        
        # Get output value
        sidechain_data.eimp = -sidechain_data.av_sasa * sidechain_data.tr_new * -0.5 + (-0.018 * sidechain_data.av_sasa * -0.5)
        eimp_tot = sum(sidechain_data.eimp) # output value

        if (i==0): previous_imp = eimp_tot

        diff_imp = eimp_tot - transfer_exp[res3] 

        
        
        if (abs(diff_imp) > abs(previous_imp - transfer_exp[res3])):
            # Change sign
            sign *= -1

        if (sign > 0): ind += i
        else:          ind -= i
        
        previous_imp = eimp_tot
        
        # if (abs(diff_imp) < 0.1):
        if (abs(diff_imp) < float(i)/N):
            #print(sidechain_data)
            print("Final diff=",abs(diff_imp))
            print("Final Eimp=",eimp_tot)
            # eimp_new_all[res3] = eimp_tot
            print("Eimp_exp=", transfer_exp[res3])
            print("--- i:{0} {1:9.3f} < {2:9.3f}".format(i, abs(diff_imp), float(i)/N))
            success=True
            break



        # Update tr_new
        # diff_ener += sign * 0.001
        
        
        #diff_ener = sign * 0.001
        diff_ener = (float(ind)/N / (1-residues_ducarmeType.tr_diff))
        
        # residues_ducarmeType.tr_new = residues_ducarmeType.tr_ref + residues_ducarmeType.norm_ref * diff_ener
        residues_ducarmeType.tr_new = residues_ducarmeType.tr_ref + diff_ener

        # is_bb_and_not_sidechain = residues_ducarmeType.index.isin(backbone_data.index)
        if (fixed_ducarmeTypes_indices.values.any()):
            # residues_ducarmeType.update(backbone_data.loc[fixed_ducarmeTypes_indices]) # Stay fixed fixed_ducarmeTypes in residues_ducarmeType
            residues_ducarmeType.tr_new.update(residues_ducarmeType.loc[fixed_ducarmeTypes_indices].tr_ref) # Stay fixed "fixed_ducarmeTypes" in residues_ducarmeType
        # residues_ducarmeType.tr_new += diff_ener

        tr_diff = (residues_ducarmeType.tr_new - residues_ducarmeType.tr_ref) / residues_ducarmeType.tr_ref
        residues_ducarmeType.tr_diff = tr_diff

        print("--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 > {3:9.3f}".format(i,N, abs(diff_imp), float(i)/N), end="\r")
        #print(str(i) + "\n" + str(round(diff_imp, 2)) + "\n" + str(round(float(ind)/N, 2)), end="\r")


        # cost_imp = (diff_imp**2)**0.5
        # # update gradient
        # residues_ducarmeType.gradient = 2*(diff_imp)*residues_ducarmeType.tr_new

        # if (prev_gradient) is not None:
        #     beta = residues_ducarmeType.gradient**2 / prev_gradient**2
        #     residues_ducarmeType.direction = -residues_ducarmeType.gradient + beta*residues_ducarmeType.direction
        # else:
        #     residues_ducarmeType.direction = -residues_ducarmeType.gradient
        # prev_gradient = residues_ducarmeType.gradient
        # residues_ducarmeType.tr_new += diff_ener * residues_ducarmeType.direction


    if (not success):
        print("="*80)
        print("adjusted transfert not found for res ", res3)
        print("="*80)

    eimp_new_all[res3] = eimp_tot
    print("-> Residue atomic types")
    print(residues_ducarmeType)
    print("-> Backbone")
    print(backbone_data)
    print("-> Sidechain")
    print(sidechain_data)

    # final_residues_ducarmeType[res3] = residues_ducarmeType.tr_new

    final_residues_ducarmeType[res3].update(residues_ducarmeType.tr_new)

    final_residues_ducarmeType.loc["eimp_new", eimp_new_all.keys()] = list(eimp_new_all.values())

    print("-> Final residues ducarmeType")
    print(final_residues_ducarmeType)

    # quit()


quit()

# pdb_count = 0

# for filename in os.listdir(res_pdbs_dir_path):
#     if(not filename.endswith(".pdb")): continue
#     res3 = filename[:3]
#     res1 = residueName3To1[res3]
#     pdb_count[res1] += 1
#     res_data = residues_data_df_grouped.get_group(residueName3To1[res3])
    
#     with open(os.path.join(res_pdbs_dir_path, filename), 'r') as f:
#         for line in f:
#             atom_name = line[12:16].strip()
#             res3 = line[17:20]
#             sasa = line[60:66] # wrote in temp factor by freesasa
#             # Add sasa of each particle
#             #if (atom_name in bb_atoms_names): continue
#             atom_type = res_data.loc[res_data['atomType'] == atom_name]
#             atom_type['av_sasa'] += float(sasa)

#     print(res_data)
    
# #     residue_data["av_sasa"] /= pdb_count

# for res1, data in residues_data_df_grouped:
#     print("-"*80)
#     print("->", res1)
#     print(data)
    
# quit()


# Read nametotype file of all residues
nametotype_path = "/Users/andrelanrezac/Dev/BS_Project/biospring/data/forcefield/nametotype.txt"
nametotype_df = pd.read_csv(nametotype_path, sep="\t", names=["name", "type"], skip_blank_lines=True, keep_default_na=False).dropna()

# Keep only name and type column
nametoype_all_df = nametotype_df[["name", "type"]]

# Add residue code 1 column
nametoype_all_df["res1"] = nametoype_all_df["name"].apply(lambda x: x[0][0])

# Group by res1
nametoype_all_groups = nametoype_all_df.groupby("res1")


for res_name3 in residueName3To1:
    print("-"*80)
    print("->", res_name3)

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
    global residue_types
    residue_types = residue_group
    
    # Filter by unique atom type, each atom type id associated by one transfer energy.
    # Different atom types can share the same transfer energy.
    # All transfert energy ajustments are made in accordance with the 7 atomic types of Ducarme1998.
    # Add references types (Ducarme1998) column -> filter by unique atomic type -> drop type name & res1 columns
    residue_types[["types_ref", "tr_ref"]] = residue_types.apply(lambda x: type_to_ducarmeTypes[x["type"]], axis=1, result_type='expand')
    residue_types = residue_types.drop_duplicates(subset=["types_ref"])
    residue_types = residue_types.drop(["name", "type", "res1"], axis=1)

    # add tr_ratio with respect to Csp3 atomtype reference
    residue_types["tr_ratio"] = residue_types["tr_ref"] / -0.1050

    # Reset index from 0
    residue_types = residue_types.reset_index(drop=True)
    
    # Initialize residue_data
    global residue_data
    residue_data = residue_group[["name", "type"]].copy()
    # Get atomname, remove res1 to particle name
    residue_data["atomname"] = residue_data["name"].apply(lambda x: x[1:])
    residue_data = residue_data.drop("name", axis=1)
    # Get types_ref from type
    residue_data[["types_ref", "tr_ref"]] = residue_data.apply(lambda x: type_to_ducarmeTypes[x["type"]], axis=1, result_type='expand')
    residue_data = residue_data.drop("type", axis=1)

    # Set atomname column to index
    residue_data = residue_data.set_index("atomname")

    # Split residue_data to backbone_data and  sidechain_data
    bb_atoms_names = ["CA","C","O","H","HA","N","H2","H3"]
    #global is_bb
    is_bb = residue_data.index.isin(bb_atoms_names) 

    # Compute the average sasa value for each particle and update/add av_sasa to sidechain_data
    # Call before update_transfert_energies to keep av_sava column in backbone_data and sidechain_data
    compute_average_sasa(res_name3)

    diff_Csp3 = 0.00
    new_Csp3 = 0.00
    sign = 1
    update_transfert_energies(diff_Csp3)
    
    

    # Adjuts transfert by modifying tr_diff
    success=False
    for i in range(10000):
        # Reset eimp
        sidechain_data["eimp"] = 0.0

        # Compute Eimp
        sidechain_data.eimp = -sidechain_data.av_sasa * sidechain_data.tr_new * -0.5 + (-0.018 * sidechain_data.av_sasa * -0.5)
        eimp_tot = sum(sidechain_data.eimp)
        abs_diff_imp = abs(eimp_tot - transfer_exp[res_name3])


        if (i > 0): # keep initial tr values for the first step
            #diff_Csp3 = np.random.uniform(-0.5, 0.5) * 2*abs(ref_diff_imp)
            diff_Csp3 += sign * 0.0001

            # Check if with the new diff_Csp3 value we move away from the ref value
            if (abs_diff_imp > abs(last_imp - transfer_exp[res_name3])):
                # Change sign
                sign *= -1
                #new_Csp3 = update_transfert_energies(-diff_Csp3)
            else:
                new_Csp3 = update_transfert_energies(diff_Csp3)

        last_imp = eimp_tot

        # print("new_Csp3: ", new_Csp3)
        #print(res_name3, " eimp_tot: ", eimp_tot, "eimp_exp: ", transfer_exp[res_name3])

        if (i == 0):
            ref_diff_imp = eimp_tot - transfer_exp[res_name3]
            print("imp ini:", eimp_tot)
            print("diff_imp:", ref_diff_imp)
        
        if (abs_diff_imp < 0.1):
            #print(sidechain_data)
            print("Final diff=",abs_diff_imp)
            print("Final Eimp=",eimp_tot)
            print("Eimp_exp=", transfer_exp[res_name3])
            all_new_Csp3_tr[res_name3] = new_Csp3
            success=True
            break
        # else:
        #     print("Diff=",abs_diff_imp)

    print("-> Residue atomic types")
    print(residue_types)
    print("-> Backbone")
    print(backbone_data)
    print("-> Sidechain")
    print(sidechain_data)
        
    
    if (not success):
        print("="*80)
        print("adjusted transfert not found for res ", res_name3)
        print("="*80)

    for res, tr in all_new_Csp3_tr.items():
        print(res, " ", tr)

    
"""
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

"""