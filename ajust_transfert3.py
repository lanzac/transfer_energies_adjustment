#!/usr/bin/env python3
import pandas as pd
pd.options.display.float_format = '{0:0.4f}'.format
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'
import os
from os import system

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
        "PRO" : "P",
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
        "P" : "PRO",
        "S" : "SER",
        "T" : "THR",
        "W" : "TRP",
        "Y" : "TYR",
        "V" : "VAL"
    }

# Table III Fauchère&Pliska x KCALTOKJ = Fig3 Ducarme1998
transfer_exp = {
        "ALA" : -1.29704,
        "ARG" : 4.22584,
        "ASN" : 2.51,
        "ASP" : 3.22168,
        "CYS" : -6.44336,
        "GLN" : 0.92048,
        "GLU" : 2.67776,
        "GLY" : 0,
        "HIS" : -0.54392,
        "ILE" : -7.53,
        "LEU" : -7.11,
        "LYS" : 4.14216,
        "MET" : -5.14632,
        "PHE" : -7.48936,
        "PRO" : -3.01248,
        "SER" : 0.16736,
        "THR" : -1.08784,
        "TRP" : -9.414,
        "TYR" : -4.01664,
        "VAL" : -5.10448,
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
final_residues_ducarmeType = pd.DataFrame(columns=residueName3To1.keys(), index= list(ducarmeType_to_transferEnergies.keys()) + ["eimp_ref", "eimp_new", "eimp_ini"])

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

new_Csp3 = residues_ducarmeType.loc["Csp3", "tr_ref"] # Initialization
# METHOD 1
# add tr_ratio with respect to Csp3 atomtype reference
method=1
# match method:
#     case 1:
residues_ducarmeType["tr_ratio"] = residues_ducarmeType["tr_ref"] / new_Csp3

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

    min_method_score = 10000000.0
    best_method = 0
    # Loop on all method
    # for m in range(1,5):
    for m in [4]:

        method = m
        # Reinit
        residues_ducarmeType.tr_new = residues_ducarmeType.tr_ref
        
        delta_ener = 0.0001
        success = False
        sign = 1
        N = 100000
        ind = 0
        for i in range(N):

            sidechain_data.update(residues_ducarmeType)
            
            # Get output value
            sidechain_data.eimp = -sidechain_data.av_sasa * sidechain_data.tr_new * -0.5 + (-0.018 * sidechain_data.av_sasa * -0.5)
            eimp_tot = sum(sidechain_data.eimp) # output value

            if (i==0):
                previous_imp = eimp_tot
                print("-> Initial Ducarme sidechain transfer energie: {}".format(eimp_tot))
                final_residues_ducarmeType.loc["eimp_ini", res3] = "{0:0.4f}".format(eimp_tot)

                print("-> Sidechain")
                print(sidechain_data)


            diff_imp = eimp_tot - transfer_exp[res3] 

            if (abs(diff_imp) > abs(previous_imp - transfer_exp[res3])):
                # Change sign
                sign *= -1

            if (sign > 0): ind += i
            else:          ind -= i
            
            previous_imp = eimp_tot
            
            # if (abs(diff_imp) < 0.1):
            match method:
                case 1 | 2 | 3:
                    break_condition = abs(diff_imp) < 0.1
                case 4:
                    break_condition = abs(diff_imp) < float(i)/N

            if (break_condition):
                #print(sidechain_data)
                # print("Final diff=",abs(diff_imp))
                # print("Final Eimp=",eimp_tot)
                # eimp_new_all[res3] = eimp_tot
                # print("Eimp_exp=", transfer_exp[res3])
                # print("--- i:{0} {1:9.3f} < {2:9.3f}".format(i, abs(diff_imp), float(i)/N))
                success=True

                # Compute method score
                method_score = sum(((residues_ducarmeType.tr_new - residues_ducarmeType.tr_ref) / residues_ducarmeType.tr_ref).abs())

                print("\nMethod {} score : {}".format(m, method_score))
                if (m==0):
                    min_method_score = method_score
                    best_tr_new = residues_ducarmeType.tr_new

                if (method_score < min_method_score):
                    min_method_score = method_score
                    best_method = method
                    best_tr_new = residues_ducarmeType.tr_new
                break


            # ADJUSTMENT ---------------------------------------------------------------------------------------------------

            match method:
                case 1:
                    new_Csp3 += sign * delta_ener
                    residues_ducarmeType.tr_new = new_Csp3 * residues_ducarmeType.tr_ratio
                    print("--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 > {3:1.2f}".format(i,N, abs(diff_imp), 0.1), end="\r")
                case 2:
                    residues_ducarmeType.tr_new += sign * delta_ener
                    print("--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 > {3:1.2f}".format(i,N, abs(diff_imp), 0.1), end="\r")
                case 3:
                    tr_new = residues_ducarmeType.tr_new
                    abs_tr_new = abs(tr_new)
                    sig = 1.0/(1.0+np.exp(tr_new * np.sign(residues_ducarmeType.tr_ref) - 0.2)**-40)
                    #residues_ducarmeType.tr_new += sign * 0.001 * sig * np.sign(residues_ducarmeType.tr_ref)
                    residues_ducarmeType.tr_new += sign * 0.01 * sig
                    print("--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 > {3:1.2f}".format(i,N, abs(diff_imp), 0.1), end="\r")
                case 4:
                    tr_deviation = residues_ducarmeType.tr_diff / residues_ducarmeType.tr_ref
                    residues_ducarmeType.tr_new = residues_ducarmeType.tr_ref + (float(ind)/N / (1 - tr_deviation))
                    print("--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 > {3:9.3f}".format(i,N, abs(diff_imp), float(i)/N), end="\r")

            
            # --------------------------------------------------------------------------------------------------------------

            # is_bb_and_not_sidechain = residues_ducarmeType.index.isin(backbone_data.index)
            if (fixed_ducarmeTypes_indices.values.any()):
                residues_ducarmeType.tr_new.update(residues_ducarmeType.loc[fixed_ducarmeTypes_indices].tr_ref) # Stay fixed "fixed_ducarmeTypes" in residues_ducarmeType

            residues_ducarmeType.tr_diff = residues_ducarmeType.tr_new - residues_ducarmeType.tr_ref


    print("Best methode for residue {} is {}".format(res3, best_method))

    if (not success):
        print("="*80)
        print("adjusted transfert not found for res ", res3)
        print("="*80)

    eimp_new_all[res3] = eimp_tot
    print("-> Residue atomic types")
    # print(residues_ducarmeType)
    # print("-> Backbone")
    # print(backbone_data)
    # print("-> Sidechain")
    # print(sidechain_data)

    # final_residues_ducarmeType[res3].update(residues_ducarmeType.tr_new)
    final_residues_ducarmeType[res3].update(best_tr_new)
    final_residues_ducarmeType.loc["eimp_new", eimp_new_all.keys()] = ["{0:0.4f}".format(e) if isinstance(e,float) else e for e in eimp_new_all.values()]

    print("-> Final residues ducarmeType")
    print(final_residues_ducarmeType.to_string(col_space=4))

    # final_residues_ducarmeType.to_excel('out_table.xlsx', float_format='%.3f')
    # system("cat out_table.xlsx")