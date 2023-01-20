#!/usr/bin/env python3
import pandas as pd
pd.options.display.float_format = '{0:0.4f}'.format
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'
import os
from os import system

import sys
sys.path.insert(1, "/Users/andrelanrezac/Dev/BS_Project/biospring/scripts")
import prepare_martini_v3beta as pm


# André Lanrezac

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

bb_atoms_names = ["CA","C","O","H","HA","N","H2","H3", "BB"] # "BB" is backbone in Martini

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

# --------------------------------------------------------------------------------------------------

# BRASSEUR
def get_ducarmeType(atomType):
    for ducarmeType, atomTypes in ducarmeType_to_atomType.items():
        if atomType in atomTypes:
            return ducarmeType
    return None

ducarmeType_to_Etr = {"Csp3" : -0.105, 
                       "Csp2" : -0.0134,
                       "Hnc"  : -0.0397,
                       "Hc"   : 0.0362,
                       "O"    : 0.0403,
                       "N"    : 0.112,
                       "S"    : -0.108,
                       #"nc"   : 0.000
}

# MARTINI
# Hexadecane/Water
HD_to_Etr = {
    "P2"  : 0.21182233258948446,
    "SP2" : 0.3029719321202046,
    "P1"  : 0.18300296761132331,
    "SP2a": 0.2632068660294278,
    "SN2d": 0.14959239148435105,
    "SQp" : 1.0982732539357418,
    "SP4" : 0.3881827880290122,
    "SQn" : 1.1115282759660008,
    "TC4" : -0.11840246631403953,
    "Qn"  : 0.6686092674933387,
    "P4"  : 0.27954784028816315,
    "TC3" : -0.22579074971514512,
    "TN3d": 0.30839712156214943,
    "TN3a": 0.30839712156214943,
    "SC1" : -0.3067590812717072,
    "SC3" : -0.20261247960538684,
    "C4"  : -0.12824617415281714,
    "TP1" : 0.43230667933265593,
    "SP1" : 0.26510044060517907,
    "SC2" : -0.27456831348393546,
}

# Octanol/Water
OCOS_to_Etr = {
    "P2"  : 0.0374651744716095,
    "SP2" : 0.115508049120828,
    "P1"  : 0.015850650737988636,
    "SP2a": 0.07763655760580243,
    "SN2d": 0.0056807237272538365,
    "SQp" : 0.41279925751377883,
    "SP4" : 0.17420886096911764,
    "SQn" : 0.41279925751377883,
    "TC4" : -0.17071983515047562,
    "Qn"  : 0.21182233258948446,
    "P4"  : 0.0922219679301157,
    "TC3" : -0.22028365825867818,
    "TN3d": 0.09362055475993823,
    "TN3a": 0.09362055475993823,
    "SC1" : -0.26510044060517907,
    "SC3" : -0.19503818130238174,
    "C4"  : -0.1844439358602314,
    "TP1" : 0.16245919796577518,
    "SP1" : 0.07953013218155372,
    "SC2" : -0.24616469484766626,
}


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
# INITIALIZATION
adjustment_representation = "brasseur"


match adjustment_representation:
    case "brasseur":
        type_to_Etr = ducarmeType_to_Etr
        ref_type_ratio_name = "Csp3"
        residues_data_path = "/Users/andrelanrezac/Dev/BS_Project/biospring/data/forcefield/nametotype.txt"
        res_pdbs_dir_path = "/Users/andrelanrezac/Dev/ims_systems/impala/transfer_energies_validation/amino_acids_pdbs"
    
    case "martini":
        type_to_Etr = HD_to_Etr # Hexadecane/Water
        #type_to_Etr = OCOS_to_Etr # Octanol/Water
        ref_type_ratio_name = "P2"
        residues_data_path = "/Users/andrelanrezac/Dev/BS_Project/biospring/data/forcefield/martini_nametotype.txt"
        res_pdbs_dir_path = "/Users/andrelanrezac/Dev/ims_systems/impala/transfer_energies_validation/amino_acids_martini_pdbs"




# Final table
# Initialize the DataFrame with the given data
final_residues_type = pd.DataFrame(columns=["tr_ref"]+list(residueName3To1.keys()), index= list(type_to_Etr.keys()) + ["eimp_ref", "eimp_new", "eimp_ini"])

# Add all experimental Ducarme sidechains transfer energies values (Fig3 Ducarme1998) to eimp_ref row
final_residues_type.loc["eimp_ref", transfer_exp.keys()] = list(transfer_exp.values())

eimp_new_all = dict.fromkeys(residueName3To1.keys(), None)


print(final_residues_type)


# --------------------------------------------------------------------------------------------------

# Configure residues_ducarmeType
residues_type = pd.DataFrame.from_dict(type_to_Etr, orient='index', columns=['tr_ref'])
residues_type = residues_type.rename_axis('particleType')
# Set 'tr_diff' column to 0.0
residues_type = residues_type.assign(tr_diff=0.0)
# Set 'tr_new' column to values in 'tr_ref' column
residues_type = residues_type.assign(tr_new=residues_type['tr_ref'])

# For method 1
ref_type_ratio = residues_type.loc[ref_type_ratio_name, "tr_ref"]
residues_type["tr_ratio"] = residues_type["tr_ref"] / ref_type_ratio

# Update tr_ref column in final_residues_type
final_residues_type["tr_ref"].update(residues_type['tr_ref'])


# --------------------------------------------------------------------------------------------------
# Read nametotype file of all residues

residues_data_df = pd.read_csv(residues_data_path, sep='\t', header=None, names=['particleName', 'particleType'])
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

# Get Ducarme type from each particleType and add it to a new column called ducarmeType
match adjustment_representation:
    case "brasseur":
        residues_data_df['particleType'] = residues_data_df['particleType'].apply(get_ducarmeType)
        residues_data_df.set_index('particleType', inplace=True)
    case "martini":
        residues_data_df.set_index('particleType', inplace=True)
        

# Merge residues_ducarmeType to residues_data_df
residues_data_df = residues_data_df.join(residues_type, how="right")

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
    res_data.set_index('particleType', inplace=True)

    is_bb = res_data.atomName.isin(bb_atoms_names) 
    backbone_data  = res_data[is_bb]
    sidechain_data = res_data[~is_bb]



    methods = [1]
    success = False
    N = 100000

    missing_type = list(set(type_to_Etr.keys()) - (set(res_data.index)))
    missing_type_index = pd.Index(missing_type)
    # Get type that belong to backbone that will not be ajusted
    # Group the data by ducarmeType
    ducarmeType_grouped = is_bb.groupby('particleType')
    # Filter the groups that have all True values
    # Here backbone_index contains all index of backbone type (fixed)
    backbone_index = ducarmeType_grouped.filter(lambda x: x.all()).index

    sidechain_index = list(set(res_data.index) - set(backbone_index))

    if (adjustment_representation=="martini"):
        print(f"Types only in backbone (adjusted): {', '.join(backbone_index.values)}")
        print(f"Missing types (=NaN): {', '.join(missing_type_index.values)}")
        fixed_types_index = missing_type_index
    else:
        # Add missing typesindex in fixed_ducarmeTypes_indices
        print(f"Types only in backbone (fixed): {', '.join(backbone_index.values)}")
        print(f"Missing types (=NaN): {', '.join(missing_type_index.values)}")
        fixed_types_index = backbone_index.union(missing_type_index)

    print(f"Types only in sidechain (adjusted) : {', '.join(sidechain_index)}")

    # Don't show missing types in final table (NaN)
    # Always show sidechain type  bg:white textcolor: black
    # Show backbone type if not in sidechain type as  bg:grey textcolor: black(fixed)
    # If value varies: bold



    # Specific cases for Martini3
    if (adjustment_representation=="martini" and sidechain_data.empty):
        enable_ajustment = True
        # We adjust the backbone
        adjusted_data = backbone_data
    elif (adjustment_representation=="martini" and not sidechain_data.empty):
        enable_ajustment = True
        # Only sidechain energy will be taken into account in the acceptance condition.
        adjusted_data = sidechain_data
    else:
        enable_ajustment = True
        adjusted_data = sidechain_data

    if (enable_ajustment):
        min_method_score = 10000000.0
        best_method = 0
        # Loop on all method
        # for m in range(1,5):
        for m in methods:

            method = m
            # Reinit
            residues_type.tr_new = residues_type.tr_ref
            
            delta_ener = 0.0001
            
            sign = 1
            
            ind = 0
            for i in range(N):

                adjusted_data.update(residues_type)
                
                # Get output value
                adjusted_data.eimp = -adjusted_data.av_sasa * adjusted_data.tr_new * -0.5 + (-0.018 * adjusted_data.av_sasa * -0.5)
                eimp_tot = sum(adjusted_data.eimp) # output value

                if (i==0):
                    previous_imp = eimp_tot
                    print("-> Initial Ducarme sidechain transfer energie: {}".format(eimp_tot))
                    final_residues_type.loc["eimp_ini", res3] = "{0:0.4f}".format(eimp_tot)

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
                    method_score = sum(((residues_type.tr_new - residues_type.tr_ref) / residues_type.tr_ref).abs())

                    print("\nMethod {} score : {}".format(m, method_score))
                    if (m==1):
                        min_method_score = method_score
                        best_tr_new = residues_type.tr_new

                    if (method_score < min_method_score):
                        min_method_score = method_score
                        best_method = method
                        best_tr_new = residues_type.tr_new
                    break


                # ADJUSTMENT ---------------------------------------------------------------------------------------------------

                match method:
                    case 1:
                        ref_type_ratio += sign * delta_ener
                        residues_type.tr_new = ref_type_ratio * residues_type.tr_ratio
                        print("--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 > {3:1.2f}".format(i,N, abs(diff_imp), 0.1), end="\r")
                    case 2:
                        residues_type.tr_new += sign * delta_ener
                        print("--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 > {3:1.2f}".format(i,N, abs(diff_imp), 0.1), end="\r")
                    case 3:
                        tr_new = residues_type.tr_new
                        abs_tr_new = abs(tr_new)
                        sig = 1.0/(1.0+np.exp(tr_new * np.sign(residues_type.tr_ref) - 0.2)**-40)
                        #residues_ducarmeType.tr_new += sign * 0.001 * sig * np.sign(residues_ducarmeType.tr_ref)
                        residues_type.tr_new += sign * 0.01 * sig
                        print("--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 > {3:1.2f}".format(i,N, abs(diff_imp), 0.1), end="\r")
                    case 4:
                        tr_deviation = residues_type.tr_diff / residues_type.tr_ref
                        residues_type.tr_new = residues_type.tr_ref + (float(ind)/N / (1 - tr_deviation))
                        print("--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 > {3:9.3f}".format(i,N, abs(diff_imp), float(i)/N), end="\r")

                
                # --------------------------------------------------------------------------------------------------------------

                if (fixed_types_index.values.any()):
                    residues_type.tr_new.update(residues_type.loc[fixed_types_index].tr_ref) # Stay fixed "fixed_ducarmeTypes" in residues_ducarmeType
                    residues_type.loc[missing_type_index, "tr_new"] = np.NaN

                residues_type.tr_diff = residues_type.tr_new - residues_type.tr_ref


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
    final_residues_type[res3].update(best_tr_new)
    final_residues_type.loc["eimp_new", eimp_new_all.keys()] = ["{0:0.4f}".format(e) if isinstance(e,float) else e for e in eimp_new_all.values()]

    print("-> Final residues ducarmeType")
    print(final_residues_type.to_string(col_space=4))

    # final_residues_ducarmeType.to_excel('out_table.xlsx', float_format='%.3f')
    # system("cat out_table.xlsx")