#!/usr/bin/env python3
import pandas as pd
pd.options.display.float_format = '{0:0.4f}'.format
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'
import os
from os import system


import argparse
import sys
import exp_data_dict as exp

# André Lanrezac


def compute_average_sasa(res3, res_data, res_pdbs_dir_path):

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

def save_text_to_file(text, filename):
    text_file = open(filename, "w")
    text_file.write(text)
    text_file.close()

def get_ducarmeType(atomType):
    for ducarmeType, atomTypes in ducarmeType_to_atomType.items():
        if atomType in atomTypes:
            return ducarmeType
    return None

# ------------------------------------------------------------------------------
# Initialization ---------------------------------------------------------------

residueName3To1 = {"ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C", "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P", "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V" }
residueName1To3 = {"R":"ARG", "N":"ASN", "A":"ALA", "D":"ASP", "C":"CYS", "E":"GLU", "Q":"GLN", "G":"GLY", "H":"HIS", "I":"ILE", "L":"LEU", "K":"LYS", "M":"MET", "F":"PHE", "P":"PRO", "S":"SER", "T":"THR", "W":"TRP", "Y":"TYR", "V":"VAL"}

# "BB" is backbone in Martini representation
bb_atoms_names = ["CA","C","O","H","HA","N","H2","H3", "BB"] 

# Brasseur structure initialization --------------------------------------------
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

# Output
output_path = "../output"
nametotransfert_ajusted_lines = []
particletypetotransfert_ajusted_lines = []
particletypetotransfert_ajusted_lines.append("#RES-ParticleType\tRefTransfer\tAdjustedTransfer\n")



def adjust():

    match adjustment_representation:
        case "brasseur":
            type_to_Etr = exp.ducarmeType_to_Etr
            ref_type_ratio_name = "Csp3"
            residues_data_path = "../input_data/nametotype/brasseur_nametotype.txt"
            res_pdbs_dir_path = "../input_data/amino_acids_pdbs"
            init_ener_values = "brasseur"
        
        case "martini":
            martini_sources = {"HexadecaneBeadSurf": exp.HD_to_Etr, "OctanolBeadSurf": exp.OCOS_to_Etr, "HexadecanePartition": exp.HD_WN, "OctanolPartition": exp.OCOS_WN}
            init_ener_values = "OctanolPartition"
            type_to_Etr = martini_sources[init_ener_values] # Hexadecane/Water
            ref_type_ratio_name = "P2"
            residues_data_path = "../input_data/nametotype/martini_nametotype.txt"
            res_pdbs_dir_path = "../input_data/amino_acids_martini_pdbs"


    # Final table
    # Initialize the DataFrame with the given data
    final_residues_type = pd.DataFrame(columns=["tr_ref"]+list(residueName3To1.keys()), index= list(type_to_Etr.keys()) + ["eimp_ref", "eimp_new", "eimp_ini"])

    # Add all experimental Ducarme sidechains transfer energies values (Fig3 Ducarme1998) to eimp_ref row
    final_residues_type.loc["eimp_ref", exp.res_transfer_exp.keys()] = list(exp.res_transfer_exp.values())

    eimp_new_all = dict.fromkeys(residueName3To1.keys(), None)

    print(final_residues_type)

    # Summary for each aa which atomic type belong to backbone or sidechain
    # Fixed or Adjusted : F or A
    # Backbone or sidechain : BB or SC
    # ->> FBB, ABB, FSC (never happens in practice), ASC
    # Missing atomic type : ---
    atomic_type_status = pd.DataFrame("---", columns=residueName3To1.keys(), index=type_to_Etr.keys())


    # --------------------------------------------------------------------------------------------------

    # Configure residues_type
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
            delta_ener = 0.0001
        case "martini":
            residues_data_df.set_index('particleType', inplace=True)
            delta_ener = 0.0001
            

    # Merge residues_ducarmeType to residues_data_df
    residues_data_df = residues_data_df.join(residues_type, how="right")

    # Group data by res1
    residues_data_df_grouped = residues_data_df.groupby("res1")


    # --------------------------------------------------------------------------------------------------
    # Compute all average sasa

    for res3, res1 in residueName3To1.items():
        print("-"*80)
        print("->", res3)
        res_data = residues_data_df_grouped.get_group(residueName3To1[res3]).reset_index()
        compute_average_sasa(res3, res_data, res_pdbs_dir_path)

        # Clean a bit
        res_data.reset_index(inplace=True)
        res_data.drop(columns=["res1"], inplace=True)
        res_data.set_index('particleType', inplace=True)

        is_bb = res_data.atomName.isin(bb_atoms_names) 
        backbone_data  = res_data[is_bb]
        sidechain_data = res_data[~is_bb]


        success = False
        

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
            atomic_type_status.loc[backbone_index, res3] = "ABB"
        else:
            # Add missing typesindex in fixed_ducarmeTypes_indices
            print(f"Types only in backbone (fixed): {', '.join(backbone_index.values)}")
            print(f"Missing types (=NaN): {', '.join(missing_type_index.values)}")
            fixed_types_index = backbone_index.union(missing_type_index)
            atomic_type_status.loc[backbone_index, res3] = "FBB"

        print(f"Types only in sidechain (adjusted) : {', '.join(sidechain_index)}")

        
        atomic_type_status.loc[sidechain_index, res3] = "ASC"

        print(atomic_type_status)

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

                    diff_imp = eimp_tot - exp.res_transfer_exp[res3] 

                    if (abs(diff_imp) > abs(previous_imp - exp.res_transfer_exp[res3])):
                        # Change sign
                        sign *= -1

                    if (sign > 0): ind += i
                    else:          ind -= i
                    
                    previous_imp = eimp_tot

                    # delta ener adjustment test for martini
                    if (adjustment_representation=="martini"):

                        if(abs(diff_imp) < 1.0):
                            delta_ener = 0.0001
                        else:
                            delta_ener = 0.001
                        
                        #delta_ener = 0.10 if abs(diff_imp) > 10.0 else 0.0001
                    
                    # if (abs(diff_imp) < 0.1):
                    match method:
                        case 1 | 2 | 3:
                            break_condition = abs(diff_imp) < 0.1
                            if(break_condition):
                                print("Break--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 < {3:9.3f}".format(i,N, abs(diff_imp), 0.1), end="\r")
                        case 4:
                            break_condition = abs(diff_imp) < float(i)/N
                            if(break_condition):
                                print("Break--- i:{0:,}/{1:,} {2:9.3f} kj.mol-1 < {3:9.3f}".format(i,N, abs(diff_imp), float(i)/N), end="\r")

                    if (break_condition):
                        #print(sidechain_data)
                        # print("Final diff=",abs(diff_imp))
                        # print("Final Eimp=",eimp_tot)
                        # eimp_new_all[res3] = eimp_tot
                        # print("Eimp_exp=", exp.res_transfer_exp[res3])
                        # print("--- i:{0} {1:9.3f} < {2:9.3f}".format(i, abs(diff_imp), float(i)/N))
                        success=True

                        # Compute method score
                        method_score = sum(((residues_type.tr_new - residues_type.tr_ref) / residues_type.tr_ref).abs())
                        print("\nMethod {} score : {}".format(m, method_score))

                        #residues_type.loc[missing_type_index, "tr_new"] = np.NaN

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

        res_data.tr_new.update(best_tr_new)

        for index, row in res_data.iterrows():
            nametotransfert_ajusted_lines.append("{0}\t{1:.4f}\n".format(res1+row['atomName'], row['tr_new']))

        # save typetotransfer
        for index, row in residues_type.loc[backbone_index.union(sidechain_index)].iterrows():
            particletypetotransfert_ajusted_lines.append("{0}\t{1:.4f}\t{2:.4f}\n".format(res3+"-"+row.name, row.tr_ref,  row.tr_new))


    nametotransfer_name = "nametotransfer_ajusted_"+adjustment_representation+"_"+init_ener_values+".txt"
    save_text_to_file(''.join(nametotransfert_ajusted_lines), os.path.join(output_path, nametotransfer_name))
    print(f"-> Output nametotransfert saved to {os.path.join(output_path, nametotransfer_name)}")

    particletypetotransfer_name = "particletypetotransfer_ajusted_"+adjustment_representation+"_"+init_ener_values+".txt"
    save_text_to_file(''.join(particletypetotransfert_ajusted_lines), os.path.join(output_path, particletypetotransfer_name))
    print(f"-> Output particletypetotransfer saved to {os.path.join(output_path, particletypetotransfer_name)}")


parser = argparse.ArgumentParser(
    description="Adjust transfert energies.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


parser.add_argument('-r', type=str, choices=['brasseur', 'martini'], help="Choose representation mode.", required=True)
parser.add_argument('-m', type=int, choices=[0, 1, 2, 3, 4], help="Choose method, 0 for all.", default=0)
parser.add_argument('--n', type=int, help="Number of iterations.", default=100000)

if __name__ == '__main__':
    args = parser.parse_args(sys.argv[1:])

    if (len(sys.argv[1:]) == 0):
        parser.print_help()
        parser.exit()

    adjustment_representation = args.r

    if (args.m == 0):
        methods = [1, 2, 3, 4]
    else:
        methods = [args.m]

    N = args.n

    adjust()