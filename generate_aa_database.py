#!/usr/bin/env python3

from os import listdir, system, remove, popen
from os.path import isfile, join
import time
import sys

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

# https://stackoverflow.com/questions/22029562/python-how-to-make-simple-animated-loading-while-process-is-running
#animation = ["10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"]
animation = ["[■□□□□□□□□□]","[■■□□□□□□□□]", "[■■■□□□□□□□]", "[■■■■□□□□□□]", "[■■■■■□□□□□]", "[■■■■■■□□□□]", "[■■■■■■■□□□]", "[■■■■■■■■□□]", "[■■■■■■■■■□]", "[■■■■■■■■■■]"]

def anim(i, length, text=""):
    anim_step_index = round((i/length)*(len(animation)-1)%(len(animation)-1))
    sys.stdout.write("\r" + animation[anim_step_index] + text)
    sys.stdout.flush()

def read_file_lines(file_path):
    with open(file_path, 'r') as f:
        return f.readlines()
    
def sort_lines_by_resname(res_lines):
    residue_to_database = {res3 : [] for res3 in residueName3To1.keys()}
    res_lines_number = len(res_lines)
    for i, rl in enumerate(res_lines):
        res_name = rl[17:20]
        if (res_name in residue_to_database.keys()):
            residue_to_database[res_name].append(rl)

    return residue_to_database

def save_text_to_file(text, filename):
    text_file = open(filename, "w")
    text_file.write(text)
    text_file.close()

def remove_file(filename):
    remove(filename)

def pdb_tidy(res_filename):
    # Call pdb_tidy and save in tidy_{res_filename}
    system('pdb_tidy.py {0} > tmp'.format(res_filename))
    # Replace the original
    system('mv tmp {0}'.format(res_filename))

def pdb_reres(res_filename):
    # Call pdb_reres and save in resres_{res_filename}
    system('pdb_reres.py {0} > tmp'.format(res_filename))
    # Replace the original
    system('mv tmp {0}'.format(res_filename))


def get_residue_number(res_filename):
    residue_number = popen(" pdb_wc.py -r {0} | grep -Po '\t\K\d+' ".format(res_filename)).read().strip()
    return int(residue_number)

def pdb_selres(res_filename, i, save_path):
    system("pdb_selres.py -{0} {1} > {2}".format(i, res_filename, save_path))

def num_files_in_dir(dir_path):
    num_files = popen(" ls {0} | wc -l | grep -Po '\K\d+' ".format(dir_path)).read().strip()
    return int(num_files)

def clean_residues_database():
    system("rm database_* amino_acids_pdbs/*")

def call_martinize(source_pdb, output_pdb):
    martini_path = "/Users/andrelanrezac/Dev/tools/martini_3.0.b.3.2"
    system("martinize -f {0} -ff {2}/martini303v.partition -x {1}  -elastic".format(source_pdb, output_pdb, martini_path))
    system("cat {0} | grep '^ATOM' > tmp".format(output_pdb))
    system("mv tmp {0}".format(output_pdb))

    pdb_tidy(output_pdb)

def main():

    database_path = "./database"
    amino_acids_path = "./amino_acids_pdbs"
    aa_martini_path = "amino_acids_martini_pdbs"

    database_file_paths = [join(database_path,f) for f in listdir(database_path) if isfile(join(database_path, f))]

    # Extract lines of all pdb in the database and store them in res_lines
    res_lines = []
    for file_path in database_file_paths: 
        res_lines += read_file_lines(file_path)

    # Sort lines by residue and store in residue_to_database dictionary
    resname_to_database = sort_lines_by_resname(res_lines)

    # Loop on all resname and clean pdb database and store in {res_filename}
    # Select residues in {res_filename} and save in amino_acids_pdbs directory
    reset = False
    if (reset) :clean_residues_database()
    for res3 in residueName3To1:
        res_filename = "database_"+res3+".pdb"
        if (reset):
            res_db = ''.join(resname_to_database[res3])
            save_text_to_file(res_db, res_filename)
            pdb_tidy(res_filename)
            pdb_reres(res_filename)
        
            residue_number = get_residue_number(res_filename)

            # print("\nCall pdb_selres...")
            print("{0} residue in {1} database".format(residue_number, res3))
            for i in range(1, residue_number+1):
                anim(i, residue_number, text=res3+"_"+str(i))
                pdb_selres(res_filename, i, join(amino_acids_path, "{0}_{1}.pdb".format(res3, i)))

    num_files_aa_path = num_files_in_dir(amino_acids_path)

    # Get Martini structures from all aa pdb
    for i,file in enumerate(listdir(amino_acids_path)):
        if file.endswith(".pdb"):
            call_martinize(join(amino_acids_path, file), join(aa_martini_path, file))
            
    quit()

    # Add SASA value in each aa pdb
    print("Add SASA : {0} files in {1}".format(num_files_aa_path, amino_acids_path))
    for i,file in enumerate(listdir(amino_acids_path)):
        if file.endswith(".pdb"):
            
            command = "freesasa {0} --shrake-rupley --format=pdb --radii=naccess --hydrogen --no-warnings --resolution=377 |\
                       grep '^ATOM' > {1}".format(join(amino_acids_path, file), join(amino_acids_path, "tmp"))
            system(command)
            system('mv {1} {0}'.format(join(amino_acids_path, file), join(amino_acids_path, "tmp")))

            anim(i, num_files_aa_path)

            break
        



if __name__ == "__main__":
    main()

