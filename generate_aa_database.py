#!/usr/bin/env python3

from os import listdir, system, remove, popen, path
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

def clean_database(path):
    system("rm {0}/*".format(path))

def call_martinize(source_pdb, output_pdb):
    martini_path = "/Users/andrelanrezac/Dev/tools/martini_3.0.b.3.2"
    output_itp = path.splitext(output_pdb)[0]+".itp"
    system("martinize -f {0} -ff {1}/martini303v.partition -x tmp -o {2} &>/dev/null".format(source_pdb, martini_path, output_itp))

    # Keep only ATOM line in output_pdb
    print(output_pdb)
    system("cat tmp | grep '^ATOM' > tmp2 ; rm tmp")

    # Get Protein itp file
    prt_itp_file = popen(" grep -Po 'Protein_.\.itp' {0} ".format(output_itp)).read().strip()
    # Get atom lines of itp file
    itp_atom_lines = popen(" awk '/\[ atoms \]/{flag=1;next} flag && !/^$/ {print;next} /^$/ {flag=0}' " + prt_itp_file).read().rstrip().splitlines(keepends=True)

    # Remove useless output_itp
    remove_file(output_itp)

    # Read output_pdb file
    output_pdb_lines = read_file_lines("tmp2")

    new_pdb_lines = []
    for i,line in enumerate(itp_atom_lines):
        bead_type= line[7:11]
        # we assume that pdb has same number of line than itp atom lines
        # Replace PDB atom name by bead type
        new_pdb_lines.append(output_pdb_lines[i][:12]+bead_type+output_pdb_lines[i][16:]) 


    save_text_to_file(''.join(new_pdb_lines), output_pdb)
    remove_file("tmp2")


def main():

    general_database_path = "./general_database"
    aa_database_path = "./aa_database"
    amino_acids_path = "./amino_acids_pdbs"
    aa_martini_path = "amino_acids_martini_pdbs"

    database_file_paths = [join(general_database_path,f) for f in listdir(general_database_path) if isfile(join(general_database_path, f))]

    # Extract lines of all pdb in the database and store them in res_lines
    res_lines = []
    for file_path in database_file_paths: 
        res_lines += read_file_lines(file_path)

    # Sort lines by residue and store in residue_to_database dictionary
    resname_to_database = sort_lines_by_resname(res_lines)

    # GENERATE AA PDBS ---------------------------------------------------------
    # Loop on all resname and clean pdb database and store in {res_filename}
    # Select residues in {res_filename} and save in amino_acids_pdbs directory
    reset = False
    if (reset) :
        clean_database(aa_database_path)
        clean_database(amino_acids_path)
        for res3 in residueName3To1:
            aa_db_path = join(aa_database_path, "database_"+res3+".pdb")
            if (reset):
                res_db = ''.join(resname_to_database[res3])
                save_text_to_file(res_db, aa_db_path)
                pdb_tidy(aa_db_path)
                pdb_reres(aa_db_path)
            
                residue_number = get_residue_number(aa_db_path)

                # print("\nCall pdb_selres...")
                print("{0} residue in {1} database".format(residue_number, res3))
                for i in range(1, residue_number+1):
                    anim(i, residue_number, text=res3+"_"+str(i))
                    pdb_selres(aa_db_path, i, join(amino_acids_path, "{0}_{1}.pdb".format(res3, i)))

    
    num_files_aa_path = num_files_in_dir(amino_acids_path)

    # GENERATE MARTINI PDBS ----------------------------------------------------
    reset = False
    if (reset) : 
        clean_database(aa_martini_path)
        # Get Martini structures from all aa pdb
        for i,file in enumerate(listdir(amino_acids_path)):
            if file.endswith(".pdb"):
                call_martinize(join(amino_acids_path, file), join(aa_martini_path, file))
                print(file)
                anim(i, num_files_aa_path, file)
        
        # Remove all Protein_x.itp
        system("rm Protein_*.itp")

    # --------------------------------------------------------------------------
    enable_all_atom_pdbs = True
    enable_martini_pdbs = False
    # --------------------------------------------------------------------------

    # ADD SASA TO AA PDBS ------------------------------------------------------
    freesasa_path = "/Users/andrelanrezac/Dev/BS_Project/freesasa/share"
    if (enable_all_atom_pdbs):
        # Add SASA value in each aa pdb
        print("Add SASA : {0} files in {1}".format(num_files_aa_path, amino_acids_path))
        for i,file in enumerate(listdir(amino_acids_path)):
            if file.endswith(".pdb"):
                
                command = "freesasa {0} --shrake-rupley --format=pdb --radii=naccess --hydrogen --no-warnings --resolution=377 |\
                        grep '^ATOM' > {1}".format(join(amino_acids_path, file), join(amino_acids_path, "tmp"))
                system(command)
                system('mv {1} {0}'.format(join(amino_acids_path, file), join(amino_acids_path, "tmp")))

                anim(i, num_files_aa_path)

    # ADD SASA TO MARTINI PDBS -------------------------------------------------
    if(enable_martini_pdbs):
        # Add SASA value in each aa pdb
        print("Add SASA : {0} files in {1}".format(num_files_aa_path, aa_martini_path))
        for i,file in enumerate(listdir(aa_martini_path)):
            if file.endswith(".pdb"):
                
                command = "freesasa {0} --shrake-rupley --format=pdb --config-file={1} --hydrogen --no-warnings --resolution=377 |\
                        grep '^ATOM' > {2}".format(join(aa_martini_path, file), join(freesasa_path, "martini3.config"), join(aa_martini_path, "tmp"))
                system(command)
                system('mv {1} {0}'.format(join(aa_martini_path, file), join(aa_martini_path, "tmp")))

                anim(i, num_files_aa_path)

        



if __name__ == "__main__":
    main()

