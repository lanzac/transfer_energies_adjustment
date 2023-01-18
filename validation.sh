#!/usr/local/bin/bash

# https://github.com/openmm/pdbfixer/tree/master/pdbfixer/templates
# https://bmrb.io/referenc/nomenclature/coordinates/aa_normal_20.pdb


PDB_PATH=amino_acids_pdbs

forcefield=/Users/andrelanrezac/Dev/BS_Project/biospring/data/forcefield

declare -A residueName3To1
declare -A residuePDBID


residueName3To1=(
        ["ALA"]="A"
        ["ARG"]="R"
        ["ASN"]="N"
        ["ASP"]="D"
        ["CYS"]="C"
        ["GLU"]="E"
        ["GLN"]="Q"
        ["GLY"]="G"
        ["HIS"]="H"
        ["ILE"]="I"
        ["LEU"]="L"
        ["LYS"]="K"
        ["MET"]="M"
        ["PHE"]="F"
        ["PRO"]="P"
        ["SER"]="S"
        ["THR"]="T"
        ["TRP"]="W"
        ["TYR"]="Y"
        ["VAL"]="V"        )

# https://www.rcsb.org/ligand/ALA
# is present as a standalone ligand in 83 entries <-
# Becareful to have OXT in ligand or have

residuePDBID=(
        ["ALA"]="1XRN"
        ["ARG"]="1BCR"
        ["ASN"]="11AS"
        # ["ASP"]="D"
        # ["CYS"]="C"
        # ["GLU"]="E"
        # ["GLN"]="Q"
        ["GLY"]="1CVI"
        # ["HIS"]="H"
        # ["ILE"]="I"
        # ["LEU"]="L"
        # ["LYS"]="K"
        # ["MET"]="M"
        # ["PHE"]="F"
        # ["PRO"]="P"
        # ["SER"]="S"
        # ["THR"]="T"
        # ["TRP"]="W"
        # ["TYR"]="Y"
        # ["VAL"]="V"        
        )

is_bb_atom() {
    [[ $1 == @(CA|C|O|H|HA|N|H2|H3) ]] && return 1 || return 0
}

# FROM PDB DATABASE  (single amino-acid)------------------------------------------------------------
# rm amino_acids_pdbs/*.pdb
# for res in "${!residuePDBID[@]}"
# do
#     echo
#     echo res:"$res"
#     pdbid=${residuePDBID[$res]}
#     if ! test -f "fetch_pdbs/$pdbid.pdb"; then
#         pdb_fetch.py -biounit $pdbid > fetch_pdbs/$pdbid.pdb
#     fi
#     pdb_splitmodel.py fetch_pdbs/"$pdbid".pdb
#     mv "$pdbid"_1.pdb fetch_pdbs/"$pdbid".pdb
#     rm "$pdbid"_*.pdb
#     pdb_selchain.py -A fetch_pdbs/$pdbid.pdb | pdb_selhetatm.py  | pdb_selresname.py -"$res" | pdb_reres.py | pdb_selres.py -1 | pdb_reatom.py | sed 's/HETATM/ATOM  /g' > amino_acids_pdbs/$res.pdb

#     if ! test -f "amino_acids_pdbs/"$res"_h.pdb"; then
#         pdb2pqr30 --ff AMBER --keep-chain amino_acids_pdbs/$res.pdb amino_acids_pdbs/"$res"_h.pdb
#         rm amino_acids_pdbs/"$res"_h.log

#         addResInfo.py amino_acids_pdbs/"$res"_h.pdb tmp.pdb
#         pdb_element.py tmp.pdb > amino_acids_pdbs/"$res"_h.pdb
#     fi 

#     echo "Final pdb for res $res"
#     cat amino_acids_pdbs/"$res".pdb
#     echo "Add hydrogens to $res"
#     cat amino_acids_pdbs/"$res"_h.pdb
# done

#FROM ALL AA PDB SEQUENCE (single amino-acid)------------------------------------------------------
#Get all amino acids from the sequence aa20_h.pdb
# rm amino_acids_pdbs/*.pdb
# # pdb_delelem.py -H aa20.pdb > aa20_noH.pdb
# # pdb2pqr30 --ff AMBER --keep-chain aa20_noH.pdb aa20_H.pdb
# for res in "${!residueName3To1[@]}"
# do
#     pdb_selresname.py -$res aa20_h.pdb | pdb_reatom.py > amino_acids_pdbs/${res}_h.pdb
#     # if no hydrogens
#     # pdb_element.py amino_acids_pdbs/${res}_h.pdb | pdb_delelem.py -H > amino_acids_pdbs/${res}.pdb
#     # rm amino_acids_pdbs/${res}_h.pdb
# done

# FROM LOCAL DATABASE (multiple occ of amino-acids)-------------------------------------------------
# rm amino_acids_pdbs/*.pdb

for res in "${!residueName3To1[@]}" # <-tmp
do # <-tmp

#res="ARG"
grep -h $res database/*.pdb > tmp.pdb
pdb_tidy.py tmp.pdb | pdb_reres.py > database_$res.pdb
rm tmp.pdb

i=1
extract=true
while $extract
do
    #if ! test -f "amino_acids_pdbs/"$res"_$i.pdb"; then continue; fi
    result=$(pdb_selres.py -$i database_${res}.pdb)
    if [[ $result == END* ]]; then break; fi
    #echo "$result" > amino_acids_pdbs/${res}_$i.pdb
    echo "$result" | 
       martinize -f ${sys_path}/${PDBID}_clean.pdb -ff $martini/martini303v.partition -x ${sys_path}/${PDBID}_pdb2spn_input.pdb  -elastic
       freesasa $pdb  --shrake-rupley --format=pdb --radii=naccess --hydrogen --no-warnings --resolution=377 | 
       grep "^ATOM" > amino_acids_pdbs/${res}_${i}.pdb
       
    ((i++))
    echo $i
done

done # <-tmp

> alltransfert.txt

# for pdb in $PDB_PATH/*.pdb
# do
#     #sasa=$(freesasa $pdb  --shrake-rupley --format=pdb --radii=naccess --hydrogen --no-warnings --resolution=377 | grep "^ATOM")
#     #echo "$sasa"
#     #echo
#     sasa=$(cat $pdb)
#     eint=0.0
#     elip=0.0
#     eimp=0.0
#     tottransfert=0.0
#     while IFS= read -r p; do
#         sasa=$(echo "${p:60:6}" | sed 's/ //g') #substring and remove space
#         res3="${p:17:3}"
#         atomname=$(echo ${p:12:4} | sed 's/ //g') #substring and remove space

#         # keep only non backbone atoms
#         is_bb_atom $atomname 
#         [[ $? == 1 ]] && continue


#         res1="${residueName3To1[$res3]}"
#         particlename=$res1$atomname
        
#         transfert=$(grep "^$particlename	" $forcefield/amber.ff | awk '{print $6}'| sed 's/\r//g')
#         if [ -z "$transfert" ]; then 
#         echo ^$particlename "not found in amber.ff try with ^$atomname"
#         transfert=$(grep ^$atomname $forcefield/amber.ff | awk '{print $6}'| sed 's/\r//g')
#         fi
#         if [ -z "$transfert" ]; then 
#             echo ^$atomname "not found in amber.ff."
#             continue
#         fi

#         # eint=$(echo "$eint - $sasa * $transfert * 0.5"|bc -l)
#         # elip=$(echo "$elip - 0.018 * $sasa * 0.5"|bc -l)
#         eint="$eint - $sasa * $transfert * -0.5" # -0.5 in hydrophobic code
#         elip="$elip - 0.018 * $sasa * -0.5"
#         echo $particlename $sasa Å2 $transfert kJ.mol-1.Å-2
        
#     done <<< "$sasa"
#     #eimp=$(echo "$eint+$elip"|bc -l)
#     #tottransfert=$(echo "$eimp*-2"|bc -l)

#     # eint=$( bc -l <<< "$eint" )
#     # elip=$( bc -l <<< "$elip" )
#     # eimp=$( bc -l <<< "$eint+$elip" )
#     # tottransfert=$(echo "$eimp*-2"|bc -l)
#     echo $res3 of pdb $pdb
#     echo eint: $( bc -l <<< "$eint" )
#     echo elip: $( bc -l <<< "$elip" )
#     echo eimp:  $( bc -l <<< "$eint+$elip" )
#     #tottransfert=$(echo "($eint+$elip)*-2"|bc -l)
#     tottransfert=$(echo "($eint+$elip)"|bc -l)
#     echo tottransfert: "$tottransfert"
#     echo "$res3 $tottransfert" >> alltransfert.txt
    
# done

# cat alltransfert.txt