```sh
cd script
```

## Generate PDBs input data
> **_NOTE:_** This script takes some time. Source PDB are found in general_database

Help
```sh
./generate_aa_database.py -h
```
---
Clean aa_database, amino_acids_martini_pdbs & amino_acids_pdbs
```sh
./generate_aa_database.py --clean 
```
---
1st step:
- Generate PDBs in aa_database & amino_acids_martini_pdbs\
- Add SASA to PDBs in amino_acids_martini_pdbs
```sh
./generate_aa_database.py
```

2nd step (Martini):
```sh
./generate_aa_database.py --martini
```
---

## Adjust transfert energies

Help
```sh
./adjust_transfert.py -h
```

Adjustment methods: option -m \
0 : try all methods and keep the best result for each amino acid \
1, 2, 3 or 4: test only one method

### Adjust transfer energies from Brasseur values
```sh
./adjust_transfert.py -r brasseur -m 0
```

### Adjust transfer energies from Martini partition 
- from Hexadecane partition / theorical bead surface : kj.mol-1.Å-2
```sh
./adjust_transfert.py -r martini --oil_water hexadecane
```

- from Octanol partition / theorical bead surface : kj.mol-1.Å-2
```sh
./adjust_transfert.py -r martini --oil_water octanol
```


> **_NOTE:_** It may be necessary to widen the terminal window to better view the printed 
tables.


### Generate ff file for pdb2spn (BioSpring tool to build spring network)
- for adjusted Brasseur values
```sh
./buildff_allatoms.py \
    -f ../output/brasseur/nametotransfer_adjusted_brasseur_brasseur.txt \
    -o ../output/brasseur/brasseur_adjusted.ff
```


- for non-adjusted Hexadecane partition (kj.mol-1.Å-2) \
-> output: martini_raw_hexadecane.ff
```sh
./buildff_martini.py -ow hexadecane
```

- for non-adjusted Octanol partition (kj.mol-1.Å-2) \
-> output: martini_raw_octanol.ff
```sh
./buildff_martini.py -ow octanol 
```

- for adjusted Hexadecane partition (kj.mol-1.Å-2)
```sh
./buildff_martini.py \
    --config_file ../output/martini/nametotransfer_adjusted_martini_hexadecane.txt \
    -o ../output/martini/martini_adjusted_hexadecane.ff
```

- for adjusted Octanol partition (kj.mol-1.Å-2) 
```sh
./buildff_martini.py \
    --config_file ../output/martini/nametotransfer_adjusted_martini_octanol.txt \
    -o ../output/martini/martini_adjusted_octanol.ff
```