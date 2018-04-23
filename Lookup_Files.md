# What's in the lookup files

## `data/protein_gene.txt`

There are 2^20=65536 genes in toyLIFE. Each of them has a binary sequence that, after applying to it the rules from the HP model, folds (or not) into a 4x4 lattice. Each line of this file corresponds to a gene (line 1 corresponds to gene 0, line 2 to gene 1, and so on). The number on each line is the identification number of the resulting protein from that gene's sequence. Id number 0 is given to degenerate structures (the protein does not fold).

For instance, gene 14 (in line 15) folds into protein 16.

## `data/proteins.txt`

In toyLIFE, proteins are identified by their folding energy and their perimeter. There are 2710 proteins. The first entry in each line is the protein identification number. The second entry records the folding energy (note that protein 0, the degenerate structure, is given a symbolic enregy of 10.0). The next four entries are the sides of the perimeter of each protein, and the next four are the sides of the reflected protein, which is also functional. Each integer corresponds to a four-aminoacid side: if one of the faces of the perimeter is 14, this means that the sequence of that side is 1110.

For instance, protein 117 has a folding energy of -6, and perimeter (0,0,1,14) in its direct form, and (0,0,7,8) in its reverse form:

Direct:

0 0 0 0

0 ()() 0

0 ()() 0

0 1 1 1

Reverse:

0 0 0 0

0 ()() 1

0 ()() 1

0 0 0 1

Note that sides are read clock-wise. Many proteins share the same perimeter, differing only in their folding energy.

## `data/prot_abundance.txt`

Contains the number of genes that code for each protein.

For instance, protein 15 is coded by two different genes.

## `data/genes_per_prot.txt`

Contains the identification numbers of the genes coding for each protein.

For instance, protein 16 is coded by genes 7 and 57344.

## `data/perimeters.txt`

Each perimeter is assigned an identification number as well.

For instance, number 6 corresponds with perimeter (0,0,1,14) in its direct form, and (0,0,7,8) in its reverse form:

## `data/proteins_with_perimeter.txt`

Contains the perimeter of each protein.

For instance, protein 117 has perimeter 6.

## `data/dimers.txt`

Dimers are formed by the union of two proteins. Because this number is so huge, we only record the different dimers obtained by binding of different perimeters, and record the binding energy of such an interaction.

The first two entries correspond to the perimeters involved in the dimer. The third entry is 0 if the order in which these perimeters are written is the order in which we shuld look at the direct folding of the dimer, and 1 otherwise. The next 6 entries correspond to the perimeter in the direct form, and the next 6 to the perimeter of the inverse folding.

For instance, perimeters 24 and 25 give rise to the following dimer:

Direct:

0 0 0 1 1 0 0 1

0 ()()()()()()() 0

0 ((25))((24)) 0

0 1 0 0 0 0 0 0

Inverse:

0 1 0 0 0 0 0 0

0 ()()()()()()() 0 

0 ((24))((25)) 0

0 0 0 1 1 0 0 1

Because the third entry is 1, we know that in the direct folding, the protein to the left is the one with perimeter 25.

## `data/metabolites.txt`

Contains the sequence of all 274 metabolites up to length 8 (after taking into account symmetries) and a 1 if they can't be catabolized.

## `data/neighbors_coding.txt` and `data/neighbors_prom.txt`

Contains the neighbors of each promoter and coding sequence.

For instance, promoter 2 (in line 3) has as neighbors promoters 0, 3, 6 and 10, obtained mutating each position in the sequence:

Promoter 2:  0010

Promoter 3:  0011

Promoter 0:  0000

Promoter 6:  0110

Promoter 10: 1010

## `promoters` folder

For each file called `dimers_bind_prom*.txt`, the file contains the identification number of every dimer that is able to bind that particular promoter, its binding energy, and a 1 if it is able to activate the polymerase.

For instance, dimer 9244 is able to bind promoter 4 with energy -2.9, and when it binds this promoter, the gene will be expressed.

The same for the `proteins_bind_prom*.txt` files, but with proteins instead of dimers.

## `dim_metabolism` folder

Each file refers to a particular metabolite. Thus, '4_7.txt' refers to a metabolite of length 4 whose binary sequence translates to 7 in decimal base. All metabolites from lengths 4 to 8 are included. In each file, the first column contains the id of a dimer that is able to bind this particular metabolite. The next number is the energy of the binding between the dimer and the metabolite. The next string can represent two things: (1) if the way the dimer has bound the metabolite makes it impossible for the metabolite to be broken, then the metabolite binary sequence is shown; (2) if the metabolite can be broken, then the string represents the binary sequence that proteins in toyLIFE will 'perceive', thus defining which interactions are possible. The fourth number represents the position in which the metabolite will be broken if catabolism takes place. And the fifth number represents the orientation of the dimer when binding the metabolite: 1 is direct form and 2 is reverse form. 

## `prot_metabolism`folder

This folder is quite similar to `dim_metabolism`. There is a file for each metabolite, but instead of dimers, here we record which proteins are able to bind the corresponding metabolite. The first number is the protein id, the second one is the energy of the binding, the third number represents the orientation of the protein when binding the metabolite (1 if direct, 2 if reverse), the fourth number represents the side of the protein that binds the metabolite (from 0 to 3), and the last number represents where in the metabolite sequence the protein is binding:

   0 1 0 1 0 1 0 1 -> metabolite of size 8
   
   0 1 0 1 -> if protein binds here, last column will show a 0
   
   ()1 0 1 0 -> last column shows a 1
   
   ()()0 1 0 1 -> last column shows a 2
   
   ()()()1 0 1 0 -> last column shows a 3     
   
   ()()()()0 1 0 1 -> last column shows a 4
           
