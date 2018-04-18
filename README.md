toyLIFE - Pablo Catalan

/******************************************************************************
 *   Copyright (C) 2018  Pablo Catalan                                        *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 ******************************************************************************/

# toyLIFE

toyLIFE is a computational model inspired by cellular biology. It was created to study the genotype-phenotype map at cellular level from first principles. If you are unfamiliar with the model, please read the file "toyLIFE. pdf" in which all details are explained.

# Compiling

toyLIFE is written in C++, so you will only need a C compiler and cmake in order to start. An example file, "example.cpp" is supplied to show how to code with this library. If you want to create your own programs, I would suggest modifying "example.cpp" and running everything within the "toylife" folder.

# Description of the library and folders

toyLIFE contains analogs of genes, proteins and metabolites that interact with each other following a small set of rules. These interactions are already computed and saved in lookup files that will be stored in memory during program execution (in the Plugin class object Model, see "example.cpp"). The lookup files can be found in the following folders:

## data

- 

toyLIFE is a collection



I should comment the following programs:
-proteins.cpp
-dimers.cpp

Parent folder: toylife
Descendants: data, promoters, dim_metabolism, prot_metabolism, results

////////////////////////////////////////////////////////////////////////
//                        BEFORE YOU START                            //
////////////////////////////////////////////////////////////////////////

All programs will make use of functions defined in "toy_functions.h".
Make sure you have it in your parent folder.

First, make sure you have read 
Arias, C. F., Catalán, P., Manrubia, S., & Cuesta, J. A. (2014). 
toyLIFE: a computational framework to study the multi-level 
organisation of the genotype-phenotype map. Scientific Reports, 4, 7549.

(Available here: http://www.nature.com/articles/srep07549).

We describe the model there, and this README file assumes you are
familiar with toyLIFE.

////////////////////////////////////////////////////////////////////////
//                      HOW TO CREATE toyLIFE                         //
////////////////////////////////////////////////////////////////////////

A) DATA
You start with the 38 possible SAWs (kept in "shapes.txt").
The interaction rules are Epp=0, Ehp=-0.3, Ehh=-2.
The energy of the toyPolymerase is -11, and its sequence is PHPH
(when computing if anything binds to the polymerase, we will read
it as if it was the inferior side of any protein, and thus the 
side is 10).

A.1) toyPROTEINS
We will use the program "proteins.cpp". This will read the file
"shapes.txt" and assign a shape to each sequence of the 2^16 possible
ones. Once we have the shape, the program will assign a perimeter and
an energy to each sequence.
It will make use of some functions defined in "toy_functions.h", so
make sure you have that file in the parent folder. 

The output of this program are several files:
a) "foldings.txt", which will save the folding associated to each
    sequence.
b) "folding_energies.txt", which will save the folding energy of
    each sequence in the order that they appear in "foldings.txt".
c) "proteins.txt", which will save the energy of each toyProtein,
    and its perimeter (normal and reversed). The first toyProtein
    (line 1) is the degenerate one, and we assign energy +10 to it.
    Each toyProtein will be labeled by the line it appears on this
    file, starting with 0. So when looking up a protein in the file,
    you will have to search for line (id+1).
d) "proteins_with_perimeters.txt", the same as above, but labeling
    each perimeter and assigning it to a protein.
e) "genes_per_protein.txt", which will save which genes code for
    a given protein. Line i states how many genes code for prot i.
    Genes that do not code for a protein (that "code" for prot 0)
    are the ones that do not appear in this file.
f) "protein_gene.txt", which will save which protein (id) corresponds
    to a given gene, whose Id is given by (Line+1).

A.2) toyDIMERS
Now we generate the dimers. Use program "dimers.cpp". This will take 
the proteins we have just created, and generate the dimer that will
live in toyLIFE. As an output, we have the file "dimers.txt", very 
similar to "proteins.txt" but with the information about dimers. Note
that dimers have six sides. Sides 0,1 and 5 correspond to the protein
whose id is smaller.

B) PROMOTERS
Now we generate the promoter binding information. The idea is to tabulate
all info in files, so that programs can run faster. We use 
"prot_promoters.cpp" and "dim_promoters.cpp" for this. Once they're done
they will generate 16 files each, telling which prots (dimers) bind to
each promoter. The file will contain the prot (dimer) id, the binding
energy and a 0 (1) if it inhibits (activates) the expression of that gene.

C) METABOLISM
Lastly, we generate the information concerning the binding to metabolites.
This is very similar to what happened with promoters, but with different
sizes. As output we have files named "8_123.txt" where 8 is the size of
the metabolite and 123 is its sequence in decimal base. 
In dimers, these files will be:
Dimer id   binding_energy   sequence	position	binding
->sequence is the sequence that the dimer leaves to be broken
->position is the number of sugras that are to the left counting from
  the center of the dimer
->binding is the binding position of the dimer:
  1) ******** -> met
     ++++++++ -> dimer (1 is left subunit aka, the protein with the 
     +1 ++2 +    smallest id)
     +  ++  +
     ++++++++
     
  2) ++++++++ 
     +1 ++2 + 
     +  ++  +
     ++++++++
     ********

  3) ******** 
     ++++++++ 
     +2 ++1 +   
     +  ++  +
     ++++++++

  4) ++++++++ 
     +2 ++1 + 
     +  ++  +
     ++++++++
     ********
 Of course, when it comes to breaking the dimer, (1) and (3) are the same
 (the same goes for (2) and (4)): they leave the same sequence open.

In proteins, it will be:
Prot id	     binding_energy	orientation	side	position
->orientation: 1 if it's normal, 2 if it's reversed
->side: from 0 to 3, which side of the protein binds to the met
->position: from 0 to met.size-3, which binding position:
  // 0 1 0 1 0 1 0 1 -> metabolite of size 8
  // 0 1 0 1 -> pos 0
  //   1 0 1 0 -> pos 1
  //     0 1 0 1 -> pos 2
  //       1 0 1 0 -> pos 3
  //         0 1 0 1 -> pos 4