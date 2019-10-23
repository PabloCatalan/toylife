# README

## ABOUT THIS FOLDER

The "patterns" folder contains all the code needed to obtain the results in the following paper:

CATALÁN, P, MANRUBIA, S. AND CUESTA, J.A. (2019) Populations of genetic circuits are unable to find the fittest
  solution in a multilevel genotype-phenotype map. In preparation.

## HOW TO RUN

First, run 'make all' to compile all executables

(1) Obtain abundances of Boolean GRNs in toyLIFE: run 'grn_abundance', results will appear in 'results/grn_abundance.txt' and other files (see 'grn_abundance.cpp' for more details). This program computes the corresponding GRN for each of the 2^40 two-gene genotypes in toyLIFE. Instead of doing it by brute-force, I wrote a specific code that made use of the unique properties of two-gene genotypes, thus saving computation time. The details are in the 'patternsplugin.cpp' file.

(2) Compute all possible patterns: run 'gpmap_morphogen_choose_diffusion 31 100 0 point 1', where 31 is tissue size, 100 is pattern length and 0 specifies aperiodic boundary conditions; 'point' is one possible option for the initial morphogen signal (for more options, see the function 'morphogen_creator' in file 'patternsplugin.cpp'), and '1' specifies the diffusion scenario (in this case, 'only_P1': for more options, see file 'gpmap_morphogen_choose_diffusion.cpp'). See all output files at the end of the C++ file.

(3) Run evolutionary simulations with pattern in Figure 3b as target: run 'simulations_patternFig3'. Output files will print endpoint abundances of GRNs that generate the pattern ('results/simulations_patternFig3_theory.txt') and their abundances in genotype space ('results/simulations_patternFig3_theory.txt'). These results are then plotted in Supp Fig S12.

(4) Run evolutionary simulations with pattern 113 as target: run 'simulations_pattern113'. Output files will print endpoint abundances of GRNs that generate the pattern ('results/simulations_pattern113.txt').

(5) Run evolutionary simulations with every pattern as target: run 'simulations_prob_of_arrival_patterns'. The output will hold one row per pattern, with the following data:

Pattern id + Prob of Arrival (Avg) + Prob of Arrival (Std) + Time of Arrival (Avg) + Time of Arrival (Std) + Pattern abundance

(6) Compute algorithmic complexity of phenotypes (Supp Fig S11): run 'phenotype_complexity'.

(7) Plot Figures in the manuscript: run 'python paper_figures.py'

## CONTACT

If you have any doubts, you can find me at pcatalan [at] math.uc3m.es
