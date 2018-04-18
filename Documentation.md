# Documentation for toyLIFE programs

## Contact

If, after reading this document, you still don't know how to make this work, or there is something you don't understand, please contact me. You can find me at `pablocatalanfdez [at] gmail.com`.

## Lookup Tables

Every toyLIFE program must start with an initialization of a ToyPlugin object, using the default constructor. This creates an instance of the class ToyPlugin that contains the relevant lookup tables and allows us to compute phenotypes from genotypes. You can check what's stoed in the lookup tables by reading the file `Lookup_Files.md`.

## Genotype Format

Genotypes in toyLIFE are binary strings, stored in the program's memory as a C++ `std::string`. String sizes must be multiples of 20, because toyLIFE genotypes are formed by a group of genes of length 20. That is why in `example.cpp` the relevant variable needed when creating a new genotype is `gene_number`, the number of genes.

## Creating and Mutating a Genotype

`std::string random_genotype(int gene_number, std::uniform_real_distribution<double>& RNG, std::default_random_engine& generator);`
  
`std::string mutation(const std::string& genotype, int pos);`


The function `random_genotype` takes as an argument the number of genes we want for our genotype, and the random number generator from C++. Mutating a genotype is as simple as calling function `mutation(genotype, pos)`, where `pos` is the position in the string that we want to mutate. 

## Regulatory Phenotype

`std::vector<int> ToyPlugin::regulatory_phenotype(const std::string& genotype);`
  
The first relevant phenotype we may want to look at in toyLIFE is the regulatory one. Given a collection of genes (a genotype), toyLIFE's rules give rise to a Boolean regulatory logic, which assigns to each regulatory state an output. In a Boolean model, genes can be ON or OFF, and thus there are 2^n possible regulatory states. For instance, the state (1,0,0) in a 3-gene genotype means that gene 1 is expressed, but gene 2 and gene 3 are not. We can write each regulatory state as an `int` if we see them as binary numbers: state (1,0,0) is number 4.

So the output of function `regulatory_phenotype` is a `std::vector<int>` containing at position *i* the regulatory output of state *i*. If this doesn't make sense, go read `toyLIFE.pdf` or any text on Boolean networks. If after that nothing makes sense, you can write me an email!

## Metabolic Phenotype

`std::vector<int> ToyPlugin::metabolic_phenotype(const std::string& genotype, const std::vector<std::string>& env);` 
 
The other relevant phenotype in toyLIFE is the metabolic one. There are 274 different metabolites of up to length 8, and 60 of them can't be catabolized (see `toyLIFE.pdf`, section 'A note on toyMetabolites'). The default setting is that we want to check which of all 214 breakable metabolites a genotype is able to catabolize. Therefore, the ToyPlugin class already initializes the vector `env` with all metabolites. But you can decide to study just a subset of all metabolites: if that is the case, just create a vector of `std::string` containing the metabolites (up to size 8) that you want to study. The function `metabolic_phenotype(genotype, env)` returns a `std::vector<int>` containing a 1 in position *i* if the metabolite in position *i* in your `env` vector is catabolized by this genotype, and 0 otherwise. 

## Robustness

`double ToyPlugin::metabolic_robustness(const std::string& genotype, const std::vector<int>& phen, const std::vector<std::string>& env);`
	      
`double ToyPlugin::regulatory_robustness(const std::string& genotype, const std::vector<int>& logic_vector);`

An important property of genotype-phenotype maps is robustness, or the fraction of neutral neighbors a sequence has. The functions `regulatory_robustness`and `metabolic_robustness` compute this number for a given genotype. Both functions take as a second argument the phenotype of our current genotype, so that you don't have to compute it again, and the metabolic case also takes as an argument the vector `env` defining the metabolic environment.

## Evolvability
		
`std::map<std::vector<int>, int> ToyPlugin::metabolic_evolvability(const std::string& genotype, const std::vector<int>& phen, const std::vector<std::string>& env);`

`std::map<std::vector<int>, int> ToyPlugin::regulatory_evolvability(const std::string& genotype, const std::vector<int>& logic_vector);`


In the same way, we can compute the evolvability of a sequence, that is, how many different phenotypes are in the neighborhood of a sequence. The output here is a `std::map<std::vector<int>, int>`. The key of this map is the phenotype of a neighboring sequence, and its value is the number of times it appears in the neighborhood. 
