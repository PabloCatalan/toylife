#include <string>
#include <vector>
#include <map>

#ifndef _PATTERNS_PLUGIN
#define _PATTERNS_PLUGIN


std::vector<int> morphogen_creator(int cells, int t, const std::string& gradient);

class PatternsPlugin{
 public:
  PatternsPlugin();
  ~PatternsPlugin();
  
  //LOOKUP TABLES
  std::vector<int> prot_gen;
  std::vector<int> prot_perimeters;
  std::vector<double> prot_energies;
  std::vector<std::vector<double> > prot_prom_energies;
  std::vector<std::vector<bool> > prot_pol;
  std::map<std::pair<int,int>, int> dim_perim;
  std::vector<double> dim_bond_energy;
  std::vector<std::vector<double> > dim_prom_energies;
  std::vector<std::vector<bool> > dim_pol;
  std::vector<double> polymerase;
  std::vector<std::vector<int> > neighbors_prom;
  std::vector<std::vector<int> > neighbors_coding;
  
  //FUNCTIONS TO INITIALIZE toyLIFE
  void init_gen(std::vector<int>& prot_gen);
  void init_prot(std::vector<int>& perimeters, std::vector<double>& energies, std::vector<std::vector<double> >& prot_prom_energies, std::vector<std::vector<bool> >& prot_pol);
  void init_dim(std::map<std::pair<int, int>, int>& dim_perim, std::vector<double>& dim_energy, std::vector<std::vector<double> >& dim_prom_energies, std::vector<std::vector<bool> >& dim_pol);
  void init_polymerase(std::vector<double>& polymerase);
  void neighbors(std::vector<std::vector<int> >& neighbors_prom, std::vector<std::vector<int> >& neighbors_coding);
  
  //FUNCTIONS TO COMPUTE PHENOTYPES WITH THREE GENES
  int dimer_from_two(int prot1, int prot2);
  int find_cycle1(const std::vector<int>& v, int i1);
  int find_cycle2(const std::vector<std::vector<int> >& p);
  std::pair<int,int> periodicity(const std::vector<std::vector<int> >& p);
  std::vector<int> logic_proms(int i, int p1, int p2);
  std::map<std::vector<int>, std::vector<int> > coarsening(int p1, int p2);
  std::vector<int> translation(const std::map<std::vector<int>, std::vector<int> >::const_iterator& p1, const std::map<std::vector<int>, std::vector<int> >::const_iterator& p2);
  std::vector<int> logic_function(int p1, int p2, int pr1, int pr2);
  void genotype_str_to_vec(const std::string& genotype, int& pr1, int& pr2, int& p1, int& p2);
  std::string genotype_vec_to_str(int pr1, int pr2, int c1, int c2);
  std::vector<int> str_logic_function(const std::string& genotype);
};

#endif
