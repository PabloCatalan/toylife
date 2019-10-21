#include <string>
#include <vector>
#include <map>
#include <random>
#include <chrono>
#include "patternsplugin.h"

#ifndef _SIMULATIONS_PLUGIN
#define _SIMULATIONS_PLUGIN

typedef std::vector<int> CA;

double rate_acceptance(double f, int N);
std::vector<std::string> wright_fisher(std::vector<std::string>& pop, const std::map<std::string, double>& popfit, double mu, std::uniform_real_distribution<double>& RNG, std::default_random_engine& generator);
std::vector<std::string> quick_wright_fisher(std::vector<std::string>& pop, const std::map<std::string, double>& popfit, double mu, std::uniform_real_distribution<double>& RNG, std::default_random_engine& generator);
std::map<int, int> phenotypes_in_pop(const std::vector<std::string>& pop, const std::map<std::string, int>& popphen);

class SimulationsPlugin{
 public:
  SimulationsPlugin();
  ~SimulationsPlugin();
  
  //LOOKUP TABLES
  std::map<int,std::vector<int> > genes_per_protein;
  std::map<std::string,int> gen180;
  std::map<std::pair<int,int>, std::pair<std::vector<int>, std::vector<int> > > gen43;
  std::map<std::string, int> caid;
  std::map<std::vector<int>, int> camap;
  std::map<int, std::vector<std::string> > phenid;
  std::map<int, std::vector<std::string> > phenmap;
  std::map<int, int> caphen;
  
  //FUNCTIONS TO INITIALIZE LOOKUP TABLES
  void init_genes_per_protein(std::map<int,std::vector<int> >& genes_per_protein);
  void init_ca_phen(std::map<std::string, int>& caid, std::map<std::vector<int>, int>& camap, std::map<int, std::vector<std::string> >& phenid, std::map<int, std::vector<std::string> >& phenmap, std::map<int, int>& caphen);
  
  //FUNCTIONS
  double fitness(const std::vector<std::string>& phen, int pattern);
  std::string gen_from_prot(int pr1, int p1, int pr2, int p2, std::default_random_engine& generator, std::uniform_real_distribution<double>& RNG);
};

#endif
