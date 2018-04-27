#include "helper_functions.h"
#include "toy_plugin.h"
#include <random>
#include <chrono>

int main(){
  //INITIALIZE toyLIFE'S LOOKUP TABLES
  ToyPlugin Model;

  //RANDOM NUMBER GENERATOR
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> RNG (0.0,1.0);
  
  //GENERATE RANDOM GENOTYPE (FUNCTION IN "helper_functions.cpp")
  int gene_number=5;
  std::string genotype=random_genotype(gene_number, RNG, generator);
  int gensize=genotype.size();

  //MUTATE THE GENOTYPE (FUNCTION IN "helper_functions.cpp")
  int mutated_position=RNG(generator)*gensize;
  std::string mutgen=mutation(genotype, mutated_position);

  //COMPUTE REGULATORY LOGIC FUNCTION (FUNCTIONS IN "toy_plugin.cpp")
  std::vector<int> logic_function=Model.regulatory_phenotype(genotype);
  //COMPUTE REGULATORY ROBUSTNESS
  double rrob=Model.regulatory_robustness(genotype, logic_function);
  //COMPUTE REGULATORY EVOLVABILITY
  std::map<std::vector<int>, int> revo=Model.regulatory_evolvability(genotype, logic_function);
  
  //COMPUTE METABOLIC PHENOTYPE
  std::vector<int> phenotype=Model.metabolic_phenotype(genotype, Model.env);
  //COMPUTE METABOLIC ROBUSTNESS
  double mrob=Model.metabolic_robustness(genotype, phenotype, Model.env);
  //COMPUTE METABOLIC EVOLVABILITY
  std::map<std::vector<int>, int> mevo=Model.metabolic_evolvability(genotype, phenotype, Model.env);  

}
