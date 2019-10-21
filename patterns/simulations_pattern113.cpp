#include "../helper_functions.h"
#include "patternsplugin.h"
#include "simulationsplugin.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>
#include <cstdlib>
#include <random>
#include <chrono>

//MAIN
int main(int argc, char* argv[]){

  //PLUGINS
  PatternsPlugin Model;
  SimulationsPlugin Data;
  
  //RANDOM NUMBER GENERATOR
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> RNG (0.0,1.0);

  //PARAMETERS
  int T=100000;//length of simulation
  int R=1000;//number of replicates
  int N=10000;//population size
  
  //LOOP
  std::map<int, int> mapL;
  std::ofstream out("results/simulations_pattern113.txt");
  int PO=113;//rare pattern
  for (int rep=0; rep<R; ++rep){
    std::stringstream srep;
    srep << rep;
    //RANDOM GENOTYPE
    std::string genotype=random_genotype(2, RNG, generator);
    int len=genotype.size();
    //COMPUTE PHENOTYPE
    std::vector<int> logic=Model.str_logic_function(genotype);
    int ca=Data.camap.at(logic);
    std::vector<std::string> phen=Data.phenmap.at(ca);
    int phenid=Data.caphen.at(ca);    
    double fit=Data.fitness(phen,PO);
    std::cout << rep << "\t" << "START\t" << genotype << "\t" << vec_to_str(logic) << "\t" << phenid << std::endl;
    //GILLESPIE ALGORITHM
    double t=0.0;
    int mutations=0;
    while (t<T){
      //FOR EACH POSSIBLE MUTATION WE RECORD THESE VARIABLES:
      std::vector<double> trV(len);//vector of transition rates
      std::vector<std::string> gV(len);//vector of genotypes
      std::vector<std::vector<int> > lV(len);//vector or GRNs
      std::vector<double> fV(len);//vector of fitnesses
      std::vector<int> pV(len);//vector or patterns
      double totalTR=0.0;
      for (int p=0; p<len; ++p){
	//FIND ALL NEIGHBORS
	std::string mut=mutation(genotype,p);
	//COMPUTE THEIR FITNESS
	std::vector<int> logicmut=Model.str_logic_function(mut);
	int camut=Data.camap.at(logicmut);
	std::vector<std::string> phenmut=Data.phenmap.at(camut);
	int phenidmut=Data.caphen.at(camut);
	double fitmut=Data.fitness(phenmut,PO);
	double df=fit/fitmut;
	double R=rate_acceptance(df,N)/(double)len;
	trV[p]=R;
	gV[p]=mut;
	lV[p]=logicmut;
	pV[p]=phenidmut;
	fV[p]=fitmut;
	totalTR+=R;
      }
      //COMPUTE TRANSITION TIME
      std::exponential_distribution<double> Exp(totalTR);
      double TT=Exp(generator);
      t+=TT;//update time
      //COMPUTE TRANSITION STATE
      std::discrete_distribution<int> D(trV.begin(), trV.end());
      int m=D(generator);
      genotype=gV[m];
      fit=fV[m];
      logic=lV[m];
      phenid=pV[m];
      mutations++;
    }
    mapL[phenid]++;
    std::cout << rep << "\t" << "END\t" << genotype << "\t" << vec_to_str(logic) << "\t" << phenid << std::endl;
  }//for rep
  for (auto it=mapL.begin(); it!=mapL.end(); ++it)
    out << it->first << "\t" << it->second/(double)R << std::endl;
  out.close();
}

  
