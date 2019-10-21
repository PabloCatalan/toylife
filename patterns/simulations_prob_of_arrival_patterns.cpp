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

//AVERAGE
template <typename T>
std::pair<double,double> average(const std::vector<T>& v){
  double a1=0.0;
  double a2=0.0;
  for (int i=0; i<v.size(); ++i){
    a1+=v[i];
    a2+=v[i]*v[i];
  }
  double avg=a1/v.size();
  double var=std::sqrt(a2/v.size()-avg*avg);
  return std::make_pair(avg,var);
}

//MAIN
int main(int argc, char* argv[]){

  //PLUGINS
  PatternsPlugin Model;
  SimulationsPlugin Data;

  //READ PATTERNS ABUNDANCE
  std::map<int,double> Pab;
  {
    std::ifstream in("results/phen_size_31_100_0_point_only_P1.txt");
    check_file(in, "results/phen_size_31_100_0_point_only_P1.txt");
    int ca;
    double ab;
    while (in >> ca >> ab)
      Pab[ca]=ab;
    in.close();
  } 
  int pnumber=Pab.size();

  //RANDOM NUMBER GENERATOR
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> RNG (0.0,1.0);

  //PARAMETERS
  int T=100;//length of simulation
  int R=100;//number of replicates
  int N=10000;//population size

  //FOR ALL PATTERNS  
  std::ofstream out("results/simulations_prob_of_arrival_patterns.txt");
  for (int PO=0; PO<pnumber; ++PO){
    std::cout << "Pattern " << PO << std::endl;
    std::vector<double> probV(R,0.0);
    std::vector<double> timeV;
    //REPETITIONS
    for (int rep=0; rep<R; ++rep){
      //RANDOM GENOTYPE
      std::string genotype=random_genotype(2, RNG, generator);
      int len=genotype.size();
      //COMPUTE PHENOTYPE
      std::vector<int> logic=Model.str_logic_function(genotype);
      int ca=Data.camap.at(logic);
      std::vector<std::string> phen=Data.phenmap.at(ca);
      int phenid=Data.caphen.at(ca);    
      double fit=Data.fitness(phen,PO);
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
	//IF WE ARRIVE, BREAK
	if (phenid==PO){
	  probV[rep]=1.0;
	  timeV.push_back(t);
	  break;
	}
      }//while T
      std::cout << PO << "\t" << rep << "\t" << t << std::endl;
    }//for rep
    auto Pres=average(probV);
    auto Tres=average(timeV);
    out << PO << "\t" << Pres.first << "\t" << Pres.second << "\t" << Tres.first << "\t" << Tres.second << "\t" << Pab.at(PO) << std::endl;
    std::cout << PO << "\t" << Pres.first << "\t" << Pres.second << "\t" << Tres.first << "\t" << Tres.second << "\t" << Pab.at(PO) << std::endl;
  }//for pattern
  out.close();
}

  
