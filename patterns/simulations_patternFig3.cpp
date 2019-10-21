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
  std::map<std::string, int> mapL;
  std::ofstream out("results/simulations_patternFig3.txt");
  int PO=108;//PATTERN WITH ALL 2s (in Fig3)
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
    double fit=Data.fitness(phen,PO);//all 2s
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
    std::string strL=vec_to_str(logic);
    strL=strL.substr(0,4);
    mapL[strL]++;
    std::cout << rep << "\t" << "END\t" << genotype << "\t" << vec_to_str(logic) << "\t" << phenid << std::endl;
  }//for rep
  for (auto it=mapL.begin(); it!=mapL.end(); ++it)
    out << it->first << "\t" << it->second/(double)R << std::endl;
  out.close();

  //COMPARISON WITH THEORETICAL ABUNDANCES:
  //CHECK HOW MANY GENOTYPES MAP INTO EACH LOGIC FUNCTION THAT GENERATES THE PATTERN "ALL 2S" IN FIG2

  //READ GRN ABUNDANCES
  std::map<std::string, double> mapL2;
  std::ifstream in("results/logic_abundance.txt");
  check_file(in, "results/logic_abundance.txt");
  std::string l1;
  double d1;
  double tot=0.0;
  while (in >> l1 >> d1){
    std::string crop=l1.substr(0,4);
    std::cout << crop << std::endl;
    mapL2[crop]+=d1;
    tot +=d1;
  }

  //WRITE OUT THE THEORETICAL ABUNDANCES
  out.open("results/simulations_patternFig3_theory.txt");
  double tot2=0.0;
  for (auto it=mapL2.begin(); it!=mapL2.end(); ++it){
    auto LF=it->first;
    auto ab=it->second;
    if (LF[0]=='2' && LF[2]=='2')
      tot2+=ab;
  }
  for (auto it=mapL2.begin(); it!=mapL2.end(); ++it){
    auto LF=it->first;
    auto ab=it->second;
    if (LF[0]=='2' && LF[2]=='2')
      out << LF << "\t" << ab << "\t" << ab/tot2 << std::endl;
  }
  out.close();
}

  
