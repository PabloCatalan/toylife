#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "../helper_functions.h"
#include "patternsplugin.h"
#include "simulationsplugin.h"

SimulationsPlugin::SimulationsPlugin(){
  //LOOKUP TABLES
  init_genes_per_protein(genes_per_protein);
  init_ca_phen(caid, camap, phenid, phenmap, caphen);
}
SimulationsPlugin::~SimulationsPlugin(){
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void SimulationsPlugin::init_genes_per_protein(std::map<int,std::vector<int> >& genes_per_protein){
  std::ifstream in("../data/protein_gene.txt");
  check_file(in, "../data/protein_gene.txt");
  int p1;
  int g1=0;
  while (in >> p1){
    if (!genes_per_protein.count(p1))
      genes_per_protein[p1]=std::vector<int>();
    genes_per_protein[p1].push_back(g1);
    ++g1;
  }
  in.close();
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//GENOTYPE FROM PROTEINS AND PROMOTERS
std::string SimulationsPlugin::gen_from_prot(int pr1, int p1, int pr2, int p2, std::default_random_engine& generator, std::uniform_real_distribution<double>& RNG){
  std::string gen;
  int c1=RNG(generator)*genes_per_protein.at(p1).size();
  c1=genes_per_protein.at(p1)[c1];
  int c2=RNG(generator)*genes_per_protein.at(p2).size();
  c2=genes_per_protein.at(p2)[c2];
  gen += dectobin(pr1, 4);
  gen += dectobin(c1, 16);
  gen += dectobin(pr2, 4);
  gen += dectobin(c2, 16);
  return gen;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//CELLULAR AUTOMATA ASSOCIATE TO EACH LOGIC FUNCTION
void SimulationsPlugin::init_ca_phen(std::map<std::string, int>& caid, std::map<std::vector<int>, int>& camap, std::map<int, std::vector<std::string> >& phenid, std::map<int, std::vector<std::string> >& phenmap, std::map<int, int>& caphen){
  //CA IDENTITY
  {    
    //CA THAT APPEAR
    std::ifstream in("results/cellular_automata.txt");
    check_file(in, "results/cellular_automata.txt");
    std::string s1;
    int id;
    while (in >> s1 >> id){
      caid[s1]=id;
    }
    in.close();

    //CA THAT APPEAR WHEN P1 DIFFUSES
    {
      std::ifstream in("results/logic_ca_1.txt");
      check_file(in, "results/logic_ca_1.txt");
      std::string s1, s2;
      while (in >> s1 >> s2){
	std::vector<int> v=str_to_vec(s1);
	int id=caid[s2];
	camap[v]=id;
      }
      in.close();
    }
  }
  //PHENOTYPES ASSOCIATED TO EACH CELLULAR AUTOMATA
  {
    //IDS ASSOCIATED TO EACH PHENOTYPE
    {
      std::ifstream in("results/phen_id_31_100_0_point_only_P1.txt");
      check_file(in, "results/phen_id_31_100_0_point_only_P1.txt");
      std::string line;
      while (std::getline(in, line)){
	std::stringstream sread;
	sread << line;
	int p1;
	double d1;
	sread >> p1 >> d1;
	std::vector<std::string> phen;
	for (int t=0; t<100; ++t){
	  std::string s1;
	  std::getline(in,s1);
	  phen.push_back(s1);
	}
	phenid[p1]=phen;
      }
      in.close();
    }
    //CA ASSOCIATED TO EACH PHEN
    std::ifstream in("results/ca_phen_31_100_0_point_only_P1.txt");
    check_file(in, "results/ca_phen_31_100_0_point_only_P1.txt");
    int ca;
    int p1;
    while (in >> ca >> p1){
      phenmap[ca]=phenid[p1];
      caphen[ca]=p1;
    }
    in.close();  
  }    
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//WRIGHT FISHER (TOO SLOW)
std::vector<std::string> wright_fisher(std::vector<std::string>& pop, const std::map<std::string, double>& popfit, double mu, std::uniform_real_distribution<double>& RNG, std::default_random_engine& generator){
  int N=pop.size();
  std::vector<std::string> newpop(N);
  double F=0.0;
  //FITNESS MAP
  std::vector<double> fitlevels(N);
  for (int i=0; i<N; ++i){
    double f1=popfit.at(pop[i]);
    fitlevels[i]=f1+F;
    F+=f1;
  }
  //IF FITNESS IS 0
  if (d_equal(F, 0.0)){//choose at randome
    for (int g=0; g<N; ++g){
      int m1=RNG(generator)*N;
      newpop[g]=pop[m1];
    }
    return newpop;
  }
  //ELSE GENERATE NEW POPULATION
  for (int i=0; i<N; ++i){
    std::string newgen;
    //CHOOSE ONE INDIVIDUAL TO GIVE BIRTH
    double m1=RNG(generator)*F;
    for (int g=0; g<N; ++g){
      if (m1<fitlevels[g]){
	newgen=pop[g];
	break;
      }
    }
    if (newgen=="")
      std::getchar();
    //MUTATE
    if (RNG(generator)<mu){
      int pos=RNG(generator)*newgen.size();
      newgen=mutation(newgen, pos);
    }
    newpop[i]=newgen;
  }
  return newpop;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//QUICK WRIGHT FISHER
std::vector<std::string> quick_wright_fisher(std::vector<std::string>& pop, const std::map<std::string, double>& popfit, double mu, std::uniform_real_distribution<double>& RNG, std::default_random_engine& generator){
  int N=pop.size();
  std::vector<std::string> newpop(N);
  double F=0.0;
  //FITNESS MAP
  std::map<double, std::vector<std::string> > fitmap;
  std::vector<double> fitlevels(N);
  for (int i=0; i<N; ++i){
    double f1=popfit.at(pop[i]);
    if (!fitmap.count(f1))
      fitmap[f1]=std::vector<std::string>();
    fitmap[f1].push_back(pop[i]);
    F+=f1;
  }
  //IF FITNESS IS 0
  if (d_equal(F, 0.0)){//choose at randome
    for (int g=0; g<N; ++g){
      int m1=RNG(generator)*N;
      newpop[g]=pop[m1];
    }
    return newpop;
  }
  //GENERATE NEW POPULATION
  for (int i=0; i<N; ++i){
    std::string newgen;
    //CHOOSE ONE INDIVIDUAL TO GIVE BIRTH
    double m1=RNG(generator)*F;
    double basic=0.0;
    for (int g=0; g<fitmap.size(); ++g){
      auto it=std::next(fitmap.begin(), g);
      double f1=it->first;
      double N1=it->second.size();
      double compare=f1*N1+basic;
      if (m1<compare){
	int mut=RNG(generator)*it->second.size();
	newgen=it->second[mut];
	break;
      }
      basic+=N1*f1;
    }      
    //MUTATE
    if (RNG(generator)<mu){
      int pos=RNG(generator)*newgen.size();
      newgen=mutation(newgen, pos);
    }
    newpop[i]=newgen;
  }
  return newpop;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//FITNESS (time 100, cells 31)
double SimulationsPlugin::fitness(const std::vector<std::string>& phen, int pattern){
  //PATTERN IS THE TARGET, WITH FITNESS 1
  //THE REMAINING PATTERNS HAVE LESS FITNESS
  //AS THEY ARE MORE DIFFERENT FROM TARGET
  std::vector<std::string> ref=phenid.at(pattern);  
  double d1=0.0;
  int p=phen.size();
  for (int t=0; t<p; ++t)
    d1 += str_hamming(phen[t], ref[t]);
  d1 /= 100.0*ref[0].size();//we average the differences over 100 time steps
  if (d1<0.0001)
    return 1.0;
  else if (d1>0.9999)
    return 1e-8;//to avoid NaNs when dividing by zero
  else
    return 1.0-d1;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//RATE OF ACCEPTANCE
double rate_acceptance(double f, int N){
  if (d_equal(f,1.0))
    return 1.0;//neutral mutations are accepted at mutation rate
  else
    return N*(f-1)/(std::pow(f,N)-1);
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//PHENOTYPES IN POPULATION
std::map<int, int> phenotypes_in_pop(const std::vector<std::string>& pop, const std::map<std::string, int>& popphen){
  std::map<int,int> phenpop;
  int N=pop.size();
  for (int p=0; p<N; ++p){
    int phen=popphen.at(pop[p]);
    phenpop[phen]++;
  }
  return phenpop;
}  
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
