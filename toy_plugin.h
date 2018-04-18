#ifndef _TOY_PLUGIN
#define _TOY_PLUGIN

#include "helper_functions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>
#include <unordered_map>
#include <cstdlib>
#include <iterator>
#include <random>
#include <chrono>

//STRUCTURES
typedef int Prot;
typedef std::string Met;

struct Dim{
  Dim();//constructor
  Dim(const int, const int, const int);//constructor with id
  bool empty() const;//checks if the Dimer is empty
  int id;
  Prot p1;
  Prot p2;
  bool operator<(const Dim& rhs) const{//in order to compare and use this in maps
    if (id < rhs.id)
      return 1;
    else if (id == rhs.id)
      if (p1<rhs.p1)
	return 1;
      else if (p1==rhs.p1)
	if (p2<rhs.p2)
	  return 1;
    return 0;
  }//operator<
  bool operator==(const Dim& rhs) const{
    return (id==rhs.id && p1==rhs.p1 && p2==rhs.p2);
  }
};

struct OWM{
  OWM();//constructor
  OWM(const Prot, const Dim&, const Met&);
  bool empty() const;
  Prot prot;
  Dim dim;
  Met met;
  bool operator<(const OWM& rhs) const{//in order to compare and use this in maps
    if (prot < rhs.prot)
      return 1;
    else if (prot == rhs.prot)
      if (dim<rhs.dim)
	return 1;
      else if (dim.id==rhs.dim.id && dim.p1==rhs.dim.p1 && dim.p2==rhs.dim.p2)
	if (met<rhs.met)
	  return 1;
    return 0;
  }//operator<
};

typedef std::map<Prot,int> mapa_prot;
typedef std::map<Dim, int> mapa_dim;
typedef std::map<Met,int> mapa_met;
typedef std::map<OWM,int> mapa_owm;

//for dimeros_metabolitos
typedef std::pair<std::string, int> pairDM;
struct Dmet{//for dimeros_metabolitos
  Dmet();//constructor
  Dmet(const double, const std::string&, const int);//constructor with id
  double eg;
  std::string seq;
  int pos;
};

template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

struct KeyHash {
  std::size_t operator()(const pairDM& k) const{
    std::size_t seed = 0;
    hash_combine(seed, k.first);
    hash_combine(seed, k.second);
    return seed;//std::hash<std::string>()(k.first)^(std::hash<int>()(k.second << 1));
  }
};

struct KeyEqual {
  bool operator()(const pairDM& lhs, const pairDM& rhs) const{
    return lhs.first == rhs.first && lhs.second == rhs.second;
  }
};

typedef std::unordered_map<pairDM, Dmet, KeyHash, KeyEqual> mapa_dmet;
typedef std::pair<Met,Prot> pairPM;//it is the same as pairDM
typedef std::unordered_map<pairPM, double, KeyHash, KeyEqual> mapa_pmet;//for proteinas_metabolitos
typedef std::map<std::pair<int, std::string>, std::pair<double, int> > mapa_pbreak;

class ToyPlugin{
 public:
  ToyPlugin();
  ~ToyPlugin();

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
  std::vector<std::string> mets;		
  std::map<std::string,int> not_broken;
  std::vector<std::string> env;
  std::vector<std::vector<int> > neighbors_prom;
  std::vector<std::vector<int> > neighbors_coding;
  mapa_dmet dim_met;
  mapa_pmet prot_met;
  mapa_pbreak prot_breaking;

  //FUNCTIONS
  std::vector<int> regulatory_phenotype(const std::string& genotype);
  std::vector<int> metabolic_phenotype(const std::string& genotype, const std::vector<std::string>& env);
  double metabolic_robustness(const std::string& genotype, const std::vector<int>& phen, const std::vector<std::string>& env);
  double regulatory_robustness(const std::string& genotype, const std::vector<int>& logic_vector);
  std::map<std::vector<int>, int> metabolic_evolvability(const std::string& genotype, const std::vector<int>& phen, const std::vector<std::string>& env);
  std::map<std::vector<int>, int> regulatory_evolvability(const std::string& genotype, const std::vector<int>& logic_vector);

 private:
  //FUNCTIONS TO INITIALIZE LOOKUP TABLES
  void init_gen(std::vector<int>& prot_gen);
  void init_prot(std::vector<int>& perimeters, std::vector<double>& energies, std::vector<std::vector<double> >& prot_prom_energies, std::vector<std::vector<bool> >& prot_pol);
  void init_dim(std::map<std::pair<int, int>, int>& dim_perim, std::vector<double>& dim_energy, std::vector<std::vector<double> >& dim_prom_energies, std::vector<std::vector<bool> >& dim_pol);
  void init_polymerase(std::vector<double>& polymerase);
  void simplified_mets(std::vector<std::string>& mets, std::map<std::string,int>& not_broken, std::vector<std::string>& env);
  void neighbors(std::vector<std::vector<int> >& neighbors_prom, std::vector<std::vector<int> >& neighbors_coding);
  void init_dim_met(mapa_dmet& dim_met);
  void init_prot_met(mapa_pmet& prot_met);
  void init_prot_breaking(mapa_pbreak& prot_breaking);
  //FUNCTIONS TO COMPUTE PHENOTYPE
  int promoter_expression(const std::vector<std::pair<int,int> >& genotype, std::pair<mapa_prot, mapa_dim>& objects);
  bool finding_cycles(std::vector<std::pair<int, std::pair<mapa_prot, mapa_dim> > >& iterations, int expression);
  void dimerization(std::pair<mapa_prot, mapa_dim>& objects);
  std::vector<std::pair<int, std::pair<mapa_prot, mapa_dim> > > boolean_network(const std::vector<std::pair<int,int> >& genotype);
  std::vector<int> complete_boolean_function(const std::vector<std::pair<int,int> >& genotype);
  void reacting(std::pair<mapa_prot, mapa_dim>& objects, mapa_met& met, mapa_owm& objects_with_mets, double& fit);
  int metabolism(const std::vector<std::pair<int,int> >& genotype, std::pair<mapa_prot, mapa_dim> objects, mapa_met& met);
  std::vector<int> vec_phenotype(const std::vector<std::pair<int,int> >& genotype, const std::vector<std::string>& env);
  std::vector<std::pair<int,int> > genotype_str_to_vec(const std::string& genotype);
};

#endif
