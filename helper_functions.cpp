#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>
#include <chrono>
#include <numeric> 
#include <functional>
#include <stdexcept>

#include "helper_functions.h"

//CHECK_FILE
void check_file(std::ifstream& file, const std::string& name){
  if (!file){
    std::cerr << "File " << name << " doesn't exist!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return;
}
//INT_POW
int int_pow(int x, int p){
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = int_pow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}
//INT_MIN
int int_min(int a,int b){
  if (a<b)
    return a;
  else 
    return b;
}
//D_EQUAL
bool d_equal(double e1, double e2){//== operator for double
  return (std::abs(e1-e2)<1e-8);
}
//D_LESS
bool d_less(double e1, double e2){//< operator for double
  return (e1-e2<-1e-8);
}
//D_MIN
double d_min(double e1, double e2){
  if (d_less(e1,e2))
    return e1;
  else
    return e2;
}
//UNION_SEQ
double union_seq(int side1, int side2){
  //OUTPUT 
  //Given to protein sides, UNION_SEQ returns their binding energy (ENERGY)
  //
  // 0 0 0 1 | 1 0 0 0   --> this should return -4
  // 0     1 | 1     0
  // 0     0 | 0     0
  // 0 0 0 0 | 0 0 0 0

  //ARGUMENTS
  //(1) (int) SIDE1: the 1st protein side
  //(2) (int) SIDE2: the 2nd protein side

  double energy=0.0;
  //We transform the sides into std::strings
  std::string side_1=dectobin(side1,4);
  std::string side_2=dectobin(side2,4);
  for (int i=0;i<side_1.size();++i){
    if ((side_1[i]=='0' && side_2[3-i]=='1') || (side_1[i]=='1' && side_2[3-i]=='0'))
      energy-=0.3;
    else if (side_1[i]=='1' && side_2[3-i]=='1')
      energy-=2;
  }
  return energy;
}
//UNION_SEQ3
double union_seq3(int side1, const std::string& met){
  //OUTPUT
  //Given one side and a metabolite, UNION_SEQ3 returns the minimal binding energy

  //ARGUMENTS
  //(1) (int) SIDE1: the protein side
  //(2) (const std::string&) MET: the metabolite sequence

  std::vector<double> energy(met.size()-3,0.0);//in a metabolite of size 8 there are 5 sides
  // 0 1 0 1 0 1 0 1 -> metabolite of size 8
  // 0 1 0 1 -> 1st side
  //   1 0 1 0 -> 2nd side
  //     0 1 0 1 -> 3rd side
  //       1 0 1 0 -> 4th side
  //         0 1 0 1 -> 5th side

  //WE LOOK FOR THE MINIMAL BINDING ENERGY
  //HERE WE DON'T CHECK FOR REPEATED MINIMA: IF A PROTEIN CAN BIND A METABOLITE IN DIFFERENT PLACES, IT STILL BINDS IT
  double e_min=0.0;
  for (int i=0; i<energy.size(); ++i){
    energy[i]=union_seq(side1, bintodec(met.substr(i,4)));//we check every possible combination
    if (d_less(energy[i],e_min))
      e_min=energy[i];
  }
  
  return e_min;
}
//HAMMING
int vec_hamming(const std::vector<int>& v1, const std::vector<int>& v2){
  int ham=0;
  int p=v1.size();
  for (int i=0; i<p; ++i)
    if (v1[i]!=v2[i])
      ++ham;
  return ham;
}
struct HammingFunc
{
    inline int operator()(char s1,char s2)
    {
        return s1 == s2?0:1;
    }
};
int str_hamming(const std::string& s1, const std::string& s2){
  int diff = std::inner_product(s1.begin(),s1.end(),s2.begin(),0,std::plus<int>(),HammingFunc());
  return diff;
}
//REVERSE
std::string reverse(const std::string& s1){
  std::string s2(s1.size(),'0');
  for (int i=0; i<s1.size(); ++i)
    s2[i]=s1[s1.size()-1-i];
  return s2;
}
//RANDOM_GENOTYPE
std::vector<std::pair<int,int> > vec_random_genotype(int gene_number, std::uniform_real_distribution<double>& RNG, std::default_random_engine& generator){
  int total_genes=65536;
  std::vector<std::pair<int,int> > genotype(gene_number);
  //GENERATE RANDOM GENOTYPE
  for (int i=0; i<gene_number; ++i){
    //First element in the pair is the promoter, second is the coding
    genotype[i].first=RNG(generator)*16;
    genotype[i].second=RNG(generator)*total_genes;
  }
  return genotype;
}
//STR_RANDOM_GENOTYPE
std::string random_genotype(int gene_number, std::uniform_real_distribution<double>& RNG, std::default_random_engine& generator){
  int n=20*gene_number;
  std::string genotype(n, '0');
  for (int g=0; g<n; ++g){
    double mu=RNG(generator);
    if (mu>0.5)
      genotype[g]='1';
  }
  return genotype;
}
//MUTATION
std::string mutation(const std::string& genotype, int pos){
  std::string mut=genotype;
  if (mut[pos]=='0')
    mut[pos]='1';
  else
    mut[pos]='0';
  return mut;
}	     
