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

std::string flatten(const std::vector<std::string>& p){
  std::string s;
  for (int i=0; i<p.size(); ++i)
    s+=p[i];
  return s;
}

std::string binary(const std::string& p){
  std::string s;
  for (int i=0; i<p.size(); ++i){
    if (p[i]=='0')
      s+="00";
    else if (p[i]=='1')
      s+="01";
    else if (p[i]=='2')
      s+="10";
    else if (p[i]=='3')
      s+="11";
  }
  return s;
}


//LEMPEL ZIV ALGORITHM
int lz77(const std::string& p){
  std::map<std::string, int> lz;
  for (int c=0; c<p.size(); ++c){
    std::string p1(1,p[c]);
    if (!lz.count(p1))
      lz[p1]=lz.size();
  }

  //   Initialize the dictionary to contain all strings of length one.
  // Find the longest string W in the dictionary that matches the current input.
  // Emit the dictionary index for W to output and remove W from the input.
  // Add W followed by the next symbol in the input to the dictionary.
  // Go to Step 2.

				      
  // lz['0']=0;
  // lz['1']=2;
  std::string comp;
  std::string pnew=p;
  while (pnew.size()){
    int c=0;
    std::string d, d2;
    d+=pnew[c];
    while (1){
      ++c;
      d2=d+pnew[c];
      if (!lz.count(d2))//not in the dictionary
	break;
      d=d2;
    }
    comp+=' '+lz[d];
    pnew=pnew.substr(c);
    lz[d2]=lz.size();
  }
  
  return lz.size()-1;
}

double complexity(const std::string& p){
  int n=p.size();
  std::string s0(n,'0');
  std::string s1(n,'1');
  if (p==s0 || p==s1)
    return std::log2(n);
  else
    return 0.5*std::log2(n)*(lz77(p)+lz77(reverse(p)));
}

//MAIN
int main(int argc, char* argv[]){
  
  // std::string n1="010101010100000011101001";
  // auto d=lz77(n1);
  // // n1="ababababa";
  // // d=lz77(n1);
  // return 0;

  //PLUGINS
  PatternsPlugin Model;
  SimulationsPlugin Data;

  //READ GRN ABUNDANCE
  std::map<std::string,double> Lab;
  {
    std::ifstream in("results/grn_abundance.txt");
    check_file(in, "results/grn_abundance.txt");
    std::string ca;
    double ab;
    while (in >> ca >> ab)
      Lab[ca]=ab;
    in.close();
  }
  int lnumber=Lab.size();

  //READ CA ABUNDANCE
  std::map<int,double> Cab;
  {
    std::ifstream in("results/ca_abundance.txt");
    check_file(in, "results/ca_abundance.txt");
    int ca;
    double ab;
    while (in >> ca >> ab)
      Cab[ca]=ab;
    in.close();
  }
  int cnumber=Cab.size();

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

  //CELLULAR AUTOMATA
  std::map<int,std::string> cmap;
  {    
    std::ifstream in("results/cellular_automata.txt");
    check_file(in, "results/cellular_automata.txt");
    std::string s1;
    int id;
    while (in >> s1 >> id){
      cmap[id]=s1;
    }
    in.close();
  }
    
  //COMPLEXITY GRNs
  std::cout << "Logic" << std::endl;
  {
    std::ofstream out("results/complexity_grn.txt");
    for (auto it=Lab.begin(); it!=Lab.end(); ++it){
      auto L=it->first;
      double ab=it->second;
      auto bL=binary(L);
      auto C=complexity(bL);
      out << C << "\t" << ab << std::endl;
    }
    out.close();
  }
  
  //COMPLEXITY CELLULAR AUTOMATA
  std::cout << "Cellular automata" << std::endl; 
  {
    std::ofstream out("results/complexity_ca.txt");
    for (auto it=Cab.begin(); it!=Cab.end(); ++it){
      int id=it->first;
      auto L=cmap.at(id);
      double ab=it->second;
      auto bL=binary(L);
      auto C=complexity(bL);
      out << C << "\t" << ab << std::endl;
    }
    out.close();
  }

  //COMPLEXITY PATTERNS
  std::cout << "Patterns" << std::endl;
  {
    std::ofstream out("results/complexity_pattern.txt");
    for (auto it=Data.phenid.begin(); it!=Data.phenid.end(); ++it){
      auto L=flatten(it->second);
      int id=it->first;
      double ab=Pab.at(id);
      auto bL=binary(L);
      auto C=complexity(bL);
      out << C << "\t" << ab << std::endl;
    }
    out.close();
  }
}
