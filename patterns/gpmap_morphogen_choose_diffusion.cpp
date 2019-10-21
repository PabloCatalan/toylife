#include "../helper_functions.h"
#include "patternsplugin.h"
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

int main(int argc, char* argv[]){

  PatternsPlugin Model;
  
  //SPATIAL DISPOSITION
  std::stringstream ss;
  ss << argv[1];
  std::stringstream ss2;
  ss2 << argv[2];
  std::stringstream sring;
  sring << argv[3];
  std::stringstream ssm;
  ssm << argv[4];
  std::stringstream ssdiff;
  ssdiff << argv[5];
  int cells;
  int time;
  bool ring;
  std::string gradient;
  ss >> cells;
  ss2 >> time;
  sring >> ring;
  ssm >> gradient;

  //DIFFUSION TEXT
  std::string sdiff="nothing";
  if (ssdiff.str()=="1")
    sdiff="only_P1";
  else if (ssdiff.str()=="2")
    sdiff="only_P2";
  else if (ssdiff.str()=="3")
    sdiff="only_P1_P2";
  else if (ssdiff.str()=="4")
    sdiff="also_dimer";  

  //LOGIC ABUNDANCE
  std::map<std::string, double> logsize;
  {
    std::ifstream in("results/logic_abundance.txt");
    check_file(in, "results/logic_abundance.txt");
    std::string s1;
    double d1;
    while (in >> s1 >> d1)
      logsize[s1]=d1;
    in.close();
  }
  //READ ALL CELLULAR AUTOMATA (ONLY THOSE THAT APPEAR WHEN P1 DIFFUSES)
  typedef std::vector<int> CA;
  std::map<CA, int> camap;
  std::map<CA, int> caid;
  std::map<int,double> casize;
  {
    //CA THAT APPEAR
    {
      std::ifstream in("results/cellular_automata.txt");
      check_file(in, "results/cellular_automata.txt");
      double d1;
      std::string s1;
      int id=0;
      while (in >> s1 >> d1){
	std::vector<int> v=str_to_vec(s1);
	caid[v]=id;
	++id;
      }
      in.close();
    }
    
    //CA THAT APPEAR WHEN P1 DIFFUSES
    std::map<std::string, int> cap2;
    std::map<std::string, double> cap2size;
    {
      std::ifstream in("results/logic_ca_"+ssdiff.str()+".txt");
      check_file(in, "results/logic_ca_"+ssdiff.str()+".txt");
      std::string s1, s2;
      while (in >> s1 >> s2){
	double logs=logsize[s1];
	std::vector<int> v=str_to_vec(s2);	
	int id=caid[v];
	camap[v]=id;
	casize[id]+=logs;
      }
      in.close();
    }
  }
  //MORPHOGEN
  std::vector<std::vector<int> > morphogen(time, std::vector<int>(cells,0));
  for (int t=0; t<time; ++t)
    morphogen[t]=morphogen_creator(cells, t, gradient);
  {
    std::ofstream out("results/morphogen_"+ss.str()+"_"+ss2.str()+"_"+sring.str()+"_"+gradient+".txt");
    for (int i=0; i<morphogen.size(); ++i){
      for (int j=0; j<morphogen[j].size(); ++j)
	out << morphogen[i][j];
      out << std::endl;
    }
    out.close();
  }
  //LOOP OVER ALL CA
  typedef std::vector<std::vector<int> > phen;
  std::map<std::pair<int,int>, double> mp;
  std::map<std::vector<std::vector<int> >, double> mapid;
  std::map<phen, int> phenid;
  std::map<int,int> caphen;
  auto it=camap.begin();
  for (; it!=camap.end(); ++it){
    //GENOTYPE
    CA ca1=it->first;
    double d1=casize[it->second];
    //PHENOTYPE
    std::vector<std::vector<int> > phenotype(time, std::vector<int>(cells,0));
    //INPUT
    std::vector<int> input(cells, 0);
    //COMPUTE PHENOTYPE
    for (int t=0; t<time; ++t){
      //ADD MORPHOGEN
      std::vector<int> morphoact=morphogen[t];
      std::vector<int> input2(cells, 0);//recompute input each time step
      for (int c=0; c<cells; ++c){
	//BUILD INPUT STATE FROM INPUT IN C-1, C AND C+1
	std::vector<int> in2(3,0);
	int prev=(c-1)%cells;
	int next=(c+1)%cells;
	in2[0]=input[prev];
	in2[1]=input[c];
	in2[2]=input[next];
	//IMPACT OF MORPHOGEN
	//0 no hay nada, 1 esta la p2, 2 esta la p1, 3 estan las dos (o el dimero), 4 estan el dimero y p1, 5 están el dimero y p2
	if (morphoact[prev]==1){//if morphoact==1 it means P1 is morphogen
	  switch(input[prev]){
	  case 0: in2[0]=2; break;
	  case 1: in2[0]=3; break;
	  }
	}
	else if (morphoact[prev]==2){//if morphoact==2 it means P2 is morphogen
	  switch(input[prev]){
	  case 0: in2[0]=1; break;
	  case 2: in2[0]=3; break;
	  }
	}
	if (morphoact[c]==1){
	  switch(input[c]){
	  case 0: in2[1]=2; break;
	  case 1: in2[1]=3; break;
	  }
	}
	else if (morphoact[c]==2){
	  switch(input[c]){
	  case 0: in2[0]=1; break;
	  case 2: in2[0]=3; break;
	  }
	}
	if (morphoact[next]==1){
	  switch(input[next]){
	  case 0: in2[2]=2; break;
	  case 1: in2[2]=3; break;
	  }
	}
	else if (morphoact[next]==2){
	  switch(input[next]){
	  case 0: in2[0]=1; break;
	  case 2: in2[0]=3; break;
	  }
	}
	if (!ring){//if boundary conditions are not periodic
	  if (c==0)
	    in2[0]=0;
	  else if (c==cells-1)
	    in2[2]=0;
	}
	int inputstate=base_n_to_10(vec_to_str(in2), 4);
	//COMPUTE OUTPUT USING CA
	input2[c]=ca1[inputstate];
      }//for all cells
      input=input2;//for next time step
      phenotype[t]=input;
    }//for all times
    
    //CHECK SPATIOTEMPORAL PERIODICITY
    // std::pair<int,int> per1=Model.periodicity(phenotype);
    // mp[per1]+=d1;
    //SAVE PHENOTYPE
    mapid[phenotype]+=d1;   
    if (!phenid.count(phenotype))
      phenid[phenotype]=phenid.size();  
    caphen[it->second]=phenid[phenotype];//add id of the CA     
  }//for ca
  
  //OUTPUT PATTERNS
  std::ofstream out("results/phen_id_"+ss.str()+"_"+ss2.str()+"_"+sring.str()+"_"+gradient+"_"+sdiff+".txt");
  for (auto it=mapid.begin(); it!=mapid.end(); ++it){
    std::vector<std::vector<int> > p1=it->first;
    out << phenid[p1] << "\t" << it->second << std::endl;
    for (int t=0; t<time; ++t){
      for (int c=0; c<cells; ++c)
	out << p1[t][c];
      out << std::endl;
    }
  }
  out.close();
  //OUTPUT PATTERN SIZE
  out.open("results/phen_size_"+ss.str()+"_"+ss2.str()+"_"+sring.str()+"_"+gradient+"_"+sdiff+".txt");
  for (auto it=mapid.begin(); it!=mapid.end(); ++it){
    std::vector<std::vector<int> > p1=it->first;
    out << phenid[p1] << "\t" << it->second << std::endl;
  }
  out.close();
  //SAVE PHENS FOR EACH CA
  out.open("results/ca_phen_"+ss.str()+"_"+ss2.str()+"_"+sring.str()+"_"+gradient+"_"+sdiff+".txt");
  for (auto it=caphen.begin(); it!=caphen.end(); ++it)
    out << it->first << "\t" << it->second << std::endl;
  out.close();
}

