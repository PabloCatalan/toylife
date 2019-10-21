#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>

#include "patternsplugin.h"
#include "../helper_functions.h"

PatternsPlugin::PatternsPlugin(){  
  //INITIALIZE ALL OF TOYLIFE'S LOOKUP TABLES
  init_gen(prot_gen);
  init_prot(prot_perimeters, prot_energies, prot_prom_energies, prot_pol);
  init_dim(dim_perim, dim_bond_energy, dim_prom_energies, dim_pol);
  init_polymerase(polymerase);  
  neighbors(neighbors_prom, neighbors_coding);
}
PatternsPlugin::~PatternsPlugin(){
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//INIT GEN
void PatternsPlugin::init_gen(std::vector<int>& prot_gen){//reads to table which protein is given by which gene
  std::string name="../data/protein_gene.txt";
  std::ifstream in(name.c_str());
  check_file(in, name);
  int total_genes=65536;
  prot_gen.resize(total_genes);
  for (int i=0; i<total_genes; ++i)
    in >> prot_gen[i];
  in.close();
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//INIT POLYMERASE
void PatternsPlugin::init_polymerase(std::vector<double>& polymerase){
  double pol_energy=-11.0;
  for (int i=0; i<16; ++i){
    if (i==5)
      polymerase.push_back(pol_energy-4.0);
    else if (i==7)
      polymerase.push_back(pol_energy-4.3);
    else if (i==11)
      polymerase.push_back(pol_energy-2.9);
    else if (i==13)
      polymerase.push_back(pol_energy-4.3);
    else if (i==14)
      polymerase.push_back(pol_energy-2.9);
    else if (i==15)
      polymerase.push_back(pol_energy-4.6);
    else
      polymerase.push_back(0.0);
  }
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//INIT PROTEIN INFORMATION
void PatternsPlugin::init_prot(std::vector<int>& perimeters, std::vector<double>& energies, std::vector<std::vector<double> >& prot_prom_energies, std::vector<std::vector<bool> >& prot_pol){//reads to table the perimeter and energies associated to each protein
  //NUMBER OF PROTEINS
  int total_proteins=0;
  std::string name="../data/proteins_with_perimeter.txt";
  std::ifstream in(name.c_str());
  check_file(in,name);
  std::string lee;
  while (std::getline(in, lee))
    ++total_proteins;
  perimeters.resize(total_proteins);
  energies.resize(total_proteins);
  in.close();
  //PERIMETERS
  name="../data/proteins_with_perimeter.txt";
  in.open(name.c_str());
  check_file(in, name);
  int id, per;
  while (in >> id >> per)
    perimeters[id]=per;
  in.close();
  perimeters[0]=-1;
  //ENERGIES
  name="../data/proteins.txt";
  in.open(name.c_str());
  check_file(in, name);
  double dirt, en;
  in >> id >> en;//open chain
  energies[id]=en;
  while (in >> id >> en){
    energies[id]=en;
    for (int j=0; j<8; j++)//each protein has 8 sides (4 for each isomer), but we don't need them now
      in >> dirt;
  }
  in.close();
      
  //LOOP OVER ALL PROMOTERS
  prot_prom_energies.resize(16);
  prot_pol.resize(16);
  for (int i=0; i<16; ++i){
    std::vector<double> energias_prom(total_proteins,0.0);      
    std::vector<bool> polim(total_proteins,0);
    std::stringstream ss;//create a stringstream
    ss << i;//add number to the stream    
    std::string name="../promoters/proteins_bind_prom"+ss.str()+".txt";
    std::ifstream in(name.c_str());
    check_file(in,name);
    int ind_prot;
    double prom_energ;
    bool prom_pol;
    while (in >> ind_prot){
      in >> prom_energ >> prom_pol;
      energias_prom[ind_prot]=prom_energ;
      polim[ind_prot]=prom_pol;
    }
    prot_prom_energies[i]=energias_prom;
    prot_pol[i]=polim;
    in.close();
  }
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//INIT DIMER INFORMATION
void PatternsPlugin::init_dim(std::map<std::pair<int, int>, int>& dim_perim, std::vector<double>& dim_energy, std::vector<std::vector<double> >& dim_prom_energies, std::vector<std::vector<bool> >& dim_pol){//reads to table the information relating to dimers
  int total_dimers=0;
  std::string name="../data/dimers.txt";
  std::ifstream in(name.c_str());
  check_file(in,name);
  int id=0;
  std::pair<int,int> perims;
  int basura;
  while (in >> perims.first >> perims.second >> basura){
    double e1;
    in >> e1;
    dim_energy.push_back(e1);
    int dirt;
    for (int j=0; j<12; ++j)
      in >> dirt;
    dim_perim[perims]=id;
    ++id;
  }
  in.close();
  total_dimers=id;
  
  //LOOP OVER ALL PROMOTERS
  dim_prom_energies.resize(16);
  dim_pol.resize(16);
  for (int i=0; i<16; ++i){
    std::vector<double> energies(total_dimers,0.0);      
    std::vector<bool> pol(total_dimers,0);
    std::stringstream ss;//create a stringstream
    ss << i;//add number to the stream
    std::string name="../promoters/dimers_bind_prom"+ss.str()+".txt";
    std::ifstream in(name.c_str());
    check_file(in,name);
    int ind_dim; 
    double prom_energ;
    bool prom_pol;     
    while (in >> ind_dim){
      in >> prom_energ >> prom_pol;
      energies[ind_dim]=prom_energ;
      pol[ind_dim]=prom_pol;
    }
    dim_prom_energies[i]=energies;
    dim_pol[i]=pol;
    in.close();
  }//for i
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//NEIGHBORS
void PatternsPlugin::neighbors(std::vector<std::vector<int> >& neighbors_prom, std::vector<std::vector<int> >& neighbors_coding){
  neighbors_prom=std::vector<std::vector<int> >(16, std::vector<int>(4));
  std::ifstream file_prom("../data/neighbors_prom.txt");
  check_file(file_prom, "../data/neighbors_prom.txt");
  for (int i=0; i<16; ++i)
    for (int g=0; g<4; ++g)
      file_prom >> neighbors_prom[i][g];
  file_prom.close();

  neighbors_coding=std::vector<std::vector<int> >(65536, std::vector<int>(16));
  std::ifstream file_coding("../data/neighbors_coding.txt");
  check_file(file_coding, "../data/neighbors_coding.txt");
  for (int i=0; i<65536; ++i)
    for (int g=0; g<16; ++g)
      file_coding >> neighbors_coding[i][g];
  file_coding.close();
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//DIMER FROM TWO
int PatternsPlugin::dimer_from_two(int prot1, int prot2){
  //OUTPUT
  //Given two proteins, PAIRS_OF_DIMS returns the dimer they form. If they don't form any dimer, it returns -1. 

  //ARGUMENTS
  //(1) (int) PROT1: first protein
  //(2) (int) PROT2: second protein
  //(3) (const std::vector<int>&) PROT_PERIMETERS: vector of size #total_proteins, lists the perimeter of every protein
  //(4) (const std::map<std::pair<int,int>,int>&) DIM_PERIM: map holding as key the pair of perimeters that form each 
  //     dimer, and as value the index of the dimers

  //WE SORT THE PROTEINS IN ASCENDING ORDER
  if (prot1>prot2){
    int prot_reserva=prot2;
    prot2=prot1;
    prot1=prot_reserva;
  }

  std::pair<int,int> perim=std::make_pair(prot_perimeters[prot1], prot_perimeters[prot2]);
  std::map<std::pair<int,int>,int>::const_iterator it=dim_perim.find(perim);
  if (it!=dim_perim.end())//if this pair forms a dimer
    return it->second;
  else//if it doesn't
    return -1;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//FIND CYCLE
int PatternsPlugin::find_cycle1(const std::vector<int>& v, int i1){
  int size=v.size();
  std::vector<int> v1=v;
  std::map<int,int> m;
  for (int i=0; i<v.size(); ++i){
    if (!m.count(v[i]))
      m[v[i]]=m.size()-1;
    v1[i]=m[v[i]];
  }
  std::string s=vec_to_str(v1);
  std::string s2=s.substr(0, size/2);
  if (i1==1){//temporal pattern
    int state=v[size-1];
    for (int i=size-2; i>-1; --i)
      if (v[i]==state)
	return (size-1-i);
  }
  else if (i1==0){//spatial pattern
    bool b1=0;
    std::string s3=reverse(s2);
    if (s2==s3)
      b1=1;
    int pos=s.find(s2, 1);
    if (b1 && pos!=1)
      return -1;//symmetrical
    else
      return pos;
  }
}
//FIND CYCLE 2
int PatternsPlugin::find_cycle2(const std::vector<std::vector<int> >& p){
  int size=p.size();
  std::vector<int> state=p.back();
  for (int i=size-2; i>-1; --i)
    if (p[i]==state)
      return (size-1-i);
}  
//PERIODICITY
std::pair<int,int> PatternsPlugin::periodicity(const std::vector<std::vector<int> >& p){
  std::pair<int,int> period;
  int cells=p[0].size();//size of spatial grid
  int time=p.size();//size of temporal grid
  int t0=time/2;//we abandon the transient state
  //1-SPATIAL PERIODICITY
  std::vector<int> space(cells);
  for (int c=0; c<cells; ++c){
    //we have to translate the grid info into an integer
    std::vector<int> grid(time-t0);
    for (int t=t0; t<time; ++t)
      grid[t-t0]=p[t][c];
    std::string s1=vec_to_str(grid);
    int state=base_n_to_10(s1, 4);
    space[c]=state;
  }
  period.first=find_cycle1(space, 0);
  //2-TEMPORAL PERIODICITY
  std::vector<int> timing(time);
  for (int t=0; t<time; ++t){
    std::vector<int> t1=p[t];
    std::string s1=vec_to_str(t1);
    int state=base_n_to_10(s1, 4);
    timing[t]=state;
  }
  period.second=find_cycle2(p);
  return period;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//LOGIC_PROMS
std::vector<int> PatternsPlugin::logic_proms(int i, int p1, int p2){
  std::vector<int> encendido(6,0);
  int dimer=dimer_from_two(p1, p2);
  double energy_dimers;
  if (dimer==-1)
    energy_dimers=0.0;
  else
    energy_dimers=dim_bond_energy[dimer]+prot_energies[p1]+prot_energies[p2];  
    
  for (int j=0; j<6; ++j){//loop over all initial states
    switch(j){
    case 0:{// (0,0) no protein is expressed
      if (d_less(polymerase[i],0.0))//if the polymerase binds the promoter, the promoter will be on
	encendido[j]=1;
      break;
    }//case 0
    case 1:{// (0,1) only p2 present
      if (d_less(prot_prom_energies[i][p2],0.0) &&//if it binds to the promoter
	  d_less(prot_prom_energies[i][p2]+prot_energies[p2],polymerase[i]) &&//AND if union is stronger than the poly
	  prot_pol[i][p2]==1)//AND if, once bound to the promoter, activates it
	encendido[j]=1;
      else if (!d_less(prot_prom_energies[i][p2],0.0) ||//if it doesn't bind to the promoter OR
	       d_less(polymerase[i],prot_prom_energies[i][p2]+prot_energies[p2])){//if the polymerase is stronger
	if (d_less(polymerase[i],0.0))//then it depends on the polymerase
	  encendido[j]=1;
      }	    
      break;
    }//case 1
    case 2:{// (1,0) only p1 present (this is the same as case 1)
      if (d_less(prot_prom_energies[i][p1],0.0) && 
	  d_less(prot_prom_energies[i][p1]+prot_energies[p1],polymerase[i]) && 
	  prot_pol[i][p1]==1)
	encendido[j]=1;
      else if (!d_less(prot_prom_energies[i][p1],0.0) ||
	       d_less(polymerase[i],prot_prom_energies[i][p1]+prot_energies[p1])){
	if (d_less(polymerase[i],0.0))
	  encendido[j]=1;
      }	    
      break;
    }//case 2
    case 3:{// (1,1) p1 and p2 are present
      if (dimer != -1){//if they form a dimer
	if (d_less(dim_prom_energies[i][dimer],0.0) &&//if the dimer binds the promoter 
	    d_less(dim_prom_energies[i][dimer]+energy_dimers,polymerase[i]) &&
	    dim_pol[i][dimer]==1)//if union is stronger than polymerase AND activates promoter
	  encendido[j]=1;
	else if (!d_less(dim_prom_energies[i][dimer],0.0) ||//if it doesn't bind to the promoter OR
		 d_less(polymerase[i],dim_prom_energies[i][dimer]+energy_dimers)){//if the polymerase is stronger
	  if (d_less(polymerase[i],0.0))//then it depends on the polymerase
	    encendido[j]=1;
	}	    
      }//if dimer!=-1
      else{//no dimer, we compare both proteins
	double b1=std::min(prot_prom_energies[i][p1]+prot_energies[p1]-polymerase[i],0.0);
	double b2=std::min(prot_prom_energies[i][p2]+prot_energies[p2]-polymerase[i],0.0);
	if (d_less(b1,b2)){
	  if (d_less(prot_prom_energies[i][p1],0.0) && 
	      d_less(b1,0.0) && 
	      prot_pol[i][p1]==1)
	    encendido[j]=1;
	  else if (!d_less(prot_prom_energies[i][p1],0.0) ||
		   d_equal(b1,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}	    
	else if (d_less(b2,b1)){
	  if (d_less(prot_prom_energies[i][p2],0.0) && 
	      d_less(b2,0.0) && 
	      prot_pol[i][p2]==1)
	    encendido[j]=1;
	  else if (!d_less(prot_prom_energies[i][p2],0.0) ||
		   d_equal(b2,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}
	else if (d_equal(b1,b2)){
	  if (d_less(polymerase[i],0.0))
	    encendido[j]=1;
	}
      }
      break;
    }//case 3
    case 4:{// dimer AND  p1
      int dimer2=dimer_from_two(p1, p1);
      double energy_dimers2=dim_bond_energy[dimer2]+prot_energies[p1]+prot_energies[p1];
      if (dimer2==-1)
	energy_dimers2=0.0;
      if (dimer==-1 && dimer2==-1)//no dimer
	encendido[j]=encendido[3];
      else if (dimer==-1 && dimer2!=-1){//homodimer! competition between homodimer and p2
	double b1=std::min(prot_prom_energies[i][p2]+prot_energies[p2]-polymerase[i],0.0);
	double b2=std::min(dim_prom_energies[i][dimer2]+energy_dimers2-polymerase[i],0.0);
	if (d_less(b1,b2)){//p1 binds
	  if (d_less(prot_prom_energies[i][p2],0.0) && 
	      d_less(b1,0.0) && 
	      prot_pol[i][p2]==1)
	    encendido[j]=1;
	  else if (!d_less(prot_prom_energies[i][p2],0.0) ||
		   d_equal(b1,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}	    
	else if (d_less(b2,b1)){//dimer binds
	  if (d_less(dim_prom_energies[i][dimer2],0.0) && 
	      d_less(b2,0.0) && 
	      dim_pol[i][dimer2]==1)
	    encendido[j]=1;
	  else if (!d_less(dim_prom_energies[i][dimer2],0.0) ||
		   d_equal(b2,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}
	else if (d_equal(b1,b2)){//both bind
	  if (d_less(polymerase[i],0.0))
	    encendido[j]=1;
	}
      }
      else if (dimer!=-1){//competition dimer vs p1
	double b1=std::min(prot_prom_energies[i][p1]+prot_energies[p1]-polymerase[i],0.0);
	double b2=std::min(dim_prom_energies[i][dimer]+energy_dimers-polymerase[i],0.0);
	if (d_less(b1,b2)){//p1 binds
	  if (d_less(prot_prom_energies[i][p1],0.0) && 
	      d_less(b1,0.0) && 
	      prot_pol[i][p1]==1)
	    encendido[j]=1;
	  else if (!d_less(prot_prom_energies[i][p1],0.0) ||
		   d_equal(b1,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}	    
	else if (d_less(b2,b1)){//dimer binds
	  if (d_less(dim_prom_energies[i][dimer],0.0) && 
	      d_less(b2,0.0) && 
	      dim_pol[i][dimer]==1)
	    encendido[j]=1;
	  else if (!d_less(dim_prom_energies[i][dimer],0.0) ||
		   d_equal(b2,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}
	else if (d_equal(b1,b2)){//both bind
	  if (d_less(polymerase[i],0.0))
	    encendido[j]=1;
	}
      }
    }//case 4
    case 5:{// dimer AND  p2
      int dimer2=dimer_from_two(p2, p2);
      double energy_dimers2=dim_bond_energy[dimer2]+prot_energies[p2]+prot_energies[p2];
      if (dimer2==-1)
	energy_dimers2=0.0;
      if (dimer==-1 && dimer2==-1)//no dimer
	encendido[j]=encendido[3];
      else if (dimer==-1 && dimer2!=-1){//homodimer! competition between homodimer and p1
	double b1=std::min(prot_prom_energies[i][p1]+prot_energies[p1]-polymerase[i],0.0);
	double b2=std::min(dim_prom_energies[i][dimer2]+energy_dimers2-polymerase[i],0.0);
	if (d_less(b1,b2)){//p1 binds
	  if (d_less(prot_prom_energies[i][p1],0.0) && 
	      d_less(b1,0.0) && 
	      prot_pol[i][p1]==1)
	    encendido[j]=1;
	  else if (!d_less(prot_prom_energies[i][p1],0.0) ||
		   d_equal(b1,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}	    
	else if (d_less(b2,b1)){//dimer binds
	  if (d_less(dim_prom_energies[i][dimer2],0.0) && 
	      d_less(b2,0.0) && 
	      dim_pol[i][dimer2]==1)
	    encendido[j]=1;
	  else if (!d_less(dim_prom_energies[i][dimer2],0.0) ||
		   d_equal(b2,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}
	else if (d_equal(b1,b2)){//both bind
	  if (d_less(polymerase[i],0.0))
	    encendido[j]=1;
	}
      }
      else if (dimer!=-1){//competition dimer vs p2
	double b1=std::min(prot_prom_energies[i][p2]+prot_energies[p2]-polymerase[i],0.0);
	double b2=std::min(dim_prom_energies[i][dimer]+energy_dimers-polymerase[i],0.0);
	if (d_less(b1,b2)){//p2 binds
	  if (d_less(prot_prom_energies[i][p2],0.0) && 
	      d_less(b1,0.0) && 
	      prot_pol[i][p2]==1)
	    encendido[j]=1;
	  else if (!d_less(prot_prom_energies[i][p2],0.0) ||
		   d_equal(b1,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}	    
	else if (d_less(b2,b1)){//dimer binds
	  if (d_less(dim_prom_energies[i][dimer],0.0) && 
	      d_less(b2,0.0) && 
	      dim_pol[i][dimer]==1)
	    encendido[j]=1;
	  else if (!d_less(dim_prom_energies[i][dimer],0.0) ||
		   d_equal(b2,0.0)){
	    if (d_less(polymerase[i],0.0))
	      encendido[j]=1;
	  }
	}
	else if (d_equal(b1,b2)){//both bind
	  if (d_less(polymerase[i],0.0))
	    encendido[j]=1;
	}
      }
    }//case 5
    }//switch
  }//for j
  return encendido;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//COARSENING
std::map<std::vector<int>, std::vector<int> > PatternsPlugin::coarsening(int p1, int p2){
  //OUTPUT
  //Given a set of two proteins, COARSENING collapses the combinations os promoters that yield the same logic function
  //It returns a std::map that holds, as a key, the logic function of the equivalence class, and as a value, the promoters that
  //are members of the equivalence class.

  //OUTPUT OBJECT
  std::map<std::vector<int>, std::vector<int> > coarse_proms;
  
  //LOOP OVER ALL PROMOTERS
  for (int i=0; i<16; ++i){
    std::vector<int> encendido=logic_proms(i, p1, p2);//this will hold what happens to the promoter i at each initial state
    coarse_proms[encendido].push_back(i);
  }
  
  return coarse_proms;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//TRANSLATION
std::vector<int> PatternsPlugin::translation(const std::map<std::vector<int>, std::vector<int> >::const_iterator& p1, const std::map<std::vector<int>, std::vector<int> >::const_iterator& p2)
{
  //OUTPUT
  //Given the logic functions of all promoters, it returns the combined logic function. It is a simple translation function.
  //For example, if from state i prom1 goes to 1, prom2 goes to 0, and prom3 goes to 1, logic_function[i]=5.

  //ARGUMENTS
  //(1) (const std::map<std::vector<int>, std::vector<int> >&) COARSE_PROMS: a map holding the logic functions of all promoters
  //     (see COARSENING function)
  //(2) (const std::vector<std::map<std::vector<int>, std::vector<int> >::const_iterator>&) PROM: iterator over COARSE_PROMS

  std::vector<int> logic_function(6,0);
  for (int i=0; i<6; ++i){
    if (p1->first[i]==0 && p2->first[i]==1)
      logic_function[i]=1;
    else if (p1->first[i]==1 && p2->first[i]==0)
      logic_function[i]=2;
    else if (p1->first[i]==1 && p2->first[i]==1)
      logic_function[i]=3;
  }
  return logic_function;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//LOGIC_FUNCTION
std::vector<int> PatternsPlugin::logic_function(int p1, int p2, int pr1, int pr2){
  std::vector<std::vector<int> > logic_vector(2);
  logic_vector[0]=logic_proms(pr1, p1, p2);//this will hold what happens to the promoter i at each initial state
  logic_vector[1]=logic_proms(pr2, p1, p2);//this will hold what happens to the promoter i at each initial state
  int dimer=dimer_from_two(p1, p2);
  std::vector<int> empty(6,0);
  if (p1==0)
    logic_vector[0]=empty;
  if (p2==0)
    logic_vector[1]=empty;
  
  std::vector<int> logic_function(6,0);
  for (int i=0; i<6; ++i){
    if (logic_vector[0][i]==0 && logic_vector[1][i]==1)
      logic_function[i]=1;
    else if (logic_vector[0][i]==1 && logic_vector[1][i]==0)
      logic_function[i]=2;
    else if (logic_vector[0][i]==1 && logic_vector[1][i]==1)
      logic_function[i]=3;
  }
  if (dimer==-1)
    logic_function.push_back(0);
  else
    logic_function.push_back(1);
  return logic_function;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//GENOTYPE_STR_TO_VEC
void PatternsPlugin::genotype_str_to_vec(const std::string& genotype, int& pr1, int& pr2, int& p1, int& p2){
  pr1=bintodec(genotype.substr(0, 4));
  p1=prot_gen[bintodec(genotype.substr(4, 16))];
  pr2=bintodec(genotype.substr(20, 4));
  p2=prot_gen[bintodec(genotype.substr(24, 16))];
  return;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//GENOTYPE_VEC_TO_STR
std::string PatternsPlugin::genotype_vec_to_str(int pr1, int pr2, int c1, int c2){
  std::string genotype;
  genotype += dectobin(pr1,4);
  genotype += dectobin(c1, 16);
  genotype += dectobin(pr2,4);
  genotype += dectobin(c2, 16);
  
  return genotype;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
std::vector<int> PatternsPlugin::str_logic_function(const std::string& genotype){
  int pr1, pr2, p1, p2;
  genotype_str_to_vec(genotype, pr1, pr2, p1, p2);
  return logic_function(p1, p2, pr1, pr2);
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//ADD MORPHOGEN TO THE CELLS
std::vector<int> morphogen_creator(int cells, int t, const std::string& gradient){
  //La idea es que el morfogeno se distribuya como un decaimiento
  //exponencial, de tal forma que la célula 0 tenga 1 de morfógeno
  // y la ceĺula input.size()-1 tenga thr (bajo)
  double thr=0.005;
  if (gradient=="exp"){
    //el morfógeno es la prot 1
    std::vector<int> m1(cells);
    int rsize=cells-1;
    double k=-std::log(thr)/(rsize);
    for (int c=0; c<cells; ++c){
      double cmorph=std::exp(-k*c);
      int period=1.0/cmorph;
      if (((t+1)%period)==0) 
	m1[c]=1;//receives morphogen
    }
    return m1;
  }
  else if (gradient=="gauss"){
    //el morfógeno es la prot 1
    std::vector<int> m1(cells);
    int rsize=cells-1;
    double k=-std::log(thr)/(rsize*rsize);
    for (int c=0; c<cells; ++c){
      double cmorph=std::exp(-k*c*c);
      int period=1.0/cmorph;
      if (((t+1)%period)==0) 
	m1[c]=1;//receives morphogen
    }
    return m1;
  }
  else if (gradient=="linear"){
    //el morfógeno es la prot 1
    std::vector<int> m1(cells);
    int rsize=cells-1;
    double k=(thr-1.0)/rsize;
    for (int c=0; c<cells; ++c){
      double cmorph=1.0+k*c;
      int period=1.0/cmorph;
      //std::cout << c << "\t" << cmorph << "\t" << period << std::endl;
      if (((t+1)%period)==0) 
	m1[c]=1;//receives morphogen
    }
    return m1;
  }
  else if (gradient=="linearP2"){
    //el morfógeno es la prot 1
    std::vector<int> m1(cells);
    int rsize=cells-1;
    double k=(thr-1.0)/rsize;
    for (int c=0; c<cells; ++c){
      double cmorph=1.0+k*c;
      int period=1.0/cmorph;
      //std::cout << c << "\t" << cmorph << "\t" << period << std::endl;
      if (((t+1)%period)==0) 
	m1[c]=2;//receives morphogen
    }
    return m1;
  }
  else if (gradient=="point"){
    //el morfógeno es la prot 1
    std::vector<int> m1(cells);
    if (t==0){
      m1[cells/2]=1;
      if (cells%2==0)//if grid is even
	m1[cells/2-1]=1;//prot 1 as input in two initial cells
    }
    return m1;
  }
  else if (gradient=="pointleft"){
    //el morfógeno es la prot 1
    std::vector<int> m1(cells);
    if (t==0){
      m1[0]=1;
    }
    return m1;
  }
  else if (gradient=="pointleftP2"){
    //el morfógeno es la prot 1
    std::vector<int> m1(cells);
    if (t==0){
      m1[0]=2;
    }
    return m1;
  }
  else if (gradient=="linear2"){
    //el morfógeno es la prot 1
    std::vector<int> m1(cells);
    int rsize=cells-1;
    double k=(thr-1.0)/rsize;
    for (int c=0; c<cells; ++c){
      double cmorph=1.0+k*c;
      int period=c+1;
      if (((t+1)%period)==0) 
	m1[c]=1;//receives morphogen
    }
    return m1;
  }
  else{
    std::vector<int> m1(cells);
    return m1;
  }
}
