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

typedef std::vector<int> CA;

int main(int argc, char* argv[]){

  //LOAD TOYLIFE DATA
  PatternsPlugin Model;
  int total_proteins=Model.prot_energies.size();
  
  //HOW MANY GENES CODE FOR A GIVEN PROTEIN
  std::vector<int> prots_gen(total_proteins,0);
  std::string name="../data/prot_abundance.txt";
  std::ifstream file_prots_gen(name.c_str());
  check_file(file_prots_gen, name);
  int b1;
  for (int i=0; i<total_proteins; i++)
    file_prots_gen >> b1 >> prots_gen[i];
  file_prots_gen.close();

  //LOOP
  std::map<CA, double> mapa;
  std::map<std::vector<int>, double> GRN_ab;
  std::map<std::vector<int>, CA> mapGRN0;//NOTHING DIFFUSES
  std::map<std::vector<int>, CA> mapGRN1;//ONLY P1 DIFFUSES
  std::map<std::vector<int>, CA> mapGRN2;//ONLY P2 DIFFUSES
  std::map<std::vector<int>, CA> mapGRN3;//BOTH DIFFUSE, DIMER DOES NOT
  std::map<std::vector<int>, CA> mapGRN4;//BOTH DIFFUSE, DIMER AS WELL
  double suma=0;

  //LOOP OVER ALL PROTEIN PAIRS (THE PAIR P1-P2 IS NOT THE SAME AS P2-P1)
  for (int p1=0; p1<total_proteins; ++p1){
    std::cout << p1 << std::endl;
    for (int p2=0; p2<total_proteins; ++p2){
      //if (!(p1==0 || p2==0)) continue;
      int dimer=Model.dimer_from_two(p1, p2);//obtain dimer from both proteins
      std::map<std::vector<int>, std::vector<int> > prom=Model.coarsening(p1, p2);//group similar promoters according to their activity
      //LOOP OVER ALL PROMOTER PAIRS
      for (auto pr1=prom.begin(); pr1!=prom.end(); ++pr1){
	auto pr2=prom.begin();
	for (; pr2!=prom.end(); ++pr2){
	  std::vector<int> GRN=Model.translation(pr1, pr2);//compute GRN function with these promoters 
	  //STATE 0 IS (0,0) (NOTHING IS EXPRESSED)
	  //STATE 1 IS (0,1) (P2 IS EXPRESSED)
	  //STATE 2 IS (1,0) (P1 IS EXPRESSED)
	  //STATE 3 IS (1,1) (P1 AND P2 ARE EXPRESSED)
	  
	  //ADAPT GRN FUNCTION TO FULFILL GENERAL CONVENTION
	  //WHEN EITHER P1 OR P2 ARE PROTEIN 0, THIS IS A NON-FUNCTIONAL PROTEIN, SO WE ASSUME NOTHING IS EXPRESSED
	  //WE THEN CHANGE THE CORRESPONDING STATES:
	  //IF P1=0, WE CHANGE STATE 2 TO 0, AND STATE 3 TO 1
	  if (p1==0 && p2!=0){
	    int a1=prots_gen[p1];
	    int a2=prots_gen[p2];
	    int b1=pr1->second.size();
	    int b2=pr2->second.size();
	    double a=a1*b1;//number of genes that can act as p1 here
	    double b=a2*b2;//number of genes that can act as p2 here
	    double d1=a*b;
	    for (int l=0; l<GRN.size(); ++l){
	      if (GRN[l]==2)
		GRN[l]=0;
	      else if (GRN[l]==3)
		GRN[l]=1;
	    }
	  }
	  //IF P2=0, WE CHANGE STATE 1 TO 1, AND STATE 3 TO 2
	  else if (p2==0 && p1!=0){
	    int a1=prots_gen[p1];
	    int a2=prots_gen[p2];
	    int b1=pr1->second.size();
	    int b2=pr2->second.size();
	    double a=a1*b1;//number of genes that can act as p1 here
	    double b=a2*b2;//number of genes that can act as p2 here
	    double d1=a*b;
	    for (int l=0; l<GRN.size(); ++l){
	      if (GRN[l]==1)
		GRN[l]=0;
	      else if (GRN[l]==3)
		GRN[l]=2;
	    }
	  }
	  //IF BOTH PROTEINS ARE 0, WE CHANGE ALL STATES TO 0
	  else if (p1==0 && p2==0){
	    for (int l=0; l<GRN.size(); ++l)
	      GRN[l]=0;
	  }
	    
	  //DISTINGUISH GRN FUNCTIONS WITH AND WITHOUT DIMER
	  if (dimer==-1)
	    GRN.push_back(0);
	  else
	    GRN.push_back(1);
	  //COUNT GRN FUNCTIONS
	  int a1=prots_gen[p1];
	  int a2=prots_gen[p2];
	  int b1=pr1->second.size();
	  int b2=pr2->second.size();
	  double a=a1*b1;//number of genes that can act as p1 here
	  double b=a2*b2;//number of genes that can act as p2 here
	  double d1=a*b;
	  suma +=d1;//update total count
	  
	  GRN_ab[GRN]+=d1;//update count for this particular GRN

	  //NOW WE CONSIDER FOUR DIFFUSION SCENARIOS IN ORDER TO TRANSLATE FROM GRN TO CELLULAR AUTOMATA	  
	  //(0) NOTHING DIFFUSES,
	  //(1) ONLY P1 DIFFUSES,
	  //(2) ONLY P2 DIFFUSES,
	  //(3) BOTH PROTEINS DIFFUSE, BUT NOT THE DIMER AND
	  //(4) EVERYTHING DIFFUSES

	  //FOR THIS PAPER WE WILL ONLY CONSIDER SCENARIO (1)

	  //THERE ARE 4^3 INPUT SCENARIOS
	  int L=int_pow(4, 3);

	  //(0) NOTHING DIFFUSES
	  CA caso0(L,0);
	  for (int k=0; k<L; ++k){
	    std::vector<int> input=base_10_to_n(k, 4);
	    int iset=input[1];
	    caso0[k]=GRN[iset];
	  }//for k
	  mapa[caso0]+=d1;
	  mapGRN0[GRN]=caso0;
	  
	  //(1) ONLY P1 DIFFUSES
	  CA caso1(L,0);
	  for (int k=0; k<L; ++k){
	    std::vector<int> input=base_10_to_n(k, 4);
	    int iset=input[1];
	    switch(iset){//0 there is nothing, 1 there is P2, 2 there is P1, 3 there are both (or the dimer), 4 dimer + p1, 5 dimer + p2
	    case 0:{//nothing in the cell
	      if (dimer==-1){//no dimer
		if (input[0]==2 || input[0]==3 || input[2]==2 || input[2]==3)//p1 somewhere
		  caso1[k]=GRN[2];//p1
		else
		  caso1[k]=GRN[0];//nothing
	      }
	      else{//dimer
		if (input[0]==2 || input[2]==2)//p1 somewhere
		  caso1[k]=GRN[2];//p1
		else
		  caso1[k]=GRN[0];//nothing
	      }
	      break;
	    }//case 0
	    case 1:{//only p2
	      if (dimer==-1){//no dimer
		if (input[0]==2 || input[0]==3 || input[2]==2 || input[2]==3)//p1 somewhere
		  caso1[k]=GRN[3];//p1 plus p2
		else
		  caso1[k]=GRN[1];//only p2
	      }
	      else{//dimer
		if (input[0]==2 || input[2]==2)//p1 somewhere
		  caso1[k]=GRN[3];//p1 plus p2
		else
		  caso1[k]=GRN[1];//only p2
	      }
	      break;
	    }//case 1
	    case 2:{//only p1
	      caso1[k]=GRN[2];//no matter what happens, nothing will diffuse into the cell
	      break;
	    }//case 2
	    case 3:{//p1 and p2
	      if (dimer==-1){//no dimer
		caso1[k]=GRN[3];//no matter what happens, nothing will diffuse into the cell
	      }
	      else{//dimer
		if (input[0]==2 || input[2]==2)//p2 somewhere
		  caso1[k]=GRN[4];//dimer and p1
		else
		  caso1[k]=GRN[3];
	      }
	      break;
	    }//case 3
	    }//switch
	  }//for k
	  mapa[caso1]+=d1;
	  mapGRN1[GRN]=caso1;

	  //(2) ONLY P2 DIFFUSES
	  CA caso2(L,0);
	  for (int k=0; k<L; ++k){
	    std::vector<int> input=base_10_to_n(k, 4);
	    int iset=input[1];
	    switch(iset){
	    case 0:{//nothing in the cell
	      if (dimer==-1){//no dimer
		if (input[0]==1 || input[0]==3 || input[2]==1 || input[2]==3)//p1 somewhere
		  caso2[k]=GRN[1];//p1
		else
		  caso2[k]=GRN[0];//nothing
	      }
	      else{//dimer
		if (input[0]==1 || input[2]==1)//p1 somewhere
		  caso2[k]=GRN[1];//p1
		else
		  caso2[k]=GRN[0];//nothing
	      }
	      break;
	    }//case 0
	    case 2:{//only p1
	      if (dimer==-1){//no dimer
		if (input[0]==1 || input[0]==3 || input[2]==1 || input[2]==3)//p1 somewhere
		  caso2[k]=GRN[3];//p1 plus p2
		else
		  caso2[k]=GRN[2];//only p2
	      }
	      else{//dimer
		if (input[0]==1 || input[2]==1)//p1 somewhere
		  caso2[k]=GRN[3];//p1 plus p2
		else
		  caso2[k]=GRN[2];//only p2
	      }
	      break;
	    }//case 1
	    case 1:{//only p2
	      caso2[k]=GRN[1];//no matter what happens, nothing will diffuse into the cell
	      break;
	    }//case 2
	    case 3:{//p1 and p2
	      if (dimer==-1){//no dimer
		caso2[k]=GRN[3];//no matter what happens, nothing will diffuse into the cell
	      }
	      else{//dimer
		if (input[0]==1 || input[2]==1)//p2 somewhere
		  caso2[k]=GRN[5];//dimer and p1
		else
		  caso2[k]=GRN[3];
	      }
	      break;
	    }//case 3
	    }//switch
	  }//for k
	  mapa[caso2]+=d1;
	  mapGRN2[GRN]=caso2;
	  
	  //(3) BOTH DIFFUSE - BUT NOT THE DIMER
	  CA caso3(L,0);
	  for (int k=0; k<L; ++k){
	    std::vector<int> input=base_10_to_n(k, 4);
	    int iset=input[1];
	    switch(iset){
	    case 0:{//nothing in the cell
	      if (dimer==-1){//no dimer
		int i1=0;
		if ((input[0]==1 && input[2]==1) ||
		    (input[0]==0 && input[2]==1) ||
		    (input[0]==1 && input[2]==0))
		  i1=1;//only p2 in the cell
		else if ((input[0]==2 && input[2]==2) ||
			 (input[0]==0 && input[2]==2) ||
			 (input[0]==2 && input[2]==0))
		  i1=2;//only p1 in the cell
		else if ((input[0]==1 && input[2]==2) || 
			 (input[0]==2 && input[2]==1) || 
			 input[0]==3  || input[2]==3)
		  i1=3;//both proteins in the cell
		caso3[k]=GRN[i1];
	      }
	      else{//dimer
		int i1=0;
		if ((input[0]==1 && input[2]==1) ||
		    (input[0]==0 && input[2]==1) ||
		    (input[0]==1 && input[2]==0) ||
		    (input[0]==1 && input[2]==3) ||
		    (input[0]==3 && input[2]==1))
		  i1=1;//only p2 in the cell
		else if ((input[0]==2 && input[2]==2) ||
			 (input[0]==0 && input[2]==2) ||
			 (input[0]==2 && input[2]==0) ||
			 (input[0]==2 && input[2]==3) ||
			 (input[0]==3 && input[2]==2))
		  i1=2;//only p1 in the cell
		else if ((input[0]==1 && input[2]==2) || 
			 (input[0]==2 && input[2]==1))
		  i1=3;//both proteins in the cell
		caso3[k]=GRN[i1];
	      }
	      break;
	    }//case 0
	    case 1:{//only p2
	      if (dimer==-1){//no dimer
		if (input[0]==2 || input[0]==3 || input[2]==2 || input[2]==3)//p1 somewhere
		  caso3[k]=GRN[3];//p1 plus p2
		else
		  caso3[k]=GRN[1];//only p2
	      }
	      else{//dimer
		if (input[0]==2 || input[2]==2)//p1 somewhere
		  caso3[k]=GRN[3];//p1 plus p2
		else
		  caso3[k]=GRN[1];//only p2
	      }
	      break;
	    }//case 1
	    case 2:{//only p1
	      if (dimer==-1){//no dimer
		if (input[0]==1 || input[0]==3 || input[2]==1 || input[2]==3)//p1 somewhere
		  caso3[k]=GRN[3];//p1 plus p2
		else
		  caso3[k]=GRN[2];//only p2
	      }
	      else{//dimer
		if (input[0]==1 || input[2]==1)//p1 somewhere
		  caso3[k]=GRN[3];//p1 plus p2
		else
		  caso3[k]=GRN[2];//only p2
	      }
	      break;
	    }//case 2
	    case 3:{//p1 and p2
	      if (dimer==-1){//no dimer
		caso3[k]=GRN[3];//no matter what happens, nothing will diffuse into the cell
	      }
	      else{//dimer
		if ((input[0]==1 && input[2]==2) ||
		    (input[0]==2 && input[2]==1))//dimer is formed again
		  caso3[k]=GRN[3];
		else if (input[0]==1 || input[2]==1)//p2 somewhere
		  caso3[k]=GRN[5];//dimer and p2
		else if (input[0]==2 || input[2]==2)//p1 somewhere
		  caso3[k]=GRN[4];//dimer and p1
	      }
	      break;
	    }//case 3
	    }//switch
	  }//for k
	  mapa[caso3]+=d1;
	  mapGRN3[GRN]=caso3;

	  //(4) BOTH DIFFUSE, DIMER AS WELL
	  CA caso4(L,0);
	  for (int k=0; k<L; ++k){
	    std::vector<int> input=base_10_to_n(k, 4);
	    int iset=input[1];
	    switch(iset){
	    case 0:{//nothing in the cell
	      if (dimer==-1){//no dimer
		int i1=0;
		if ((input[0]==1 && input[2]==1) ||
		    (input[0]==0 && input[2]==1) ||
		    (input[0]==1 && input[2]==0))
		  i1=1;//only p2 in the cell
		else if ((input[0]==2 && input[2]==2) ||
			 (input[0]==0 && input[2]==2) ||
			 (input[0]==2 && input[2]==0))
		  i1=2;//only p1 in the cell
		else if ((input[0]==1 && input[2]==2) || 
			 (input[0]==2 && input[2]==1) || 
			 input[0]==3  || input[2]==3)
		  i1=3;//both proteins in the cell
		caso4[k]=GRN[i1];
	      }
	      else{//dimer
		int i1=0;
		if ((input[0]==1 && input[2]==1) ||
		    (input[0]==0 && input[2]==1) ||
		    (input[0]==1 && input[2]==0))
		  i1=1;//only p2 in the cell
		else if ((input[0]==2 && input[2]==2) ||
			 (input[0]==0 && input[2]==2) ||
			 (input[0]==2 && input[2]==0))
		  i1=2;//only p1 in the cell
		else if ((input[0]==1 && input[2]==2) || 
			 (input[0]==2 && input[2]==1) ||
			 (input[0]==0 && input[2]==3) ||
			 (input[0]==3 && input[2]==0) ||
			 (input[0]==3 && input[2]==3))
		  i1=3;//both proteins in the cell
		else if ((input[0]==1 && input[2]==3) ||
			 (input[0]==3 && input[2]==1))
		  i1=5;
		else if ((input[0]==2 && input[2]==3) ||
			 (input[0]==3 && input[2]==2))
		  i1=4;
		caso4[k]=GRN[i1];
	      }
	      break;
	    }//case 0
	    case 1:{//only p2
	      if (dimer==-1){//no dimer
		if (input[0]==2 || input[0]==3 || input[2]==2 || input[2]==3)//p1 somewhere
		  caso4[k]=GRN[3];//p1 plus p2
		else
		  caso4[k]=GRN[1];//only p2
	      }
	      else{//dimer
		if ((input[0]==2 && input[2]<3) ||
		    (input[0]<3 && input[2]==2))//p1 somewhere
		  caso4[k]=GRN[3];//p1 plus p2
		else if ((input[0]!=2 && input[2]==3) ||
			 (input[0]==3 && input[2]!=2))//p2 and dimer
		  caso4[k]=GRN[5];
		else if ((input[0]==2 && input[2]==3) ||
			 (input[0]==3 && input[2]==2))//p1 and dimer
		  caso4[k]=GRN[3];//dimer is formed
		else
		  caso4[k]=GRN[1];//only p2
	      }
	      break;
	    }//case 1
	    case 2:{//only p1
	      if (dimer==-1){//no dimer
		if (input[0]==1 || input[0]==3 || input[2]==1 || input[2]==3)//p1 somewhere
		  caso4[k]=GRN[3];//p1 plus p2
		else
		  caso4[k]=GRN[2];//only p2
	      }
	      else{//dimer
		if ((input[0]==1 && input[2]<3) ||
		    (input[0]<3 && input[2]==1))//p2 somewhere
		  caso4[k]=GRN[3];//p1 plus p2
		else if ((input[0]!=1 && input[2]==3) ||
			 (input[0]==3 && input[2]!=1))//p1 and dimer
		  caso4[k]=GRN[4];
		else if ((input[0]==1 && input[2]==3) ||
			 (input[0]==3 && input[2]==1))//p2 and dimer
		  caso4[k]=GRN[3];//dimer is formed
		else
		  caso4[k]=GRN[2];//only p2
	      }
	      break;
	    }//case 2
	    case 3:{//p1 and p2
	      if (dimer==-1){//no dimer
		caso4[k]=GRN[3];//no matter what happens, nothing will change the state
	      }
	      else{//dimer
		if ((input[0]==1 && input[2]==2) ||
		    (input[0]==2 && input[2]==1))//dimer is formed again
		  caso4[k]=GRN[3];
		else if (input[0]==1 || input[2]==1)//p2 somewhere
		  caso4[k]=GRN[5];//dimer and p2
		else if (input[0]==2 || input[2]==2)//p1 somewhere
		  caso4[k]=GRN[4];//dimer and p1
	      }
	      break;
	    }//case 3
	    }//switch
	  }//for k
	  mapa[caso4]+=d1;
	  mapGRN4[GRN]=caso4;
	}//pr2
      }//pr1
    }//p2
  }//p1

  //CHECK COUNT IS CORRECT (we have considered every case)
  std::cout << suma << "\t" << std::pow(2, 40.0) << std::endl;

  //OUTPUT

  //WRITE CELLULAR AUTOMATA
  std::ofstream out("results/cellular_automata.txt");
  int id1=0;
  std::map<std::vector<int>, int> caid;
  for (auto it=mapa.begin(); it!=mapa.end(); ++it){
    std::string s1=vec_to_str(it->first);
    out << s1 << "\t" << id1 << std::endl;
    caid[it->first]=id1;
    id1++;
  }
  out.close();
  //WRITE ABUNDANCE OF GRNs
  out.open("results/GRN_abundance.txt");
  for (auto it=GRN_ab.begin(); it!=GRN_ab.end(); ++it)
    out << vec_to_str(it->first) << "\t" << it->second << std::endl;
  out.close();
  //WRITE THE CELLULAR AUTOMATA THAT CORRESPOND TO EACH GRN IN EACH DIFFUSION SCENARIO
  // out.open("results/GRN_ca_0.txt");
  // for (auto it=mapGRN0.begin(); it!=mapGRN0.end(); ++it)
  //   out << vec_to_str(it->first) << "\t" << vec_to_str(it->second) << std::endl;
  // out.close(); 
  out.open("results/GRN_ca_1.txt");
  for (auto it=mapGRN1.begin(); it!=mapGRN1.end(); ++it)
    out << vec_to_str(it->first) << "\t" << vec_to_str(it->second) << std::endl;
  out.close();
  // out.open("results/GRN_ca_2.txt");
  // for (auto it=mapGRN2.begin(); it!=mapGRN2.end(); ++it)
  //   out << vec_to_str(it->first) << "\t" << vec_to_str(it->second) << std::endl;
  // out.close();
  // out.open("results/GRN_ca_3.txt");
  // for (auto it=mapGRN3.begin(); it!=mapGRN3.end(); ++it)
  //   out << vec_to_str(it->first) << "\t" << vec_to_str(it->second) << std::endl;
  // out.close();
  // out.open("results/GRN_ca_4.txt");
  // for (auto it=mapGRN4.begin(); it!=mapGRN4.end(); ++it)
  //   out << vec_to_str(it->first) << "\t" << vec_to_str(it->second) << std::endl;
  // out.close();

  //ABUNDANCES OF CELLULAR AUTOMATA (SCENARIO 1)
  std::map<int,double> Cab;
  for (auto it=mapGRN1.begin(); it!=mapGRN1.end(); ++it){
    auto GRN=it->first;
    auto castr=it->second;
    auto id=caid.at(castr);
    Cab[id]+=GRN_ab.at(GRN);
  }
  out.open("results/ca_abundance.txt");
  for (auto it=Cab.begin(); it!=Cab.end(); ++it)
    out << it->first << "\t" << it->second << std::endl;
  out.close();
}
