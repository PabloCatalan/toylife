#include "toy_plugin.h"
#include "helper_functions.h"

//STRUCTURE FUNCTIONS
Dim::Dim():
  id(-1), p1(0), p2(0){}
Dim::Dim(const int Id, const int P1, const int P2):
  id(Id), p1(P1), p2(P2){}
bool Dim::empty() const{
  return (id==-1);
    }
OWM::OWM():
  prot(-1), dim(), met(){}
OWM::OWM(const Prot p1, const Dim& d1, const Met& m1):
  prot(p1), dim(d1), met(m1){}
bool OWM::empty() const{
  return (prot==-1 && dim.empty() && met.empty());
}
Dmet::Dmet():
  eg(0.0), seq(), pos(0){}
Dmet::Dmet(const double E, const std::string& S, const int P):
  eg(E), seq(S), pos(P){}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
ToyPlugin::ToyPlugin(){  
  //INITIALIZE ALL OF TOYLIFE'S LOOKUP TABLES
  init_gen(prot_gen);
  init_prot(prot_perimeters, prot_energies, prot_prom_energies, prot_pol);
  init_dim(dim_perim, dim_bond_energy, dim_prom_energies, dim_pol);
  init_polymerase(polymerase);  
  simplified_mets(mets, not_broken, env);
  neighbors(neighbors_prom, neighbors_coding);
  init_dim_met(dim_met);
  init_prot_met(prot_met);
  init_prot_breaking(prot_breaking);  
}
ToyPlugin::~ToyPlugin(){
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//INIT GEN
void ToyPlugin::init_gen(std::vector<int>& prot_gen){//reads to table which protein is given by which gene
  std::string name="data/protein_gene.txt";
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
void ToyPlugin::init_polymerase(std::vector<double>& polymerase){
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
void ToyPlugin::init_prot(std::vector<int>& perimeters, std::vector<double>& energies, std::vector<std::vector<double> >& prot_prom_energies, std::vector<std::vector<bool> >& prot_pol){//reads to table the perimeter and energies associated to each protein
  //NUMBER OF PROTEINS
  int total_proteins=0;
  std::string name="data/proteins_with_perimeter.txt";
  std::ifstream in(name.c_str());
  check_file(in,name);
  std::string lee;
  while (std::getline(in, lee))
    ++total_proteins;
  perimeters.resize(total_proteins);
  energies.resize(total_proteins);
  in.close();
  //PERIMETERS
  name="data/proteins_with_perimeter.txt";
  in.open(name.c_str());
  check_file(in, name);
  int id, per;
  while (in >> id >> per)
    perimeters[id]=per;
  in.close();
  perimeters[0]=-1;
  //ENERGIES
  name="data/proteins.txt";
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
    std::string name="promoters/proteins_bind_prom"+ss.str()+".txt";
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
void ToyPlugin::init_dim(std::map<std::pair<int, int>, int>& dim_perim, std::vector<double>& dim_energy, std::vector<std::vector<double> >& dim_prom_energies, std::vector<std::vector<bool> >& dim_pol){//reads to table the information relating to dimers
  int total_dimers=0;
  std::string name="data/dimers.txt";
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
    std::string name="promoters/dimers_bind_prom"+ss.str()+".txt";
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
//SIMPLIFIED METS
void ToyPlugin::simplified_mets(std::vector<std::string>& mets, std::map<std::string,int>& not_broken, std::vector<std::string>& env){
  std::ifstream in("data/metabolites.txt");
  check_file(in, "data/metabolites.txt");
  std::string s1;
  int nb;
  while (in >> s1 >> nb){
    mets.push_back(s1);
    if (nb)
      not_broken[s1]=1;
  }
  in.close();
  
  for (int j=0; j<mets.size(); ++j)
    if (!not_broken.count(mets[j]))
      env.push_back(mets[j]);

  return;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//NEIGHBORS
void ToyPlugin::neighbors(std::vector<std::vector<int> >& neighbors_prom, std::vector<std::vector<int> >& neighbors_coding){
  neighbors_prom=std::vector<std::vector<int> >(16, std::vector<int>(4));
  std::ifstream file_prom("data/neighbors_prom.txt");
  check_file(file_prom, "data/neighbors_prom.txt");
  for (int i=0; i<16; ++i)
    for (int g=0; g<4; ++g)
      file_prom >> neighbors_prom[i][g];
  file_prom.close();

  neighbors_coding=std::vector<std::vector<int> >(65536, std::vector<int>(16));
  std::ifstream file_coding("data/neighbors_coding.txt");
  check_file(file_coding, "data/neighbors_coding.txt");
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
//DIM_MET
void ToyPlugin::init_dim_met(mapa_dmet& dim_met){//key: Met; value: map -> key (Id of the dimer), value: Dmet, structure with the following members: (eg) binding energy, (seq) sequence left for proteins to bind to, (pos) how many sugars are to the left of the center 
  for (int i=4; i<9; ++i){
    std::map<int, std::vector<int> > mapa;
    for (int j=0; j<pow(2,i); ++j){
      std::string s1=dectobin(j,i);
      std::string s2=reverse(s1);
      int h=bintodec(s2);
      if (j<h || j==h){
	mapa[j].push_back(j);
	mapa[j].push_back(h);
      }
    }
    std::stringstream ss;
    ss << i;
    for (std::map<int, std::vector<int> >::const_iterator it=mapa.begin(); it!=mapa.end(); ++it){//Para todas las secuencias
      std::stringstream ss_met;
      ss_met << it->second[0];
      std::string name="dim_metabolism/"+ss.str()+"_"+ss_met.str()+".txt";
      std::ifstream file_dimeros_mets(name.c_str());
      check_file(file_dimeros_mets, name);
      int index;
      std::string met_id=dectobin(it->second[0],i);
      while (file_dimeros_mets >> index){
	double energy;
	std::string sequence;
	int binding;
	int position;
	file_dimeros_mets >> energy >> sequence >> position >> binding;
	dim_met[std::make_pair(met_id, index)]=Dmet(energy, sequence, position);
      }
      file_dimeros_mets.close();
    }//for mets
  }//for sizes
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//PROTEINS_METABOLITES
void ToyPlugin::init_prot_met(mapa_pmet& prot_met){
  for (int i=4; i<9; ++i){
    std::map<int, std::vector<int> > mapa;
    for (int j=0; j<pow(2,i); ++j){
      std::string s1=dectobin(j,i);
      std::string s2=reverse(s1);
      int h=bintodec(s2);
      if (j<h || j==h){
	mapa[j].push_back(j);
	mapa[j].push_back(h);
      }
    }
    std::stringstream ss;
    ss << i;
    for (std::map<int, std::vector<int> >::const_iterator it=mapa.begin(); it!=mapa.end(); ++it){//Para todas las secuencias
      std::stringstream ss_prot;
      ss_prot << it->second[0];
      std::string name="prot_metabolism/"+ss.str()+"_"+ss_prot.str()+".txt";
      std::ifstream file_dimeros_prots(name.c_str());
      check_file(file_dimeros_prots, name);
      int index;
      std::string met_id=dectobin(it->second[0],i);
      while (file_dimeros_prots >> index){
	double energy;
	double basura;
	file_dimeros_prots >> energy >> basura >> basura >> basura;
	prot_met[std::make_pair(met_id, index)]=energy;
      }
      file_dimeros_prots.close();
    }//for mets
  }//for sizes
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//PROT_BREAKING
void ToyPlugin::init_prot_breaking(mapa_pbreak& prot_breaking){
  std::ifstream in("prot_metabolism/prot_breaking.txt");
  check_file(in, "prot_breaking.txt");
  int id1;
  std::string s1;
  while (in >> id1 >> s1){
    std::pair<int,std::string> p1(id1,s1);
    double d1;
    int id2;
    in >> d1 >> id2;
    prot_breaking[p1]=std::make_pair(d1,id2);
  }
  in.close();
  return;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//PROMOTER_EXPRESSION
int ToyPlugin::promoter_expression(const std::vector<std::pair<int,int> >& genotype, std::pair<mapa_prot, mapa_dim>& objects){
  //LOOK UP THE BINDING ENERGY OF EVERY PROTEIN TO THE PROM
  std::vector<bool> promoter_active(genotype.size(),0);

  //LOOP OVER ALL PROMOTERS
  for (int p=0; p<genotype.size(); ++p){
    int prom=genotype[p].first;
    std::pair<std::map<Prot,double>, std::map<Dim,double> > binding_energies;
    for (mapa_prot::const_iterator it_prot=objects.first.begin(); it_prot != objects.first.end(); ++it_prot)//for every protein in OBJECTS.FIRST
      if (d_less(prot_prom_energies[prom][it_prot->first],0.0) &&//if the protein can bind to PROM AND
	  d_less(prot_prom_energies[prom][it_prot->first]+prot_energies[it_prot->first]-polymerase[prom],0.0))//it can bind it more strongly than POLY
	binding_energies.first[it_prot->first]=prot_prom_energies[prom][it_prot->first]+prot_energies[it_prot->first]-polymerase[prom];//we substract the binding energy of the polymerase; if it is greater than zero, binding energy is zero.
    for (mapa_dim::const_iterator it_dim=objects.second.begin(); it_dim != objects.second.end(); ++it_dim)//for every dimer in OBJECTS.SECOND
      if (d_less(dim_prom_energies[prom][it_dim->first.id],0.0) &&//if the dimer can bind to PROM AND
	  d_less(dim_prom_energies[prom][it_dim->first.id]+dim_bond_energy[it_dim->first.id]+prot_energies[it_dim->first.p1]+prot_energies[it_dim->first.p2]-polymerase[prom], 0.0))//it can bind it more strongly than POLY
	binding_energies.second[it_dim->first]=dim_prom_energies[prom][it_dim->first.id]+dim_bond_energy[it_dim->first.id]+prot_energies[it_dim->first.p1]+prot_energies[it_dim->first.p2]-polymerase[prom];//we substract the binding energy of the polymerase; if it is greater than zero, binding energy is zero 
 
    //WE LOOK FOR THE MINIMUM
    while (1){
      double min=0.0;
      bool i_min=0;//"0" means the minimum is a protein; "1" means it is a dimer
      std::map<Prot,double>::const_iterator prot_min;
      std::map<Dim, double>::const_iterator dim_min;
      for (std::map<Prot,double>::const_iterator it1=binding_energies.first.begin(); it1 != binding_energies.first.end(); ++it1){//loop over all proteins, we save its id (i_min) and its energy (min)
	if (d_less(it1->second,min)){
	  min=it1->second;
	  prot_min=it1;
	}//if binding<min
      }//for i 
      for (std::map<Dim, double>::const_iterator it2=binding_energies.second.begin(); it2 != binding_energies.second.end(); ++it2){//loop over all dimers, we save its id (i_min) and its energy (min)
	if (d_less(it2->second,min)){
	  min=it2->second;
	  i_min=1;
	  dim_min=it2;
	}//if binding<min
      }//for i
      
      if (d_equal(min,0.0)){//If the minimal is 0, no protein or dimer binds, and it all depends on the polymerase
	promoter_active[p]=d_less(polymerase[prom],0.0);
	break;//we exit the loop
      }
      
      //We check if there is more than one minimum, and eliminate the repeated ones
      int count_mins=0;
      if (!i_min){//if the minimum is a protein
	for (std::map<Prot,double>::const_iterator it1=std::next(prot_min); it1 != binding_energies.first.end(); ++it1)//we start the loop on the object following the minimal (the others have greater binding energies)
	  if (d_equal(it1->second,min))//if the binding energies are the same
	    count_mins++;//we sum to the counter
	for (std::map<Dim, double>::const_iterator it2=binding_energies.second.begin(); it2 != binding_energies.second.end(); ++it2)//we compare the minimum protein with the dimers
	  if (d_equal(it2->second, min))
	    count_mins++;
      }//if minimum is a protein
      else{//minimum is a dimer (then we don't need to check repetitions in proteins, because we examined the proteins first)
	for (std::map<Dim, double>::const_iterator it2=std::next(dim_min); it2 != binding_energies.second.end(); ++it2)//we start the loop on the object following the minimal (the others have greater binding energies)
	  if (d_equal(it2->second,min))//if the binding energies are the same
	    count_mins++;//we sum to the counter
      }//minimum is a dimer
      
      if (count_mins){//if there are repeated energies
	for (std::map<Prot,double>::iterator it1=binding_energies.first.begin(); it1 != binding_energies.first.end(); ++it1)
	  if (d_equal(it1->second,min))//if the binding energies are the same
	    binding_energies.first.erase(it1);//we sum to the counter
	for (std::map<Dim, double>::iterator it2=binding_energies.second.begin(); it2 != binding_energies.second.end(); ++it2)
	  if (d_equal(it2->second, min))
	    binding_energies.second.erase(it2);
      }
      else{//if there is only one minimum, we check if PROM is activated
	if (!i_min)
	  promoter_active[p]=prot_pol[prom][prot_min->first];//return the minimum
	else
	  promoter_active[p]=dim_pol[prom][dim_min->first.id];//return the minimum
	break;//we exit the loop
      }
    }//while TRUE
  }//for promoters
  
  //EXPRESSION
  objects.first.clear();//now we know which promoters are expressed, all objects disappear
  objects.second.clear();
  for (int i=0; i<promoter_active.size(); ++i)
    if (promoter_active[i])
      objects.first[prot_gen[genotype[i].second]]++;//we add the proteins that are expressed

  //OUTPUT
  int translated=0;//as every gene has a place in the genotype, we can express every state as an integer: (1,1,0) = 2^2+2^1 = 6
  for (int i=0; i<promoter_active.size(); ++i)
    translated += promoter_active[i]*int_pow(2,promoter_active.size()-i-1);
  return translated;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//FINDING_CYCLES
bool ToyPlugin::finding_cycles(std::vector<std::pair<int, std::pair<mapa_prot, mapa_dim> > >& iterations, int expression){
  //OUTPUT
  //Given a vector of all visited states and the last state visited, TERNA_FINDING_CYCLES checks if this state has been visited
  //before. If it has, then it transforms the ITERATIONS vector and leaves only the states belonging to the cycle; it returns 
  //TRUE. If it hasn't, it returns FALSE.

  //ARGUMENTS
  //(1) (std::vector<int>& ITERATIONS: the vector of visited states
  //(2) (int) EXPRESSION: the last state visited (it is not part of the ITERATIONS vector)

  //LOOP OVER ITERATIONS VECTOR
  for (int i=0; i!=iterations.size(); ++i)
    if (iterations[i].first==expression){//if, at any point, this state has been visited before, we have a cycle
      iterations.erase(iterations.begin(),iterations.begin()+i);//we erase all the previous states, as they are transient
      return 1;
    }
  
  return 0;//no cycle is found
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//DIMERIZATION
void ToyPlugin::dimerization(std::pair<mapa_prot, mapa_dim>& objects){
  //WHICH DIMERS FORM
  std::map<double, mapa_dim> binding_energies;//creates a map of energies, so that dimers with the same final energy are classified in the same group. Each energy has a map of dimers, consisting of the dimer and the number of times it can be formed  
  for (mapa_prot::const_iterator it1=objects.first.begin(); it1!=objects.first.end(); ++it1){//we only want to explore the combinations of two proteins
    for (mapa_prot::const_iterator it2=it1; it2!=objects.first.end(); ++it2){//remember that proteins are sorted (and that they can form dimers with themselves)
      if (it1!=it2 ||
	  (it1==it2 && it1->second>1)){//Prots will only dimerize with themselves if there is more than one copy
	std::pair<int,int> perim=std::make_pair(prot_perimeters[it1->first], prot_perimeters[it2->first]);
	std::map<std::pair<int,int>,int>::const_iterator it_dim=dim_perim.find(perim);
	if (it_dim!=dim_perim.end()){//if this pair forms a dimer, we form the tuple
	  double b_energy=dim_bond_energy[it_dim->second]+prot_energies[it1->first]+prot_energies[it2->first];
	  binding_energies[b_energy][Dim(it_dim->second,it1->first,it2->first)]+=int_min(it1->second,it2->second);
	}//if a dimer can be formed
      }//if it1!=it2
    }//for it2
  }//for it1
  
  //NOW WE HAVE A MAP OF ENERGIES
  //LOOP OVER EACH ENERGY
  for (std::map<double, mapa_dim>::iterator it=binding_energies.begin(); it!= binding_energies.end(); ++it){
    std::map<Prot,int> number;//will hold how many complexes each protein is involved in (if it's greater than objects.first[i], then no dimers is formed
    for (mapa_dim::iterator it_dim=it->second.begin(); it_dim!=it->second.end(); ++it_dim){//for each possible dimer
      number[it_dim->first.p1]+=it_dim->second;
      number[it_dim->first.p2]+=it_dim->second;
    }
    //Now, for every possible complex, we check if the proteins involved are involved in too many complexes
    for (mapa_dim::iterator it_dim=it->second.begin(); it_dim!=it->second.end(); ++it_dim)//for each dimer
      if (objects.first[it_dim->first.p1]<number[it_dim->first.p1] ||//if either of the two prots take part
	  objects.first[it_dim->first.p2]<number[it_dim->first.p2])//in more dimers that they can
	it->second.erase(it_dim);//this dimer will not form
    //Now we form the remaining dimers
    for (mapa_dim::iterator it_dim=it->second.begin(); it_dim!=it->second.end(); ++it_dim){//for each dimer
      objects.second[it_dim->first]+=it_dim->second;
      objects.first[it_dim->first.p1]-=it_dim->second;
      objects.first[it_dim->first.p2]-=it_dim->second;
    }
  }//for each energy

  for (mapa_prot::iterator it1=objects.first.begin(); it1!=objects.first.end(); ++it1)
    if (it1->second==0)
      objects.first.erase(it1);//eliminate the proteins that have 0 copies
}  
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//BOOLEAN NETWORK
std::vector<std::pair<int, std::pair<mapa_prot, mapa_dim> > > ToyPlugin::boolean_network(const std::vector<std::pair<int,int> >& genotype){

  //COMPUTES THE ATTRACTOR OF THE DYNAMICS OF A toyLIFE ORGANISM
  std::pair<mapa_prot, mapa_dim> objects;
  std::vector<std::pair<int, std::pair<mapa_prot, mapa_dim> > > cycle(1,std::make_pair(0, objects));//first state is 0
  int nu_state=promoter_expression(genotype, objects);//compute next state
  if (nu_state==0)//the new state is also 0
    return std::vector<std::pair<int, std::pair<mapa_prot, mapa_dim> > >(); 
  dimerization(objects);//dimerize
  cycle.push_back(std::make_pair(nu_state, objects));
    
  while (1){
    int new_state=promoter_expression(genotype, objects);
    if (finding_cycles(cycle, new_state))//if the state has already been visited
      break;//this is the only way out of the while loop
    dimerization(objects);//dimerize
    cycle.push_back(std::make_pair(new_state, objects));//if the state hasn't been visited, we keep looking
  }
  return cycle;//this is the vector containing the attractor states
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//COMPLETE BOOLEAN FUNCTION
std::vector<int> ToyPlugin::complete_boolean_function(const std::vector<std::pair<int,int> >& genotype){

  //COMPUTES THE ATTRACTOR OF THE DYNAMICS OF A toyLIFE ORGANISM
  int g=genotype.size();
  int state_space=int_pow(2, g);
  std::vector<int> logic_function(state_space, 0);
  std::vector<bool> visited(state_space, 0);
  for (int i=0; i<state_space; ++i){
    if (visited[i]) continue;
    std::vector<int> who=str_to_vec(dectobin(i, g));
    std::pair<mapa_prot, mapa_dim> objects;
    for (int p=0; p<g; ++p)
      if (who[p])
	objects.first[prot_gen[genotype[p].second]]++;
    dimerization(objects);//dimerize
    int nu_state=promoter_expression(genotype, objects);//compute next state   
    logic_function[i]=nu_state;
  }
  return logic_function;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//REACTING
void ToyPlugin::reacting(std::pair<mapa_prot, mapa_dim>& objects, mapa_met& met, mapa_owm& objects_with_mets, double& fit){
  //GIVEN A SET OF OBJECTS, IT TELLS US WHAT BINDS WHAT AND WHAT BREAKS WHAT
  std::map<double, mapa_owm> binding_energies;//creates a map of energies, so that complexes with the same final energy are classified in the same group. Each energy has a map of OWM, which has three members: prot, dim and met (if there is no Prot in the complex, it will be -1; if there is no Dim, it will be (-1,0,0), and if there is no Met it will be the empty string, "") and [VALUE] an int counting how many possible molecules can be formed of this type
  //First we check if any Dimer forms
  for (mapa_prot::const_iterator it1=objects.first.begin(); it1!=objects.first.end(); ++it1)//we only want to explore the combinations of two proteins
    for (mapa_prot::const_iterator it2=it1; it2!=objects.first.end(); ++it2)//remember that proteins are sorted (and that they can form dimers with themselves)
      if (it1!=it2 ||
	  (it1==it2 && it1->second>1)){//Prots will only dimerize with themselves if there is more than one copy
	std::pair<int,int> perim=std::make_pair(prot_perimeters[it1->first], prot_perimeters[it2->first]);
	std::map<std::pair<int,int>,int>::const_iterator it_dim=dim_perim.find(perim);
	if (it_dim!=dim_perim.end()){//if this pair forms a dimer, we form the tuple
	  double b_energy=dim_bond_energy[it_dim->second]+prot_energies[it1->first]+prot_energies[it2->first];//energy of the complex
	  Dim new_dim(it_dim->second, it1->first, it2->first);
	  OWM new_complex(-1, new_dim, Met());
	  binding_energies[b_energy][new_complex]=int_min(it1->second, it2->second);
	}//if a dimer can be formed
      }//if it1!=it2
  
  //Then we check if any Protein binds any Met
  for (mapa_prot::const_iterator it_prot=objects.first.begin(); it_prot!=objects.first.end(); ++it_prot)
    for (mapa_met::const_iterator it_met=met.begin(); it_met!=met.end(); ++it_met){
      pairPM par1=std::make_pair(it_met->first, it_prot->first);
      try{//see if the prot binds the met
	// std::map<Met, std::map<Prot,double> >::const_iterator iter_tuplita=prot_met.find(it_met->first);
	// std::map<Prot,double> tuplita=iter_tuplita->second;
	double e1=prot_met.at(par1);
	if (d_less(e1,0.0)){//if the protein can bind the metabolite
	  double b_energy=e1+prot_energies[it_prot->first];//energy of the complex
	  OWM new_complex(it_prot->first, Dim(), it_met->first);
	  binding_energies[b_energy][new_complex]=int_min(it_prot->second, it_met->second);
	}
      }
      catch (std::out_of_range){}//if it doesn't, don't do anything
    }
							
  //Then we check if any Dimer binds any Met
  for (mapa_dim::const_iterator it_dim=objects.second.begin(); it_dim!=objects.second.end(); ++it_dim)
    for (mapa_met::const_iterator it_met=met.begin(); it_met!=met.end(); ++it_met){
      pairDM par1=std::make_pair(it_met->first, it_dim->first.id);
      try{//look and see if the dimer binds to the met
	Dmet tupla=dim_met.at(par1);
	if (d_less(tupla.eg,0.0)){//if the dimer can bind the metabolite
	  double b_energy=tupla.eg+dim_bond_energy[it_dim->first.id]+prot_energies[it_dim->first.p1]+prot_energies[it_dim->first.p2];//energy of the complex
	  OWM new_complex(-1,it_dim->first, it_met->first);
	  binding_energies[b_energy][new_complex]=int_min(it_dim->second, it_met->second);
	}
      }
      catch (std::out_of_range){}//if it doesn't don't do anything
    }
  
  //Finally, we check if any D+M breaks
  for (mapa_owm::const_iterator it_dm=objects_with_mets.begin(); it_dm!=objects_with_mets.end(); ++it_dm){
    if (it_dm->first.dim.id==-1) continue;//this means it's a P+M complex
    for (mapa_prot::const_iterator it_prot=objects.first.begin(); it_prot!=objects.first.end(); ++it_prot){
      pairDM par1=std::make_pair(it_dm->first.met, it_dm->first.dim.id);
      Dmet tupla=dim_met.at(par1);
      if (tupla.seq.size()!=16) continue;//only check this if the D+M can actually be broken!
      std::pair<int, std::string> par2(it_prot->first, tupla.seq);
      try{//look and see if the protein can break this dimer
	double b_energy=prot_breaking.at(par2).first+prot_energies[it_prot->first];
	int c_min=prot_breaking.at(par2).second;
	if (c_min<2)//the sequences should be ordered in such a way that the first half is always the first protein
	  b_energy += prot_energies[it_dm->first.dim.p1];//sum energy of first subunit
	// if ((c_min<2 && tupla.bind<3) ||//if the prot binds to the first subunit when the dimer is not reversed 
	//     (c_min>1 && tupla.bind>2))//or if it binds to the second subunit when the dimer is reversed
	//   b_energy += prot_energies[it_dm->first.dim.p1];//sum energy of first subunit
	else
	  b_energy += prot_energies[it_dm->first.dim.p2];//sum energy of second subunit
	double met_energy=dim_bond_energy[it_dm->first.dim.id]+prot_energies[it_dm->first.dim.p1]+prot_energies[it_dm->first.dim.p2];//E of the existing complex
	if (d_less(b_energy, met_energy)){//then it could break
	  OWM new_complex(it_prot->first,it_dm->first.dim, it_dm->first.met);
	  binding_energies[b_energy][new_complex]=int_min(it_prot->second, it_dm->second);
	}
      }//try
      catch (std::out_of_range){}//it doesn't bind to this sequence
    }//for prots
  }//for d+m complexes
  
  //NOW WE HAVE A MAP OF ENERGIES
  //LOOP OVER EACH ENERGY AND SEE WHAT HAPPENS
  for (std::map<double, mapa_owm>::iterator it=binding_energies.begin(); it!= binding_energies.end(); ++it){
    //std::cout << "Energy! " << it->first << std::endl;
    mapa_prot number_prots;//will hold how many complexes each protein is involved in (if it's greater than objects.first[i], then no dimers is formed
    mapa_dim number_dims;//will hold how many complexes each protein is involved in (if it's greater than objects.first[i], then no dimers is formed
    mapa_met number_mets;//will hold how many complexes each protein is involved in (if it's greater than objects.first[i], then no dimers is formed
    mapa_owm number_dm;//because std::map doesn't work well with boost::tuple, this vector will serve as an index for the elements in objects_with_mets;
    for (mapa_owm::iterator it_tuple=it->second.begin(); it_tuple != it->second.end(); ++it_tuple)//for each complex
      if (it_tuple->first.met.empty()){//protein-protein complex
	number_prots[it_tuple->first.dim.p1]+=it_tuple->second;
	number_prots[it_tuple->first.dim.p2]+=it_tuple->second;
      }
      else if (it_tuple->first.prot!=-1){//there are proteins and mets
	if (it_tuple->first.dim.empty()){//protein-met complex (no dimers)
	  number_prots[it_tuple->first.prot]+=it_tuple->second;
	  number_mets[it_tuple->first.met]+=it_tuple->second;
	}
	else{//protein breaking a d+m complex
	  number_prots[it_tuple->first.prot]+=it_tuple->second;
	  number_dm[OWM (-1,it_tuple->first.dim,it_tuple->first.met)]+=it_tuple->second;
	}
      }//proteins and mets are involved
      else{//dimer-met complex
	number_dims[it_tuple->first.dim]+=it_tuple->second;
	number_mets[it_tuple->first.met]+=it_tuple->second;
      }
    //NOW, FOR EVERY POSSIBLE COMPLEX, WE CHECK IF THE UNITS ARE INVOLVED IN TOO MANY COMPLEXES AND DISCARD THOSE
    for (mapa_owm::iterator it_tuple=it->second.begin(); it_tuple != it->second.end(); ++it_tuple)//for each complex
      if (it_tuple->first.met.empty()){//protein-protein complex
	if (objects.first[it_tuple->first.dim.p1]<number_prots[it_tuple->first.dim.p1] ||//prot1 more involved
	    objects.first[it_tuple->first.dim.p2]<number_prots[it_tuple->first.dim.p2])//prot2 more involved
	  it->second.erase(it_tuple);//this dimer will not form
      }
      else if (it_tuple->first.prot!=-1){//proteins and mets are involved
	if (it_tuple->first.dim.empty()){//protein-met complex
	  if (objects.first[it_tuple->first.prot]<number_prots[it_tuple->first.prot] ||//proteins more involved
	      met[it_tuple->first.met]<number_mets[it_tuple->first.met])//met more involved
	    it->second.erase(it_tuple);//this P+M will not form
	}
	else{//protein breaking a d+m complex
	  OWM new_tuple(-1,it_tuple->first.dim,it_tuple->first.met);
	  if (objects.first[it_tuple->first.prot]<number_prots[it_tuple->first.prot] ||//protein more involved
	      objects_with_mets[new_tuple]<number_dm[new_tuple])//d+m more involved
	    it->second.erase(it_tuple);//this D+M will not be broken
	}
      }//proteins and mets are involved
      else{//dimer-met complex 
	if (objects.second[it_tuple->first.dim]<number_dims[it_tuple->first.dim] ||//dimer more involved
	    met[it_tuple->first.met]<number_mets[it_tuple->first.met])//met more involved
	  it->second.erase(it_tuple);//this D+M will not form
      }
    //NOW WE FORM THE COMPLEXES THAT REMAIN
    for (mapa_owm::const_iterator it_tuple=it->second.begin(); it_tuple != it->second.end(); ++it_tuple)//for each complex
      if (it_tuple->first.met.empty()){//protein-protein complex
	objects.second[it_tuple->first.dim]+=it_tuple->second;//we add as many copies as possible
	objects.first[it_tuple->first.dim.p1]-=it_tuple->second;//we eliminate the Prot copies
	objects.first[it_tuple->first.dim.p2]-=it_tuple->second;
      }
      else if (it_tuple->first.prot!=-1){//there are proteins and mets
	if (it_tuple->first.dim.empty()){//protein-met complex
	  objects_with_mets[it_tuple->first]+=it_tuple->second;
	  objects.first[it_tuple->first.prot]-=it_tuple->second;
	  met[it_tuple->first.met]-=it_tuple->second;
	}
	else{//protein breaking a d+m complex
	  objects_with_mets[OWM (-1,it_tuple->first.dim,it_tuple->first.met)]-=it_tuple->second;
	  objects.first[it_tuple->first.prot]-=it_tuple->second;
	  //Now we add the rests of the mets to Met
	  pairDM par1=std::make_pair(it_tuple->first.met, it_tuple->first.dim.id);
	  Dmet tupla=dim_met.at(par1);
	  std::string pmet1=it_tuple->first.met.substr(0,tupla.pos);//we split the met in the two parts defined by the Dimer when it binds
	  std::string pmet2=it_tuple->first.met.substr(pmet1.size());
	  met[pmet1]+=it_tuple->second;
	  met[pmet2]+=it_tuple->second;
	  fit += it_tuple->second;//this means it's broken!
	}
      }//proteins and mets
      else{//dimer-met complex (it should check this: if (boost::get<0>(boost::get<1>(*it_tuple))!=-1)
	objects_with_mets[it_tuple->first]+=it_tuple->second;
	objects.second[it_tuple->first.dim]-=it_tuple->second;
	met[it_tuple->first.met]-=it_tuple->second;
      }
  }//for all energies
  
  //NOW WE ELIMINATE THE ELEMENTS THAT HAVE BEEN USED UP
  for (mapa_prot::iterator it_prot=objects.first.begin(); it_prot!=objects.first.end(); ++it_prot)
    if (it_prot->second==0)
      objects.first.erase(it_prot);//eliminate the proteins that have 0 copies
  
  for (mapa_dim::iterator it_dim=objects.second.begin(); it_dim!=objects.second.end(); ++it_dim)
    if (it_dim->second==0)
      objects.second.erase(it_dim);//eliminate dimers that have 0 copies
  
  for (mapa_met::iterator it_met=met.begin(); it_met!=met.end(); ++it_met)
    if (it_met->second==0)
      met.erase(it_met);//eliminate mets that have 0 copies

  for (mapa_owm::iterator it_owm=objects_with_mets.begin(); it_owm!=objects_with_mets.end(); ++it_owm)
    if (it_owm->second==0)
      objects_with_mets.erase(it_owm);
  
  return;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//METABOLISM
int ToyPlugin::metabolism(const std::vector<std::pair<int,int> >& genotype, std::pair<mapa_prot, mapa_dim> objects, mapa_met& met){
  int gene_number=genotype.size();
  mapa_owm objects_with_mets;
  int cuenta=0;//if it hasn't broken the Met in 2^n steps, then it won't be able to do it
  double bas=0.0;
  while (1){//we loop over the cycle until the metabolite is broken or it's let go
    //START THE CYCLE: BINDING
    reacting(objects, met, objects_with_mets, bas);
    if (objects_with_mets.empty())//that means no one binds
      return 0;
    //REGULATION
    int new_state=promoter_expression(genotype, objects);
    //NOW THE NEW PROTEINS WILL DIMERIZE IN THE PRESENCE OF THE METABOLITE-COMPOUND
    mapa_owm::iterator it_owm1=objects_with_mets.begin();
    reacting(objects, met, objects_with_mets, bas);
    mapa_met::iterator it_met1=met.begin();
    if (it_met1!=met.end() && it_met1->first.size()<8)//that means something is broken
      return 1;
    met.clear();//at the end of the cycle, all metabolites that haven't been bound by anything disappear
    //Now the P+M and D+M will become part of Met too (the Ps and Ds will disappear)
    for (mapa_owm::iterator it=objects_with_mets.begin(); it!=objects_with_mets.end(); ++it)
      met[it->first.met]+=it->second;
    objects_with_mets.clear();//and all objects with mets will disappear of course
    ++cuenta;
    if (cuenta>int_pow(2,gene_number))//if it hasn't broken the Met in 2^n steps, then it won't be able to do it
      return 0;
  }//while
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//PHENOTYPE
std::vector<int> ToyPlugin::vec_phenotype(const std::vector<std::pair<int,int> >& genotype, const std::vector<std::string>& env){
  //ALGORITHM
  //COMPUTE THE CYCLE FIRST (WITHOUT METABOLITES)
  std::vector<std::pair<int,std::pair<mapa_prot, mapa_dim> > > cycle=boolean_network(genotype);
 
  //NOW WE CHECK IF IN ANY MOMENT OF THE CYCLE, IT CAN BREAK ANYTHING
  std::vector<int> phenotype(env.size());
  for (int m=0; m<env.size(); ++m){//loop over all metabolites
    std::map<Met,int> met;
    for (int i=0; i<cycle.size(); ++i){
      met[env[m]]++;//this will allow us to generalize for more than one metabolite
      phenotype[m]=metabolism(genotype, cycle[i].second, met);
      met.clear();
      if (phenotype[m])
	break;
    }
  }
  
  return phenotype;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//GENOTYPE_STR_TO_VEC
std::vector<std::pair<int,int> > ToyPlugin::genotype_str_to_vec(const std::string& genotype){
  int gene_number=genotype.size()/20;
  std::vector<std::pair<int,int> > vec_genotype(gene_number);
  for (int g=0; g<gene_number; ++g){
    vec_genotype[g].first=bintodec(genotype.substr(20*g, 4));
    vec_genotype[g].second=bintodec(genotype.substr(20*g+4, 16));
  }
  return vec_genotype;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//STR LOGIC FUNCTION
std::vector<int> ToyPlugin::regulatory_phenotype(const std::string& genotype){
  std::vector<std::pair<int,int> > vec_genotype=genotype_str_to_vec(genotype);
  return complete_boolean_function(vec_genotype);
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//STR PHENOTYPE
std::vector<int> ToyPlugin::metabolic_phenotype(const std::string& genotype, const std::vector<std::string>& env){
  std::vector<std::pair<int,int> > vec_genotype=genotype_str_to_vec(genotype);
  return vec_phenotype(vec_genotype, env);
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//METABOLIC ROBUSTNESS
double ToyPlugin::metabolic_robustness(const std::string& genotype, const std::vector<int>& phen, const std::vector<std::string>& env){
  double rob=0.0;
  int gensize=genotype.size(); 
  for (int g=0; g<gensize; ++g){
    std::string mutgen=mutation(genotype, g);
    std::vector<int> mutphen=metabolic_phenotype(mutgen, env);
    if (mutphen==phen)
      ++rob;
  }
  return rob/gensize;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//REGULATORY ROBUSTNESS
double ToyPlugin::regulatory_robustness(const std::string& genotype, const std::vector<int>& logic_vector){
  double rob=0.0;
  int gensize=genotype.size(); 
  for (int g=0; g<gensize; ++g){
    std::string mutgen=mutation(genotype, g);
    std::vector<int> mutphen=regulatory_phenotype(mutgen);
    if (mutphen==logic_vector)
      ++rob;
  }
  return rob/gensize;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//METABOLIC EVOLVABILITY
std::map<std::vector<int>, int> ToyPlugin::metabolic_evolvability(const std::string& genotype, const std::vector<int>& phen, const std::vector<std::string>& env){
  std::map<std::vector<int>,int> evolvability;
  int gensize=genotype.size(); 
  for (int g=0; g<gensize; ++g){
    std::string mutgen=mutation(genotype, g);
    std::vector<int> mutphen=metabolic_phenotype(mutgen, env);
    if (mutphen!=phen)
      evolvability[mutphen]++;
  }
  return evolvability;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//REGULATORY EVOLVABILITY
std::map<std::vector<int>, int> ToyPlugin::regulatory_evolvability(const std::string& genotype, const std::vector<int>& logic_vector){
  std::map<std::vector<int>,int> evolvability;
  int gensize=genotype.size(); 
  for (int g=0; g<gensize; ++g){
    std::string mutgen=mutation(genotype, g);
    std::vector<int> mutphen=regulatory_phenotype(mutgen);
    if (mutphen!=logic_vector)
      evolvability[mutphen]++;
  }
  return evolvability;
}
