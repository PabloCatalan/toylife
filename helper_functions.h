#ifndef _HELPER_FUNCTIONS
#define _HELPER_FUNCTIONS

//Header for toyLIFE's general functions
#include <vector> //std::vector
#include <string> //std::string
#include <fstream> //std::ofstream
#include <map> //std::map
#include <cmath>//std::pow
#include <sstream>//std::stringstream
#include <chrono>
#include <random>

//BINOMIAL COEFFICIENT
template <class T = unsigned long>
T binomial_coefficient(unsigned long n, unsigned long k) {
    unsigned long i;
    T b;
    if (0 == k || n == k) {
        return 1;
    }
    if (k > n) {
        return 0;
    }
    if (k > (n - k)) {
        k = n - k;
    }
    if (1 == k) {
        return n;
    }
    b = 1;
    for (i = 1; i <= k; ++i) {
        b *= (n - (k - i));
        if (b < 0) return -1; /* Overflow */
        b /= i;
    }
    return b;
}
//CHECK_FILE
void check_file(std::ifstream& file, const std::string& name);
//INT_POW
int int_pow(int x, int p);
//INT_MIN
int int_min(int a, int b);
//EQUAL_VECTORS
bool equal_vectors(std::vector<int> vec1, std::vector<int> vec2);
//D_EQUAL
bool d_equal(double e1, double e2);
//D_LESS
bool d_less(double e1, double e2);
//D_MIN
double d_min(double e1, double e2);
//UNION_SEQ
double union_seq(int side1, int side2);
double union_seq3(int side1, const std::string& met);
//DECIMAL TO BINARY
inline std::string dectobin(int decimal,const int size){//converts decimal integers into binary std::strings
  std::string binary(size,'0');
  for (int i=size-1; i>-1; --i)
    if (decimal/std::pow(2,i)>=1){
      binary[size-i-1]='1';
      decimal=decimal-pow(2,i);
    }
  return binary;
};
//BINARY TO DECIMAL
inline int bintodec(const std::string& binary){//converts binary strings into decimal integers
  int decimal=0;
  int binary_size=binary.size();
  for(int i=0; i<binary_size; ++i)
    decimal += std::pow(2,binary_size-i-1)*(binary[i]-'0');
  return decimal;
};
//BASE 10 TO BASE n
inline std::vector<int> base_10_to_n(int decimal, int n){//converts decimal integers into (4,4,4)
  int size=3;
  std::vector<int> seq(size,0);
  for (int i=size-1; i>-1; --i){
    int r=decimal%n;
    for (int j=0; j<n; ++j)
      if (r==j)
	seq[i]=j;
    decimal/=n;
    if (decimal==0)
      break;
  }
  return seq;
};
//BASE n TO BASE 10
inline int base_n_to_10(const std::string& binary, int n){//converts binary strings into decimal integers
  int decimal=0;
  int binary_size=binary.size();
  for(int i=0; i<binary_size; ++i)
    decimal += std::pow(n,binary_size-i-1)*(binary[i]-'0');
  return decimal;
};
//VECTOR TO STRING
inline std::string vec_to_str(const std::vector<int>& vec){
  std::string str(vec.size(),'0');
  for (int i=0; i<str.size(); ++i){
    std::stringstream ss;
    ss << vec[i];
    ss >> str[i];
  }
  return str;
};
//STRING TO VECTOR
inline std::vector<int> str_to_vec(const std::string& str){
  std::vector<int> vec(str.size(),0);
  for (int i=0; i<vec.size(); ++i)
    vec[i]=str[i]-'0';
  return vec;
};
//HAMMING
int vec_hamming(const std::vector<int>& v1, const std::vector<int>& v2);
int str_hamming(const std::string& v1, const std::string& v2);
//REVERSE
std::string reverse(const std::string& s1);
//RANDOM_GENOTYPE
std::vector<std::pair<int,int> > vec_random_genotype(int gene_number, std::uniform_real_distribution<double>& RNG, std::default_random_engine& generator);
std::string random_genotype(int gene_number, std::uniform_real_distribution<double>& RNG, std::default_random_engine& generator);
//MUTATION
std::string mutation(const std::string& genotype, int pos);

#endif
