// December 07, 2023: experiment with vector hashing. Arguments are
// supplied without keywords: first, the number of entries per vector
// next, the maximum value of an entry, third, the number of vectors
// in a reference set, fourth, the number in a test set. The program
// will probe for and identify any cases in which the hash of a vector
// in the test set (which is discarded) matches a the hash of a vector
// in the comparison set

#include "solarsyst_dyn_geo01.h"
#include "cmath"
  
int main(int argc, char *argv[])
{
  vector <long> lvec;
  long_index lindex = long_index(0,0);
  vector <long_index> lindvec;
  long i=0;
  long j=0;
  long k=0;
  long num_entry=stol(argv[1]);  
  long maxval=stol(argv[2]);
  long numref=stol(argv[3]);
  long lnum=0;
  long lhash=0;
  string seedstring = argv[4];
  vector <vector <long>> refset;
  
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  // Load reference set
  refset={};
  for(i=0;i<numref;i++) {
    lvec={};
    for(j=0;j<num_entry;j++) {
      lnum = maxval*unitvar(generator);
      lvec.push_back(lnum);
    }
    sort(lvec.begin(), lvec.end());
    // Check for duplicate entries in the vector
    for(j=1;j<num_entry;j++) {
      while(lvec[j]==lvec[j-1]) lvec[j]+=1;
    }
    // Load into reference vector set
    refset.push_back(lvec);
    // Calculate hash
    lhash = blend_vector(lvec);
    //cout << "Hash for reference vector " << i << " is " << refhash[i] << "\n";
    lindex = long_index(lhash,i);
    lindvec.push_back(lindex);
  }
  // Sort lindvec by hash value
  sort(lindvec.begin(), lindvec.end(), lower_long_index());
  for(i=1;i<numref;i++) {
    if(lindvec[i].lelem==lindvec[i-1].lelem) {
      cout << "Hashes match for vectors " << lindvec[i].index << " and " << lindvec[i-1].index << "\n";
      for(k=0; k<num_entry; k++) cout << refset[lindvec[i].index][k] << " ";
      cout << "\n";
      for(k=0; k<num_entry; k++) cout << refset[lindvec[i-1].index][k] << " ";
      cout << "\n";
      for(k=0; k<num_entry; k++) {
	if(refset[lindvec[i].index][k] != refset[lindvec[i-1].index][k]) {
	  cout << "Mismatch at vector element " << k << ": " << refset[lindvec[i].index][k] << " != " <<  refset[lindvec[i-1].index][k] << "\n";
	  return(1);
	}
      }
    }
  }
 
  return(0);
}
