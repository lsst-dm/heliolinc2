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
  long i=0;
  long j=0;
  long k=0;
  long num_entry=stol(argv[1]);  
  long maxval=stol(argv[2]);  
  long numref=stol(argv[3]);  
  long numtest=stol(argv[4]);
  long lnum=0;
  string seedstring = argv[5];
  vector <long> refhash;
  vector <vector <long>> refset;
  
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  // Load reference set
  refhash={};
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
    lnum = blend_vector(lvec);
    refhash.push_back(lnum);
    cout << "Hash for reference vector " << i << " is " << refhash[i] << "\n";
  }
  // Create test set
  for(i=0; i<numtest; i++) {
    lvec={};
    for(j=0;j<num_entry;j++) {
      lnum = maxval*unitvar(generator);
      lvec.push_back(lnum);
    }
    lnum = blend_vector(lvec);
    // Search for matches
    for(j=0;j<numref;j++) {
      if(lnum==refhash[j]) {
	cout << "Matching hash: test " << i << " matches ref " << j << ": " << lnum << "==" << refhash[j] << "\n";
	for(k=0; k<num_entry; k++) cout << lvec[k] << " ";
	cout << "\n";	
	for(k=0; k<num_entry; k++) cout << refset[j][k] << " ";
	cout << "\n";	
	for(k=0;k<num_entry;k++) {
	  if(lvec[k]!=refset[j][k]) {
	    cout << "Mismatch at element " << k << ": " << lvec[k] << "!=" << refset[j][k] << "\n";
	    return(1);
	  }
	}
      }
    }
    if(i%1000==0) cout << "Test number " << i << ", hash = " << lnum << "\n";
  }
  return(0);
}
