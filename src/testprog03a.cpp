// December 05, 2023: experiment with vector.insert functionality

#include "std_lib_facilities.h"
#include "cmath"
    
int main(int argc, char *argv[])
{
  vector <long> lvec;
  long i=0;
  
  cout << "Confirm lvec is empty: size = " << lvec.size() << "\n";
  lvec.insert(lvec.begin(),2);

  for(i=0; i<long(lvec.size()); i++) {
    cout << lvec[i] << " ";
  }
  cout << "\n";

  lvec.insert(lvec.begin(),0);
  for(i=0; i<long(lvec.size()); i++) {
    cout << lvec[i] << " ";
  }
  cout << "\n";
  lvec.insert(lvec.begin()+1,1);
  for(i=0; i<long(lvec.size()); i++) {
    cout << lvec[i] << " ";
  }
  cout << "\n";
  lvec.insert(lvec.begin()+3,3);
  for(i=0; i<long(lvec.size()); i++) {
    cout << lvec[i] << " ";
  }
  cout << "\n";
  lvec.insert(lvec.end(),5);
  for(i=0; i<long(lvec.size()); i++) {
    cout << lvec[i] << " ";
  }
  cout << "\n";
  lvec.insert(lvec.end()-1,4);
  for(i=0; i<long(lvec.size()); i++) {
    cout << lvec[i] << " ";
  }
  cout << "\n";
  
  return(0);
}
