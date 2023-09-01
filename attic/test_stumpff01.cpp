// August 28, 2023.
// Test routine for converting Keplerian orbital elements
// to Cartesian state vectors.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
  
int main()
{
  double x=0.0l;
  int testnum=1;
  int testct=0;
  int status=0;
  double c0,c1,c2,c3;
  c0=c1=c2=c3=0.0l;
  double d0,d1,d2,d3;
  d0=d1=d2=d3=0.0l;
  int runtype=2;
  double xrange;
  
  cout << "Enter the number of tests to run\n";
  cin >> testnum;
  cout << testnum << " tests will be run\n";
  cout << "Enter 0 to run the analytical version, 1 to run the\n";
  cout << "continuing-fraction version, or 2 to run both\n";
  cin >> runtype;
  cout << runtype << "\n";
  cout << "Enter the range for possible values of x\n";
  cin >> xrange;
  cout << xrange << "\n";
 
  for(testct=0; testct<testnum; testct++) {
    x = xrange*(2.0l*rand()/RAND_MAX - 1.0l);
    if(runtype==0 || runtype==2) {
      status = Stumpff_func(x, &c0, &c1, &c2, &c3);
      if(status!=0) {
	cerr << "ERROR: Stumpff_func returned error status " << status << "\n";
	return(status);
      }
    }
    if(runtype==2) cout << c0 << " " << c1 << " " << c2 << " " << c3 << "\n";
    if(runtype==1 || runtype==2) {
      Stumpff_func_cf(x, &d0, &d1, &d2, &d3);
      if(status!=0) {
	cerr << "ERROR: Stumpff_func_cf returned error status " << status << "\n";
	return(status);
      }
    }
    if(runtype==2) {
      cout << c0 << " " << c1 << " " << c2 << " " << c3 << "\n";
      cout << x << ": " << std::scientific << (c0-d0)/c0 << " " << (c1-d1)/c1 <<" " << (c2-d2)/c2 << " " << (c3-d3)/c3 << "\n";
    }
  }
  return(0);
}
