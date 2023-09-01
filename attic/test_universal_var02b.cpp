// August 28, 2023.
// Test routine for converting Keplerian orbital elements
// to Cartesian state vectors.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
  
int main()
{
  point3d startpos = point3d(0.0L,0.0L,0.0L);
  point3d startvel = point3d(0.0L,0.0L,0.0L);
  double vesc=0.0l;
  double heliodist=AU_KM;
  double v=0.0l;
  double MGsun=GMSUN_KM3_SEC2;
  point3d outpos = point3d(0.0L,0.0L,0.0L);
  point3d outvel = point3d(0.0L,0.0L,0.0L);
  double mjdstart = 60000.0l;
  double mjdend = 60020.0l;
  int testnum=1;
  int testct=0;
  
  cout << "Enter the number of tests to run\n";
  cin >> testnum;
  cout << testnum << " tests will be run\n";
  for(testct=0; testct<testnum; testct++) {
    startpos.x = (2.0l*rand()/RAND_MAX - 1.0l)*3.0l*AU_KM;
    startpos.y = (2.0l*rand()/RAND_MAX - 1.0l)*3.0l*AU_KM;
    startpos.z = (2.0l*rand()/RAND_MAX - 1.0l)*3.0l*AU_KM;

    heliodist = vecabs3d(startpos);
    vesc = sqrt(2.0*MGsun/heliodist);
    
    startvel.x = (4.0l*rand()/RAND_MAX - 1.0l)*vesc;
    startvel.y = (4.0l*rand()/RAND_MAX - 1.0l)*vesc;
    startvel.z = (4.0l*rand()/RAND_MAX - 1.0l)*vesc;

    v = vecabs3d(startvel);

    while(v<=vesc) {
      startvel.x = (4.0l*rand()/RAND_MAX - 1.0l)*vesc;
      startvel.y = (4.0l*rand()/RAND_MAX - 1.0l)*vesc;
      startvel.z = (4.0l*rand()/RAND_MAX - 1.0l)*vesc;

      v = vecabs3d(startvel);
    }

    Kepler_univ_int(MGsun, mjdstart, startpos, startvel, mjdend, outpos, outvel);
  }

  return(0);
}
