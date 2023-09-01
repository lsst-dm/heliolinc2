// August 28, 2023.
// Test routine for converting Keplerian orbital elements
// to Cartesian state vectors.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
  
int main()
{
  point3d startpos = point3d(0.0L,0.0L,0.0L);
  point3d endpos = point3d(0.0L,0.0L,0.0L);
  point3d startvel = point3d(0.0L,0.0L,0.0L);
  double MGsun=GMSUN_KM3_SEC2;
  double a = 0.0l;
  int testnum=1;
  int testct=0;
  string seedstring;
  //double a;

  cout << "Enter a random number seed\n";
  cin >> seedstring;
  cout << seedstring << "\n";
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  cout << "Enter the number of tests to run\n";
  cin >> testnum;
  cout << testnum << " tests will be run\n";
  for(testct=0; testct<testnum; testct++) {
    startpos.x = (unitvar(generator)-0.5l)*AU_KM*4.0l;
    startpos.y = (unitvar(generator)-0.5l)*AU_KM*4.0l;
    startpos.z = (unitvar(generator)-0.5l)*AU_KM*4.0l;
    endpos.x = startpos.x + (unitvar(generator)-0.5l)*AU_KM*0.4l;
    endpos.y = startpos.y + (unitvar(generator)-0.5l)*AU_KM*0.4l;
    endpos.z = startpos.z + (unitvar(generator)-0.5l)*AU_KM*0.4l;
    double r1 = vecabs3d(startpos);
    double r2 = vecabs3d(endpos);
    while(r1<0.5l*AU_KM || r2<0.5l*AU_KM) {
      // Too close to the sun, select new random positions.
      startpos.x = (unitvar(generator)-0.5l)*AU_KM*4.0l;
      startpos.y = (unitvar(generator)-0.5l)*AU_KM*4.0l;
      startpos.z = (unitvar(generator)-0.5l)*AU_KM*4.0l;
      endpos.x = startpos.x + (unitvar(generator)-0.5l)*AU_KM*0.4l;
      endpos.y = startpos.y + (unitvar(generator)-0.5l)*AU_KM*0.4l;
      endpos.z = startpos.z + (unitvar(generator)-0.5l)*AU_KM*0.4l;
      r1 = vecabs3d(startpos);
      r2 = vecabs3d(endpos);
    }
    point3d pdiff = point3d(endpos.x - startpos.x, endpos.y - startpos.y, endpos.z - startpos.z);
    double c = vecabs3d(pdiff);
    double dtp = ((r1+r2+c)*sqrt(r1+r2+c) - (r1+r2-c)*sqrt(r1+r2-c))/6.0l/sqrt(MGsun);
    dtp/=SOLARDAY;
    double timediff = unitvar(generator)*30.0 + dtp;
    startvel = Twopoint_Kepler_v1(MGsun, startpos, endpos, timediff, 1.0, &a, KEPTRANSITMAX, 0);
  }
  return(0);
}
