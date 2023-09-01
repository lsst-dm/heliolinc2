// August 28, 2023.
// Test routine for converting Keplerian orbital elements
// to Cartesian state vectors.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
  
int main()
{
  point3d startpos = point3d(0.0L,0.0L,0.0L);
  point3d endpos = point3d(0.0L,0.0L,0.0L);
  point3d outpos = point3d(0.0L,0.0L,0.0L);
  point3d startvel = point3d(0.0L,0.0L,0.0L);
  point3d endvel = point3d(0.0L,0.0L,0.0L);
  double MGsun=GMSUN_KM3_SEC2;
  int testnum=1;
  int testct=0;
  string seedstring;
  //double a;

  cout << "Enter a random number seed\n";
  cin >> seedstring;
  cout << seedstring << "\n";
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine
  ofstream outstream1 {"diffjunk01a"};

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
    cout << "startpos: " << startpos.x << " " << startpos.y << " " << startpos.z << "\n";
    cout << "endpos: " << endpos.x << " " << endpos.y << " " << endpos.z << "\n";
    point3d pdiff = point3d(endpos.x - startpos.x, endpos.y - startpos.y, endpos.z - startpos.z);
    double c = vecabs3d(pdiff);
    double dtp = ((r1+r2+c)*sqrt(r1+r2+c) - (r1+r2-c)*sqrt(r1+r2-c))/6.0l/sqrt(MGsun);
    dtp/=SOLARDAY;
    //double timediff = unitvar(generator)*2.0l*dtp;
    double timediff = 0.5l*dtp + unitvar(generator)*dtp;
    cout << "dtp = " << dtp << ", timediff = " << timediff << "\n";
    int status = Twopoint_Kepler_vstar(MGsun, startpos, endpos, timediff, startvel, KEPTRANSITMAX);
    if(status!=0) {
      cerr << "Warning: Twopoint_Kepler_vstar() returned error status " << status << "\n";
    }
    cout << "Twopoint_Kepler_vstar() finds initial velocity " << startvel.x << " " << startvel.y << " " << startvel.z << "\n";
    double u = dotprod3d(startpos,startvel);
    double v1 = vecabs3d(startvel);
    double a = r1*MGsun/(2.0l*MGsun - v1*v1*r1);
    double n,e;
    if(a>0.0) {
      n = sqrt(MGsun/a/a/a);
      e = sqrt(DSQUARE(1.0l-r1/a) + DSQUARE(u/n/a/a));
    } else {
      e = sqrt(DSQUARE(1.0l-r1/a) - DSQUARE(u/sqrt(-MGsun*a)));
    }
    cout << "a = " << a/AU_KM << " AU, e = " << e << "\n";
    Kepler_univ_int(MGsun, 0.0l, startpos, startvel, timediff, outpos, endvel);
    cout << "Integrating this orbit produces an ending position of:\n";
    cout << outpos.x << " " << outpos.y << " " << outpos.z << ", which may be compared with:\n";
    cout << endpos.x << " " << endpos.y << " " << endpos.z << ": differences are:\n";
    cout << outpos.x-endpos.x << " " << outpos.y-endpos.y << " " << outpos.z-endpos.z << "\n";
    outstream1 << endpos.x << " " << endpos.y << " " << endpos.z << " " << outpos.x << " " << outpos.y << " " << outpos.z << " " << a/AU_KM << " " << e << " " << endpos.x-outpos.x << " " << endpos.y-outpos.y << " " << endpos.z-outpos.z << "\n";
  }
  outstream1.close();
  return(0);
}
