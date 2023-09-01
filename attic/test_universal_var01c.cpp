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
  point3d outpos1 = point3d(0.0L,0.0L,0.0L);
  point3d outvel1 = point3d(0.0L,0.0L,0.0L);
  point3d outpos2 = point3d(0.0L,0.0L,0.0L);
  point3d outvel2 = point3d(0.0L,0.0L,0.0L);
  double mjdstart = 60000.0l;
  double mjdend = 60020.0l;
  int testnum=1;
  int testct=0;
  point3d posdiff = point3d(0.0L,0.0L,0.0L);
  point3d veldiff = point3d(0.0L,0.0L,0.0L);
  ofstream outstream1 {"diffjunk01a"};

  cout << "Enter the number of tests to run\n";
  cin >> testnum;
  cout << testnum << " tests will be run\n";
  for(testct=0; testct<testnum; testct++) {
    startpos.x = (2.0l*rand()/RAND_MAX - 1.0l)*3.0l*AU_KM;
    startpos.y = (2.0l*rand()/RAND_MAX - 1.0l)*3.0l*AU_KM;
    startpos.z = (2.0l*rand()/RAND_MAX - 1.0l)*3.0l*AU_KM;

    heliodist = vecabs3d(startpos);
    vesc = sqrt(2.0*MGsun/heliodist);
    
    startvel.x = (2.0l*rand()/RAND_MAX - 1.0l)*vesc;
    startvel.y = (2.0l*rand()/RAND_MAX - 1.0l)*vesc;
    startvel.z = (2.0l*rand()/RAND_MAX - 1.0l)*vesc;

    v = vecabs3d(startvel);

    while(v>=vesc) {
      startvel.x = (2.0l*rand()/RAND_MAX - 1.0l)*vesc;
      startvel.y = (2.0l*rand()/RAND_MAX - 1.0l)*vesc;
      startvel.z = (2.0l*rand()/RAND_MAX - 1.0l)*vesc;

      v = vecabs3d(startvel);
    }

    //cout << fixed << setprecision(6) << mjdstart << " " << mjdend << " " << startpos.x << " " << startpos.y  << " " << startpos.z  << " " << startvel.x << " " << startvel.y  << " " << startvel.z  << "\n";
    Kepler_fg_func_int(MGsun, mjdstart, startpos, startvel, mjdend, outpos1, outvel1);
    //cout << fixed << setprecision(6) << mjdstart << " " << mjdend << " " << startpos.x << " " << startpos.y  << " " << startpos.z  << " " << startvel.x << " " << startvel.y  << " " << startvel.z  << "\n";
    Kepler_univ_int(MGsun, mjdstart, startpos, startvel, mjdend, outpos2, outvel2);

    cout << fixed << setprecision(6) << outpos1.x << " " << outpos1.y  << " " << outpos1.z  << " " << outvel1.x << " " << outvel1.y  << " " << outvel1.z  << "\n";
    cout << fixed << setprecision(6) << outpos2.x << " " << outpos2.y  << " " << outpos2.z  << " " << outvel2.x << " " << outvel2.y  << " " << outvel2.z  << "\n";
    double u = dotprod3d(startpos,startvel);
    double v1 = vecabs3d(startvel);
    double a = heliodist*MGsun/(2.0l*MGsun - v1*v1*heliodist);
    double n,e;
    if(a>0.0) {
      n = sqrt(MGsun/a/a/a);
      e = sqrt(DSQUARE(1.0l-heliodist/a) + DSQUARE(u/n/a/a));
    } else {
      e = sqrt(DSQUARE(1.0l-heliodist/a) - DSQUARE(u/sqrt(-MGsun*a)));
    }
    posdiff.x = outpos1.x - outpos2.x;
    posdiff.y = outpos1.y - outpos2.y;
    posdiff.z = outpos1.z - outpos2.z;
    veldiff.x = outvel1.x - outvel2.x;
    veldiff.y = outvel1.y - outvel2.y;
    veldiff.z = outvel1.z - outvel2.z;
    outstream1 << vecabs3d(posdiff) << " " << vecabs3d(veldiff)  << " " << a/AU_KM << " " << e << " " << startpos.x << " " << startpos.y << " " << startpos.z << " " << startvel.x << " " << startvel.y << " " << startvel.z << "\n";
  }
  outstream1.close();

  return(0);
}
