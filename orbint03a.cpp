// December 07, 2021: orbint03a.cpp: 
// Working towards fast 2-body orbit propagation based on
// Goldstein and the Kepler Equation.

#include "std_lib_facilities.h"
#include "solarsyst_dyn_geo01.h"
#include "cmath"

int main()
{
  string stest,planetfile;
  int planetnum=0;
  int planetct=0;
  int polyorder=1;
  int i=0;
  long double GMsun;
  vector <point3LD> sunpos;
  vector <point3LD> sunvel;
  point3LD sunposnow = point3LD(0,0,0);
  point3LD sunvelnow = point3LD(0,0,0);
  point3LD targpos1 = point3LD(0,0,0);
  point3LD tp1 = point3LD(0,0,0);
  point3LD targpos2 = point3LD(0,0,0);
  point3LD targvel1 = point3LD(0,0,0);
  point3LD tv1 = point3LD(0,0,0);
  point3LD targvel2 = point3LD(0,0,0);
  vector <long double> mjd;
  long double mjdstart=0.0;
  long double mjdend=0.0;
  long runnum,runct;
  runnum = runct=1;

  cout << "How many orbit integrations do you want to run?\n";
  cin >> stest;
  runnum=stol(stest);
  cout << "Integration will be run " << runnum << " times.\n";

  cout << "Enter the product GM for the Sun, in km^3/sec^2\n";
  cin >> stest;
  GMsun = stold(stest);

  cout << "Enter the ephemeris file for the Sun\n";
  cin >> planetfile;
  sunpos={};
  sunvel={};
  read_horizons_fileLD(planetfile,mjd,sunpos,sunvel);
     
  // Get starting time, position, and velocity for the target object.
  cout << "Enter the starting mjd, 3-D barycentric position,\nand 3-D barycentric velocity for your target object\n";
  cin >> stest;
  mjdstart = stold(stest);
  cin >> stest;
  targpos1.x = stold(stest);
  cin >> stest;
  targpos1.y = stold(stest);
  cin >> stest;
  targpos1.z = stold(stest);
  cin >> stest;
  targvel1.x = stold(stest);
  cin >> stest;
  targvel1.y = stold(stest);
  cin >> stest;
  targvel1.z = stold(stest);

  cout << "Enter the time to which we will integrate the orbit of\n";
  cout << "your target object.\n";
  cin >> stest;
  mjdend = stold(stest);
  
  for(runct=0;runct<runnum;runct++) {
    planetposvel01LD(mjdstart,3, mjd, sunpos, sunvel, sunposnow, sunvelnow);
    tp1.x = targpos1.x - sunposnow.x;
    tp1.y = targpos1.y - sunposnow.y;
    tp1.z = targpos1.z - sunposnow.z;
    tv1.x = targvel1.x - sunvelnow.x;
    tv1.y = targvel1.y - sunvelnow.y;
    tv1.z = targvel1.z - sunvelnow.z;
    Keplerint(GMsun,mjdstart,tp1,tv1,mjdend,targpos2,targvel2);
    planetposvel01LD(mjdend,3, mjd, sunpos, sunvel, sunposnow, sunvelnow);
    targpos2.x += sunposnow.x;
    targpos2.y += sunposnow.y;
    targpos2.z += sunposnow.z;
    targvel2.x += sunvelnow.x;
    targvel2.y += sunvelnow.y;
    targvel2.z += sunvelnow.z;
  }
  cout << fixed << setprecision(6) << "Target position: " << targpos2.x << " "  << targpos2.y << " "  << targpos2.z << "\n";
  cout << fixed << setprecision(6) << "Target velocity: " << targvel2.x << " "  << targvel2.y << " "  << targvel2.z << "\n";

  return(0);
}
