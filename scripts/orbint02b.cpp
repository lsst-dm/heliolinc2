// December 07, 2021: orbint02b.cpp: 
// Test different routines for integrating orbits in the Solar System,
// using the new library solarsyst_dyn_geo01.h.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

int main()
{
  string stest,planetfile;
  int planetnum=0;
  int planetct=0;
  int polyorder=1;
  int i=0;
  vector <long double> planetmasses;
  vector <point3LD> planetpos;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  point3LD targpos1 = point3LD(0,0,0);
  point3LD targpos2 = point3LD(0,0,0);
  point3LD targvel1 = point3LD(0,0,0);
  point3LD targvel2 = point3LD(0,0,0);
  vector <long double> mjd;
  vector <long double> mjdtest;
  long double mjdstart=0.0;
  long double mjdend=0.0;
  long runnum,runct;
  runnum = runct=1;

  cout << "How many orbit integrations do you want to run?\n";
  cin >> stest;
  runnum=stol(stest);
  cout << "Integration will be run " << runnum << " times.\n";

  cout << "What polynomial order do you want to use for dynamical integrations?\n";
  cin >> stest;
  polyorder=stoi(stest);
  cout << "Input poly order is " << polyorder << "\n";
  
  // Find out how many massive perturbers we are going to use.
  cout << "How many massive bodies will you use?\n";
  cin >> stest;
  planetnum = stoi(stest);
  cout << "You have requested " << planetnum << " perturbing bodies.\n";
    
  // Read masses and positions for perturbing bodies from the user.
  for(planetct=0; planetct<planetnum; planetct++) {
    cout << "Enter the product GM for planet number " << planetct+1 << ", in km^3/sec^2\n";
    cin >> stest;
    planetmasses.push_back(stold(stest));
  }
  for(planetct=0; planetct<planetnum; planetct++) {
    cout << "Enter the ephemeris file for planet " << planetct+1 << "\n";
    cin >> planetfile;
    mjdtest={};
    temppos={};
    tempvel={};
    read_horizons_fileLD(planetfile,mjdtest,temppos,tempvel);
    if(planetct==0) mjd=mjdtest;
    else {
      for(i=0;i<mjd.size();i++) {
	if(mjdtest[i]!=mjd[i]) {
	  cout << "ERROR: time vectors do not match for input planet files\n";
	  cout << planetct+1 << " and 1!\n";
	  return(1);
	}
      }
    }
    for(i=0;i<temppos.size();i++) {
      planetpos.push_back(temppos[i]);
    }
  }
     
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
  if(mjdend>=mjdstart) {
    // We are integrating forward in time
    if(polyorder<1) {
      for(runct=0;runct<runnum;runct++) {
	integrate_orbit_constac(planetnum, mjd, planetmasses, planetpos, mjdstart, targpos1, targvel1, mjdend, targpos2, targvel2);
      }
    } else if(polyorder<2) {
      for(runct=0;runct<runnum;runct++) {
	integrate_orbit01LD(planetnum, mjd, planetmasses, planetpos, mjdstart, targpos1, targvel1, mjdend, targpos2, targvel2);
      }
    }
    else {
      for(runct=0;runct<runnum;runct++) {
	integrate_orbit02LD(polyorder,planetnum, mjd, planetmasses, planetpos, mjdstart, targpos1, targvel1, mjdend, targpos2, targvel2);
      }
    }
  }
  else {
    // We are going backward in time.
    // Reverse the order and the sign of the mjd vector.
    for(i=0;i<mjd.size();i++) mjdtest[mjd.size()-i-1] = -mjd[i];
    for(i=0;i<mjd.size();i++) mjd[i] = mjdtest[i];
    // Reverse the order of the planet vectors
    for(planetct=0; planetct<planetnum; planetct++) {
      for(i=0;i<mjd.size();i++) {
	temppos[mjd.size()-i-1] = planetpos[planetct*mjd.size() + i];
      }
      for(i=0;i<mjd.size();i++) {
	planetpos[planetct*mjd.size() + i] = temppos[i];
      }
    }
    // Reverse the direction of the input velocity
    targvel1.x = -targvel1.x;
    targvel1.y = -targvel1.y;
    targvel1.z = -targvel1.z;
    // Reverse the sign of the start and end times
    mjdstart = -mjdstart;
    mjdend = -mjdend;
    if(polyorder<1) {
      for(runct=0;runct<runnum;runct++) {
	integrate_orbit_constac(planetnum, mjd, planetmasses, planetpos, mjdstart, targpos1, targvel1, mjdend, targpos2, targvel2);
      }
    } else if(polyorder==1) {
      for(runct=0;runct<runnum;runct++) {
	integrate_orbit01LD(planetnum, mjd, planetmasses, planetpos, mjdstart, targpos1, targvel1, mjdend, targpos2, targvel2);
      }
    } else {
      for(runct=0;runct<runnum;runct++) {
	integrate_orbit02LD(polyorder,planetnum, mjd, planetmasses, planetpos, mjdstart, targpos1, targvel1, mjdend, targpos2, targvel2);
      }
    }
  }
  
  cout << fixed << setprecision(6) << "Target position: " << targpos2.x << " "  << targpos2.y << " "  << targpos2.z << "\n";
  cout << fixed << setprecision(6) << "Target velocity: " << targvel2.x << " "  << targvel2.y << " "  << targvel2.z << "\n";
  
  return(0);
}
