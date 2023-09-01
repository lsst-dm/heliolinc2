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
  point3d outpos = point3d(0.0L,0.0L,0.0L);
  point3d outvel = point3d(0.0L,0.0L,0.0L);
  double timediff = 0.0l;
  int testnum=1;
  int testct=0;
  //double a;
  
  cout << "Enter the number of tests to run\n";
  cin >> testnum;
  cout << testnum << " tests will be run\n";
  for(testct=0; testct<testnum; testct++) {
    cout << "Enter the starting x,y,z position in km\n";
    cin >> startpos.x;
    cin >> startpos.y;
    cin >> startpos.z;
    cout << startpos.x << " " << startpos.y << " " << startpos.z << "\n";
    cout << "Enter the ending x,y,z position in km\n";
    cin >> endpos.x;
    cin >> endpos.y;
    cin >> endpos.z;
    cout << endpos.x << " " << endpos.y << " " << endpos.z << "\n";
    cout << "Enter the time-difference in days\n";
    cin >> timediff;
    cout << timediff << "\n";

    //    int status = Twopoint_Kepler_vel(MGsun, startpos, endpos, timediff, startvel, &a, KEPTRANSITMAX);
    int status = Twopoint_Kepler_vstar(MGsun, startpos, endpos, timediff, startvel, KEPTRANSITMAX);
    if(status!=0) {
      cerr << "Warning: Twopoint_Kepler_vstar() returned error status " << status << "\n";
    }
    cout << "Twopoint_Kepler_vstar() finds initial velocity " << startvel.x << " " << startvel.y << " " << startvel.z << "\n";
    Kepler_univ_int(MGsun, 0.0l, startpos, startvel, timediff, outpos, outvel);
    cout << "Integrating this orbit produces an ending position of:\n";
    cout << outpos.x << " " << outpos.y << " " << outpos.z << ", which may be compared with:\n";
    cout << endpos.x << " " << endpos.y << " " << endpos.z << ": differences are:\n";
    cout << outpos.x-endpos.x << " " << outpos.y-endpos.y << " " << outpos.z-endpos.z << "\n";
  }

  return(0);
}
