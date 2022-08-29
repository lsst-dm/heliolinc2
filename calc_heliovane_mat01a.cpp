// August 03, 2022: calc_heliovane_mat01a.cpp: 
// Create a matrix of possible values of heliocentric ecliptic longitude
// and its first two time derivatives (angular velocity along the ecliptic,
// and the corresponding angular accleration). The ecliptic longitude
// itself is set to evenly spaced values from 0 to 360 degrees. The
// velocity is set within a user-defined range in units of degrees per day,
// but if this range includes points that could not be reached in
// a bound orbit (require a tangential velocity exceeding the solar
// escape velocity at the specified minimum heliocentric distance),
// the program will exit with an error. The user also sets the step size
// for acceleration, in deg/day^2. There is an expectation that the 
// range of acceleration will be +/- 2 Omega*sqrt(2GMsun/r^3 - Omega^2),
// where Omega is the angular velocity. This acceleration forms a very
// generous envelope enclosing nearly all physically possible cases until
// inclination exceeds 70 degrees, and a large fraction of them even
// at higher inclinations. If the user wishes to use a narrower range,
// the parameter accelscale can be set to less than 1.0. In the unlikely
// event that the user wishes a wider range in acceleration, the
// parameter may of course be set to more than 1.0. Regardless of
// the acceleration step size and range, an acceleration of 0.0
// will be probed. Large acceleration ranges will not produce any
// errors, since there are rare cases for which the acceleration can
// in fact become arbitrarily large.
//
//
// Descriptions of related programs:
// January 17, 2022: calc_accel_mat01b.cpp:
// Like calc_accel_mat01a.cpp, but instead of printing a specific
// set of acceleration quantiles, prints out specified min and max
// quantiles (e.g. 5th to 95th percentiles) plus evenly sampled
// steps between them. Note that the steps are even in the actual
// acceleration, not even quantiles of the distribution (which is
// what calc_accel_mat01a gives). The units of the steps are GMsun/r^2.
// The even sampling in physical accleration rather than in quantiles
// of the distribution ensures we can configure the sampling of
// distance, velocity, acceleration quantiles so that we capture
// even rare objects with unusual orbits, rather than being sensitive
// only to objects in more common orbits.
//
// Description of ancestor program calc_accel_mat01a.cpp:
// Probe the search space of the heliolinc algorithm (that is,
// heliocentric distance and radial vecloity), and for each grid
// point in this space, calculate quantiles in acceleration for
// asteroids that appear at that point (that is, will at some
// point in their orbits have the corresponding heliocentric
// distance and radial velocity). Write a file providing the
// distance/velocity grid points and their associated heliocentric
// radial accelerations, for use as input in heliolinc runs using
// projectpairs04c.cpp.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MINDIST 0.59L // Default minimum distance in AU.
#define VELSTART -3.0 // Default minimum angular velocity in deg/day 
#define VELSTEP 0.5 // Default angular velocity step size in deg/day
#define VELNUM 13 // Default number of angular velocities to probe
#define ACCSCALE 1.0 // Default scaling for acceleration range, in terms
                     // of the nominal upper bound on acceleration.
#define ACCELSTEP 0.1 // Default acceleration step in units of deg/day^2
                      // Note: at 0.59 AU, the nominal upper bound on
                      // acceleration is 0.23 deg/day^2




static void show_usage()
{
  cerr << "Usage: calc_heliovane_mat01a -outfile output_file -mindist mindist -velstep velstart velstep velnum -accelstep accelstep accelscale\n";
}
    
int main(int argc, char *argv[])
{
  string outfile;
  ofstream outstream1;
  long double MGsun = GMSUN_KM3_SEC2;
  long double mindist = MINDIST;
  long double velstart = VELSTART;
  long double velstep = VELSTEP;
  int velnum = VELNUM;
  int i=0;
  int velct;
  int act,accelnum;
  long double accelscale = ACCSCALE;
  long double accelstep = ACCELSTEP;
  long double maxvel,maxangvel;
  long double accelrange;
  long double angvel,angacc;
  
  if(argc<3) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-outfile" || string(argv[i]) == "-of" || string(argv[i]) == "-outf" || string(argv[i]) == "--outputfile" || string(argv[i]) == "--ofile" || string(argv[i]) == "--fileout") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "output filename keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mindist" || string(argv[i]) == "-md" || string(argv[i]) == "-distmin" || string(argv[i]) == "-dmin" || string(argv[i]) == "--mindist" || string(argv[i]) == "--distmin" || string(argv[i]) == "--distancemin") {
      if(i+1 < argc) {
	//There is still something to read;
	mindist=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Minimum distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-velstep" || string(argv[i]) == "-vs" || string(argv[i]) == "-vstep" || string(argv[i]) == "-vel" || string(argv[i]) == "--velstep" || string(argv[i]) == "--velocitystep" || string(argv[i]) == "--velocity") {
      if(i+1 < argc) {
	//There is still something to read;
	velstart=stold(argv[++i]);
      }
      else {
	cerr << "Velocity step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	velstep=stold(argv[++i]);
      }
      else {
	cerr << "Velocity step keyword supplied with only one of the three required arguments\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	velnum=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Velocity step keyword supplied with only two of the three required arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-accstep" || string(argv[i]) == "-acc"  || string(argv[i]) == "-accelstep" || string(argv[i]) == "-as" || string(argv[i]) == "--acceleration" || string(argv[i]) == "--accelstep" || string(argv[i]) == "--accelerationstep") {
      if(i+1 < argc) {
	//There is still something to read;
	accelstep=stold(argv[++i]);
      }
      else {
	cerr << "output acceleration step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	accelscale=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Acceleration step keyword supplied with only one of the two required arguments\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
    
  cout << "output file name " << outfile << "\n";
  cout << "minimum distance " << mindist << " AU\n";
  cout << "minimum and maximum angular velocities " << velstart << " and " << velstart + velstep*double(velnum-1) << " deg/day\n";
  cout << "angular velocity step size " << velstep << " deg/day\n";
  cout << "number of angular velocity steps " << velnum << "\n";
  cout << "acceleration scale factor " << accelscale << "\n";
  cout << "acceleration step size " << accelstep << " deg/day^2\n";

  maxvel = sqrt(2.0L*MGsun/mindist/AU_KM);
  cout << "Escape velocity at minimum distance is " << maxvel << " km/sec\n";
  maxangvel = maxvel/mindist/AU_KM*SOLARDAY*DEGPRAD;
  cout << "This corresponds to an angular velocity of " << maxangvel << " deg/day\n";
  if(velstart<-maxangvel) {
    cout << "ERROR: starting angular velocity " << velstart << " is out of the physical range\n";
    return(1);
  }
  if(velstart + velstep*double(velnum-1) > maxangvel) {
    cout << "ERROR: ending angular velocity " << velstart + velstep*double(velnum-1) << " is out of the physical range\n";
    return(1);
  }

  accelrange = 2*maxangvel*sqrt(2.0L*MGsun/intpowLD(mindist*AU_KM,3))*SOLARDAY;
  cout << "Maximum angular acceleration should be " << accelrange << " deg/day^2\n";
  
  outstream1.open(outfile,ios_base::out);

  outstream1 << "#longitude_vel(deg/day) longitude_acc(deg/day^2)\n";
  for(velct=0;velct<velnum;velct++) {
    angvel = velstart + velstep*double(velct);
    accelrange = accelscale*2*angvel*sqrt(2.0L*MGsun/intpowLD(mindist*AU_KM,3))*SOLARDAY;
    accelnum = accelrange/accelstep;
    for(act = -accelnum ; act <= accelnum ; act++) {
      angacc = accelstep*double(act);
      outstream1 << fixed << setprecision(6) << angvel << " " << angacc << "\n";
    }
  }
  outstream1.close();
  return(0);
}
