// October 13, 2023: calc_vanemat01:
// Like calc_accel_mat02b, but creates hypothesis matrices for
// heliovane, rather than heliolinc. That is, the fundamental variable
// is not heliocentric distance, but heliocentric ecliptic longitude.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define FIXRANGE 0.10l // Distance from Earth, in AU, at which sampling
                      // intervals are no longer scaled with geocentric distance
#define STEPMAX 10000000
#define ESCSCALE 0.99l // Velocities nominally equal to the escape velocity will
                      // be scaled by this factor, to ensure that the orbits
                      // remain formally bound.

static void show_usage()
{
  cerr << "Usage: calc_vanemat01 -minheliodist minimum distance from sun (AU) -mingeodist min distance from Earth (AU) -distfracstep \ndistances sampling interval as a fraction of geodist -timescale characteristic timescale in days (typically 1/4 the total time spanned by the input data) -minsunelong minimum solar elongation (deg) -outfile output file\n";
  cerr << "\n\nDespite lack of valid inputs, we still run with reasonable defaults for everything:\n\n";
}

int main(int argc, char *argv[])
{
  double minheliodist = 0.5l;
  double mingeodist = FIXRANGE;
  double distfracstep = 0.1l;
  double diststep = 0.01;
  double timescale = 4.0l; // Expected to be 1/4 the total time spanned by the input data.
  double minsunelong = 40.0l; // Not many surveys are effective within 40 degrees of the sun.
  string outfile = "testvanemat01.csv";
  long i=0;
  int verbose=0;
  double finerange=0.0l;
  double lambda,lambda_dot,lambda_ddot,lamstep,dot_step,ddot_step;
  lambda = lambda_dot = lambda_ddot = lamstep = dot_step = ddot_step = 0.0l;
  double v_esc=0.0l;
  double max_angvel=0.0l;
  double max_angacc=0.0l;
  int velnum=0;
  int velct=0;
  int accelnum=0;
  int accelct=0;
  double coarselimit = 135.0l;
  double dist = 0.0l;
  ofstream outstream1;  
  
  if(argc<3) show_usage();
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-mingeodist" || string(argv[i]) == "-mingeo" || string(argv[i]) == "-mingeodistance" || string(argv[i]) == "-minimumgeodist" || string(argv[i]) == "--minimum_geocentric_distance" || string(argv[i]) == "--mingeodist" || string(argv[i]) == "--mingeo" || string(argv[i]) == "--mingeocendist") {
      if(i+1 < argc) {
	//There is still something to read;
	mingeodist = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input minimum geocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minheliodist" || string(argv[i]) == "-minhelio" || string(argv[i]) == "-minheliodistance" || string(argv[i]) == "-minimumheliodist" || string(argv[i]) == "--minimum_heliocentric_distance" || string(argv[i]) == "--minheliodist" || string(argv[i]) == "--minhelio" || string(argv[i]) == "--minheliocendist") {
      if(i+1 < argc) {
	//There is still something to read;
	minheliodist = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input minimum geocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-distfracstep" || string(argv[i]) == "-fracstep" || string(argv[i]) == "-diststep" || string(argv[i]) == "-dfs" || string(argv[i]) == "--distfracstep" ) {
      if(i+1 < argc) {
	//There is still something to read;
	distfracstep = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Fractional distance step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timescale" || string(argv[i]) == "-chartime") {
      if(i+1 < argc) {
	//There is still something to read;
	timescale = stod(argv[++i]); // Expected to be 1/4 the total time spanned by the data.
	i++;
      }
      else {
	cerr << "Characteristic timescale keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minsunelong" || string(argv[i]) == "-minelong") {
      if(i+1 < argc) {
	//There is still something to read;
	minsunelong = stod(argv[++i]); 
	i++;
      }
      else {
	cerr << "Characteristic timescale keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--out" || string(argv[i]) == "--outfile") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Verbosity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  cout << "Minimum geocentric distance will be " << mingeodist << "AU\n";
  cout << "Minimum heliocentric distance will be " << minheliodist << "AU\n";
  cout << "With a sampling equal to " << distfracstep << " of the minimum inferred geocentric distance\n";
  cout << "Characteristic timescale will be " << timescale << " days\n";
  cout << "Output file to be written is called " << outfile << "\n";

  // Catch required parameters if missing
  if(mingeodist<=0.0l) {
    cout << "\nERROR: minimum geocentric distance must be strictly positive\n";
    show_usage();
    return(1);
  } else if(minheliodist<=0.0l) {
    cout << "\nERROR: minimum heliocentric distance must be strictly positive\n";
    show_usage();
    return(1);
  } else if(distfracstep<=0.0l) {
    cout << "\nERROR: Fractional distance step must be strictly positive \n";
    show_usage();
    return(1);
  } else if(timescale<=0.0l) {
    cout << "\nERROR: Characteristic timescale must be strictly positive \n";
    show_usage();
    return(1);
  } else if(outfile.size()<=0) {
    cout << "\nERROR: output filename is required\n";
    show_usage();
    return(1);
  }
  if(sin(minsunelong/DEGPRAD)>minheliodist) {
    cerr << "WARNING: at specified minimum solar elongation of " << minsunelong << " degrees,\n";
    cerr << "the specified minimum heliocentric distance (" << minheliodist << " AU) will never be accessible\n";
    minheliodist = sin(minsunelong/DEGPRAD);
    cerr << "Resetting miniumum heliocentric distance to " << minheliodist << " AU\n";
  }
  v_esc = sqrt(2.0l*GMSUN_KM3_SEC2/(minheliodist*AU_KM)); // km/sec
  max_angvel = (v_esc/(minheliodist*AU_KM))*DEGPRAD*SOLARDAY; // degrees/day
  cout << "Maximum escape velocity is " << v_esc << " km/sec,\nleading to angular velocity " << max_angvel << " deg/day.\n";
  
  outstream1.open(outfile);
  outstream1 << "#lambda(deg) lambda_dot(deg/day) lambda_ddot(deg/day^2)\n";
  outstream1 << fixed << setprecision(6);
  // We take the zeropoint of heliocentric ecliptic longitude to be the longitude
  // of Earth at the reference time, and we approximate Earth's orbital motion at
  // 1 degree per day (for purposes of crudely setting a limit on the range of
  // fine sampling).
  finerange = 2.0l*timescale + mingeodist*DEGPRAD;
  coarselimit = 180.0l - minsunelong;
  if(verbose>=1) cout << "coarse range +/-" << coarselimit << ", fine range +/-" << finerange << "\n";
  lambda = -coarselimit;
  while(lambda<=coarselimit) {
    if(fabs(lambda) <= finerange) {
      // We use a constant, maximally fine sampling inside finerange
      diststep = distfracstep*mingeodist; // in AU, nominally
      lamstep = diststep*DEGPRAD;
      cout << "Inside finerange\n";
    } else {
      // What is the minimum distance the object could be from Earth?
      dist = 2.0l*sin(fabs(lambda)/2.0l/DEGPRAD);
      if(dist<mingeodist) dist=mingeodist;
      diststep = distfracstep*dist; // in AU, nominally
      lamstep = diststep*DEGPRAD;
      cout << "outside finerange\n";
    }
    lambda_dot = -max_angvel;
    dot_step = 2.0l*lamstep/timescale;
    velnum = ceil(2.0l*max_angvel/dot_step);
    if(verbose>=1) cout << "lambda = " << lambda << ", lamstep = " << lamstep << ", max_angvel = " << max_angvel << "\n";
    if(velnum<=1) {
      // Angular velocity step is so large, we'll only probe angular velocity=0.0;
      // This also means angular acceleration zero.
      lambda_dot = lambda_ddot = 0.0;
      outstream1 << lambda << " " << lambda_dot << " " << lambda_ddot << "\n";
    } else {
      // Customize the angular velocity sampling to have approximately the nominal
      // step size, but optimally sample the full range with the specified number of steps.
      // Now for the loop, probing everything except the endpoints
      for(velct=0;velct<velnum;velct++) {
	lambda_dot = -max_angvel + (double(velct)+0.5)*2.0l*max_angvel/double(velnum);
	max_angacc = fabs(lambda_dot)*sqrt(2.0l*GMSUN_KM3_SEC2/intpowD(minheliodist*AU_KM,3.0)); // deg/day/sec
	max_angacc *= SOLARDAY; // deg/day^2
	ddot_step = 4.0l*dot_step/timescale;
	accelnum = ceil(2.0l*max_angacc/ddot_step);
      
	if(accelnum<=1) {
	  // Acceleration step is so large, we'll only use zero acceleration.
	  lambda_ddot = 0.0;
	  outstream1 << lambda << " " << lambda_dot << " " << lambda_ddot << "\n";
	} else {
	  // Customize angular acceleration sampling to have approximately the nominal
	  // step size, but optimally sample the full range with the specified
	  // number of steps.
	  for(accelct=0;accelct<accelnum;accelct++) {
	    lambda_ddot = -max_angacc + (double(accelct)+0.5)*2.0l*max_angacc/double(accelnum);
	    outstream1 << lambda << " " << lambda_dot << " " << lambda_ddot << "\n";
	  }
	}
      }
    }
    lambda += lamstep;
  }
  outstream1.close();	
  return(0);
}
