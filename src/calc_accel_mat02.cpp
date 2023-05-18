// May 08, 2023: calc_accel_mat02

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
  cerr << "Usage: calc_accel_mat02 -mindist minimum distance (AU) -maxdist maximum distance (AU) -distfracstep \ndistances sampling interval as a fraction of geodist\n-velstep velocity sampling interval in km/sec at 1AU geodist -timescale characteristic timescale in days (used for acceleration step) -outfile output file\n";
}

int main(int argc, char *argv[])
{
  double mindist = 1.0l + FIXRANGE;
  double maxdist = 100.0l;
  double distfracstep = 0.1l;
  double velstep = 20.0l;
  double timescale = 10.0l;
  string outfile = "camout02.csv";
  long i=0;
  long stepct=0;
  int verbose=0;
  ofstream outstream1;
  double dist,velkm,velAU,vesc,accelk,accelnorm,g0;
  dist = velkm = velAU = vesc = accelk = accelnorm = g0 = 0.0l;
  double distkm=0.0l;
  double geodist=0.0l;
  double nvelstep,ndiststep,accelstep,minacc,maxacc,vtan,xi;
  nvelstep = ndiststep = accelstep = minacc = maxacc = vtan = xi = 0.0l;
  long accelnum,accelct,velnum,velct;
  
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-mindist" || string(argv[i]) == "-mind" || string(argv[i]) == "-mindistance" || string(argv[i]) == "-minimumdist" || string(argv[i]) == "--minimum_distance" || string(argv[i]) == "--mindist" || string(argv[i]) == "--mind" || string(argv[i]) == "--mindist") {
      if(i+1 < argc) {
	//There is still something to read;
	mindist = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input minimum heliocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxdist" || string(argv[i]) == "-maxd" || string(argv[i]) == "-maxdistance" || string(argv[i]) == "-maximumdist" || string(argv[i]) == "--maximum_distance" || string(argv[i]) == "--maxdist" || string(argv[i]) == "--maxd" || string(argv[i]) == "--maxdist") {
      if(i+1 < argc) {
	//There is still something to read;
	maxdist = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input maximum heliocentric distance keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-velstep" || string(argv[i]) == "-vstep") {
      if(i+1 < argc) {
	//There is still something to read;
	velstep = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Velocity sampling interval keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timescale" || string(argv[i]) == "-chartime") {
      if(i+1 < argc) {
	//There is still something to read;
	timescale = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Characteristic timescale keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "--outfile" || string(argv[i]) == "--out" || string(argv[i]) == "--outfile") {
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

  cout << "Distances will be probed from " << mindist << " to " << maxdist << "AU\n";
  cout << "With a sampling equal to " << distfracstep << " of the minimum inferred geocentric distance\n";
  cout << "Velocity sampling will be " << velstep << " km/sec, scaled by the minimum inferred geocentric distance in AU\n";
  cout << "Characteristic timescale will be " << timescale << " days\n";
  cout << "Output file to be written is called " << outfile << "\n";

  // Catch required parameters if missing
  if(mindist<=0.0l) {
    cout << "\nERROR: minimum distance must be strictly positive\n";
    show_usage();
    return(1);
  } else if(maxdist<=0.0l || maxdist<=mindist) {
    cout << "\nERROR: maximum distnace must be strictly greater than minimum distance\n";
    show_usage();
    return(1);
  } else if(distfracstep<=0.0l) {
    cout << "\nERROR: Fractional distance step must be strictly positive \n";
    show_usage();
    return(1);
  } else if(velstep<=0.0l) {
    cout << "\nERROR: Velocity step must be strictly positive \n";
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

  outstream1.open(outfile);
  outstream1 << "#r(AU) rdot(AU/day) mean_accel\n";
  dist = mindist;
  while(dist<=maxdist && stepct<=STEPMAX) {
    distkm = dist*AU_KM;
    vesc = sqrt(2.0l*GMSUN_KM3_SEC2/distkm); // km/sec
    geodist = fabs(dist-1.0l);
    if(geodist<FIXRANGE) geodist = FIXRANGE;
    nvelstep = velstep*geodist; // km/sec
    g0 = GMSUN_KM3_SEC2/distkm/distkm; // Solar gravity in km/sec^2
    accelstep = nvelstep/timescale/SOLARDAY; // Acceleration required to change velocity
                                             // by one velocity step over the characteristic timescale
    if(verbose>=1) cout << "Distance is " << dist << " AU = " << distkm << " km, vesc = " << vesc << " km/sec, velstep = " << nvelstep << "km/sec, gt = " << g0*timescale*SOLARDAY << "\n";
    if(nvelstep > 2.0l*vesc) {
      // Velocity step is so large, we'll only probe radial velocity=0
      // Since radial velocity is zero, tangential velocity can be anywhere from 0 to vesc.
      // WE TAKE POSITIVE ACCELERATION AS INWARD TOWARD THE SUN
      maxacc = g0;
      minacc = -g0;
      velkm = 0.0;
      velAU = velkm*SOLARDAY/AU_KM;
      if(verbose>=1) cout << "vel = " << velAU << " accelstep = " << accelstep << ", range is " << minacc << " to " << maxacc << "\n";
      if(accelstep > (maxacc-minacc)) {
	// Acceleration step is so large, we'll only use average acceleration
	outstream1 << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	stepct++;
      } else {
	// Customize the acceleration sampling to have approximately the nominal
	// step size, but span the full range with an integer number of steps.
	accelnum = ceil((maxacc-minacc)/accelstep);
	for(accelct=0;accelct<=accelnum;accelct++) {
	  accelk = minacc + double(accelct)*(maxacc-minacc)/double(accelnum);
	  accelnorm = accelk/g0;
	  outstream1 << dist << " " << velAU << " " << accelnorm << "\n";
	  if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << accelnorm << "\n";
	  stepct++;
	}
      }
    } else {
      // Customize the velocity sampling to have approximately the nominal
      // step size, but span the full range with an integer number of steps.
      velnum = ceil(2.0l*vesc/nvelstep);
      // Deal with the full escape velocity outside the loop.
      // There is only one possible acceleration in this case: full GMsun/r^2
      velAU = -vesc*ESCSCALE*SOLARDAY/AU_KM; // ESCSCALE=0.99 prevents roundoff error from producing formally unbound orbits.
      outstream1 << dist << " " << velAU << " " << 1.0 << "\n";
      if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << 1.0 << "\n";
      stepct++;
      // Now for the loop, probing everything except the endpoints
      for(velct=1;velct<velnum;velct++) {
	velkm = -vesc + double(velct)*2.0l*vesc/double(velnum);
	velAU = velkm*SOLARDAY/AU_KM;
	xi = velkm/vesc;
	maxacc = g0;
	minacc = g0*(1.0l - 2.0l*sqrt(1.0l - xi*xi));
	if(verbose>=1) cout << "vel = " << velAU << " accelstep = " << accelstep << ", range is " << minacc << " to " << maxacc << "\n";
	if(accelstep > (maxacc-minacc)) {
	  // Acceleration step is so large, we'll only use average acceleration
	  outstream1 << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	  if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	  stepct++;
	} else {
	  // Customize the acceleration sampling to have approximately the nominal
	  // step size, but to span the full range with an integer number of steps.
	  accelnum = ceil((maxacc-minacc)/accelstep);
	  for(accelct=0;accelct<=accelnum;accelct++) {
	    accelk = minacc + double(accelct)*(maxacc-minacc)/double(accelnum);
	    accelnorm = accelk/g0;
	    outstream1 << dist << " " << velAU << " " << accelnorm << "\n";
	    if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << accelnorm << "\n";
	    stepct++;
	  }
	}
      }
      // Deal with the full escape velocity outside the loop.
      // There is only one possible acceleration in this case: full GMsun/r^2
      velAU = vesc*ESCSCALE*SOLARDAY/AU_KM; // ESCSCALE=0.99 prevents roundoff error from producing formally unbound orbits.
      outstream1 << dist << " " << velAU << " " << 1.0 << "\n";
      if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << 1.0 << "\n";
      stepct++;
    }
    dist += distfracstep*geodist;
  }
	
  return(0);
}
