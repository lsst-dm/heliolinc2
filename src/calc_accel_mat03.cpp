// October 27, 2023: calc_accel_mat03:
// Uses the insight of calc_accel_mat02b that acceleration sampling
// can be inferred from velocity sampling, but does not use the analogous
// inference from distance to velocity. Velocity sampling is therefore
// decoupled from distance sampling, and now instead is coupled to the
// clustering radius, which is now a parameter.
//Also takes the total time spanned
// by the data as the timescale input, as opposed to the use in calc_accel_mat02b
// of a 'characteristic time' somewhat confusingly defined as 1/4 the total timespan.

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
  cerr << "Usage: calc_accel_mat03 -mindist minimum distance (AU) -maxdist maximum distance (AU) -distfracstep \ndistances sampling interval as a fraction of geodist -clustrad clustering radius (km) -velscale undersampling factor for velocity -totaltimespan total time spanned by the data (days) -outfile output file\n";
  cerr << "\n\nDespite lack of valid inputs, we still run with reasonable defaults for everything:\n\n";
}

int main(int argc, char *argv[])
{
  double mindist = 1.0l + FIXRANGE;
  double maxdist = 100.0l;
  double distfracstep = 0.1l;
  double diststep = 0.01;
  double velstep = 20.0l;
  double vel_undersamp = 1.0l;
  double clustrad = 100000.0l;
  double totaltime = 16.0l;
  string outfile = "camout02.csv";
  long i=0;
  long stepct=0;
  int verbose=0;
  ofstream outstream1;
  double dist,velkm,velAU,vesc,accelk,accelnorm,g0;
  dist = velkm = velAU = vesc = accelk = accelnorm = g0 = 0.0l;
  double distkm=0.0l;
  double geodist=0.0l;
  double accelstep,minacc,maxacc,vtan,xi;
  accelstep = minacc = maxacc = vtan = xi = 0.0l;
  long accelnum,accelct,velnum,velct;

  if(argc<3) show_usage();
  
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
    } else if(string(argv[i]) == "-clustrad" || string(argv[i]) == "-cr" || string(argv[i]) == "-crad" || string(argv[i]) == "-cluster" || string(argv[i]) == "--clustrad" || string(argv[i]) == "--clusterradius" || string(argv[i]) == "--clusterrad") {
      if(i+1 < argc) {
	//There is still something to read;
	clustrad=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Clustering radius keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-velscale" || string(argv[i]) == "-velunder" || string(argv[i]) == "-velfac" || string(argv[i]) == "-vs" || string(argv[i]) == "-vf" ) {
      if(i+1 < argc) {
	//There is still something to read;
	vel_undersamp = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Velocity undersampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-totaltimespan" || string(argv[i]) == "-totaltime" || string(argv[i]) == "-timespan") {
      if(i+1 < argc) {
	//There is still something to read;
	totaltime = stod(argv[++i]); // The total time spanned by the data.
	i++;
      }
      else {
	cerr << "Total time span keyword supplied with no corresponding argument\n";
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

  cout << "Distances will be probed from " << mindist << " to " << maxdist << "AU\n";
  cout << "With a sampling equal to " << distfracstep << " of the minimum inferred geocentric distance\n";
  cout << "Total time span is " << totaltime << " days\n";
  cout << "Clustering radius is " << clustrad << " km\n";
  cout << "Velocity undersampling factor is " << vel_undersamp << "\n";
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
  } else if(totaltime<=0.0l) {
    cout << "\nERROR: Total time span must be strictly positive \n";
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
    geodist = fabs(dist-1.0l);
    if(geodist<FIXRANGE) geodist = FIXRANGE;
    diststep = distfracstep*geodist; // AU
    distkm = dist*AU_KM;
    vesc = sqrt(2.0l*GMSUN_KM3_SEC2/distkm); // km/sec
    velstep = 2.0l*vel_undersamp*clustrad*geodist/totaltime/SOLARDAY; // km/sec
    g0 = GMSUN_KM3_SEC2/distkm/distkm; // Solar gravity in km/sec^2
    accelstep = 4.0l*velstep/totaltime/SOLARDAY; 
    if(verbose>=1) cout << "Distance is " << dist << " AU = " << distkm << " km, vesc = " << vesc << " km/sec, velstep = " << velstep << "km/sec, gt = " << g0*totaltime*SOLARDAY << "\n";
    velnum = ceil(2.0l*vesc/velstep);
    if(velnum<=1) {
      // Velocity step is so large, we'll only probe radial velocity=0
      // Since radial velocity is zero, tangential velocity can be anywhere from 0 to vesc.
      // WE TAKE POSITIVE ACCELERATION AS INWARD TOWARD THE SUN
      maxacc = g0;
      minacc = -g0;
      velkm = 0.0;
      velAU = velkm*SOLARDAY/AU_KM;
      if(verbose>=1) cout << "vel = " << velAU << " accelstep = " << accelstep << ", range is " << minacc << " to " << maxacc << "\n";
      accelnum = ceil((maxacc-minacc)/accelstep);
      if(accelnum<=1) {
	// Acceleration step is so large, we'll only use average acceleration
	outstream1 << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	stepct++;
      } else {
	// Customize the acceleration sampling to have approximately the nominal
	// step size, but optimally sample the full range with the specified number of steps
	for(accelct=0;accelct<accelnum;accelct++) {
	  accelk = minacc + (double(accelct)+0.5l)*(maxacc-minacc)/double(accelnum);
	  accelnorm = accelk/g0;
	  outstream1 << dist << " " << velAU << " " << accelnorm << "\n";
	  if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << accelnorm << "\n";
	  stepct++;
	}
      }
    } else {
      // Customize the velocity sampling to have approximately the nominal
      // step size, but optimally sample the full range with the specified number of steps.
      // Now for the loop, probing everything except the endpoints
      for(velct=0;velct<velnum;velct++) {
	velkm = -vesc + (double(velct)+0.5)*2.0l*vesc/double(velnum);
	velAU = velkm*SOLARDAY/AU_KM;
	xi = velkm/vesc;
	maxacc = g0;
	minacc = g0*(1.0l - 2.0l*sqrt(1.0l - xi*xi));
	if(verbose>=1) cout << "vel = " << velAU << " accelstep = " << accelstep << ", range is " << minacc << " to " << maxacc << "\n";
	accelnum = ceil((maxacc-minacc)/accelstep);
	if(accelnum<=1) {
	  // Acceleration step is so large, we'll only use average acceleration
	  outstream1 << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	  if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	  stepct++;
	} else {
	  // Customize the acceleration sampling to have approximately the nominal
	  // step size, but optimally sample the full range with the specified number of steps
	  for(accelct=0;accelct<accelnum;accelct++) {
	    accelk = minacc + (double(accelct)+0.5l)*(maxacc-minacc)/double(accelnum);
	    accelnorm = accelk/g0;
	    outstream1 << dist << " " << velAU << " " << accelnorm << "\n";
	    if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << accelnorm << "\n";
	    stepct++;
	  }
	}
      }
    }
    dist += diststep;
  }
	
  return(0);
}
