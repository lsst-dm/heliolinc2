// March 07, 2024: calc_heliohypmat.cpp
// Like the calc_accel_mat programs, but uses a very different approach
// to the problem of calculating a hypothesis matrix for heliolinc.
// Works from empirical power laws.

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
  cerr << "Usage: calc_heliohypmat -mindist minimum distance (AU) -maxdist maximum distance (AU) -distconst -distpwr -velconst -velpwr -accconst -accpwr -accchangerad -clustchangerad -outfile output file\n";

}

int main(int argc, char *argv[])
{
  double mindist = 1.0l + FIXRANGE;
  double maxdist = 100.0l;
  double diststep = 1.0l;
  double velstep = 1.0e-2;
  double distconst = 0.258336;
  double distpwr = 1.369798;
  double velconst = 2.713272e-3;
  double velpwr = 0.704204;
  double accconst = 7.180432e-7;
  double accpwr = 0.0778692;
  double accchangerad = 1.5; // geocentric distance beyond which the
                             // acceleration sampling no longer varies.
  double clustchangerad = 0.5; // geocentric distance within which the
                              // clustering radius no longer changes.
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
  string outfile;
  
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
    } else if(string(argv[i]) == "-distconst" || string(argv[i]) == "-distconst" || string(argv[i]) == "-distconst" || string(argv[i]) == "-distconst" || string(argv[i]) == "--distconst" ) {
      if(i+1 < argc) {
	//There is still something to read;
	distconst = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Constant coefficient for distance power law keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-distpwr" || string(argv[i]) == "-distpwr" || string(argv[i]) == "-distpwr" || string(argv[i]) == "-distpwr" || string(argv[i]) == "--distpwr" || string(argv[i]) == "--distpwr" || string(argv[i]) == "--distpwr") {
      if(i+1 < argc) {
	//There is still something to read;
	distpwr=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Exponent for distance power law keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-velconst" || string(argv[i]) == "-velconst" || string(argv[i]) == "-velconst" || string(argv[i]) == "-velconst" || string(argv[i]) == "-velconst" ) {
      if(i+1 < argc) {
	//There is still something to read;
	velconst = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Constant coefficient for velocity power law keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-velpwr" || string(argv[i]) == "-velpwr" || string(argv[i]) == "-velpwr") {
      if(i+1 < argc) {
	//There is still something to read;
	velpwr = stod(argv[++i]); // The total time spanned by the data.
	i++;
      }
      else {
	cerr << "Exponent for velocity power law keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-accconst" || string(argv[i]) == "-accconst" || string(argv[i]) == "-accconst" || string(argv[i]) == "-accconst" || string(argv[i]) == "-accconst" ) {
      if(i+1 < argc) {
	//There is still something to read;
	accconst = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Constant coefficient for acceleration power law keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-accpwr" || string(argv[i]) == "-accpwr" || string(argv[i]) == "-accpwr") {
      if(i+1 < argc) {
	//There is still something to read;
	accpwr = stod(argv[++i]); 
	i++;
      }
      else {
	cerr << "Exponent for acceleration power law keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-accchangerad" || string(argv[i]) == "-accchangerad" || string(argv[i]) == "-accchangerad") {
      if(i+1 < argc) {
	//There is still something to read;
	accchangerad = stod(argv[++i]); 
	i++;
      }
      else {
	cerr << "Transition distance for acceleration power law keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clustchangerad" || string(argv[i]) == "-clustchangerad" || string(argv[i]) == "-clustchangerad") {
      if(i+1 < argc) {
	//There is still something to read;
	clustchangerad = stod(argv[++i]); 
	i++;
      }
      else {
	cerr << "Transition distance for cluster scaling keyword supplied with no corresponding argument\n";
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
    diststep = distconst*pow(geodist,distpwr); // AU
    distkm = dist*AU_KM;
    vesc = sqrt(2.0l*GMSUN_KM3_SEC2/distkm)*SOLARDAY/AU_KM; // AU/day
    velstep = velconst*pow(geodist,velpwr); // AU/day
    g0 = GMSUN_KM3_SEC2/distkm/distkm; // Solar gravity in km/sec^2
    if(geodist<accchangerad) {
      accelstep = accconst*pow(geodist,accpwr); // km/sec^2
    } else {
      accelstep = accconst*pow(accchangerad,accpwr);
    }
    // Scale all the steps according to the clustering radius
    if(geodist<clustchangerad) {
      diststep *= clustchangerad;
      velstep *= clustchangerad;
      accelstep *= clustchangerad;
    } else {
      diststep *= geodist;
      velstep *= geodist;
      accelstep *= geodist;
    }
      
    if(verbose>=1) cout << "Distance is " << dist << " AU = " << distkm << " km, vesc = " << vesc*AU_KM/SOLARDAY << " km/sec, velstep = " << velstep*AU_KM/SOLARDAY << "km/sec, accelstep = " << accelstep << " km/sec^2, GMsun = " << g0 << " km/sec^2\n";
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
	velAU = -vesc + (double(velct)+0.5)*2.0l*vesc/double(velnum);
	velkm = velAU*AU_KM/SOLARDAY;
	xi = velAU/vesc;
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
