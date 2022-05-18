// May 02, 2022: interstellar_sim01a.cpp
// Fit an input observations file using the downhill simplex method.
// The observation file must contain MJD, RA, Dec,
// astrometric uncertainty, obscode.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define NUMPOS 3
#define MJDCOL 3
#define RACOL 6
#define DECCOL 8
#define MINOBSINTERVAL 1.0 // Minimum time-between-images in seconds
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME 1.5 // Default max inter-image time interval
                    // for tracklets, in hours (will be converted
                    // to days before use).
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)
#define DEBUG 0
#define PHASE_G 0.15
#define IDNUMLEN 7

// Note: configfile contains the masses and ephemerides for all of the planets.
// configfile also contains the name of an observatory code file. The input
// units of the state vector are AU for positions and km/sec for velocities.
// Internally, the program will use AU for positions and AU/timescale for
// velocities. The default value of timescale is one year, but it can be
// set to other values as needed. The observation file must contain
// MJD, RA, Dec, magnitude, astrometric uncertainty, magnitude uncertainty,
// and obscode. The RA and Dec must be in decimal degrees. The astrometric
// uncertainty must be in arcseconds.
static void show_usage()
{
  cerr << "Usage: calc_accel_interstellar01a -cfg configfile -ranseed random_number_seed -simnum simnum -diststep dist0 diststep distnum -velstep vel0 velstep velnum -accelstep minquant maxquant accelstep -timestep timestep -outfile outfile \n or, at minimum: \ncalc_accel_interstellar01a -cfg configfile -ranseed random_number_seed -simnum simnum -outfile outfile \n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  int status=0;
  ifstream instream1;
  ofstream outstream1;
  vector <heliogridpoint> gridvec;
  heliogridpoint hgp = heliogridpoint(0L,0L,{});
  string stest;
  string configfile;
  string outfile;
  long double x,y,z,vx,vy,vz;
  x=y=z=vx=vy=vz=0.0L;
  double usigma = 0.0l;
  double vsigma = 0.0l;
  double wsigma = 0.0l;
  double umean = 0.0l;
  double vmean = 0.0l;
  double wmean = 0.0l;
  double uvel = 0.0l;
  double vvel = 0.0l;
  double wvel = 0.0l;
  double GalRA = 0.0l;
  double GalDec = 0.0l;
  double bmax = 0.0l;
  double acoef,xmin,xmax,vinf,encounter_dist;
  long simnum=0;
  long simct=0;
  int configread=0;
  long double ldval=0.0L;
  double dval=0.0L;
  int reachedeof=0;
  double impactpar;
  long double dist;
  long double E,lscalar,a,e,coshH,H0,omega,t0;
  long double random_number;
  long double timestep = 86400.0L;
  long double tnow = 0.0L;
  vector <long double> hvel;
  vector <long double> hrad;
  vector <long double> hacc;
  long double drdt,dvdt;
  long tnum;
  string seedstring;
  int goodsimct=0;
  int evergood=0;
  long double GMsun = GMSUN_KM3_SEC2;
  double minquant=0.0;
  double maxquant=1.0;
  long double mindist = 1.0L;
  long double maxdist = 1.0L;
  long double minvel = -0.02L;
  long double diststep = 0.1L;
  long double velstep = -0.001L;
  int distnum = 91;
  int distct = 0;
  int velnum = 41;
  int velct = 0;
  double accelstep = 0.1;
  int accelstepnum=0;
  double minaccel=0.0l;
  double maxaccel=0.0l;
  vector <long double> accelvec;
  long double r1,vrad,ldc,lvc,accel,fullgrav;
 
  if(argc<9) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-c" || string(argv[i]) == "-cfg" || string(argv[i]) == "-config" || string(argv[i]) == "--config" || string(argv[i]) == "--configfile" || string(argv[i]) == "--configuration" || string(argv[i]) == "--ConfigFile") {
      if(i+1 < argc) {
	//There is still something to read;
	configfile=argv[++i];
	i++;
	// Read the configuration file. This must happen here, so that
	// default values supplied in this file can be overwritten later
	// if the user desires.
	// Read configuration file.
	ifstream instream1 {configfile};
	if(!instream1) {
	  cerr << "ERROR: can't open input config file " << configfile << "\n";
	  return(1);
	}
	// Read Gaussian sigma on the distribution of the U component of the velocity
	status=readconfigd(instream1,&usigma);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&usigma);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Gaussian sigma on U velocity component read as " << usigma << " km/sec\n";
	// Read Gaussian sigma on the distribution of the V component of the velocity
	status=readconfigd(instream1,&vsigma);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&vsigma);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Gaussian sigma on V velocity component read as " << vsigma << " km/sec\n";
	// Read Gaussian sigma on the distribution of the W component of the velocity
	status=readconfigd(instream1,&wsigma);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&wsigma);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Gaussian sigma on W velocity component read as " << wsigma << " km/sec\n";
	// Read the mean U componenent of the velocity
	status=readconfigd(instream1,&umean);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&umean);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Mean U component of the velocity component read as " << umean << " km/sec\n";
	// Read the mean V component of the velocity
	status=readconfigd(instream1,&vmean);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&vmean);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Mean V component of the velocity component read as " << vmean << " km/sec\n";
	// Read the mean W component of the velocity
	status=readconfigd(instream1,&wmean);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&wmean);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Mean W component of the velocity component read as " << wmean << " km/sec\n";
	// Read maximum impact parameter in AU
	status=readconfigd(instream1,&bmax);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&bmax);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Maximum impact parameter read as " << bmax << " AU\n";
	// Read default minimum distance
	status=readconfigLD(instream1,&mindist);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mindist);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default minimum distance read as " << mindist << " AU\n";
	// Read default step size for distance
	status=readconfigLD(instream1,&diststep);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&diststep);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default step size for distance read as " << diststep << " AU\n";
	// Read default number of distance samples
	status=readconfigint(instream1,&distnum);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&distnum);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default number of distance samples read as " << distnum << "\n";
	// Read default minimum velocity
	status=readconfigLD(instream1,&minvel);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&minvel);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default minimum velocity read as " << minvel << " AU/day\n";
	// Read default step size for velocity
	status=readconfigLD(instream1,&velstep);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&velstep);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default step size for velocity read as " << velstep << " AU/day\n";
	// Read default number of velocity samples
	status=readconfigint(instream1,&velnum);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&velnum);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default number of velocity samples read as " << velnum << "\n";
	// Read default minimum quantile for acceleration
	status=readconfigd(instream1,&minquant);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&minquant);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default minimum quantile for acceleration read as " << minquant << "\n";
	// Read default maximum quantile for acceleration
	status=readconfigd(instream1,&maxquant);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&maxquant);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default maximum quantile for acceleration read as " << maxquant << "\n";
	// Read default acceleration step size
	status=readconfigd(instream1,&accelstep);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&accelstep);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default acceleration step read as " << accelstep << " GMsun/r^2\n";

	// Close input stream that was reading the config file.
	instream1.close();
	configread=1;
      } else {
	cerr << "Configuration file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }
    // Stop execution if we haven't read a config file successfully by now.
    if(configread!=1) {
      cerr << "ERROR: configuration file must be supplied and successfully";
      cerr << "read before any of the other parameters. This is because\n";
      cerr << "it contains defaults that can optionally be overridden by\n";
      cerr << "user-supplied values in later arguments.\n";
      return(2);
    } else cout << "Configuration file read successfully\n";
    if(string(argv[i]) == "-rs" || string(argv[i]) == "-ranseed" || string(argv[i]) == "-rseed" || string(argv[i]) == "-seed" || string(argv[i]) == "--ranseed" || string(argv[i]) == "--randomseed" || string(argv[i]) == "--randomnumberseed" || string(argv[i]) == "--random_number_seed") {
      if(i+1 < argc) {
	// There is still something to read;
	seedstring=argv[++i];
	i++;
      } else {
	cerr << "Random number seed keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-diststep" || string(argv[i]) == "-ds" || string(argv[i]) == "-dstep" || string(argv[i]) == "-dist" || string(argv[i]) == "--diststep" || string(argv[i]) == "--distancestep" || string(argv[i]) == "--distance") {
      if(i+1 < argc) {
	//There is still something to read;
	mindist=stold(argv[++i]);
      }
      else {
	cerr << "Distance step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	diststep=stold(argv[++i]);
      }
      else {
	cerr << "Distance step keyword supplied with only one of the three required arguments\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	distnum=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Distance step keyword supplied with only two of the three required arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-velstep" || string(argv[i]) == "-vs" || string(argv[i]) == "-vstep" || string(argv[i]) == "-vel" || string(argv[i]) == "--velstep" || string(argv[i]) == "--velocitystep" || string(argv[i]) == "--velocity") {
      if(i+1 < argc) {
	//There is still something to read;
	minvel=stold(argv[++i]);
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
    }  else if(string(argv[i]) == "-accstep" || string(argv[i]) == "-acc"  || string(argv[i]) == "-accelstep" || string(argv[i]) == "-as" || string(argv[i]) == "--acceleration" || string(argv[i]) == "--accelstep" || string(argv[i]) == "--accelerationstep") {
      if(i+1 < argc) {
	//There is still something to read;
	minquant=stold(argv[++i]);
      }
      else {
	cerr << "output acceleration step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	maxquant=stold(argv[++i]);
      }
      else {
	cerr << "Acceleration step keyword supplied with only one of the three required arguments\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	accelstep=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Acceleration step keyword supplied with only two of the three required arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timestep" || string(argv[i]) == "-ts" || string(argv[i]) == "-tstep" || string(argv[i]) == "-time" || string(argv[i]) == "--timestep" || string(argv[i]) == "--time" || string(argv[i]) == "--tempstep") {
      if(i+1 < argc) {
	//There is still something to read;
	timestep=stold(argv[++i]);
      }
      else {
	cerr << "Time step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-simnum" || string(argv[i]) == "-sn" || string(argv[i]) == "-snum" || string(argv[i]) == "--simnum") {
      if(i+1 < argc) {
	//There is still something to read;
	simnum=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Simulation number keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outorb" || string(argv[i]) == "--outorbits") {
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
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
    
  maxdist = mindist + ((long double)(distnum-1))*diststep;
    
  cout.precision(17);  
  cout << "randum number seed string " << seedstring << "\n";
  cout << "input configuration file " << configfile << "\n";
  cout << "input gaussian sigma values for U, V, W velocities " << usigma << " " << vsigma << " " << wsigma << " km/sec\n";
  cout << "input mean values for U, V, W velocities " << umean << " " << vmean << " " << wmean << " km/sec\n";
  cout << "number of encounters to simulate: " << simnum << "\n";
  cout << "output file " << outfile << "\n";

  
  // Allocate grid vector with empty acceleration vectors
  gridvec={};
  for(distct=0; distct<distnum; distct++) {
    for(velct=0; velct<velnum; velct++) {
      hgp = heliogridpoint(mindist+diststep*(long double)distct, minvel+velstep*(long double)velct,{});
      gridvec.push_back(hgp);
      //cout << "Loaded point " << distct*velnum+velct << " = " << gridvec.size()-1 << "\n";
    }
  }
  
  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  // Loop on simulated objects.
  for(simct=0; simct<simnum; simct++) {

    // 3. Assign v_infinity and b_infinity (the initial impact parameter).
    // Assign U, V, W velocities from a Gaussian approximation
    // to the distribution of local disk stars
    uvel = umean + usigma*gaussian_deviate_mt(generator);
    vvel = vmean + vsigma*gaussian_deviate_mt(generator);
    wvel = wmean + wsigma*gaussian_deviate_mt(generator);
    // Now we have a realistic sun-relative velocity for an interstellar object.
    // Calculate total velocity (vinf means infinite distance
    // from the sun: i.e., not yet altered by the sun's gravity).
    vinf = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);
    //cout << "vinf = " << vinf << " km/sec\n";
    // assign the initial impact parameter.
    impactpar = sqrt(bmax*bmax*unitvar(generator));
    // This generates a differential distribution with probability
    // density proportional to the impact parameter, which is
    // appropriate because the differential area element is r*dr
    // Previously I had used this one:
    // impactpar = pow(bmax*bmax*bmax*unitvar(generator),0.333333);
    // It gives probability density proportional to the square of the
    // impact parameter, but is not correct because it is the cumulative
    // distribution that goes as impactpar^2. The differential
    // distribution, which is the relevant one here, is linearly
    // proportional to impactpar.
  
    // 4. Calculate encounter distance
    encounter_dist = impactpar*impactpar*AU_KM*AU_KM*vinf*vinf/(GMsun + sqrt(GMsun*GMsun + vinf*vinf*vinf*vinf*impactpar*impactpar*AU_KM*AU_KM));
    //cout << "encounter_dist = " << encounter_dist/AU_KM << " AU\n";
     
   // 5. Reject the object if the encounter distance is greater than the maximum distance
    // at which the object could be detected: i.e., there is no chance of a detection.
    if(encounter_dist/AU_KM<=maxdist) {
      // The object will get close enough to the sun to be of interest.
      E = 0.5L*vinf*vinf;
      lscalar = vinf*impactpar*AU_KM;
      //cout << "Energy = " << E << ", lscalar = " << lscalar << "\n";

      a = -GMsun*0.5L/E; // By convention, semimajor axis a is negative
                     // for hyperbolic orbits.
      e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/GMsun/GMsun);
      cout << "Simulated object number " << simct << " of " << simnum << " a = " << a/AU_KM << ", e = " << e << "\n";

      // Determine the value of H0 for which the heliocentric distance
      // is equal to maxdist
      omega = sqrt(-GMsun/(a*a*a));
      coshH = (a-maxdist*AU_KM)/(a*e);
      if(coshH>=1.0L) H0 = -acosh(coshH);
      else {
	cerr << "WARNING: cosh(H) = " << coshH << "\n";
	return(1);
      }
      t0 = (e*sinh(H0) - H0)/omega;
      //cout << "time since perihelion = " << t0 << " seconds = " << t0/86400.0L << " days = " << t0/31556952.0l << " years\n";
      // Setup for Keplerian integration spanning the period when the object
      // is within a distance maxdist of the Sun.
      tnum=0;
      dist = maxdist*AU_KM;
      hrad = hvel = hacc = {};
      //cout << "maxdist = " << maxdist << "\n";
      while(tnum==0 || dist/AU_KM < maxdist) {
	hrad.push_back(dist);
	drdt = -a*e*omega*sinh(H0)/(e*cosh(H0)-1.0L);
	dvdt = (-a*e*e*omega*omega + a*e*omega*omega*cosh(H0))/intpowLD(e*cosh(H0)-1.0L,3);
	hvel.push_back(drdt);
	hacc.push_back(dvdt);
	tnum++;
	tnow = t0+((long double)tnum)*timestep;
	H0=hyp_transcendental(tnow*omega,e,KEPTRANSTOL);
	dist = a*(1.0L - e*cosh(H0));
      }
      //cout << "Loaded ephemeris vectors with " << hrad.size() << " = " << hvel.size() << " " << hacc.size() << " elements\n";
      
      for(i=0;i<hrad.size();i++) {
	if(i>0 && i<hrad.size()-1) {
	  //cout << "point " << i << " dist = " << hrad[i]/AU_KM << " vel = " << hvel[i] << " = " << 0.5*(hrad[i+1]-hrad[i-1])/timestep << " acc = " << hacc[i] << " = " << 0.5*(hvel[i+1]-hvel[i-1])/timestep << " = " << ((hrad[i+1]-hrad[i])/timestep - (hrad[i]-hrad[i-1])/timestep)/timestep << " fullgrav " << GMsun/hrad[i]/hrad[i] << "\n";
	}
	// Convert heliocentric distance to AU
	r1 = hrad[i]/AU_KM;
	// Convert heliocentric velocity to AU/day
	vrad = hvel[i]*SOLARDAY/AU_KM;
	// Find which gridpoint we're looking at.
	ldc = (r1 - mindist+0.5*diststep)/diststep;
	lvc = (vrad - minvel+0.5*velstep)/velstep;
	distct = ldc;
	velct = lvc;
	if(distct>=0 && distct<distnum && velct>=0 && velct<velnum) {
	  // Scale acceleration relative to GMsun/r^2
	  fullgrav = -GMsun/hrad[i]/hrad[i];
	  //cout << "Accel = " << hacc[i] << " fullgrav = " << fullgrav << " ratio = ";
	  accel = hacc[i] / fullgrav;
	  //cout << accel << "\n";
	  gridvec[distct*velnum+velct].acc.push_back(accel);
	}
      }
    }
  }
  // Open output file and write header
  outstream1.open(outfile);
  outstream1 << "#r(AU) rdot(AU/day) norm mean_accel mean_eccen mean_semimajor mean_perihelion\n";
  for(distct=0; distct<distnum; distct++) {
    for(velct=0; velct<velnum; velct++) {
      //cout << "distct = " << distct << " velct = " << velct << " accelsize = " << gridvec[distct*velnum+velct].acc.size() << "\n";
      if(gridvec[distct*velnum+velct].acc.size() > 0) {
	// Load vector of accelerations
	accelvec={};
	for(i=0; i<gridvec[distct*velnum+velct].acc.size(); i++) {
	  accelvec.push_back(gridvec[distct*velnum+velct].acc[i]);
	  //cout << "accel = " << gridvec[distct*velnum+velct].acc[i] << "\n";
	}
	// Sort vector of accelerations.
	sort(accelvec.begin(), accelvec.end());
	// Extract requested quantiles;
	j=minquant*accelvec.size();
	minaccel = accelvec[j];
	j=maxquant*accelvec.size();
	maxaccel = accelvec[j];
	// Round to nearest step.
	minaccel = accelstep*round(minaccel/accelstep);
	maxaccel = accelstep*round(maxaccel/accelstep);
	accelstepnum = round((maxaccel-minaccel)/accelstep);
	//cout << "minaccel = " << minaccel << " maxaccel = " << maxaccel << " accelstepnum = " << accelstepnum << "\n";
	for(j=0;j<=accelstepnum;j++) {
	  outstream1 << fixed << setprecision(6) << gridvec[distct*velnum+velct].dist << " " << gridvec[distct*velnum+velct].vel << " " << accelvec.size() <<  " " << minaccel+double(j)*accelstep << " -99.9 -99.9 -99.9\n";
	}
      }
    }
  }
  outstream1.close();

  return(0);
}




