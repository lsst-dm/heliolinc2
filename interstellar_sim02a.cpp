// July 19, 2022: interstellar_sim02a.cpp
// Like interstellar_sim01a.cpp, but instead of going straight for the
// production of simulated observations, creates a catalog of orbital
// elements for simulated objects of interest, in the format MJD,
// barycentric state vectors (X, Y, Z, VX, VY, VZ), absolute magnitude.
//
// Description of ancestor program:
// May 02, 2022: interstellar_sim01a.cpp
// Creates simulated observations of interstellar objects, given an input
// image catalog with MJD and boresight RA, Dec for each image. Models the
// population of interstellar objects as a power law in absolute magnitude,
// with v_infinity relative to the sun distributed according to a Gaussian
// ellipsoid in terms of the Galactic U, V, and W velocities. The central
// coordinate and the widths (Gaussian sigmas) of this ellipsoid in the
// U, V, and W velocity coordinates are set by the user. Given a U, V, W
// velocity selected randomly from this distribution, the program also
// randomly selects an impact parameter (b_infinity) distributed according
// to the area of concentric rings P(b) = C0 b db. The actual encounter
// distance is calculated from b_infinity and the total velocity, and if
// the resulting value is close enough to the sun that the object might
// be detectable given its absolute magnitude, the U, V, and W velocity
// and the impact parameter are converted into a hyperbolic Keplerian orbit
// starting from a set (very large, e.g., 2000 AU) distance from the sun.
// Start time is adjusted so the perihelion will occur within some prescribed range
// of time. The resulting hyperbolic Keplerian orbit is then re-evaluated
// at a prescribed start time, presumed to be just before the earliest
// possible perihelion time, so that the object will then be much closer
// to the Sun. The evaluation of the hyperbolic Keplerian orbit yields new
// barycentric state vectors, which are then propagated forward with
// full n-body integration to produce a finely sampled ephemeris spanning
// all of the image times. This ephemeris is interpolated to predict exact
// positions and magnitudes (accounting for phase) on each image where the
// object is detected. Detection depends on the object's RA and Dec lying
// within a specified angular arc (image radius) of the image boresight,
// and the magnitude being brighter than a specified constant limiting
// magnitude.
//
// As an example of reasonable values for the U, V, W distribution,
// Charles Francis and Erik Anderson 2009, New Astronomy, 14, 615,
// find that most stars within 300 pc of the Sun form a local
// velocity ellipsoid (an ellipsoidal Gaussian distribution
// of velocities), with sigma_U = 32.6 km/sec, sigma_V = 22.4 km/sec,
// and sigma_W = 16.5 km/sec. The mean velocities relative to the Sun
// (that is, the center of the elliposoid) are -10.2, -18.3, and -7.4 km/sec,
// for U, V, and W respectively.
//

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define PHASE_G 0.15
#define IDNUMLEN 7
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds

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
  cerr << "Usage: interstellar_sim02a -cfg configfile -ranseed random_number_seed -mjdstart mjdstart -mjdend mjdend -simnum simnum -outfile outfile \n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  int status=0;
  ifstream instream1;
  string stest;
  string configfile;
  string planetfile;
  string outfile;
  long double x,y,z,vx,vy,vz;
  x=y=z=vx=vy=vz=0.0L;
  long double mjdstart = 0.0L;
  long double mjdend = 0.0L;
  long double mjd_perihelion_min = 0.0L;
  long double mjd_perihelion_max = 0.0L;
  long double mjd_perihelion = 0.0L;
  int planetfile_startpoint=0;
  int planetfile_endpoint=0;
  long double GMsun = 0.0L;
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
  double startdist = 0.0l;
  double bmax = 0.0l;
  double Hmin = 0.0l;
  double Hmax = 0.0l;
  double Hslope = 0.0l;
  double limiting_mag = 0.0l;
  double acoef,xmin,xmax,vinf,encounter_dist,absmag,maxdist;
  long simnum=0;
  long simct=0;
  int configread=0;
  int polyorder=0;
  int planetnum=0;
  int planetct=0;
  int pctSun=0;
  int pctEarth=0;
  long double ldval=0.0L;
  double dval=0.0L;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> Earthpos;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  string lnfromfile;
  vector <string> linestringvec;
  int badread=0;
  int reachedeof=0;
  point3LD outpos = point3LD(0,0,0);
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  point3LD hyppos = point3LD(0,0,0);
  point3LD hypvel = point3LD(0,0,0);
  point3LD unitbary = point3LD(0,0,0);
  point3LD lvec = point3LD(0L,0L,0L);
  double newRA,newDec,impactpar,maxdetdist;
  long double outRA,outDec;
  double racenter,deccenter,dist,pa,posRA,posDec,imRA,imDec;
  long double ldRA,ldDec,vesc;
  long double E,lscalar,a,e,coshH,H0,radvel,omega,t0,mjd_infall;
  int fieldnum=0;
  int MJDcol = 1;
  int RAcol = 2;
  int Deccol = 3;
  ofstream outstream1;
  long double random_number;
  long double astromsigma=0.1L;
  string seedstring;
  int goodsimct=0;
  int evergood=0;
  int idnl = IDNUMLEN;
  int compfac = 1;
  string idnumstring;
  
  if(argc<13) {
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
	// Read polyorder
	status=readconfigint(instream1,&polyorder);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&polyorder);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Polynomial order for orbit integration read as " << polyorder << "\n";
	// Read the number of planets.
	status=readconfigint(instream1,&planetnum);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&planetnum);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Number of planets read as " << planetnum << "\n";
	// Read the index of the Sun within the planet vectors
	status=readconfigint(instream1,&pctSun);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&pctSun);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "The Sun is planet number " << pctSun << "\n";
	// Read the index of Earth within the planet vectors
	status=readconfigint(instream1,&pctEarth);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&pctEarth);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Earth is planet number " << pctEarth << "\n";
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the mass for planet number planetct
	  status=readconfigLD(instream1,&ldval);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigLD(instream1,&ldval);
	  }
	  planetmasses.push_back(ldval);
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else {
	    cout << "MG for planet " << planetct << " read as " << planetmasses[planetct] << " km^3/sec^2\n";
	    if(planetct==pctSun) GMsun = planetmasses[planetct];
	  }
	}
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the ephemeris file for planet number planetct
	  status=readconfigstring(instream1,planetfile);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigstring(instream1,planetfile);
	  }
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else cout << "Ephemeris file for planet " << planetct << " is named " << planetfile << "\n";
	  mjdtest={};
	  temppos={};
	  tempvel={};
	  read_horizons_fileLD(planetfile,mjdtest,temppos,tempvel);
	  if(planetct==0) planetmjd=mjdtest;
	  else {
	    for(j=0;j<planetmjd.size();j++) {
	      if(mjdtest[j]!=planetmjd[j]) {
		cout << "ERROR: time vectors do not match for input planet files\n";
		cout << planetct+1 << " and 1!\n";
		return(1);
	      }
	    }
	  }
	  for(j=0;j<temppos.size();j++) {
	    planetpos.push_back(temppos[j]);
	  }
	  if(planetct == pctEarth) Earthpos = temppos;
	  if(planetct == pctSun) {
	    Sunpos = temppos;
	    Sunvel = tempvel;
	  }
	  cout << "Finished reading ephemeris file " << planetfile << "\n";
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

	// Read starting distance from the sun in AU
	status=readconfigd(instream1,&startdist);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&startdist);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Starting distance from the sun read as " << startdist << " AU\n";
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
	// Read default minimum H magnitude
	status=readconfigd(instream1,&Hmin);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&Hmin);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Minimum H magnitude read as " << Hmin << "\n";
	// Read default maximum H magnitude
	status=readconfigd(instream1,&Hmax);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&Hmax);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Maximum H magnitude read as " << Hmax << "\n";
	// Read default power law slope.
	status=readconfigd(instream1,&Hslope);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&Hslope);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "H magnitude power law slope read as " << Hslope << "\n";
	// Read default limiting magnitude
	status=readconfigd(instream1,&limiting_mag);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigd(instream1,&limiting_mag);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Limiting magnitude read as " << limiting_mag << "\n";
	// Read default starting MJD
	status=readconfigLD(instream1,&mjdstart);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjdstart);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default starting MJD read as " << mjdstart << "\n";
	// Read default ending MJD
	status=readconfigLD(instream1,&mjdend);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjdend);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default ending MJD read as " << mjdend << "\n";
	// Read earliest allowed time of perihelion
	status=readconfigLD(instream1,&mjd_perihelion_min);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjd_perihelion_min);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Earliest allowed time of perihelion read as MJD " << mjd_perihelion_min << "\n";
	// Read latest allowed time of perihelion
	status=readconfigLD(instream1,&mjd_perihelion_max);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjd_perihelion_max);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Latest allowed time of perihelion read as MJD " << mjd_perihelion_max << "\n";
	// Read astrometric error in arcseconds.
	// Note: this is vestigial from interstellar_sim01a, since astromsigma
	// is not used in the current program, but it is retained so that the
	// configuration files will be interchangeable.
	status=readconfigLD(instream1,&astromsigma);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&astromsigma);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Latest allowed time of perihelion read as MJD " << mjd_perihelion_max << "\n";
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
    } else if(string(argv[i]) == "-mjdstart" || string(argv[i]) == "-ms" || string(argv[i]) == "-MJDstart" || string(argv[i]) == "--mjdstart" || string(argv[i]) == "--MJDstart" || string(argv[i]) == "--ModifiedJulianDaystart" || string(argv[i]) == "--modifiedjuliandaystart") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdstart=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD start keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdend" || string(argv[i]) == "-me" || string(argv[i]) == "-MJDend" || string(argv[i]) == "--mjdend" || string(argv[i]) == "--MJDend" || string(argv[i]) == "--ModifiedJulianDayend" || string(argv[i]) == "--modifiedjuliandayend") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdend=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD end keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }   else if(string(argv[i]) == "-mjdend" || string(argv[i]) == "-me" || string(argv[i]) == "-MJDend" || string(argv[i]) == "--mjdend" || string(argv[i]) == "--MJDend" || string(argv[i]) == "--ModifiedJulianDayend" || string(argv[i]) == "--modifiedjuliandayend") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdend=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD end keyword supplied with no corresponding argument\n";
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

  cout.precision(17);  
  cout << "randum number seed string " << seedstring << "\n";
  cout << "input configuration file " << configfile << "\n";
  cout << "input starting MJD " << mjdstart << "\n";
  cout << "input ending MJD " << mjdend << "\n";
  cout << "Perihelion must occur between MJD " << mjd_perihelion_min << " and " << mjd_perihelion_max << "\n";
  cout << "input gaussian sigma values for U, V, W velocities " << usigma << " " << vsigma << " " << wsigma << " km/sec\n";
  cout << "input mean values for U, V, W velocities " << umean << " " << vmean << " " << wmean << " km/sec\n";
  cout << "starting heliocentric distance: " << startdist << " AU\n";
  cout << "min and max H magnitudes: " << Hmin << " " << Hmax << "\n";
  cout << "power law slope for H magnitude: " << Hslope << "\n";
  cout << "number of encounters to simulate: " << simnum << "\n";
  cout << "output file " << outfile << "\n";

  // Match mjdstart and mjdend to planet file.
  
  planetfile_startpoint = planetfile_endpoint = -99;
  for(j=0;j<planetmjd.size();j++) {
    if(fabs(planetmjd[j]-mjdstart) < IMAGETIMETOL/SOLARDAY) planetfile_startpoint = j;
    if(fabs(planetmjd[j]-mjdend) < IMAGETIMETOL/SOLARDAY) planetfile_endpoint = j;
  }
  if(planetfile_startpoint<0 || planetfile_endpoint<0) {
    cerr << "ERROR: mjdstart " << mjdstart << " and/or mjdend " << mjdend << " could not be matched\nto any timestep in the planet files\n";
    return(1);
  } else {
    cout << "mjdstart " << mjdstart << " corresponds to timestep " << planetfile_startpoint << " in the planet files\n";
    cout << "mjdend " << mjdend << " corresponds to timestep " << planetfile_endpoint << " in the planet files\n";
  }
  
  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  // Open output file and write header
  outstream1.open(outfile);
  outstream1 << "stringID,MJD,ISO-barycenter(J2000x)(km),ISO-barycenter(J2000y)(km),ISO-barycenter(J2000z)(km),ISO-barycenter(J2000vx)(km/s),ISO-barycenter(J2000vy)(km/s),ISO-barycenter(J2000vz)(km/s),absmag,uvel,vvel,wvel,vinf,impactpar,a,e,encounter_dist,mjd_perihelion\n";

  // Keep track of leading zeros for string IDs
  idnl = IDNUMLEN-1;
  compfac = 10;
  idnumstring={};
  for(i=0; i<idnl; i++) idnumstring.push_back('0');

  // Loop on simulated objects.
  for(simct=0; simct<simnum; simct++) {
    evergood=0; // Was object ever on an image? We assume no until proven otherwise.
    // 1. Randomly assign an absolute magnitude H
    acoef = 1.0l/(Hslope*log(10.0));
    xmin = exp(Hmin/acoef);
    xmax = exp(Hmax/acoef);
    x = xmin + (xmax-xmin)*unitvar(generator);
    absmag = acoef*log(x);
    outstream1.precision(10);  
    
    // 2. Calculate maximum distance at which this object could be
    // detected, under ideal circumstances.
    maxdist = sqrt(0.25l + sqrt(exp((limiting_mag-absmag)*EFOLDS_PER_MAG))) + 0.5l;
    cout << "maxdist = " << maxdist << "\n";

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
    cout << "vinf = " << vinf << " km/sec\n";
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
    
    if(impactpar >= startdist) {
      cerr << "ERROR: you have enabled impact parameters up to " << bmax << " AU, and\n";
      cerr << "one was actually chosen (" << impactpar << " AU) that is larger than\n";
      cerr << "the starting distance " << startdist << " AU.\nThis creates a logical and geometrical contradiction.\n";
      return(4);
    }
    // 4. Calculate encounter distance
    encounter_dist = impactpar*impactpar*AU_KM*AU_KM*vinf*vinf/(GMsun + sqrt(GMsun*GMsun + vinf*vinf*vinf*vinf*impactpar*impactpar*AU_KM*AU_KM));
     cout << "encounter_dist = " << encounter_dist/AU_KM << " AU\n";
     
   // 5. Reject the object if the encounter distance is greater than the maximum distance
    // at which the object could be detected: i.e., there is no chance of a detection.
    if(encounter_dist/AU_KM<=maxdist) {
      evergood=1;
      // The object has at least some chance of becoming detectable.
      // 6. Calculate the RA, Dec corresponding to the initial velocity (i.e., the radiant).
      // Turn the sun-relative velocities in the Galactic U,V,W system into an equivalent
      // unit vector in Galactic coordinates.
      uvw_to_galcoord(uvel,vvel,wvel,GalRA,GalDec);
      // Convert to celestial coordinates
      poleswitch02(GalRA,GalDec,NCPGAL_LON,NGPDEC,NGPRA,newRA,newDec);
      cout << "radiant coords: Galactic " << GalRA << " " << GalDec << " celestial " << newRA << " " << newDec << "\n";
      // Starting position should be near the radiant
      racenter = newRA-180.0l;
      while(racenter<0.0l) racenter+=360.0l;
      deccenter = -newDec;
      // But offset by an amount determined by the impact parameter
      // 7. Calculate angular offset arcsine(b_infinity/startdist).
      dist = DEGPRAD*asin(impactpar/startdist);
      
      // 8. Assign a random position angle, and calculate initial heliocentric RA, Dec.
      pa = unitvar(generator)*360.0l;
      arc2cel01(racenter,deccenter,dist,pa,posRA,posDec);
      cout << "dist = " << dist << ", pa = " << pa << ", position coords = " << posRA << " " << posDec << "\n";

      // 9. Scale up v_infinity slightly to account for solar gravity at startdist.
      vesc = sqrt(2.0l*GMsun/startdist/AU_KM);
      vinf = sqrt(vinf*vinf + vesc*vesc);
      
      // 10. Calculate starting state vectors r and v.
      // position
      ldRA = posRA;
      ldDec = posDec;
      celestial_to_stateunitLD(ldRA,ldDec,unitbary);
      startpos.x = unitbary.x*startdist*AU_KM;
      startpos.y = unitbary.y*startdist*AU_KM;
      startpos.z = unitbary.z*startdist*AU_KM;
      // velocity
      ldRA = newRA;
      ldDec = newDec;
      celestial_to_stateunitLD(ldRA,ldDec,unitbary);
      startvel.x = unitbary.x*vinf;
      startvel.y = unitbary.y*vinf;
      startvel.z = unitbary.z*vinf;
 
      // 11. Solve hyperbolic Keplerian orbit.

      E = 0.5L*vinf*vinf - GMsun/startdist/AU_KM;
      lvec = crossprod3LD(startpos,startvel);
      lscalar = sqrt(dotprod3LD(lvec,lvec));
      cout << "Energy = " << E << ", lscalar = " << lscalar << " = " << vinf*impactpar*AU_KM << "\n";

      a = -GMsun*0.5L/E; // By convention, semimajor axis a is negative
                     // for hyperbolic orbits.
      e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/GMsun/GMsun);
      cout << "Simulated object number " << simct << " of " << simnum << " a = " << a/AU_KM << ", e = " << e << "\n";

      // 12. Find time-to-perihelion.
      // Calculate the value of the hyperbolic anomaly H at the starting time
      if(e>1.0L) {
	coshH = (a-startdist*AU_KM)/(a*e);
	if(coshH>=1.0L) H0 = acosh(coshH);
	else {
	  cerr << "ERROR: cosh(H) = " << coshH << "\n";
	  return(1);
	}
      }
      radvel = dotprod3LD(startpos,startvel)/(startdist*AU_KM);
      cout << "H0 = " << H0 << ", radial velocity = " << radvel << " km/sec\n";
  
      if(radvel>=0) {
	// We are moving outward from perihelion: H0 will be correct
	;
      } else {
	// We are moving inward towards perihelion: H0 needs adjustment.
	H0 = -H0;
      }
      cout << "H0 = " << H0 << ", radial velocity = " << radvel << " km/sec\n";
      // Calculate time since perihelion using H0.
      omega = sqrt(-GMsun/(a*a*a));
      cout << "omega = " << omega << "\n";
      t0 = (e*sinh(H0) - H0)/omega;
      cout << "time since perihelion = " << t0 << " seconds = " << t0/86400.0L << " days = " << t0/31556952.0l << " years\n";
      
      // 13. Define the time corresponding to the starting state vectors
      //     such that the perihelion will occur within a specified interval.
      // Assign random time-of-perihelion within specified limits
      mjd_perihelion = mjd_perihelion_min + (mjd_perihelion_max - mjd_perihelion_min)*unitvar(generator);
      // Calculate start time in MJD given this value for the time of perihelion
      mjd_infall = mjd_perihelion + t0/86400.0L;

      // 14. Integrate forward to a time which, if it corresponds to mjdstart,
      //     will produce a perihelion at the required instant.
      // Integrate orbit from the MJD corresponding to the start of infall,
      // down to the mjdstart, the requested start time for the n-body integration.
      Hyper_Kepint(GMsun, mjd_infall, startpos, startvel, mjdstart, hyppos, hypvel);

      // Add back in the Solar coordinates at mjdstart
      hyppos.x += Sunpos[planetfile_startpoint].x;
      hyppos.y += Sunpos[planetfile_startpoint].y;
      hyppos.z += Sunpos[planetfile_startpoint].z;
      hypvel.x += Sunvel[planetfile_startpoint].x;
      hypvel.y += Sunvel[planetfile_startpoint].y;
      hypvel.z += Sunvel[planetfile_startpoint].z;

      // stringID,MJD,ISO-barycenter(J2000x)(km),ISO-barycenter(J2000y)(km),ISO-barycenter(J2000z)(km),ISO-barycenter(J2000vx)(km/s),ISO-barycenter(J2000vy)(km/s),ISO-barycenter(J2000vz)(km/s),absmag,uvel,vvel,wvel,vinf,impactpar,a,e,encounter_dist,mjd_perihelion\n";

      outstream1 << "iso" << idnumstring << goodsimct << "," << mjdstart << "," << hyppos.x << "," << hyppos.y << "," << hyppos.z << "," << hypvel.x << "," << hypvel.y << "," << hypvel.z << "," << absmag << "," << uvel << "," << vvel << "," << wvel << "," << vinf << "," << impactpar << ",";
      outstream1 << a/AU_KM << "," << e << "," << encounter_dist/AU_KM << "," << mjd_perihelion << "\n";
      
      cout << "Hyppos: " << hyppos.x/AU_KM << " " << hyppos.y/AU_KM << " " << hyppos.z/AU_KM << "\n";
      
      // Close if-statement checking if the object ever comes close
      // enough to the Sun to be plausibly detectable.
    } else {
      cout << "Rejected: encounter distance was " << encounter_dist/AU_KM << ": too large\n";
    }
    if(evergood==1) {
      goodsimct++;
      if(goodsimct>=compfac) {
	compfac *= 10;
	if(idnl>0) idnl--;
	idnumstring={};
	for(i=0; i<idnl; i++) idnumstring.push_back('0');
      }
    }
    // Close loop on count of simulated objects
  }
  outstream1.close();
  
  return(0);
}
