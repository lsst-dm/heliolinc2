/// July 12, 2023: interstellar_sim03a.cpp
//
// October 16, 2023: changed to not output an ephemeris file unless
// requested, and to enable different prefixes for object string IDs.
//
// Creates simulated orbits for of interstellar objects. Models the
// population of interstellar objects as a power law in absolute magnitude,
// with v_infinity relative to the sun distributed according to a Gaussian
// ellipsoid in terms of the Galactic U, V, and W velocities. The central
// coordinate and the widths (Gaussian sigmas) of this ellipsoid in the
// U, V, and W velocity coordinates are set by the user. Given a U, V, W
// velocity selected randomly from this distribution, the program also
// randomly selects an impact parameter (b_infinity) distributed according
// to the area of concentric rings P(b) ~ b db. The actual encounter
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
#define DEBUG 0
#define PHASE_G 0.15
#define IDNUMLEN 7

#define TIME_BEFORE_PERIHELION 1000.0 // State vectors will be calculated for this far before perihelion, in days.
#define PRINTEVERY 10

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
  cerr << "Usage: interstellar_sim03a -cfg configfile -ranseed random_number_seed -mjdstart mjdstart -mjdend mjdend -prefix prefix -simnum simnum -outfile1 state vector file -outfile2 ephemeris file (optional) \n";
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
  string outstate,outephem;
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
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> targpos;
  vector <point3LD> targvel;
  point3LD targ_to_sun = point3LD(0,0,0);
  point3LD targvel_to_sunvel = point3LD(0,0,0);
  point3LD targ_to_obs = point3LD(0,0,0);
  point3LD obs_to_sun = point3LD(0,0,0);
  long double cosphase=0L;
  long double phaseang=0L;
  long double sunelong=0L;
  long double phi1,phi2;
  long double obsdist=1.0L;
  long double sundist=1.0L;
  long double obsfromsun=1.0L;
  long double phaseslope = PHASE_G;
  long double light_travel_time;
  double obsmag;
  vector <point3LD> observerpos;
  vector <long double> targMJD;
  string lnfromfile;
  vector <string> linestringvec;
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  point3LD hyppos = point3LD(0,0,0);
  point3LD hypvel = point3LD(0,0,0);
  point3LD barypos = point3LD(0,0,0);
  point3LD baryvel = point3LD(0,0,0);
  point3LD unitbary = point3LD(0,0,0);
  point3LD lvec = point3LD(0L,0L,0L);
  double newRA,newDec,impactpar;
  long double outRA,outDec;
  double racenter,deccenter,dist,pa,posRA,posDec;
  long double ldRA,ldDec,vesc;
  long double E,lscalar,a,e,coshH,H0,radvel,omega,t0,mjd_infall,mjd_epoch;
  ofstream outstream1;
  ofstream outstream2;
  string seedstring;
  int goodsimct=0;
  int idnl = IDNUMLEN;
  int compfac = 1;
  string idnumstring;
  string prefix = "iso";
  
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
	    for(j=0;j<long(planetmjd.size());j++) {
	      if(mjdtest[j]!=planetmjd[j]) {
		cout << "ERROR: time vectors do not match for input planet files\n";
		cout << planetct+1 << " and 1!\n";
		return(1);
	      }
	    }
	  }
	  for(j=0;j<long(temppos.size());j++) {
	    planetpos.push_back(temppos[j]);
	  }
	  if(planetct == pctEarth) {
	    Earthpos = temppos;
	    Earthvel = temppos;
	  } else if(planetct == pctSun) {
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
    } else if(string(argv[i]) == "-prefix" || string(argv[i]) == "-pf") {
      if(i+1 < argc) {
	//There is still something to read;
	prefix=argv[++i];
	i++;
      }
      else {
	cerr << "Prefix keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out1" || string(argv[i]) == "-outfile1") {
      if(i+1 < argc) {
	//There is still something to read;
	outstate=argv[++i];
	i++;
      }
      else {
	cerr << "Output state vector file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out2" || string(argv[i]) == "-outfile2") {
      if(i+1 < argc) {
	//There is still something to read;
	outephem=argv[++i];
	i++;
      }
      else {
	cerr << "Output state vector file keyword supplied with no corresponding argument\n";
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
  cout << "prefix for object identifiers: " << prefix << "\n";
  cout << "output state vector file " << outstate;
  if(outephem.size()>0) cout << ", and ephemeris file " << outephem << "\n";
  else cout << ", and no output ephemeris file\n";

  // Initialize random number generator
  seed_seq seed (seedstring.begin(),seedstring.end());
  mt19937_64 generator (seed);   // mt19937 is a standard mersenne_twister_engine

  // Open output file and write header
  if(outephem.size()>0) outstream1.open(outephem);
  outstream2.open(outstate);
  if(outephem.size()>0) outstream1 << "stringID absmag uvel vvel wvel vinf impactpar a e encounter_dist mjd_perihelion sundist obsdist sunelong phaseang MJD outRA outDec obsmag Ast-Sun(J2000x)(km) Ast-Sun(J2000y)(km) Ast-Sun(J2000z)(km) Ast-Sun(J2000vx)(km/s) Ast-Sun(J2000vy)(km/s) Ast-Sun(J2000vz)(km/s) Ast-solarbary_x(km) Ast-solarbary_y(km) Ast-solarbary_z(km) Ast-solarbary_vx(km/s) Ast-solarbary_vy(km/s) Ast-solarbary_vz(km/s) radiantRA radiantDec\n";
  outstream2 << "stringID absmag uvel vvel wvel vinf impactpar a e encounter_dist mjd_perihelion mjd_epoch sundist Ast-Sun(J2000x)(km) Ast-Sun(J2000y)(km) Ast-Sun(J2000z)(km) Ast-Sun(J2000vx)(km/s) Ast-Sun(J2000vy)(km/s) Ast-Sun(J2000vz)(km/s) Ast-solarbary_x(km) Ast-solarbary_y(km) Ast-solarbary_z(km) Ast-solarbary_vx(km/s) Ast-solarbary_vy(km/s) Ast-solarbary_vz(km/s) radiantRA radiantDec\n";

  // Keep track of leading zeros for string IDs
  idnl = IDNUMLEN-1;
  compfac = 10;
  idnumstring={};
  for(i=0; i<idnl; i++) idnumstring.push_back('0');

  // Loop on simulated objects.
  for(simct=0; simct<simnum; simct++) {
    // 1. Randomly assign an absolute magnitude H
    acoef = 1.0l/(Hslope*log(10.0));
    xmin = exp(Hmin/acoef);
    xmax = exp(Hmax/acoef);
    x = xmin + (xmax-xmin)*unitvar(generator);
    absmag = acoef*log(x);
    if(outephem.size()>0) outstream1.precision(10);  
    outstream2.precision(10);  
    
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
      mjd_infall = mjd_perihelion + t0/86400.0L; // Note that t0 is expected to be negative
      // Primary epoch should be a long time before perihelion, so we start n-body
      // integration before most of the planetary encounters
      mjd_epoch = mjd_perihelion - TIME_BEFORE_PERIHELION;
      // Revise mjd_epoch so that it lands exactly on a timestep from the planet files
      double min_err = LARGERR2;
      planetfile_startpoint = -99;
      for(j=0;j<long(planetmjd.size());j++) {
	if(fabs(planetmjd[j]-mjd_epoch) < min_err) {
	  min_err = fabs(planetmjd[j]-mjd_epoch);
	  planetfile_startpoint = j;
	}
      }
      planetfile_endpoint = planetmjd.size() - polyorder-2;
      if(planetfile_startpoint<0 || planetfile_startpoint>=long(planetmjd.size())) {
	cerr << "ERROR: mjd_epoch " << mjd_epoch << " could not be matched\nto any timestep in the planet files\n";
	return(1);
      } else {
	cout << "mjd_epoch " << mjd_epoch << " corresponds to timestep " << planetfile_startpoint << " in the planet files\n";
      }
      mjd_epoch = planetmjd[planetfile_startpoint];
  
      // Integrate orbit from the MJD corresponding to the start of infall,
      // down to the selected mjd_epoch, the requested start time for the n-body integration.
      Hyper_Kepint(GMsun, mjd_infall, startpos, startvel, mjd_epoch, hyppos, hypvel);
      sundist = sqrt(hyppos.x*hyppos.x + hyppos.y*hyppos.y + hyppos.z*hyppos.z);
      // Add back in the Solar coordinates at mjd_epoch
      barypos.x = hyppos.x + Sunpos[planetfile_startpoint].x;
      barypos.y = hyppos.y + Sunpos[planetfile_startpoint].y;
      barypos.z = hyppos.z + Sunpos[planetfile_startpoint].z;
      baryvel.x = hypvel.x + Sunvel[planetfile_startpoint].x;
      baryvel.y = hypvel.y + Sunvel[planetfile_startpoint].y;
      baryvel.z = hypvel.z + Sunvel[planetfile_startpoint].z;
 
      cout << "Hyppos: " << hyppos.x/AU_KM << " " << hyppos.y/AU_KM << " " << hyppos.z/AU_KM << " " << sundist/AU_KM << "\n";
      outstream2 << prefix << idnumstring << goodsimct << " " << absmag << " " << uvel << " " << vvel << " " << wvel << " " << vinf << " " << impactpar << " ";
      outstream2 << a/AU_KM << " " << e << " " << encounter_dist/AU_KM << " " << mjd_perihelion << " ";
      outstream2 << mjd_epoch << " " << sundist/AU_KM << " ";
      outstream2 << hyppos.x << " " << hyppos.y << " " << hyppos.z << " ";
      outstream2 << hypvel.x << " " << hypvel.y << " " << hypvel.z << " ";
      outstream2 << barypos.x << " " << barypos.y << " " << barypos.z << " ";
      outstream2 << baryvel.x << " " << baryvel.y << " " << baryvel.z << " ";
      outstream2 << racenter << " " << deccenter << "\n";
	  
      // 15. Perform n-body integration beginning with the state vectors
      //     produced by the hyperbolic Keplerian integration at mjdstart.
      integrate_orbit04LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, barypos, baryvel, planetfile_startpoint, planetfile_endpoint, targMJD, targpos, targvel);
      // Print out ephemeris every 10 days
      for(long ptct=planetfile_startpoint; ptct<planetfile_endpoint; ptct++) {
	if(ptct%PRINTEVERY==0) {
	  targMJD[ptct-planetfile_startpoint];
 
	  targ_to_obs.x = Earthpos[ptct].x - targpos[ptct-planetfile_startpoint].x;
	  targ_to_obs.y = Earthpos[ptct].y - targpos[ptct-planetfile_startpoint].y;
	  targ_to_obs.z = Earthpos[ptct].z - targpos[ptct-planetfile_startpoint].z;
	  targ_to_sun.x = Sunpos[ptct].x - targpos[ptct-planetfile_startpoint].x;
	  targ_to_sun.y = Sunpos[ptct].y - targpos[ptct-planetfile_startpoint].y;
	  targ_to_sun.z = Sunpos[ptct].z - targpos[ptct-planetfile_startpoint].z;
	  targvel_to_sunvel.x = Sunvel[ptct].x - targvel[ptct-planetfile_startpoint].x;
	  targvel_to_sunvel.y = Sunvel[ptct].y - targvel[ptct-planetfile_startpoint].y;
	  targvel_to_sunvel.z = Sunvel[ptct].z - targvel[ptct-planetfile_startpoint].z;

	  obsdist = sqrt(targ_to_obs.x*targ_to_obs.x + targ_to_obs.y*targ_to_obs.y + targ_to_obs.z*targ_to_obs.z);
	  sundist = sqrt(targ_to_sun.x*targ_to_sun.x + targ_to_sun.y*targ_to_sun.y + targ_to_sun.z*targ_to_sun.z);
	  // Calculate observed magnitude not accounting for phase
	  obsmag = absmag + 2.5*log10(obsdist*obsdist*sundist*sundist/AU_KM/AU_KM/AU_KM/AU_KM);

	  // Object might be visible, calculate phase angle
	  cosphase = dotprod3LD(targ_to_obs,targ_to_sun)/obsdist/sundist;
	  if(cosphase>1.0L) {
	    cout << "WARNING: trying to take arccos of 1.0 + " << cosphase-1.0L << "\n";
	    phaseang=0.0L;
	  } else if(cosphase<-1.0L) {
	    cout << "WARNING: trying to take arccos of -1.0 - " << cosphase+1.0L << "\n";
	    phaseang=M_PI;
	  } else phaseang = acos(cosphase);
	  // Correct the magnitude calculated above for phase effects
	  phi1 = exp(-PHASECONST_A1*pow(tan(phaseang/2.0),PHASECONST_B1));
	  phi2 = exp(-PHASECONST_A2*pow(tan(phaseang/2.0),PHASECONST_B2));
	  obsmag -= 2.5*log10((1.0-phaseslope)*phi1 + phaseslope*phi2);
	  // Also calculate solar elongation, while we have the values handy
	  obs_to_sun.x = Sunpos[ptct].x - Earthpos[ptct].x;
	  obs_to_sun.y = Sunpos[ptct].y - Earthpos[ptct].y;
	  obs_to_sun.z = Sunpos[ptct].z - Earthpos[ptct].z;
	  targ_to_obs.x = targpos[ptct-planetfile_startpoint].x - Earthpos[ptct].x; // Now it is really obs_to_targ
	  targ_to_obs.y = targpos[ptct-planetfile_startpoint].y - Earthpos[ptct].y;
	  targ_to_obs.z = targpos[ptct-planetfile_startpoint].z - Earthpos[ptct].z;
	  obsfromsun = sqrt(obs_to_sun.x*obs_to_sun.x + obs_to_sun.y*obs_to_sun.y + obs_to_sun.z*obs_to_sun.z);
	  cosphase = dotprod3LD(targ_to_obs,obs_to_sun)/obsfromsun/obsdist;
	  if(cosphase>1.0L) {
	    cout << "WARNING: trying to take arccos of 1.0 + " << cosphase-1.0L << "\n";
	    sunelong=0.0L;
	  } else if(cosphase<-1.0L) {
	    cout << "WARNING: trying to take arccos of -1.0 - " << cosphase+1.0L << "\n";
	    sunelong=M_PI;
	  } else sunelong = acos(cosphase);
	  // Calculate the object's precise position, accounting for light travel time
	  light_travel_time = obsdist*1000.0/CLIGHT; // Factor of 1000 converts obsdist to meters
	  // Light-travel-time corrected version of coordinates relative to the observer
	  targ_to_obs.x -= light_travel_time*targvel[ptct-planetfile_startpoint].x;
	  targ_to_obs.y -= light_travel_time*targvel[ptct-planetfile_startpoint].y;
	  targ_to_obs.z -= light_travel_time*targvel[ptct-planetfile_startpoint].z;
	  // Light-travel-time corrected observer-target distance
	  obsdist = sqrt(targ_to_obs.x*targ_to_obs.x + targ_to_obs.y*targ_to_obs.y + targ_to_obs.z*targ_to_obs.z);
	  // Calculate unit vector
	  targ_to_obs.x /= obsdist;
	  targ_to_obs.y /= obsdist;
	  targ_to_obs.z /= obsdist;
	  // Project onto the celestial sphere.
	  stateunitLD_to_celestial(targ_to_obs, outRA, outDec);
	  // Write a bunch of data to the output file.
	  if(outephem.size()>0) {
	    outstream1 << "iso" << idnumstring << goodsimct << " " << absmag << " " << uvel << " " << vvel << " " << wvel << " " << vinf << " " << impactpar << " ";
	    outstream1 << a/AU_KM << " " << e << " " << encounter_dist/AU_KM << " " << mjd_perihelion << " ";
	    outstream1 << sundist/AU_KM << " " << obsdist/AU_KM << " " << sunelong*DEGPRAD << " " << phaseang*DEGPRAD << " ";
	    outstream1 << planetmjd[ptct] << " " << outRA << " " << outDec << " " << obsmag << " ";
	    // Note that we reverse the sign here because the previously calculated
	    // quantities are the negatives of the true state vectors: that is, their
	    // vector origin is at the object, pointing toward the sun. The negatives
	    // we apply here solve this problem and produce the correct sun-to-object state vectors.
	    outstream1 << -targ_to_sun.x << " " << -targ_to_sun.y << " " << -targ_to_sun.z << " ";
	    outstream1 << -targvel_to_sunvel.x << " " << -targvel_to_sunvel.y << " " << -targvel_to_sunvel.z << " ";
	    outstream1 << targpos[ptct-planetfile_startpoint].x << " " << targpos[ptct-planetfile_startpoint].y << " " << targpos[ptct-planetfile_startpoint].z << " ";
	    outstream1 << targvel[ptct-planetfile_startpoint].x << " " << targvel[ptct-planetfile_startpoint].y << " " << targvel[ptct-planetfile_startpoint].z << " ";
	    outstream1 << racenter << " " << deccenter << "\n";
	  }
	  // Close conditional for printing
	}
	// Close loop over all valid times
      }
      // Close if-statement checking if the object ever comes close
      // enough to the Sun to be plausibly detectable.
      goodsimct++;
      if(goodsimct>=compfac) {
	compfac *= 10;
	if(idnl>0) idnl--;
	idnumstring={};
	for(i=0; i<idnl; i++) idnumstring.push_back('0');
      }
    } else {
      cout << "Rejected: encounter distance was " << encounter_dist/AU_KM << ": too large\n";
    }
    // Close loop on count of simulated objects
  }
  if(outephem.size()>0) outstream1.close();
  
  return(0);
}
