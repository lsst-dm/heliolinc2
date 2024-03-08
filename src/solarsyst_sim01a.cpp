/// June 08, 2022: solarsyst_sim01a.cpp
// Given an input distribution in terms of Keplerian orbital parameters,
// simulate a population of solar system objects, creating an output
// catalog of simulated detections that is suitable for testing heliolinc.

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
  cerr << "Usage: solarsyst_sim01a -cfg configfile -ranseed random_number_seed -images imfile -imrad image radius(deg) -mjdstart mjdstart -mjdend mjdend -astromerr 1-D astrometric error (arcsec) -simnum simnum -outfile outfile \n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  int status=0;
  ifstream instream1;
  string stest;
  string configfile;
  string imfile;
  string planetfile;
  string outfile;
  long double minperidist,maxperidist,minapdist,maxapdist;
  long double incmean,incsigma,maxdetdist;
  long double x,y,z,vx,vy,vz;
  x=y=z=vx=vy=vz=0.0L;
  long double mjdstart = 0.0L;
  long double mjdend = 0.0L;
  double imrad = IMAGERAD; // radius from image center to most distant corner (deg).
  int planetfile_startpoint=0;
  int planetfile_endpoint=0;
  long double GMsun = 0.0L;
  double Hmin = 0.0l;
  double Hmax = 0.0l;
  double Hslope = 0.0l;
  double limiting_mag = 0.0l;
  double acoef,xmin,xmax,absmag;
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
  vector <point3LD> targpos;
  vector <point3LD> targvel;
  point3LD targ_to_sun = point3LD(0,0,0);
  point3LD targvel_to_sunvel = point3LD(0,0,0);
  point3LD targ_to_obs = point3LD(0,0,0);
  point3LD obs_to_sun = point3LD(0,0,0);
  point3LD Sunposnow = point3LD(0,0,0);
  point3LD Sunvelnow = point3LD(0,0,0);
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
  vector <double> imageRA;
  vector <double> imageDec;
  vector <long double> imageMJD;
  point3LD imagepos = point3LD(0,0,0);
  point3LD imagevel = point3LD(0,0,0);
  int imnum=0;
  int imct=0;
  string lnfromfile;
  vector <string> linestringvec;
  int badread=0;
  int reachedeof=0;
  //  double obslon = 289.26345L;
  //  double plxcos = 0.865020L;
  //  double plxsin = -0.500901L;
  double obslon = 289.25058l; // Changed to match X05 obscode
  double plxcos = 0.864981l;
  double plxsin = -0.500958l;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  point3LD startpos = point3LD(0,0,0);
  point3LD startvel = point3LD(0,0,0);
  point3LD Keppos = point3LD(0,0,0);
  point3LD Kepvel = point3LD(0,0,0);
  point3LD unitbary = point3LD(0,0,0);
  point3LD lvec = point3LD(0L,0L,0L);
  double newRA,newDec;
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
  keplerian_orbit keporb=keplerian_orbit(0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L);
  long double peridist,apdist;
  
  // keporb.semimaj_axis        in AU
  // keporb.eccentricity        unitless
  // keporb.inclination         in degrees
  // keporb.long_ascend_node    Longitude of the ascending node, in degrees
  // keporb.arg_perihelion      Argument of perihelion, in degrees
  // keporb.mean_anom           Mean anomaly at the epoch, in degrees
  // keporb.mjd_epoch           Epoch for the orbit in MJD
  // keporb.mean_daily_motion   in degrees/day

  if(argc<15) {
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
	// Read minimum perihelion distance
	status=readconfigLD(instream1,&minperidist);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&minperidist);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Minimum perihelion distance read as " << minperidist << " AU\n";
	// Read maximum perihelion distance
	status=readconfigLD(instream1,&maxperidist);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&maxperidist);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Maximum perihelion distance read as " << maxperidist << " AU\n";
	// Read minimum aphelion distance
	status=readconfigLD(instream1,&minapdist);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&minapdist);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Minimum aphelion distance read as " << minapdist << " AU\n";
	// Read maximum aphelion distance
	status=readconfigLD(instream1,&maxapdist);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&maxapdist);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Maximum aphelion distance read as " << maxapdist << " AU\n";
	// Read mean value for the inclination
	status=readconfigLD(instream1,&incmean);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&incmean);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Mean value for orbital inclination read as " << incmean << " degrees\n";
	// Read Gaussian sigma for the inclination
	status=readconfigLD(instream1,&incsigma);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&incsigma);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Gaussian sigma for the inclination read as " << incsigma << " degrees\n";
	// Read the maximum detection distance in AU
	status=readconfigLD(instream1,&maxdetdist);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&maxdetdist);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Maximum detection distance read as " << maxdetdist << " AU\n";
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
	// Read astrometric error in arcseconds
	status=readconfigLD(instream1,&astromsigma);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&astromsigma);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Astrometric error read as " << astromsigma << " arcsec\n";
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
    } else if(string(argv[i]) == "-i" || string(argv[i]) == "-im" || string(argv[i]) == "-imf" || string(argv[i]) == "-images" || string(argv[i]) == "-imfile" || string(argv[i]) == "--imagefile" || string(argv[i]) == "--imfile" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
	// There is still something to read;
	imfile=argv[++i];
	i++;
      } else {
	cerr << "Image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-imrad") {
      if(i+1 < argc) {
	//There is still something to read;
        imrad=stod(argv[++i]);
	i++;
	if(!isnormal(imrad) || imrad<=0.0) {
	  cerr << "Error: invalid image radius (" << imrad << " deg) supplied.\n";
	  cerr << "Image radius must be strictly positive!\n";
	  return(2);
	}
      }
      else {
	cerr << "Output image radius keyword supplied with no corresponding argument\n";
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
    }  else if(string(argv[i]) == "-astromerr" || string(argv[i]) == "-astromsig" || string(argv[i]) == "-as"  || string(argv[i]) == "-sigast"  || string(argv[i]) == "-sigastrom" || string(argv[i]) == "--astromerror" || string(argv[i]) == "--astromsigma") {
      if(i+1 < argc) {
	//There is still something to read;
	astromsigma=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Astrometric error keyword supplied with too few corresponding arguments\n";
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
  cout << "input image file " << imfile << "\n";
  cout << "input starting MJD " << mjdstart << "\n";
  cout << "input ending MJD " << mjdend << "\n";
  cout << "Perihelion distance must be between " << minperidist << " and " << maxperidist << "\n";
  cout << "Aphelion distance must be between " << minapdist << " and " << maxapdist << "\n";
  cout << "Inclination follows a Gaussian distribution with mean = " << incmean << " and sigma = " << incsigma << "\n";
  cout << "min and max H magnitudes: " << Hmin << " " << Hmax << "\n";
  cout << "power law slope for H magnitude: " << Hslope << "\n";
  cout << "1-D Gaussian error added to output astrometry: " << astromsigma << " arcsec\n";
  cout << "number of encounters to simulate: " << simnum << "\n";
  cout << "output file " << outfile << "\n";

  if(minperidist>maxperidist) {
    cerr << "ERROR: minimum perihelion distance is greater than maximum\n";
    return(1);
  }
  if(minapdist>maxapdist) {
    cerr << "ERROR: minimum aphelion distance is greater than maximum\n";
    return(1);
  }
  if(maxperidist>maxapdist) {
    cerr << "ERROR: maximum perihelion distance is greater than\n";
    cerr << "maximum aphelion distances, meaning that there are\n";
    cerr << "allowed values for the perihelion distance for which\n";
    cerr << "no corresponding allowed values for the aphelion distance exist\n";
    return(1);
  }

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

  // Read input image file.
  instream1.open(imfile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << imfile << "\n";
    return(1);
  }
  badread=0;
  while(reachedeof==0 && !instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(lnfromfile.size()>10 && (reachedeof==0 || reachedeof==1)) {
      fieldnum = get_col_vector01(lnfromfile, linestringvec);
      if(fieldnum>=MJDcol) {
	ldval = stold(linestringvec[MJDcol-1]);
	imageMJD.push_back(ldval);
      } else badread=1;
      if(fieldnum>=RAcol) {
	dval = stod(linestringvec[RAcol-1]);
	imageRA.push_back(dval);
      } else badread=1;
      if(fieldnum>=Deccol) {
	dval = stod(linestringvec[Deccol-1]);
	imageDec.push_back(dval);
      } else badread=1;
    }
    cout << imageMJD.size() << " " << imageMJD[imageMJD.size()-1] << " " << imageRA[imageMJD.size()-1] << " " << imageDec[imageMJD.size()-1] << "\n";
  }
  imnum=imageMJD.size();
  if(imnum<2 || imageRA.size()!=imnum || imageDec.size()!=imnum) {
    cerr  << fixed << setprecision(6) << "Error: image time vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << imnum << " " << imageRA.size() << " " << imageDec.size() << "\n";
    return(1);
  }
  if(badread>=1){
    cerr << "Error reading " << imfile << ": apparent short line\n";
    return(1);
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  instream1.close();
  
  // Calculate the exact position of the observer at the time of each image.
  observerpos={};
  for(imct=0;imct<imnum;imct++) {
    observer_barycoords01LD(imageMJD[imct], 5, obslon, plxcos, plxsin, planetmjd, Earthpos, outpos);
    observerpos.push_back(outpos);
  }

  // Open output file and write header
  outstream1.open(outfile);
  outstream1 << "stringID absmag semimaj_axis eccentricity inclination long_ascend_node arg_perihelion mean_anom mjd_epoch mean_daily_motion sundist obsdist sunelong phaseang imageMJD outRA outDec obsmag Ast-Sun(J2000x)(km) Ast-Sun(J2000y)(km) Ast-Sun(J2000z)(km) Ast-Sun(J2000vx)(km/s) Ast-Sun(J2000vy)(km/s) Ast-Sun(J2000vz)(km/s) perfectRA perfectDec\n";
  // keporb.semimaj_axis        in AU
  // keporb.eccentricity        unitless
  // keporb.inclination         in degrees
  // keporb.long_ascend_node    Longitude of the ascending node, in degrees
  // keporb.arg_perihelion      Argument of perihelion, in degrees
  // keporb.mean_anom           Mean anomaly at the epoch, in degrees
  // keporb.mjd_epoch           Epoch for the orbit in MJD
  // keporb.mean_daily_motion   in degrees/day

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
    
    // 3. Assign the orbital inclination from a Gaussian distribution
    keporb.inclination = incmean + incsigma*gaussian_deviate_mt(generator);

    // 4. Assign the perihelion and aphelion distance from uniform distributions.
    peridist = minperidist + (maxperidist - minperidist)*unitvar(generator);
    if(peridist<minapdist) {
      apdist = minapdist + (maxapdist - minapdist)*unitvar(generator);
    } else {
      apdist = peridist + (maxapdist - peridist)*unitvar(generator);
    }

    // 5. Solve for semimajor axis and eccentricity given these
    //    perihelion and aphelion distances.
    keporb.semimaj_axis = 0.5L*peridist + 0.5L*apdist;
    keporb.eccentricity = (apdist-peridist)/(apdist+peridist);

    // 6. Assign the three dynamically unimportant Keplerian
    //    orbital elements from a uniform distribution on [0,360)
    keporb.long_ascend_node = 360.0L*unitvar(generator);
    keporb.arg_perihelion = 360.0L*unitvar(generator);
    keporb.mean_anom = 360.0L*unitvar(generator);
    if(keporb.long_ascend_node >= 360.0L) keporb.long_ascend_node -= 360.0L;
    if(keporb.arg_perihelion >= 360.0L) keporb.arg_perihelion -= 360.0L;
    if(keporb.mean_anom >= 360.0L) keporb.mean_anom -= 360.0L;
    keporb.mean_daily_motion = sqrt(LDSQUARE(KCONST)/intpowLD(keporb.semimaj_axis,3.0))*DEGPRAD;
    keporb.mjd_epoch = mjdstart;
 
    // 7. Convert Keplerian orbital elements to dynamical state vectors.
    Kepler2dyn(mjdstart, keporb, outpos, outvel);
   // Convert from AU and AU/day to km and km/sec
    outpos.x*=AU_KM;
    outpos.y*=AU_KM;
    outpos.z*=AU_KM;
    outvel.x*=AU_KM/SOLARDAY;
    outvel.y*=AU_KM/SOLARDAY;
    outvel.z*=AU_KM/SOLARDAY;
    cout << "simulated object " << simct << ", state vectors in km and km/sec:\n";
    cout << outpos.x << " " << outpos.y << " " << outpos.z << "\n";
    cout << outvel.x << " " << outvel.y << " " << outvel.z << "\n";

    // Add back in the Solar coordinates at mjdstart
    planetposvel01LD(mjdend,polyorder,planetmjd,Sunpos,Sunvel,Sunposnow,Sunvelnow);
    outpos.x += Sunpos[planetfile_startpoint].x;
    outpos.y += Sunpos[planetfile_startpoint].y;
    outpos.z += Sunpos[planetfile_startpoint].z;
    outvel.x += Sunvel[planetfile_startpoint].x;
    outvel.y += Sunvel[planetfile_startpoint].y;
    outvel.z += Sunvel[planetfile_startpoint].z;
      
    // 8. Perform n-body integration beginning with the state vectors
    // produced by Kepler2dyn
    integrate_orbit04LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, outpos, outvel, planetfile_startpoint, planetfile_endpoint, targMJD, targpos, targvel);
    for(imct=0;imct<imnum;imct++)
      {
	planetposvel01LD(imageMJD[imct],polyorder,targMJD,targpos,targvel,imagepos,imagevel);
	planetposvel01LD(imageMJD[imct],polyorder,planetmjd,Sunpos,Sunvel,Sunposnow,Sunvelnow);
 
	targ_to_obs.x = observerpos[imct].x - imagepos.x;
	targ_to_obs.y = observerpos[imct].y - imagepos.y;
	targ_to_obs.z = observerpos[imct].z - imagepos.z;
	targ_to_sun.x = Sunposnow.x - imagepos.x;
	targ_to_sun.y = Sunposnow.y - imagepos.y;
	targ_to_sun.z = Sunposnow.z - imagepos.z;
	targvel_to_sunvel.x = Sunvelnow.x - imagevel.x;
	targvel_to_sunvel.y = Sunvelnow.y - imagevel.y;
	targvel_to_sunvel.z = Sunvelnow.z - imagevel.z;

	obsdist = sqrt(targ_to_obs.x*targ_to_obs.x + targ_to_obs.y*targ_to_obs.y + targ_to_obs.z*targ_to_obs.z);
	sundist = sqrt(targ_to_sun.x*targ_to_sun.x + targ_to_sun.y*targ_to_sun.y + targ_to_sun.z*targ_to_sun.z);
	// Calculate observed magnitude not accounting for phase
	obsmag = absmag + 2.5*log10(obsdist*obsdist*sundist*sundist/AU_KM/AU_KM/AU_KM/AU_KM);
	if(obsmag<limiting_mag  && (sundist/AU_KM)<=maxdetdist) {
	  //cout << simct << ": obsmag = " << obsmag << ", sundist = " << sundist/AU_KM << ", obsdist = " << obsdist/AU_KM << "\n";
	  // Note: maxdetdisk is designed to enable the artificial
	  // imposition of a maximum heliocentric distance, beyond which
	  // objects are simply not considered. The original planned use
	  // case is to examine the detectability of NEOs less than 1AU
	  // from the sun, regardless of whether their orbits are confined
	  // to heliocentric distance less than 1AU or not.
	  
	  // Object might be visible, calculate phase angle.
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
	  obs_to_sun.x = Sunposnow.x - observerpos[imct].x;
	  obs_to_sun.y = Sunposnow.y - observerpos[imct].y;
	  obs_to_sun.z = Sunposnow.z - observerpos[imct].z;
	  targ_to_obs.x = imagepos.x - observerpos[imct].x; // Now it is really obs_to_targ
	  targ_to_obs.y = imagepos.y - observerpos[imct].y;
	  targ_to_obs.z = imagepos.z - observerpos[imct].z;
	  obsfromsun = sqrt(obs_to_sun.x*obs_to_sun.x + obs_to_sun.y*obs_to_sun.y + obs_to_sun.z*obs_to_sun.z);
	  cosphase = dotprod3LD(targ_to_obs,obs_to_sun)/obsfromsun/obsdist;
	  if(cosphase>1.0L) {
	    cout << "WARNING: trying to take arccos of 1.0 + " << cosphase-1.0L << "\n";
	    sunelong=0.0L;
	  } else if(cosphase<-1.0L) {
	    cout << "WARNING: trying to take arccos of -1.0 - " << cosphase+1.0L << "\n";
	    sunelong=M_PI;
	  } else sunelong = acos(cosphase);
	    
	  if(obsmag<limiting_mag) {
	    // Object might still be visible, calculate its position
	    // Instantaneous velocity
	    planetpos01LD(imageMJD[imct],polyorder,targMJD,targvel,imagevel);
	    // Initial approximation of the coordinates relative to the observer
	    imagepos.x -= observerpos[imct].x;
	    imagepos.y -= observerpos[imct].y;
	    imagepos.z -= observerpos[imct].z;
	    light_travel_time = obsdist*1000.0/CLIGHT; // Factor of 1000 converts obsdist to meters
	    // Light-travel-time corrected version of coordinates relative to the observer
	    imagepos.x -= light_travel_time*imagevel.x;
	    imagepos.y -= light_travel_time*imagevel.y;
	    imagepos.z -= light_travel_time*imagevel.z;
	    // Light-travel-time corrected observer-target distance
	    obsdist = sqrt(imagepos.x*imagepos.x + imagepos.y*imagepos.y + imagepos.z*imagepos.z);
	    // Calculate unit vector
	    imagepos.x /= obsdist;
	    imagepos.y /= obsdist;
	    imagepos.z /= obsdist;
	    // Project onto the celestial sphere.
	    stateunitLD_to_celestial(imagepos, outRA, outDec);
	    dist = distradec01(imageRA[imct],imageDec[imct],outRA,outDec);
	    // See if it's on the image.
	    if(dist<imrad) {
	      evergood=1;
	      // It is on the image. Write a bunch of data to the output file.
	      outstream1 << "ieo" << idnumstring << goodsimct << " " << absmag << " ";
	      outstream1 << keporb.semimaj_axis << " " << keporb.eccentricity << " " << keporb.inclination << " ";
	      outstream1 << keporb.long_ascend_node << " " << keporb.arg_perihelion << " " << keporb.mean_anom << " ";
	      outstream1 << keporb.mjd_epoch << " " << keporb.mean_daily_motion << " ";
	      outstream1 << sundist/AU_KM << " " << obsdist/AU_KM << " " << sunelong*DEGPRAD << " " << phaseang*DEGPRAD << " ";
	      outstream1 << imageMJD[imct] << " " << outRA + astromsigma*gaussian_deviate_mt(generator)/3600.0L/cos(outDec/DEGPRAD) << " " << outDec + astromsigma*gaussian_deviate_mt(generator)/3600.0L << " " << obsmag << " ";
	      // Note that we reverse the sign here because the previously calculated
	      // quantities are the negatives of the true state vectors: that is, their
	      // vector origin is at the object, pointing toward the sun. The negatives
	      // we apply here solve this problem and produce the correct sun-to-object state vectors.
	      outstream1 << -targ_to_sun.x << " " << -targ_to_sun.y << " " << -targ_to_sun.z << " ";
	      outstream1 << -targvel_to_sunvel.x << " " << -targvel_to_sunvel.y << " " << -targvel_to_sunvel.z << " ";
	      outstream1 << outRA << " " << outDec << "\n";
	    } else {
	      ; // Object would have been bright enough to detect, but wasn't on the image
	    }
	  } else {
	    ; // Object was too faint due to phase effects
	  }
	} else {
	  ; // Object was too faint due to mere distance.
	}
	// Close loop over all images
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
