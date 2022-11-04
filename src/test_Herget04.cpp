// November 03, 2022: test_Herget04
// Test method of Herget routine from J. M. A. Danby,
// Fundamentals of Celestial Mechanics.
// This version tests a self-contained function that does
// a complete downhill simplex Method of Herget fit.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MINHERGETDIST 0.0001L
#define HERGET_DOWNSCALE 0.9L
#define DISTPOWMAX 7
#define DISTPOWSCALE 2.0L

long double Hergetfit01(long double geodist1, long double geodist2, long double simplex_side, long double ftol, int point1, int point2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid, vector <long double> &orbit);

static void show_usage()
{
  cerr << "Usage: test_Herget04 -cfg configfile -observations obsfile -fitpoints point1 point2 -geodist dist1 dist2 -maxrchi max_reduced_chisq -simpstep simplex_step -ftol ftol -outfile outfile\n";
  cerr << "\n\nor, at minimum,\n";
  cerr << "test_Herget03 -cfg configfile -observations obsfile -outfile outfile\n";
}


int main(int argc, char *argv[])
{
  int i,j,status,polyorder,configread,fieldnum,badread;
  int MJDcol, RAcol, Deccol, magcol, obscodecol, sigastromcol, sigmagcol;
  int reachedeof,obsnum,obsct;
  string configfile, observatory_code_file, Earthfile, obsfile;
  string outfile, lnfromfile, stest;
  vector <string> linestringvec;
  ifstream instream1;
  ofstream outstream1;
  point3LD outpos = point3LD(0.0L,0.0L,0.0L);
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> observerpos;
  vector <long double> earthmjd;
  long double geodist1, geodist2, simplex_side, ldval, chisq;
  long double maxchisq, maxrchi = 1.0L;
  
  int point1, point2, distpow;
  long double ftol = 1.0e-5;
  vector <long double> obsMJD, obsRA, obsDec, magnitude, sigastrom, sigmag;
  vector <long double> fitRA, fitDec, resid;
  vector <string> obscodevec;
  char obscode[MINSTRINGLEN];
  double obslon = 289.26345L;
  double plxcos = 0.865020L;
  double plxsin = -0.500901L;
  observatory obs1 = observatory("I11",0l,0l,0l);
  vector <observatory> observatory_list = {};
  vector <long double> orbit;
  
  i = j = status = polyorder = configread = fieldnum = badread = reachedeof = 0;
  MJDcol = RAcol = Deccol = magcol = obscodecol = sigastromcol = sigmagcol = 0;
  point1 = point2 = -1;
  
  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(i<argc && (string(argv[i]) == "-c" || string(argv[i]) == "-cfg" || string(argv[i]) == "-config" || string(argv[i]) == "--config" || string(argv[i]) == "--configfile" || string(argv[i]) == "--configuration" || string(argv[i]) == "--ConfigFile")) {
      if(i+1 < argc) {
	//There is still something to read;
	configfile=argv[++i];
	i++;
	// Read the configuration file. This must happen here, so that
	// default values supplied in this file can be overwritten later
	// if the user desires.
	// Read configuration file.
	instream1.open(configfile);
	if(!instream1) {
	  cerr << "ERROR: can't open input config file " << configfile << "\n";
	  return(1);
	}
	// Read the name of the observatory code file
	status=readconfigstring(instream1,observatory_code_file);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigstring(instream1,observatory_code_file);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observatory code file read as " << observatory_code_file << "\n";      
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
	} else cout << "Polynomial order for interpolation read as " << polyorder << "\n";
	// Read the ephemeris file for the Earth
	status=readconfigstring(instream1,Earthfile);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigstring(instream1,Earthfile);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Ephemeris file for the Earth is named " << Earthfile << "\n";
	earthmjd={};
	Earthpos={};
	Earthvel={};
	read_horizons_fileLD(Earthfile,earthmjd,Earthpos,Earthvel);
	// Read default geocentric distance at point 1 from the ephemeris file
	status=readconfigLD(instream1,&geodist1);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&geodist1);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default geocentric distance at point 1 read as " << geodist1 << " AU\n";
	// Read default geocentric distance at point 2 from the ephemeris file
	status=readconfigLD(instream1,&geodist2);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&geodist2);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default geocentric distance at point 2 read as " << geodist2 << " AU\n";
	// Read initial side length for the simplex
	status=readconfigLD(instream1,&simplex_side);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&simplex_side);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Initial side length for simplex read as " << simplex_side << "\n";
	// Read observation file column holding the MJD;
	status=readconfigint(instream1,&MJDcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&MJDcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding MJD is " << MJDcol << " \n";
	// Read observation file column holding the RA;
	status=readconfigint(instream1,&RAcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&RAcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding RA is " << RAcol << " \n";
	// Read observation file column holding the Dec;
	status=readconfigint(instream1,&Deccol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&Deccol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding Dec is " << Deccol << " \n";
	// Read observation file column holding the magnitude;
	status=readconfigint(instream1,&magcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&magcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding magnitude is " << magcol << " \n";
	// Read observation file column holding the observatory code;
	status=readconfigint(instream1,&obscodecol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&obscodecol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding observatory code is " << obscodecol << " \n";
	// Read observation file column holding the astrometric uncertainty;
	status=readconfigint(instream1,&sigastromcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&sigastromcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding astrometric uncertainty is " << sigastromcol << " \n";
	// Read observation file column holding the photometric uncertainty;
	status=readconfigint(instream1,&sigmagcol);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&sigmagcol);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Observation file column holding photometric uncertainty is " << sigmagcol << " \n";
	// Close input stream that was reading the config file.
	instream1.close();
	configread=1;
	 cout << "Configuration file read successfully\n";
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
    }
    if(i<argc && (string(argv[i]) == "-o" || string(argv[i]) == "-obs" || string(argv[i]) == "-obsfile" || string(argv[i]) == "-observations" || string(argv[i]) == "--obsfile" || string(argv[i]) == "--observation" || string(argv[i]) == "--observationfile" || string(argv[i]) == "--observations")) {
      if(i+1 < argc) {
	// There is still something to read;
	obsfile=argv[++i];
	i++;
      } else {
	cerr << "Observation file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(i<argc && (string(argv[i]) == "-fitpoints" || string(argv[i]) == "-fp" || string(argv[i]) == "-fpoints" || string(argv[i]) == "--fitpoints")) {
      if(i+2 < argc) {
	//There is still something to read;
	point1=stoi(argv[++i]);
	point2=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Fit points keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(i<argc && (string(argv[i]) == "-geodist" || string(argv[i]) == "-gd" || string(argv[i]) == "-gdist" || string(argv[i]) == "--geodist")) {
      if(i+2 < argc) {
	//There is still something to read;
	geodist1=stold(argv[++i]);
	geodist2=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Geocentric distances keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    }  else if(i<argc && (string(argv[i]) == "-maxrchi" || string(argv[i]) == "-mrc" || string(argv[i]) == "-maxreducedchi" || string(argv[i]) == "--maxrchi")) {
      if(i+1 < argc) {
	//There is still something to read;
	maxrchi=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Geocentric distances keyword supplied with too few corresponding arguments\n";
	show_usage();
	return(1);
      }
    } else if(i<argc && (string(argv[i]) == "-simpstep" || string(argv[i]) == "-ss" || string(argv[i]) == "-simplexstep" || string(argv[i]) == "-simpside" || string(argv[i]) == "-simplexside" || string(argv[i]) == "--initialsimplexstep" || string(argv[i]) == "--startingsimplexside")) {
      if(i+1 < argc) {
	//There is still something to read;
	simplex_side=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Initial simplex side length keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }else if(i<argc && (string(argv[i]) == "-ftol")) {
      if(i+1 < argc) {
	//There is still something to read;
	ftol=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "ftol keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(i<argc && (string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outorb" || string(argv[i]) == "--outorbits")) {
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
  
  cout << "Reading input observation file " << obsfile << " " << instream1.rdstate() << "\n";
  // Read input observation file.
  instream1.open(obsfile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << obsfile << "\n";
    return(1);
  }
  badread=0;
  while(reachedeof==0 && !instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(lnfromfile.size()>20 && (reachedeof==0 || reachedeof==1)) {
      fieldnum = get_col_vector01(lnfromfile, linestringvec);
      if(fieldnum>=MJDcol) {
	ldval = stold(linestringvec[MJDcol-1]);
	obsMJD.push_back(ldval);
      } else badread=1;
      if(fieldnum>=RAcol) {
	ldval = stold(linestringvec[RAcol-1]);
	obsRA.push_back(ldval);
      } else badread=1;
      if(fieldnum>=Deccol) {
	ldval = stold(linestringvec[Deccol-1]);
	obsDec.push_back(ldval);
      } else badread=1;
      if(fieldnum>=magcol) {
	ldval = stold(linestringvec[magcol-1]);
	magnitude.push_back(ldval);
      } else {
	magnitude.push_back(-99.9);
	// This isn't a bad read because we can use dummy magnitudes if necessary
      }
      if(fieldnum>=sigastromcol) {
	ldval = stold(linestringvec[sigastromcol-1]);
	sigastrom.push_back(ldval);
      } else {
	sigastrom.push_back(1.0L);
	// This isn't a bad read because we can use a default
	// astrometric uncertainty of 1 arcsec if necessary
      }
      if(fieldnum>=sigmagcol) {
	ldval = stold(linestringvec[sigmagcol-1]);
	sigmag.push_back(ldval);
      } else {
	sigmag.push_back(-99.0L);
	// This isn't a bad read because we can use dummy magnitude uncertainties if necessary
      }
       if(fieldnum>=obscodecol) {
	 obscodevec.push_back(linestringvec[obscodecol-1]);
      } else badread=1;
    }
    cout << obsMJD.size() << " " << obsMJD[obsMJD.size()-1] << " " << obsRA[obsMJD.size()-1] << " " << obsDec[obsMJD.size()-1] << "\n";
  }
  obsnum=obsMJD.size();
  if(obsnum<2 || obsRA.size()!=obsnum || obsDec.size()!=obsnum || magnitude.size()!=obsnum || sigastrom.size()!=obsnum || sigmag.size()!=obsnum || obscodevec.size()!=obsnum) {
    cerr  << fixed << setprecision(6) << "Error: observation vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << obsnum << " " << obsRA.size() << " " << obsDec.size() << " " << magnitude.size() << " " << sigastrom.size() << " " << sigmag.size() << " " << obscodevec.size() << "\n";
    return(1);
  }
  if(badread>=1){
    cerr << "Error reading " << obsfile << ": apparent short line\n";
    return(1);
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  instream1.close();

  if(point1 == -1) point1 = 1;
  if(point2 == -1) point2 = obsMJD.size();
  if(point1<1) {
    cout << "WARNING: invalid starting point " << point1 << " out of " << obsMJD.size() << " was specified!\n";
    point1 = 1;
    cout << "It has been reset to a value of " << point1 << "\n";
  }
  if(point2 <= point1) {
    cout << "WARNING: specified endpoint " << point2 << " is before start point " << point1 << "\n"; 
    point2 = obsMJD.size();
    cout << "It has been reset to a value of " << point2 << "\n";
  }
  if(point2>obsMJD.size()) {
    cout << "WARNING: invalid endpoint " << point2 << " out of " << obsMJD.size() << " was specified!\n";
    point2 = obsMJD.size();
    cout << "It has been reset to a value of " << point2 << "\n";
  }
  
  cout.precision(17);  
  cout << "input configuration file " << configfile << "\n";
  cout << "input observation file " << obsfile << "\n";
  cout << "Method of Herget will be used between points " << point1 << " and " << point2 << "\n";
  cout << "Initial-guess geocentric distances are " << geodist1 << " and " << geodist2 << " AU\n";
  cout << "side length for initial simplex " << simplex_side << " AU\n";
  cout << "MJD, RA, and Dec are in columns " << MJDcol << ", " << RAcol << ", " << Deccol << ", respectively\n";
  cout << "magnitude and astrometric uncertainty (arcsec) are in columns " << magcol << " and " << sigastromcol << "\n";
  cout << "magnitude uncertainty and observatory code are in columns " << sigmagcol << " and " << obscodecol << "\n";
  cout << "ftol is " << ftol << "\n";
  cout << "output file " << outfile << "\n";

  // Read observatory code file
  instream1.open(observatory_code_file);
  if(!instream1) {
    cerr << "can't open input file " << observatory_code_file << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  while (!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    instream1 >> stest;
    stringncopy01(obscode,stest,MINSTRINGLEN);
    instream1 >> obslon;
    instream1 >> plxcos;
    instream1 >> plxsin;
    obs1 = observatory(obscode,obslon,plxcos,plxsin);
    observatory_list.push_back(obs1);
    // Skip the rest of the line
    getline(instream1,lnfromfile);
  }
  instream1.close();
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << observatory_code_file << ":\n";
  
  // Calculate the exact position of the observer at the time of each observation.
  observerpos={};
  for(obsct=0;obsct<obsnum;obsct++) {
    if(obsct==0 || (obsct>0 && obscodevec[obsct]!=obscodevec[obsct-1])) {
      // Observatory has changed: get observatory coordinates for this image.
      stringncopy01(obscode, obscodevec[obsct],MINSTRINGLEN);
      status = obscode_lookup(observatory_list,obscode,obslon,plxcos,plxsin);
      if(status>0) {
	cerr << "ERROR: obscode_lookup failed for observatory code " << obscode << "\n";
	return(3);
      }
    }
    observer_barycoords01LD(obsMJD[obsct], 5, obslon, plxcos, plxsin, earthmjd, Earthpos, outpos);
    observerpos.push_back(outpos);
  }

  maxchisq = maxrchi*(long double)obsnum;
  chisq = Hergetfit01(geodist1, geodist2, simplex_side, ftol, point1, point2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit);
  distpow=0;
  while(chisq>maxchisq && distpow<DISTPOWMAX) {
    // Try making the distances smaller, seeing if we can find a fit
    distpow++;
    geodist1/=DISTPOWSCALE;
    geodist2/=DISTPOWSCALE;
    simplex_side/=DISTPOWSCALE;
    cout << "RE-TRYING FIT WITH SMALLER DISTANCES, TRY " << distpow << ": dist = " << geodist1 << " " << geodist2 << "\n";
    chisq = Hergetfit01(geodist1, geodist2, simplex_side, ftol, point1, point2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit);
  }
  if(chisq>maxchisq) {
    // Re-set the distance parameters
    geodist1*=intpowLD(DISTPOWSCALE,DISTPOWMAX);
    geodist2*=intpowLD(DISTPOWSCALE,DISTPOWMAX);
    simplex_side*=intpowLD(DISTPOWSCALE,DISTPOWMAX);
    distpow=0;
  }
  while(chisq>maxchisq && distpow<DISTPOWMAX) {
    // Try making the distances bigger, seeing if we can find a fit
    distpow++;
    geodist1*=DISTPOWSCALE;
    geodist2*=DISTPOWSCALE;
    simplex_side*=DISTPOWSCALE;
    cout << "RE-TRYING FIT WITH SMALLER DISTANCES, TRY " << distpow << ": dist = " << geodist1 << " " << geodist2 << "\n";
    chisq = Hergetfit01(geodist1, geodist2, simplex_side, ftol, point1, point2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit);
  }
 
  cout << fixed << setprecision(8) << "Orbit a, e, MJD, pos, vel: " << orbit[0]/AU_KM << " " << orbit[1] << " " << orbit[2]  << " " << orbit[3]/AU_KM  << " " << orbit[4]/AU_KM  << " " << orbit[5]/AU_KM  << " " << orbit[6]  << " " << orbit[7]  << " " << orbit[8] << "\n";
  
  outstream1.open(outfile);
  for(obsct=0;obsct<obsnum;obsct++) {
    outstream1 << fixed << setprecision(8) << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " " << fitRA[obsct] << " " << fitDec[obsct] << " " << resid[obsct] << "\n";
  }
  outstream1.close();
  return(0);
}

long double Hergetfit01(long double geodist1, long double geodist2, long double simplex_side, long double ftol, int point1, int point2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid, vector <long double> &orbit)
{
  int Hergetpoint1, Hergetpoint2;
  long double simprange;
  long double simplex[3][2];
  long double simpchi[3];
  long double refdist[2],trialdist[2],bestdist[2];
  long double chisq, bestchi, worstchi, newchi;
  int i,j,worstpoint, bestpoint;
  int unboundsimplex[3];

  // Input points are indexed from 1; apply offset
  Hergetpoint1 = point1-1;
  Hergetpoint2 = point2-1;

  // SETUP FOR DOWNHILL SIMPLEX SEARCH
  // Define initial simplex
  simplex[0][0] = geodist1;
  simplex[0][1] = geodist2;
  simplex[1][0] = geodist1-simplex_side;
  while(simplex[1][0]<=0.0) simplex[1][0] = 0.5L*(simplex[1][0] + geodist1);
  simplex[1][1] = geodist2;
  simplex[2][0] = geodist1;
  simplex[2][1] = geodist2-simplex_side;
  while(simplex[2][1]<=0.0) simplex[2][1] = 0.5L*(simplex[2][1] + geodist2);

  // See if the simplex leads to hyperbolic orbits
  for(i=0;i<3;i++) {
    unboundsimplex[i] = Herget_unboundcheck01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec);
    if(unboundsimplex[i]!=0 && unboundsimplex[i]!=1) {
      // Input is fatally flawed, exit.
      cerr << "ERROR: fatally flawed input to downhill simplex, dists " << simplex[i][0] << " " << simplex[i][1] << "\n";
      cerr << "points " << Hergetpoint1 << " and " << Hergetpoint2 << " out of allowed range 0 to " << obsMJD.size() << "\n";
      return(1);
    }
  }
  while(unboundsimplex[0]==1 && unboundsimplex[1]==1 && unboundsimplex[2]==1) {
    // All the points are bad, shrink all the distances.
    cout << "All points are hyperbolic with simplex:\n";
    for(i=0;i<3;i++) {
      cout << simplex[i][0] << " " << simplex[i][1] << "\n";
      simplex[i][0]*=HERGET_DOWNSCALE;
      simplex[i][1]*=HERGET_DOWNSCALE;
      unboundsimplex[i] = Herget_unboundcheck01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec);
    }
  }
  // If we get here, there must be at least one good point.
  bestpoint = worstpoint = -1;
  for(i=2;i>=0;i--) {
    if(unboundsimplex[i]==0) bestpoint=i;
    if(unboundsimplex[i]==1) worstpoint=i;
  }
  if(bestpoint<0) {
    cerr << "Logically impossible case involving hyperbolic simplex points\n";
    return(1);
  }
  if(worstpoint>=0) {
    cout << "Good simplex point " << bestpoint << ": " << simplex[bestpoint][0] << " " << simplex[bestpoint][1] << "\n";
    // There is at least one bad point.
    for(i=0;i<3;i++) {
      while(unboundsimplex[i]==1 && sqrt(LDSQUARE(simplex[i][0]-simplex[bestpoint][0]) + LDSQUARE(simplex[i][1]-simplex[bestpoint][1])) > MINHERGETDIST) {
	// Bring the bad point closer to a good point until it stops being bad.
	cout << "Modifying bad simplex point " << i << ": " << simplex[i][0] << " " << simplex[i][1] << "\n";
	simplex[i][0] = 0.5L*(simplex[i][0]+simplex[bestpoint][0]);
	simplex[i][1] = 0.5L*(simplex[i][1]+simplex[bestpoint][1]);
	unboundsimplex[i] = Herget_unboundcheck01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec);
      }
    }
  }
  if(unboundsimplex[0]!=0 || unboundsimplex[1]!=0 || unboundsimplex[2]!=0) {
    // We tried everything and still couldn't come up with a good simplex.
    cerr << "Attempt to get viable simplex failed:\n";
    for(i=0;i<3;i++) {
      cerr << simplex[i][0] << " " << simplex[i][1] << " unbound = " << unboundsimplex[i] << "\n";
    }
    cerr << "Aborting\n";
    return(1);
  } else {
    cout << "Good input simplex:\n";
    for(i=0;i<3;i++) {
      cerr << simplex[i][0] << " " << simplex[i][1] << " unbound = " << unboundsimplex[i] << "\n";
    }
  }
  
  for(i=0;i<3;i++) simpchi[i]=LARGERR;
  // Normally, the while loop immediately below should execute only once,
  // but if the input distances don't lead to a bound orbit, it will
  // reduce them and try again.
  while((simpchi[0] == LARGERR || simpchi[1] == LARGERR || simpchi[2] == LARGERR) && simplex[0][0]>MINHERGETDIST) {
    // Calculate chi-square values for each point in the initial simplex
    // Note that the output vectors fitRA, fitDec, and resid are null-wiped
    // internally, so it isn't necessary to wipe them here.
    for(i=0;i<3;i++) {
      cout << "Calling Hergetchi01 with distances " << simplex[i][0] << " " << simplex[i][1] << "\n";
      simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit);
    }
    if(simpchi[0] == LARGERR || simpchi[1] == LARGERR || simpchi[2] == LARGERR) {
      cout << "Hergetchi01 returned failure code with simplex:\n";
      for(i=0;i<3;i++) {
	cout << simplex[i][0] << " " << simplex[i][1] << "\n";
	simplex[i][0]*=HERGET_DOWNSCALE;
	simplex[i][1]*=HERGET_DOWNSCALE;
      }
    }
  }
  if(simplex[0][0]<=MINHERGETDIST) {
    cerr << "ERROR: no acceptable solutions found for the Kepler two-point boundary value problem:\n";
    cerr << "Method of Herget cannot proceed with these data\n";
    return(0);
  }
  cout << "Chi-square value for input distances is " << simpchi[0] << "\n";
  
  // Find best and worst points
  worstpoint=bestpoint=0;
  bestchi = worstchi = simpchi[0];
  for(i=1;i<3;i++) {
    if(simpchi[i]<bestchi) {
      bestchi = simpchi[i];
      bestpoint=i;
    }
    if(simpchi[i]>worstchi) {
      worstchi = simpchi[i];
      worstpoint=i;
    }
  }
  simprange = (worstchi-bestchi)/bestchi;

  // LAUNCH DOWNHILL SIMPLEX SEARCH
  while(simprange>ftol) {
    cout << fixed << setprecision(6) << "Best chi-square value is " << bestchi << ", range is " << simprange << ", vector is " << simplex[bestpoint][0] << " "  << simplex[bestpoint][1] << "\n";
      
    // Try to reflect away from worst point
    // Find mean over all the points except the worst one
    refdist[0] = refdist[1] = 0.0L;
    for(i=0;i<3;i++) {
      if(i!=worstpoint) {
	refdist[0] += simplex[i][0]/2.0L;
	refdist[1] += simplex[i][1]/2.0L;
      }
    }
    // Calculate new trial point
    trialdist[0] = refdist[0] - (simplex[worstpoint][0] - refdist[0]);
    trialdist[1] = refdist[1] - (simplex[worstpoint][1] - refdist[1]);
    // Calculate chi-square value at this new point
    chisq = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit);
    if(chisq<bestchi) {
      // Very good result. Let this point replace worstpoint in the simplex
      for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
      simpchi[worstpoint]=chisq;
      // Extrapolate further in this direction: maybe we can do even better
      trialdist[0] = refdist[0] - 2.0L*(simplex[worstpoint][0] - refdist[0]);
      trialdist[1] = refdist[1] - 2.0L*(simplex[worstpoint][1] - refdist[1]);
      newchi = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit);
      if(newchi<chisq) {
	// Let this even better point replace worstpoint in the simplex
	for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
	simpchi[worstpoint]=chisq;
      }
      // This closes the case where reflecting away from the
      // worst point was a big success.
    } else {
      // Reflecting away from the worst point wasn't great, but
      // we'll see what we can manage.
      if(chisq<worstchi) { 
	// The new point was at least better than the previous worst.
	// Add it to the simplex in place of the worst point
	for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
	simpchi[worstpoint]=chisq;
      } else {
	// The new point was really no good.
	// This is the part of the story where we give up on
	// reflecting away from the worst point, and we try
	// something else.
	// First, try contracting away from the bad point,
	// instead of reflecting away from it.
	trialdist[0] = 0.5L*(simplex[worstpoint][0] + refdist[0]);
	trialdist[1] = 0.5L*(simplex[worstpoint][1] + refdist[1]);
	// Calculate chi-square value at this new point
	chisq = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit);
	if(chisq<worstchi) {
	  // The new point is better than the previous worst point
	  // Add it to the simplex in place of the worst point
	  for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
	  simpchi[worstpoint]=chisq;
	} else {
	  // Even contracting away from the bad point didn't help.
	  // Only one thing left to try: contract toward the best point.
	  // This means each point will become an average of the best
	  // point and its former self.
	  for(i=0;i<3;i++) {
	    if(i!=bestpoint) {
	      simplex[i][0] = 0.5L*(simplex[i][0] + simplex[bestpoint][0]);
	      simplex[i][1] = 0.5L*(simplex[i][1] + simplex[bestpoint][1]);
	      simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit);
	    }
	  }
	  // Close case where nothing worked but contracting around the best point.
	}
	// Close case where reflecting away from the worst point did not work. 
      }
      // Close case where reflecting away from the worst point was not a big success.
    }
    // Identify best and worst points for next iteration.
    worstpoint=bestpoint=0;
    bestchi = worstchi = simpchi[0];
    for(i=1;i<3;i++) {
      if(simpchi[i]<bestchi) {
	bestchi = simpchi[i];
	bestpoint=i;
      }
      if(simpchi[i]>worstchi) {
	worstchi = simpchi[i];
	worstpoint=i;
      }
    }
    simprange = (worstchi-bestchi)/bestchi;
  }
  cout << fixed << setprecision(6) << "Best chi-square value was " << bestchi << ", with geocentric distances " << simplex[bestpoint][0] << " and " << simplex[bestpoint][1] << "\n";
  
  // Perform fit with final best parameters
  chisq = Hergetchi01(simplex[bestpoint][0], simplex[bestpoint][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit);
  return(chisq);
}
