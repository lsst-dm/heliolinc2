// November 01, 2022: test_Herget03
// Test method of Herget routine from J. M. A. Danby,
// Fundamentals of Celestial Mechanics.
// This is the first version where we actually attempt
// a Herget fit.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define DEBUG_2PTBVP 1

int Keplerint_multipoint01(const long double MGsun, const long double mjdstart, const vector <long double> &obsMJD, const point3LD &startpos, const point3LD &startvel, vector <point3LD> &obspos, vector <point3LD> &obsvel);

long double orbitchi01(const point3LD &objectpos, const point3LD &objectvel, const long double mjdstart, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid);

long double TwopointF(long double a, long double k, long double lambda1, long double lambda2, long double deltat, long double Xsign, long double Ysign);

long double TwopointFprime(long double a, long double k, long double lambda1, long double lambda2, long double deltat, long double Xsign, long double Ysign);

int eccen_calc_fast(long double a, point3LD rvec1, point3LD rvec2, long double *e, long double *theta, long double Xsign, long double Ysign);

point3LD Twopoint_Kepler_v1(const long double GMsun, const point3LD startpos, const point3LD endpos, const long double timediff, long double Y, long double *a, long double *e, int itmax);

point3LD geodist_to_3Dpos01(long double RA, long double Dec, point3LD observerpos, long double geodist);

long double Hergetchi01(long double geodist1, long double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid);

static void show_usage()
{
  cerr << "Usage: test_Herget03 -cfg configfile -observations obsfile -fitpoints point1 point2 -geodist dist1 dist2 -simpstep simplex_step -ftol ftol -outfile outfile\n";
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
  point3LD outpos = point3LD(0.0L,0.0L,0.0L);
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> observerpos;
  vector <long double> earthmjd;
  long double geodist1, geodist2, simplex_side, ldval, simprange;
  int point1, point2, Hergetpoint1, Hergetpoint2;
  long double ftol = 1.0e-10;
  vector <long double> obsMJD, obsRA, obsDec, magnitude, sigastrom, sigmag;
  vector <long double> fitRA, fitDec, resid;
  vector <string> obscodevec;
  char obscode[MINSTRINGLEN];
  double obslon = 289.26345L;
  double plxcos = 0.865020L;
  double plxsin = -0.500901L;
  observatory obs1 = observatory("I11",0l,0l,0l);
  vector <observatory> observatory_list = {};
  long double simplex[3][2];
  long double simpchi[3];
  long double refdist[2],trialdist[2],bestdist[2];
  long double chisq, bestchi, worstchi, newchi;
  int worstpoint, bestpoint;
  
  i = j = status = polyorder = configread = fieldnum = badread = 0;
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

  // SETUP FOR DOWNHILL SIMPLEX SEARCH
  // Define initial simplex
  simplex[0][0] = geodist1;
  simplex[0][1] = geodist2;
  simplex[1][0] = geodist1+simplex_side;
  simplex[1][1] = geodist2;
  simplex[2][0] = geodist1;
  simplex[2][1] = geodist2+simplex_side;

  // Input points are indexed from 1; apply offset
  Hergetpoint1 = point1-1;
  Hergetpoint2 = point2-1;

  // Calculate chi-square values for each point in the initial simplex
  // Note that the output vectors fitRA, fitDec, and resid are null-wiped
  // internally, so it isn't necessary to wipe them here.
  for(i=0;i<3;i++) {
    simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid);
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
    chisq = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid);
    if(chisq<bestchi) {
      // Very good result. Let this point replace worstpoint in the simplex
      for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
      simpchi[worstpoint]=chisq;
      // Extrapolate further in this direction: maybe we can do even better
      trialdist[0] = refdist[0] - 2.0L*(simplex[worstpoint][0] - refdist[0]);
      trialdist[1] = refdist[1] - 2.0L*(simplex[worstpoint][1] - refdist[1]);
      newchi = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid);
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
	chisq = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid);
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
	      simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid);
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
  cout << "Best chi-square value was " << bestchi << ", with geocentric distances " << simplex[bestpoint][0] << " and " << simplex[bestpoint][1] << "\n";
  // Perform fit with final best parameters
  chisq = Hergetchi01(simplex[bestpoint][0], simplex[bestpoint][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid);
  for(obsct=0;obsct<obsnum;obsct++) {
    cout << obsMJD[obsct] << " " << obsRA[obsct] << " " << obsDec[obsct] << " " << fitRA[obsct] << " " << fitDec[obsct] << " " << resid[obsct] << "\n";
  }

  return(0);
}

// Keplerint_multipoint01: November 02, 2022: Like Keplerint, but does the
// calculation for a bunch of points simultaneously. Note that
// we assume the observation times and mjdstart are in UT1, which
// means that JPL Horizons state vectors cannot be used directly
// for mjdstart: one would have to correct the nominal value of
// mjdstart corresponding to the JPL Horizons state vectors.
// Description of ancestor program Keplerint:
// Integrate an orbit assuming we have a Keplerian 2-body problem
// with all the mass in the Sun, and the input position and velocity
// are relative to the Sun.
int Keplerint_multipoint01(const long double MGsun, const long double mjdstart, const vector <long double> &obsMJD, const point3LD &startpos, const point3LD &startvel, vector <point3LD> &obspos, vector <point3LD> &obsvel)
{
  long double e,E,a,lscalar,r0,v0,r1,v1;
  point3LD lvec = point3LD(0L,0L,0L);
  point3LD lunit = point3LD(0L,0L,0L);
  point3LD r0unit = point3LD(0L,0L,0L);
  point3LD r1unit = point3LD(0L,0L,0L);
  point3LD v1unit = point3LD(0L,0L,0L);
  point3LD targpos = point3LD(0L,0L,0L);
  point3LD targvel = point3LD(0L,0L,0L);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0L;
  long double costheta,theta0,theta1,radvel,cospsi,psi;
  costheta = theta0 = theta1 = radvel = cospsi = psi = 0L;
  long double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0L;
  long double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0L;
  int status = 0;
  long double junkra,junkdec,sinev,thetav,v1ra,v1dec;
  junkra = junkdec = sinev = thetav = v1ra = v1dec = 0L;
  int obsct=0;
  int obsnum = obsMJD.size();
 
  // Calculate scalar input position
  r0 = vecabs3LD(startpos);
  v0 = vecabs3LD(startvel);
  
  // Calculate specific energy and angular momentum
  E = 0.5L*v0*v0 - MGsun/r0;
  lvec = crossprod3LD(startpos,startvel);
  lscalar = sqrt(dotprod3LD(lvec,lvec));
  if(E>=0L) {
    //cerr << "ERROR: Keplerint finds positive total energy: " << E << "\n";
    return(1);
  }
		 
  // Calculate a and e: orbital semimajor axis and eccentricity.
  cout.precision(17);
  a = -MGsun*0.5L/E;
  e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/MGsun/MGsun);
  if(e<0L || e>=1.0L) {
    cerr << "ERROR: Keplerint finds eccentricity out of range: " << e << "\n";
    return(1);
  }

  //cout << "semimajor axis = " << a/AU_KM << " and eccentricity = " << e << "\n";
	       
  // Calculate angle theta0 from perihelion using simple ellipse geometry.
  if(e>0L) {
    costheta = ((a-a*e*e)/r0 - 1.0L)/e;
    if(costheta>=-1.0L && costheta<=1.0L) theta0 = acos(costheta);
    else {
      cerr << "ERROR: Keplerint finds costheta = " << costheta << "\n";
      return(1);
    }
  }
  radvel = dotprod3LD(startpos,startvel)/r0;
  //cout << "Radial velocity = " << radvel << " km/sec\n";
  
  if(radvel>=0) {
    // We are moving outward from perihelion: theta will be correct.
    //cout << "Moving outward from perihelion.\n";
    ;
  } else {
    // We are moving inward towards perihelion: theta needs adjustment.
    theta0 = 2.0L*M_PI - theta0;
    //cout << "Moving inward towards perihelion.\n";
  }
  //cout << "theta0 = " << theta0*DEGPRAD << "\n";
  
  // Calculate Goldstein's psi variable from theta.
  cospsi = (costheta + e)/(1.0L + costheta*e);
  if(cospsi>=-1.0L && cospsi<=1.0L) psi = acos(cospsi);
  else {
    cerr << "ERROR: Keplerint finds cospsi = " << cospsi << "\n";
    return(1);
  }
  if(radvel<0) {
    // We are moving inward towards perihelion: psi needs adjustment.
    psi = 2.0L*M_PI - psi;
  }
  //cout << "psi = " << psi*DEGPRAD << "\n";
 
  // Calculate time since perihelion using psi.
  omega = sqrt(MGsun/(a*a*a));
  //cout << "Period = " << 2.0L*M_PI/omega/SOLARDAY/365.25 << " years\n";
  t0omega = psi - e*sin(psi);
  //cout << "t0omega = " << t0omega;

  // Loop on all times-of-observation, and calculate the target position at those times
  obspos = obsvel = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    // The new time t1 for which we want to re-evaluate psi is
    // given by t0 + obsMJD[obsct]-mjdstart.
    t1omega = t0omega + (obsMJD[obsct]-mjdstart)*SOLARDAY*omega;
    while(t1omega > 2.0L*M_PI) t1omega -= 2.0L*M_PI;
    while(t1omega < 0.0L) t1omega += 2.0L*M_PI;
    // Solve Kepler's equation for psi(t1)
    psi = kep_transcendental(t1omega,e,KEPTRANSTOL);
    cospsi = cos(psi);
    // Calculate theta(t1) from psi(t1)
    if(1.0L - e*cospsi != 0.0L) {
      costheta = (cospsi - e)/(1.0L - e*cospsi);
      if(costheta >= -1.0L && costheta <= 1.0L) theta1 = acos(costheta);
      else if(costheta < -1.0L) {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = M_PI;
      } else {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = 0.0L;
      }
      if(psi>M_PI && theta1<=M_PI) theta1 = 2.0L*M_PI - theta1;
    } else {
      cerr << "Warning: e*cos(psi) = " << e*cospsi << " so 1 - e*cos(psi) = " << 1.0L - e*cospsi << "\n";
      theta1 = 0.0L;
    }
    while(theta1<0.0L) theta1 += 2.0L*M_PI;
    while(theta1>=2.0L*M_PI) theta1 -= 2.0L*M_PI;

    // Calculate r(t1) from psi(t1)
    r1 = a*(1.0L - e*cospsi);
    // Calculate v1 from r1 and the known energy
    v1 = sqrt((E +  MGsun/r1)*2.0L);
  
    // Use vector algebra to find the full vector r(t1).
    // This vector is perpendicular to lvec, and is angled by theta1-theta0
    // relative to startpos.
    // Convert angular momentum vector to spherical coordinates
    celedeproj01LD(lvec,&lra,&ldec); // Note that output is in degrees.
    celedeproj01LD(startpos,&r0ra,&r0dec); // Note that output is in degrees.
    // Transform the starting unit vector into a coordinate system with
    // the angular momentum vector at the pole, and the old pole at RA=0
    poleswitch01LD(r0ra/DEGPRAD,r0dec/DEGPRAD,lra/DEGPRAD,ldec/DEGPRAD,0.0L,newra,newdec); // Output is radians
    // Rotate starting unit vector around the angular momentum axis by
    // the calculated angle.
    newra += theta1-theta0;
    // The unit vector for the new position r1 is on the equator at this RA,
    // in the coordinate system that has the angular momentum vector at the pole.
    // Convert back to the original coordinate system.
    poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,r1ra,r1dec); // Output is radians
    // Now for the velocity. If the velocity is at right angle to the vector r1,
    // the product v1*r1 is the angular momentum. Otherwise, l/(v1*r1) is the sine
    // of the angle between v1 and r1.

    sinev = lscalar/v1/r1;
    if(sinev>=1.0L) thetav = 0.5L*M_PI;
    else if(sinev<0.0L) {
      cerr << "ERROR: negative angular momentum?\nv1,r1,v1*r1,lscalar,sinev = " << v1 << ", " << r1 << ", " << v1*r1 << ", " << lscalar << ", " << sinev << "\n";
      thetav = 0.0L;
    }
    else thetav = asin(sinev);
    if(theta1<=M_PI) {
      // Outward bound from perihelion.
      newra += thetav;
    } else {
      // Inward bound to perihelion
      newra += (M_PI - thetav);
    }
    poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,v1ra,v1dec); // Output is radians

    r1unit = celeproj01LD(r1ra*DEGPRAD,r1dec*DEGPRAD);
    v1unit =celeproj01LD(v1ra*DEGPRAD,v1dec*DEGPRAD);
  
    targpos.x = r1unit.x*r1;
    targpos.y = r1unit.y*r1;
    targpos.z = r1unit.z*r1;
    targvel.x = v1unit.x*v1;
    targvel.y = v1unit.y*v1;
    targvel.z = v1unit.z*v1;
    obspos.push_back(targpos);
    obsvel.push_back(targvel);
  }
    
  return(0);
}


// orbitchi01: November 02, 2022:
// Get chi-square value based on input state vectors,
// using 2-body Keplerian integration, rather than n-body.
// Input state vectors are expected to be in km and km/sec.
// Note that the output vectors fitRA, fitDec, and resid are
// null-wiped inside orbitchi01, so it isn't necessary for the
// calling function to wipe them.
long double orbitchi01(const point3LD &objectpos, const point3LD &objectvel, const long double mjdstart, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid)
{
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  int obsct;
  int obsnum = obsMJD.size();
  long double light_travel_time;
  point3LD outpos = point3LD(0,0,0);
  long double outRA=0L;
  long double outDec=0L;
  long double ldval=0L;
  long double chisq=0L;
  double dval;
  resid = fitRA = fitDec = {};
  int status=0;
  
  // Integrate orbit.
  status=0;
  status = Keplerint_multipoint01(GMSUN_KM3_SEC2,mjdstart,obsMJD,objectpos,objectvel,obspos,obsvel);
  if(status!=0) {
    // Keplerint_multipoint01 failed, likely because input state vectors lead
    // to an unbound orbit.
    return(LARGERR);
  }
  for(obsct=0;obsct<obsnum;obsct++) {
    // Initial approximation of the coordinates relative to the observer
    outpos.x = obspos[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - observerpos[obsct].z;
    // Initial approximation of the observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Convert to meters and divide by the speed of light to get the light travel time.
    light_travel_time = ldval*1000.0/CLIGHT;
    // Light-travel-time corrected version of coordinates relative to the observer
    outpos.x = obspos[obsct].x - light_travel_time*obsvel[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - light_travel_time*obsvel[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - light_travel_time*obsvel[obsct].z - observerpos[obsct].z;
    // Light-travel-time corrected observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    // Project onto the celestial sphere.
    stateunitLD_to_celestial(outpos, outRA, outDec);
    dval = distradec01(obsRA[obsct],obsDec[obsct],outRA,outDec);
    dval *= 3600.0L; // Convert to arcsec
    fitRA.push_back(outRA);
    fitDec.push_back(outDec);
    resid.push_back(dval);
  }
  chisq=0.0L;
  for(obsct=0;obsct<obsnum;obsct++) {
    chisq += LDSQUARE(resid[obsct]/sigastrom[obsct]);
  }
  return(chisq);
}

// TwopointF: October 26, 2022:
// Given input values for k = sqrt(GMsun), delta-t, lambda1 = sqrt(r1+r2+c),
// lambda2 = sqrt(r1+r2-c), and the semimajor axis a in km, evaluate a function
// of the semimajor axis a that, if minimized, solves the Kepler two-point
// boundary value problem. Distances are in units of km, but delta-t is in units
// of days.
long double TwopointF(long double a, long double k, long double lambda1, long double lambda2, long double deltat, long double Xsign, long double Ysign)
{
  long double kterm, asin1, sq1, asin2, sq2, Xterm;
  //  long double dkterm, dasin1, dsq1, dasin2, dsq2,a2;
  // Evaluate term-by-term, for clarity
  kterm = k*deltat*SOLARDAY/a/sqrt(a);
  asin1 = Xsign*2.0L*asin(lambda1/2.0L/sqrt(a));
  sq1 = Xsign*lambda1*sqrt(4.0L*a-lambda1*lambda1)/2.0L/a;
  asin2 = Ysign*2.0L*asin(lambda2/2.0L/sqrt(a));
  sq2 = Ysign*lambda2*sqrt(4.0L*a-lambda2*lambda2)/2.0L/a;
  if(Xsign==-1.0L) Xterm = -2.0L*M_PI;
  else Xterm = 0.0L;
  
  //a2 = a+1000.0L;
  //dkterm = k*deltat*SOLARDAY/a2/sqrt(a2) - k*deltat*SOLARDAY/a/sqrt(a);
  //dasin1 = 2.0L*asin(lambda1/2.0L/sqrt(a2)) - 2.0L*asin(lambda1/2.0L/sqrt(a)); 
  //dsq1 = lambda1*sqrt(4.0L*a2-lambda1*lambda1)/2.0L/a2 - lambda1*sqrt(4.0L*a-lambda1*lambda1)/2.0L/a;
  //dasin2 = 2.0L*asin(lambda2/2.0L/sqrt(a2)) - 2.0L*asin(lambda2/2.0L/sqrt(a));
  //dsq2 = lambda2*sqrt(4.0L*a2-lambda2*lambda2)/2.0L/a2 - lambda2*sqrt(4.0L*a-lambda2*lambda2)/2.0L/a;

  //dkterm /= 1000.0L;
  //dasin1 /= 1000.0L;
  //dsq1 /= 1000.0L;
  //dasin2 /= 1000.0L;
  //dsq2 /= 1000.0L;

  //cout << "Numder by term: " << dkterm << " " << dasin1 << " " << dsq1 << " " << dasin2 << " " << dsq2 << "\n";
  
  return(kterm - asin1 + sq1 + asin2 - sq2 +  Xterm);
}

 // TwopointFprime: October 26, 2022:
// Given input values for k = sqrt(GMsun), delta-t, lambda1 = sqrt(r1+r2+c),
// lambda2 = sqrt(r1+r2-c), and the semimajor axis a in km, evaluate the
// derivative with respect to the semimajor axis 'a' a function that,
// if minimized, solves the Kepler two-point
// boundary value problem. Distances are in units of km, but delta-t is in units
// of days.
long double TwopointFprime(long double a, long double k, long double lambda1, long double lambda2, long double deltat, long double Xsign, long double Ysign)
{
  long double kterm, asin1, sq1, asin2, sq2, foural1, foural2;
  // Terms that appear a lot, pre-calculated for simplicity
  foural1 = 4.0L*a - lambda1*lambda1;
  foural2 = 4.0L*a - lambda2*lambda2;
  
  // Evaluate term-by-term, for clarity
  kterm = -3.0L/2.0L*k*deltat*SOLARDAY/a/a/sqrt(a);
  asin1 = Xsign*lambda1/a/sqrt(foural1);
  sq1 = Xsign*(lambda1*lambda1*lambda1 - 2.0L*a*lambda1)/(2.0L*a*a*sqrt(foural1));
  asin2 = -Ysign*lambda2/a/sqrt(foural2);
  sq2 = Ysign*(2.0L*a*lambda2 - lambda2*lambda2*lambda2)/(2.0L*a*a*sqrt(foural2));
  
  //cout << "Fprime by terms: " << kterm << " " << asin1 << " " << sq1 << " " << asin2 << " " << sq2 << " " << kterm + asin1 + sq1 + asin2 + sq2 << "\n";
  
  return(kterm + asin1 + sq1 + asin2 + sq2);
}

  
// eccen_calc_fast: November 01, 2022:
// Given the semimajor axis of an ellipse, and two 3-D points
// on the ellipse, calculate the eccentricity of
// the ellipse, or -1.0 if there is no valid elliptical solution
// This version uses only the eccentric anomaly (psi), and not
// the true anomaly (theta), resulting in a single solution).
// The eccentric anomaly solution seemed on development tests
// to be slightly less accurate, but the differences appeared
// far too small to matter, and not worth the reduction in
// speed.
int eccen_calc_fast(long double a, point3LD rvec1, point3LD rvec2, long double *e, long double *theta, long double Xsign, long double Ysign)
{
  long double r1 = vecabs3LD(rvec1);
  long double r2 = vecabs3LD(rvec2);
  long double costheta,thetatest,r2test;
  long double cos_delta_theta = dotprod3LD(rvec1,rvec2)/r1/r2;
  long double delta_theta = acos(cos_delta_theta);
 
  // Parameters for solution using the eccentric anomaly
  long double c = sqrt(LDSQUARE(rvec2.x-rvec1.x) + LDSQUARE(rvec2.y-rvec1.y) + LDSQUARE(rvec2.z-rvec1.z));
  long double lambda1 = r1+r2+c;
  long double lambda2 = r1+r2-c;
  long double eps_star = 2.0L*asin(sqrt(lambda1/4.0L/a));
  long double delt_star = 2.0L*asin(sqrt(lambda2/4.0L/a));
  long double eps = Ysign*eps_star;
  long double delt = M_PI*(1.0L - Xsign) + Xsign*delt_star;
  long double delta_psi = eps - delt;
  long double amr1 = a-r1;
  long double amr2 = a-r2;
  long double cosdpsi = cos(delta_psi);
  long double cos2dpsi = LDSQUARE(cosdpsi);
  long double sindpsi = sin(delta_psi);
  long double sin2dpsi = LDSQUARE(sindpsi);
  
  long double eea = sqrt(cos2dpsi*amr1*amr1 - 2.0L*cosdpsi*amr1*amr2 + amr2*amr2 + sin2dpsi*amr1*amr1)/sindpsi/a;
  costheta = (a*(1.0L-eea*eea) - r1)/eea/r1;
  thetatest = acos(costheta);
  r2test = a*(1-eea*eea)/(1.0L +  eea*cos(thetatest + delta_theta));
  long double delta1 = fabs(r2test-r2);
  r2test = a*(1-eea*eea)/(1.0L +  eea*cos(2.0L*M_PI - thetatest + delta_theta));
  long double delta2 = fabs(r2test-r2);
  //cout << "eea deltas: " << delta1 << " " << delta2 << "\n";
  if(delta2<delta1) thetatest = 2.0L*M_PI - thetatest;
  //cout << "eea theta = " << thetatest*DEGPRAD << "\n";
  *e = eea;
  *theta = thetatest;

  return(0);
}

// Twopoint_Kepler_v1: November 01, 2022:
// Given two points in an object's orbit (as 3-D Cartesian
// vectors relative to the sun), and the time it takes to
// move from the first point to the second, solve for the 
// semimajor axis a of the object's Keplerian orbit, and
// from these derive the initial velocity in units of km/sec.
// Input positions are in units of km,
// timediff is in units of days.
// This code closely follows the derivation in Section 6.11 of
// J. M. A. Danby's Foundations of Celestial Mechanics.
point3LD Twopoint_Kepler_v1(const long double GMsun, const point3LD startpos, const point3LD endpos, const long double timediff, const long double Ysign, long double *a, long double *e, int itmax)
{
  // General quantities, and those having to do with solving for the semimajor axis.
  long double r1 = vecabs3LD(startpos);
  long double r2 = vecabs3LD(endpos);
  point3LD pdiff = point3LD(endpos.x - startpos.x, endpos.y - startpos.y, endpos.z - startpos.z);
  point3LD v1 = point3LD(0.0L,0.0L,0.0L);
  long double c = vecabs3LD(pdiff);
  long double lambda1 = sqrt(r1+r2+c);
  long double lambda2 = sqrt(r1+r2-c);
  long double k = sqrt(GMsun);
  
  // Determine the sign-specifier X
  long double ac = (r1+r2+c)/4.0L;
  long double nc = k/ac/sqrt(ac); // This is in radians per second.
  long double dc = 2.0L*asin(sqrt(lambda2/lambda1));
  long double dtc = (M_PI - dc + sin(dc))/nc;
  long double X=1.0L;
  long double Y=Ysign;
  long double aorb;
  int ai;
  long double delta_aorb;
  long double f;
  long double fprime;
  long double f2;
  long double fprime2;
  long double eccen,thetaperi;
  
  if(timediff*SOLARDAY>dtc) X=-1.0L;

  if(DEBUG_2PTBVP>1) cout << "Initial setup stuff:\n";
  if(DEBUG_2PTBVP>1) cout << "r1 = " << r1/AU_KM << ", r2 = " << r2/AU_KM << ", c = " << c/AU_KM << ", lambdas = " << lambda1 << ", " << lambda2 << ", k = " << k << ", X = " << X << ", Y = " << Y << "\n";

  //for(ai=10;ai<=200;ai++) {
  //  aorb = (long double)ai*0.01L*AU_KM;
  //  if(aorb > lambda1*lambda1/4.0L) {
  //    cout << ai << " " << aorb/AU_KM << " " << TwopointF(aorb, k, lambda1, lambda2, timediff) << "\n";
  //  }
  //}

  aorb = 0.5L*r1 + 0.5L*r2;
  f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
  fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);

  int itnum=0;
  delta_aorb = -f/fprime;
  while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
  while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
    // Don't let it drop into the NAN range
    delta_aorb /= 2.0L;
  }
 
  if(DEBUG_2PTBVP>1) cout << "0th iteration:\n";
  if(DEBUG_2PTBVP>1) cout << "a = " << aorb << " = " << aorb/AU_KM << " AU, f = " << f << ", fprime = " << fprime << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
  
  while((fabs(f) > KEPTRANSTOL || fabs(delta_aorb) > KEPTRANSTOL*AU_KM) && aorb<1000000.0L*AU_KM && itnum<itmax) {
    aorb += delta_aorb;
    f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    f2 = TwopointF(aorb+1000.0L, k, lambda1, lambda2, timediff);
    fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    fprime2 = (f2-f)/1000.0L;
    delta_aorb = -f/fprime;
    while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
    while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
      // Don't let it drop into the NAN range
      delta_aorb /= 2.0L;
    }
    // Calculate the eccentricity
    // eccen_calc_precise(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    // cout << "iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = " << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    eccen_calc_fast(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    if(DEBUG_2PTBVP>1) cout << "iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = " << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    itnum++;
  }
  if(fabs(f) > KEPTRANSTOL || fabs(delta_aorb > KEPTRANSTOL*AU_KM)) {
    // We never found a solution. Re-start near the minimum allowable:
    aorb = lambda1*lambda1/3.99L;
    f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
    fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);
    itnum=0;
    delta_aorb = -f/fprime;
    while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
    while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
      // Don't let it drop into the NAN range
      delta_aorb /= 2.0L;
    }
 
    if(DEBUG_2PTBVP>1) cout << "2nd try, 0th iteration, :\n";
    if(DEBUG_2PTBVP>1) cout << "a = " << aorb << " = " << aorb/AU_KM << " AU, f = " << f << ", fprime = " << fprime << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    while((fabs(f) > KEPTRANSTOL || fabs(delta_aorb) > KEPTRANSTOL*AU_KM)  && aorb<1000000.0L*AU_KM && itnum<itmax) {
    aorb += delta_aorb;
    f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    f2 = TwopointF(aorb+1000.0L, k, lambda1, lambda2, timediff);
    fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    fprime2 = (f2-f)/1000.0L;
    delta_aorb = -f/fprime;
    while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
    while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
      // Don't let it drop into the NAN range
      delta_aorb /= 2.0L;
    }
    // Calculate the eccentricity
    // eccen_calc_precise(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    // cout << "2nd try, iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = " << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    eccen_calc_fast(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    if(DEBUG_2PTBVP>1) cout << "2nd try, iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = " << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    itnum++;
    }
  }
  *a = aorb;
  *e = eccen;

  // Quantities having to do with solving for the velocity
  long double alpha,beta,gamma,ffunc,gfunc;
  
  // Calculate the velocity
  alpha = 0.5L*lambda1;
  beta = 0.5L*lambda2;
  gamma = alpha*sqrt(1.0L - lambda2*lambda2/4.0L/aorb) -  X*Y*beta*sqrt(1.0L - lambda1*lambda1/4.0L/aorb);
  ffunc = 1.0L - 2*gamma*gamma/r1;
  gfunc = (4.0L/k)*Y*alpha*beta*gamma;

  v1.x = (1.0L/gfunc)*(endpos.x - ffunc*startpos.x);
  v1.y = (1.0L/gfunc)*(endpos.y - ffunc*startpos.y);
  v1.z = (1.0L/gfunc)*(endpos.z - ffunc*startpos.z);
  
  return(v1);
}


// geodist_to_3Dpos01: November 02, 2022:
// Given and input RA, Dec position, observer's position, and
// geocentric distance, output a 3D position for the object.
// The geocentric distance is expected in AU, but the output
// position vector has units of km.
point3LD geodist_to_3Dpos01(long double RA, long double Dec, point3LD observerpos, long double geodist)
{
  point3LD baryvec = point3LD(0.0L,0.0L,0.0L);
  
  celestial_to_stateunitLD(RA, Dec, baryvec);
  baryvec.x = observerpos.x + baryvec.x*geodist*AU_KM;
  baryvec.y = observerpos.y + baryvec.y*geodist*AU_KM;
  baryvec.z = observerpos.z + baryvec.z*geodist*AU_KM;
  return(baryvec);
}

// Hergetchi01: November 02, 2022:
// Calculate the chi-square value between input observations and
// the result of orbit determination using the Method of Herget,
// where the Method of Herget is performed between Hergetpoint1 and
// Hergetpoint2 of the input observations, under the assumption that
// the object was a distance geodist1 from the observer when
// the observation corresponding to Hergetpoint1 was made, and
// a distance geodist2 from the observer when the observation
// corresponding to Hergetpoint2 was made. The vectors are all
// in units of km and km/sec, but geodist1 and geodist2 are
// in units of AU.
// Note that the output vectors fitRA, fitDec, and resid are
// null-wiped inside Hergetchi01, so it isn't necessary for the
// calling function to wipe them.
long double Hergetchi01(long double geodist1, long double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid)
{
  int numobs = obsMJD.size();
  if(obsRA.size() != numobs || obsDec.size() != numobs || sigastrom.size() != numobs || observerpos.size() != numobs) {
    cerr << "ERROR: Hergetchi01 finds unequal lenths among input vectors:\n";
    cerr << "observed MJD, RA, Dec, sigastrom, and observerpos have lengths " << numobs << " " << obsRA.size() << " " << obsDec.size() << " " <<  sigastrom.size() << " " << observerpos.size() << "\n";
    return(-99.0L);
  }
  if(Hergetpoint2<=Hergetpoint1 || Hergetpoint1<0 || Hergetpoint2>=numobs) {
    cerr << "ERROR: Hergetchi01 has invalid input reference points:\n";
    cerr << "Starting point " << Hergetpoint1 << " and ending point " << Hergetpoint2 << ", where allowed range is 0 to " << numobs-1 << "\n";
    return(-99.0L);
  }
  
  point3LD startpos = geodist_to_3Dpos01(obsRA[Hergetpoint1], obsDec[Hergetpoint1], observerpos[Hergetpoint1], geodist1);
  point3LD endpos = geodist_to_3Dpos01(obsRA[Hergetpoint2], obsDec[Hergetpoint2], observerpos[Hergetpoint2], geodist2);
  // Time difference should include a light-travel-time correction. The sign is determined
  // by the fact that if the object gets further away, the object time moves backward
  // relative to the observer time. Hence, if the object gets further away (i.e., geodist2>geodist1),
  // the object experiences less time than the observer beween the two observations, because
  // the observer is looking further back in time at the second observation.
  long double deltat = obsMJD[Hergetpoint2] - obsMJD[Hergetpoint1] - (geodist2-geodist1)/CLIGHT_AUDAY;
  long double a,e;
  point3LD startvel = Twopoint_Kepler_v1(GMSUN_KM3_SEC2, startpos, endpos, deltat, 1.0L, &a, &e, 100);
  long double orbchi;
  // Note that the output vectors fitRA, fitDec, and resid are null-wiped
  // internally in orbitchi01, so it isn't necessary to wipe them here.
  orbchi = orbitchi01(startpos, startvel, obsMJD[Hergetpoint1]-geodist1/CLIGHT_AUDAY, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid);
  return(orbchi);
}
