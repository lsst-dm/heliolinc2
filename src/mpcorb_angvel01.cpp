// November 20, 2023
// Read a file with lines in the format of the Minor Planet Center's MPCORB.DAT,
// and calculate the angular velocities of the asteroids 

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define VERBOSE 0

static void show_usage()
{
  cerr << "Usage: mpcorb_angvel01 -config configuration file -mpcorb input MPCORB file -outfile output file\n";
}
  
int main(int argc, char *argv[])
{
  string mpcorb,outfile;
  string lnfromfile;
  string desig;
  string configfile;
  string planetfile;
  ofstream outstream1;
  ifstream instream1;
  string epoch;
  int ondata=0;
  int polyorder=0;
  long double semimaj_axis = 0l; // in AU
  long double eccentricity = 0l; // unitless
  long double inclination = 0l;  // in degrees
  long double long_ascend_node = 0l; // Longitude of the ascending node, in degrees
  long double arg_perihelion = 0l;   // Argument of perihelion, in degrees
  long double mean_anom = 0l;        // Mean anomaly at the epoch, in degrees
  long double mjd_epoch = 0l;        // Epoch for the orbit in MJD
  long double mean_daily_motion = 0l; // in degrees/day
  long double H = 0l; 
  long double G = 0l;
  int i=0;
  long j=0;
  int status=0;
  int planetnum=0;
  int planetct=0;
  int pctSun=0;
  int pctEarth=0;
  long double ldval=0.0L;
  point3LD outpos = point3LD(0L,0L,0L);
  point3LD outvel = point3LD(0L,0L,0L);
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  point3LD Earth2targ = point3LD(0,0,0);
  point3LD Earth2targvel = point3LD(0,0,0);
  point3LD Earth2sun = point3LD(0,0,0);
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> targpos;
  vector <point3LD> targvel;
  vector <long double> targMJD;
  int planetfile_startpoint=0;
  int planetfile_endpoint=0;
  long double timespan = 0L;
  long double mjdstart = 0L;
  long double mjdend = 0L;
  long double minperidist=0L;
  long double maxperidist=1e30L;
  long double minapdist=0L;
  long double maxapdist=1e30L;
  long double incmax=1e30L;
  int configread=0;
  asteroid_orbitLD oneorb = asteroid_orbitLD(desig,semimaj_axis, eccentricity, inclination, long_ascend_node, arg_perihelion, mean_anom, mjd_epoch, mean_daily_motion,H,G);
  long double targlambda1, targbeta1, targlambda2, targbeta2, sunlambda, sunbeta;
  targlambda1 = targbeta1 = targlambda2 = targbeta2 = sunlambda = sunbeta = 0L;
  long double deltat = 0.1L;
  double lambdavel,betavel;
  lambdavel = betavel = 0.0l;
  int elongstep,latstep;
  elongstep=10;
  latstep=5;
  int latmax = 30;
  vector <int> elongcen;
  vector <int> latcen;
  point2d angvelpt = point2d(0,0);
  vector <point2d> velvec;
  vector <vector <point2d>> velmat;
  int lat,elong,latnum,elongnum,latct,elongct;
  double dlatct,delongct;
  char nsfx[64];
  string numsuffix,outname;
  
  if(argc<5) {
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
	instream1.open(configfile);
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
	    Earthvel = tempvel;
	  } else if(planetct == pctSun) {
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
	status=readconfigLD(instream1,&incmax);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&incmax);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Max value for orbital inclination read as " << incmax << " degrees\n";
	// Read default integration span
	status=readconfigLD(instream1,&timespan);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&timespan);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default timespan read as " << timespan << " days\n";
	// Read range +/- for ecliptic latitude
	status=readconfigint(instream1,&latmax);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&latmax);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Range +/- for ecliptic latitude read as " << latmax << "\n";
	// Read step size for ecliptic latitude
	status=readconfigint(instream1,&latstep);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&latstep);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Step size for ecliptic latitude read as " << latstep << "\n";
	// Read step size for solar elongation
	status=readconfigint(instream1,&elongstep);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&elongstep);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Step size for solar elongation read as " << elongstep << "\n";
	
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
    if(string(argv[i]) == "-mpcorb") {
      if(i+1 < argc) {
	//There is still something to read;
	mpcorb=argv[++i];
	i++;
      }
      else {
	cerr << "Input MPCORB file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile") {
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
  cout << "input configuration file " << configfile << "\n";
  cout << "Perihelion distance must be between " << minperidist << " and " << maxperidist << "\n";
  cout << "Aphelion distance must be between " << minapdist << " and " << maxapdist << "\n";
  cout << "Inclination must be less than = " << incmax << "\n";
  cout << "Time range for integration is " << timespan << " days\n";
  cout << "Range for ecliptic latitude is +/-" << latmax << " degrees\n";
  cout << "Step size for ecliptic latitude is " << latstep << " degrees\n";
  cout << "Step size for solar elongation is " << elongstep << " degrees\n";
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

  i=elongstep/2;
  if(2*i != elongstep) {
    cerr << "ERROR: step size in solar elongation must be an even number of degrees.\n";
    cerr << "Supplied value of " << elongstep << " is invalid\n";
    return(1);
  }
  
  // Load vectors holding central ecliptic latitudes and solar elongations
  latnum = 2*latmax/latstep + 1;
  if((latnum-1)*latstep != 2*latmax) {
    cerr << "ERROR: latitude step size " << latstep << " is not evenly divisible\n";
    cerr << "into the range +/- " << latmax << "\n";
    return(1);
  }
  elongnum = 360/elongstep;
  if(elongnum*elongstep != 360) {
    cerr << "ERROR: solar elongation step size " << elongstep << " is not evenly divisible into 360\n";
    return(1);
  }
  cout << "There will be " << latnum << " steps in ecliptic latitude, and " << elongnum << " steps in solar elongation,\n";
  cout << "for a total of " << latnum*elongnum << " steps\n";
  lat=-latmax;
  while(lat<=latmax) {
    elong = elongstep/2;
    while(elong<=(360-elongstep/2)) {
      latcen.push_back(lat);
      elongcen.push_back(elong);
      velmat.push_back(velvec);
      elong += elongstep;
    }
    lat += latstep;
  }

  
  // Open MPCORB file
  instream1.open(mpcorb);
  // Find start of data
  ondata=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad() && !ondata) {
    getline(instream1,lnfromfile);
    ondata=1;
    for(i=0;i<4;i++) if(lnfromfile[i]!='-') ondata=0;
  }
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    status = read_orbline(instream1, oneorb);

    cout << oneorb.desig << " " << oneorb.semimaj_axis << " " << oneorb.eccentricity << " " << oneorb.inclination << " " << oneorb.long_ascend_node << " " << oneorb.arg_perihelion << " " << oneorb.mean_anom << " " << oneorb.mjd_epoch << " " << oneorb.mean_daily_motion << "\n";
    
    if(status==0 && oneorb.semimaj_axis*(1.0-oneorb.eccentricity) >= minperidist && oneorb.semimaj_axis*(1.0-oneorb.eccentricity) <= maxperidist && oneorb.semimaj_axis*(1.0+oneorb.eccentricity) >= minapdist && oneorb.semimaj_axis*(1.0+oneorb.eccentricity) <= maxapdist && oneorb.inclination <= incmax) {
     
      // We will integrate this orbit from the epoch to the epoch+timespan
      mjdstart = oneorb.mjd_epoch;
      mjdend = mjdstart+timespan;
      
      // Match mjdstart and mjdend to planet file.
  
      planetfile_startpoint = planetfile_endpoint = -99;
      for(j=0;j<long(planetmjd.size());j++) {
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

      Kepler2dyn(planetmjd[planetfile_startpoint], oneorb, outpos,  outvel);
      outpos.x *= AU_KM;
      outpos.y *= AU_KM;
      outpos.z *= AU_KM;
      outvel.x *= AU_KM/SOLARDAY;
      outvel.y *= AU_KM/SOLARDAY;
      outvel.z *= AU_KM/SOLARDAY;
		 
      // Add back in the Solar coordinates at mjdstart
      outpos.x += Sunpos[planetfile_startpoint].x;
      outpos.y += Sunpos[planetfile_startpoint].y;
      outpos.z += Sunpos[planetfile_startpoint].z;
      outvel.x += Sunvel[planetfile_startpoint].x;
      outvel.y += Sunvel[planetfile_startpoint].y;
      outvel.z += Sunvel[planetfile_startpoint].z;
      
      // Perform n-body integration beginning with the state vectors
      // produced by Kepler2dyn
      integrate_orbit04LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, outpos, outvel, planetfile_startpoint, planetfile_endpoint, targMJD, targpos, targvel);
      for(j=planetfile_startpoint; j<=planetfile_endpoint; j++) {
	Earth2targ.x = targpos[j-planetfile_startpoint].x - Earthpos[j].x;
	Earth2targ.y = targpos[j-planetfile_startpoint].y - Earthpos[j].y;
	Earth2targ.z = targpos[j-planetfile_startpoint].z - Earthpos[j].z;
	Earth2targvel.x = targvel[j-planetfile_startpoint].x - Earthvel[j].x;
	Earth2targvel.y = targvel[j-planetfile_startpoint].y - Earthvel[j].y;
	Earth2targvel.z = targvel[j-planetfile_startpoint].z - Earthvel[j].z;
	Earth2sun.x = Sunpos[j].x - Earthpos[j].x;
	Earth2sun.y = Sunpos[j].y - Earthpos[j].y;
	Earth2sun.z = Sunpos[j].z - Earthpos[j].z;
	// Calculate geocentric ecliptic latitude and longitude, and the corresponding angular velocities.
	celedeproj01LD(Earth2targ, &targlambda1, &targbeta1);
	celedeproj01LD(Earth2sun, &sunlambda, &sunbeta);
	Earth2targ.x += deltat*SOLARDAY*Earth2targvel.x;
	Earth2targ.y += deltat*SOLARDAY*Earth2targvel.y;
	Earth2targ.z += deltat*SOLARDAY*Earth2targvel.z;
	celedeproj01LD(Earth2targ, &targlambda2, &targbeta2);
	lambdavel = (targlambda2 - targlambda1)/deltat;
	betavel = (targbeta2 - targbeta1)/deltat;
	angvelpt = point2d(lambdavel,betavel);
	targlambda1 -= sunlambda;
	while(targlambda1>=360.0L) targlambda1-=360.0L;
	while(targlambda1<0.0L) targlambda1+=360.0L;
	if(VERBOSE>=1) cout << "Asteroid " << oneorb.desig << " solar elongation " << targlambda1 << ", ecliptic latitude " << targbeta1 << "\n";
	dlatct = (targbeta1+double(latmax))/double(latstep) + 0.5l;
	latct = dlatct;
	delongct = double(targlambda1)/double(elongstep);
	elongct = delongct;
	if(VERBOSE>=1) cout << "angvel " << lambdavel << " " << betavel << " intcoords: " << latct << " " << elongct << " " << latct*elongnum + elongct << "\n"; 
	if(latct>=0 && latct<latnum && elongct>=0 && elongct<elongnum) {
	  velmat[latct*elongnum + elongct].push_back(angvelpt);
	}
      }
    }
  }
  instream1.close();
  for(i=0;i<latnum*elongnum;i++) {
    // Construct the name of the output file
    if(latcen[i]>=0) {
      int dnum=sprintf(nsfx,"%d",latcen[i]);
      numsuffix = nsfx;
      if(dnum<=1) numsuffix = "n0" + numsuffix;
      else numsuffix = "n" + numsuffix;
    } else if(latcen[i]<0) {
      int dnum=sprintf(nsfx,"%d",-latcen[i]);
      numsuffix = nsfx;
      if(dnum<=1) numsuffix = "s0" + numsuffix;
      else numsuffix = "s" + numsuffix;
    }
    outname = outfile + numsuffix + "_";
    int dnum=sprintf(nsfx,"%d",elongcen[i]);
    numsuffix = nsfx;
    if(dnum<=1) numsuffix = "00" + numsuffix;
    else if(dnum==2) numsuffix = "0" + numsuffix;
    outname = outname + numsuffix + ".txt";
    cout << "Writing output file " << i << " of " << velmat.size() << ", named " << outname << "\n";
    outstream1.open(outname);
    for(j=0;j<long(velmat[i].size());j++) {
      outstream1 << velmat[i][j].x << " " << velmat[i][j].y << "\n";
    }
    outstream1.close();
  }
  return(0);
}
