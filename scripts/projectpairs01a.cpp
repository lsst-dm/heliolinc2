// December 07, 2021: projectpairs01a.cpp
// Read pair files output by maketrack01a.cpp, and perform the basic
// analysis that will eventually be used by heliolinc3D, in an
// 'open book' way so that the characteristics of the resulting
// clusters can be investigated.
// Analysis performed:
// 1. Calculate the observer's precise barycentric
//    position at the moment of each observation.
// 2. Use an input (assumed accurate) distance for the target from
//    the Solar System barycenter at each observation to calculate
//    the target's barycentric coordinates from the observed RA, Dec.
// 3. Difference pairs of observations to obtain the target's
//    barycentric velocity.
// 4. Integrate the position and velocity information from each
//    pair of observations to calculate the object's location
//    at an input reference time. Hence, each pair of observations
//    is used to obtain an independent estimate of the object's
//    barycentric position and velocity at a fixed reference time
//    that is the same for all of the pairs.
// 5. The whole point of this exercise is to investigate the distribution
//    of these independently derived estimates in the 6-D space of
//    barycentric position and velocity. The tighter they cluster, the
//    easier the linking problem.

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

int iswhitespace(int c)
{
  if(c==' ' || c=='\t' || c=='\r' || c=='\n' || c=='\v' || c=='\f') return(1);
  else return(0);
}

int readconfigLD(ifstream &instream1, long double *ldval)
{
  string lnfromfile;
  string stest;
  int i=0;
  int c = '0';
  *ldval = 0L;
  getline(instream1,lnfromfile);
  if(instream1.eof()) {
    // EOF: should not happen because the file should adhere to the
    // specifications given in the calling function, and hence the
    // calling function should stop before it hits the end of the file.
    return(-1);
  } else if(instream1.fail()) {
    return(-2); // Worse problem.
  } else if(instream1.bad()) {
    return(-3); // Even worse.
  } else if(lnfromfile[0]=='#') {
    // It was a comment line. Not an error, but signal
    // the detection of a comment in the return.
    return(1);
  } else {
    // Apparently a valid line.
    c='0';
    i=0;
    while(i<lnfromfile.size() && !iswhitespace(c) && c!=EOF) {
      c=lnfromfile[i];
      if(!iswhitespace(c) && c!=EOF) stest.push_back(c);
      i++;
    }
    *ldval = stold(stest);
    return(0);
  }
}

int readconfigd(ifstream &instream1, double *dval)
{
  string lnfromfile;
  string stest;
  int i=0;
  int c = '0';
  *dval = 0L;
  getline(instream1,lnfromfile);
  if(instream1.eof()) {
    // EOF: should not happen because the file should adhere to the
    // specifications given in the calling function, and hence the
    // calling function should stop before it hits the end of the file.
    return(-1);
  } else if(instream1.fail()) {
    return(-2); // Worse problem.
  } else if(instream1.bad()) {
    return(-3); // Even worse.
  } else if(lnfromfile[0]=='#') {
    // It was a comment line. Not an error, but signal
    // the detection of a comment in the return.
    return(1);
  } else {
    // Apparently a valid line.
    c='0';
    i=0;
    while(i<lnfromfile.size() && !iswhitespace(c) && c!=EOF) {
      c=lnfromfile[i];
      if(!iswhitespace(c) && c!=EOF) stest.push_back(c);
      i++;
    }
    *dval = stod(stest);
    return(0);
  }
}

int readconfigint(ifstream &instream1, int *ival)
{
  string lnfromfile;
  string stest;
  int i=0;
  int c = '0';
  *ival = 0;
  getline(instream1,lnfromfile);
  if(instream1.eof()) {
    // EOF: should not happen because the file should adhere to the
    // specifications given in the calling function, and hence the
    // calling function should stop before it hits the end of the file.
    return(-1);
  } else if(instream1.fail()) {
    return(-2); // Worse problem.
  } else if(instream1.bad()) {
    return(-3); // Even worse.
  } else if(lnfromfile[0]=='#') {
    // It was a comment line. Not an error, but signal
    // the detection of a comment in the return.
    return(1);
  } else {
    // Apparently a valid line.
    c='0';
    i=0;
    while(i<lnfromfile.size() && !iswhitespace(c) && c!=EOF) {
      c=lnfromfile[i];
      if(!iswhitespace(c) && c!=EOF) stest.push_back(c);
      i++;
    }
    *ival = stoi(stest);
    return(0);
  }
}

int readconfigstring(ifstream &instream1, string &sval)
{
  string lnfromfile;
  string stest;
  int i=0;
  int c = '0';
  sval={};
  getline(instream1,lnfromfile);
  if(instream1.eof()) {
    // EOF: should not happen because the file should adhere to the
    // specifications given in the calling function, and hence the
    // calling function should stop before it hits the end of the file.
    return(-1);
  } else if(instream1.fail()) {
    return(-2); // Worse problem.
  } else if(instream1.bad()) {
    return(-3); // Even worse.
  } else if(lnfromfile[0]=='#') {
    // It was a comment line. Not an error, but signal
    // the detection of a comment in the return.
    return(1);
  } else {
    // Apparently a valid line.
    c='0';
    i=0;
    while(i<lnfromfile.size() && !iswhitespace(c) && c!=EOF) {
      c=lnfromfile[i];
      if(!iswhitespace(c) && c!=EOF) stest.push_back(c);
      i++;
    }
    sval = stest;
    return(0);
  }
}

static void show_usage()
{
  cerr << "Usage: projectpairs01a -cfg configfile -dets detfile -pairs pairfile -dist heliodistfile -mjd mjdref -out outfile \n";
}
    
int main(int argc, char *argv[])
{
  det_bsc o1 = det_bsc(0,0,0);
  vector <det_bsc> detvec = {};
  longpair onepair = longpair(0,0);
  vector <longpair> pairvec ={};
  string configfile;
  string indetfile;
  string inpairfile;
  string heliodistfile;
  string outfile;
  string lnfromfile;
  long double MJD,RA,Dec;
  int reachedeof=0;
  long double heliodist=0L;
  long double mjdref=0L;
  long double ldval = 0L;
  string stest;
  long double obslon, plxcos, plxsin;
  obslon = plxcos = plxsin = 0L;
  string planetfile;
  int planetnum=0;
  int planetct=0;
  int polyorder=1;
  int i=0;
  int j=0;
  int c='0';
  int status=0;
  int pctEarth=0;
  vector <long double> planetmasses;
  vector <point3LD> planetpos;
  vector <point3LD> planetpos_reverse;
  vector <point3LD> Earthpos;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  point3LD targpos1 = point3LD(0,0,0);
  point3LD targpos2 = point3LD(0,0,0);
  point3LD targvel1 = point3LD(0,0,0);
  point3LD targvel2 = point3LD(0,0,0);
  point3LD targvel_reverse = point3LD(0,0,0);
  point3LD outpos = point3LD(0,0,0);
  point3LD unitbary = point3LD(0,0,0);
  point3LD barypos1 = point3LD(0,0,0);
  point3LD barypos2 = point3LD(0,0,0);
  point3LD baryvel = point3LD(0,0,0);
  vector <point3LD> observer_barypos;
  vector <long double> mjd;
  vector <long double> mjd_reverse;
  vector <long double> mjdtest;
  long double mjdstart=0.0;
  long double mjdend=0.0;
  long double mjdavg=0.0;
  long runnum,runct;
  runnum = runct=1;
  vector <long double> heliodistvec;
  long double dist=0.0;
  long double pa=0.0;
  long double timediff=0.0;
  int i1=0; int i2=0;
  long double delta1,delta2;

  if(argc!=13)
    {
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
      }
      else {
	cerr << "Configuration file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-d" || string(argv[i]) == "-dets" || string(argv[i]) == "-det" || string(argv[i]) == "--dets" || string(argv[i]) == "--det" || string(argv[i]) == "--detection" || string(argv[i]) == "--detections") {
      if(i+1 < argc) {
	//There is still something to read;
	indetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-p" || string(argv[i]) == "-pair" || string(argv[i]) == "-pairs" || string(argv[i]) == "--pairs" || string(argv[i]) == "--pair" || string(argv[i]) == "--pairfile" || string(argv[i]) == "--pairsfile") {
      if(i+1 < argc) {
	//There is still something to read;
	inpairfile=argv[++i];
	i++;
      }
      else {
	cerr << "Pair file keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-dist" || string(argv[i]) == "-hd" || string(argv[i]) == "-helio" || string(argv[i]) == "--heliodist" || string(argv[i]) == "--dist" || string(argv[i]) == "--helio" || string(argv[i]) == "--hdist") {
      if(i+1 < argc) {
	//There is still something to read;
	heliodistfile=argv[++i];
	i++;
      }
      else {
	cerr << "Heliocentric distance keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-m" || string(argv[i]) == "-mjd" || string(argv[i]) == "-MJD" || string(argv[i]) == "--mjd" || string(argv[i]) == "--MJD" || string(argv[i]) == "--modifiedjulianday" || string(argv[i]) == "--ModifiedJulianDay") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdref=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Pair file keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outpair" || string(argv[i]) == "--outpairs") {
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
  cout << "input detection file " << indetfile << "\n";
  cout << "input pair file " << inpairfile << "\n";
  cout << "input heliocentric distance file " << heliodistfile << "\n";
  cout << "input reference MJD " << mjdref << "\n";
  cout << "output file " << outfile << "\n";

  // Read configuration file.
  ifstream instream1 {configfile};
  if(!instream1) {
    cerr << "ERROR: can't open input config file " << configfile << "\n";
    return(1);
  }
  // Read observatory longitude
  status=readconfigLD(instream1,&obslon);
  while(status==1) {
    // The line we have just read is a pure comment line,
    // so we just want to skip to the next one.
    status=readconfigLD(instream1,&obslon);
  }
  if(status<0) {
    cerr << "Error reading config file\n";
    return(1);
  } else cout << "Observatory longitude read as " << obslon << "\n";
  // Read parallax cosine for observatory latitude.
  status=readconfigLD(instream1,&plxcos);
  while(status==1) {
    // The line we have just read is a pure comment line,
    // so we just want to skip to the next one.
    status=readconfigLD(instream1,&plxcos);
  }
  if(status<0) {
    cerr << "Error reading config file\n";
    return(1);
  } else cout << "Parallax cosine for observatory latitude read as " << plxcos << "\n";
  // Read parallax sine for observatory latitude.
  status=readconfigLD(instream1,&plxsin);
  while(status==1) {
    // The line we have just read is a pure comment line,
    // so we just want to skip to the next one.
    status=readconfigLD(instream1,&plxsin);
  }
  if(status<0) {
    cerr << "Error reading config file\n";
    return(1);
  } else cout << "Parallax sine for observatory latitude read as " << plxsin << "\n";  
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
    } else cout << "MG for planet " << planetct << " read as " << planetmasses[planetct] << " km^3/sec^2\n";
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
    if(planetct==0) mjd=mjdtest;
    else {
      for(i=0;i<mjd.size();i++) {
	if(mjdtest[i]!=mjd[i]) {
	  cout << "ERROR: time vectors do not match for input planet files\n";
	  cout << planetct+1 << " and 1!\n";
	  return(1);
	}
      }
    }
    for(i=0;i<temppos.size();i++) {
      planetpos.push_back(temppos[i]);
    }
    if(planetct == pctEarth) Earthpos = temppos;
    cout << "Finished reading ephemeris file " << planetfile << "\n";
  }
  // Close input stream that was reading the config file.
  instream1.close();
  
  // Read input detection file.
  instream1.open(indetfile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << indetfile << "\n";
    return(1);
  }
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    i=0;
    j = 0;
    c='0';
    while(i<lnfromfile.size() && reachedeof == 0) {
      string stest;
      c='0';
      while(i<lnfromfile.size() && c!=' ' && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=' ' && c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==1) MJD=stod(stest);
      else if(j==2) RA=stod(stest);
      else if(j==3) Dec=stod(stest);
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0) {
      // cout<<"MJD, RA, Dec: " << MJD-floor(MJD) << " " << RA << " " << Dec << "\n";
      o1=det_bsc(MJD,RA,Dec);
      detvec.push_back(o1);
    }
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  instream1.close();

  // Read input image pair file
  instream1.open(inpairfile);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << inpairfile << "\n";
    return(1);
  }
  reachedeof=0;
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    i=0;
    j = 0;
    c='0';
    long i1=0;
    long i2=0;
    while(i<lnfromfile.size() && reachedeof == 0) {
      string stest;
      c='0';
      while(i<lnfromfile.size() && c!=' ' && c!=',' && c!='\n' && c!=EOF) {
	// We allow the file to be delimited by comma or space.
	c=lnfromfile[i];
	if(c!=',' && c!=' ' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==1) i1=stol(stest);
      else if(j==2) i2=stol(stest);
    }
    if((reachedeof == 0 || reachedeof == 1) && i2>0) {
      onepair = longpair(i1,i2);
      pairvec.push_back(onepair);
    }
  }
  if(reachedeof==1) {
      cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  instream1.close();
  
  // Read heliocentric distance file. For now, this is a hack
  // in which the distances supplied are actually barycentric
  // distances, and they are supplied for every detection
  // individually. In heliolinc3D proper, they will be based
  // on hypothesis heliocentric distances and velocities over
  // which we will scan with fine sampling.
  
  instream1.open(heliodistfile);
  while(instream1 >> heliodist) heliodistvec.push_back(heliodist);
  instream1.close();

  for(i=0;i<heliodistvec.size();i++)
    {
      cout << "Barycentric distance = " << heliodistvec[i] << "\n";
    }
  if(heliodistvec.size()!=detvec.size()) {
    cerr << "ERROR: number of heliocentric distance values does\n";
    cerr << "not match the number of input detections!\n";
    return(1);
  }

  // Calculate the exact position of the observer at the time of each observation.
  for(i=0;i<detvec.size();i++)
    {
      observer_barycoords01LD(detvec[i].MJD, 5, obslon, plxcos, plxsin, mjd, Earthpos, outpos);
      observer_barypos.push_back(outpos);
    }
  
  // Test: print out time-sorted detection table.
  ofstream outstream1 {"testinitpos01.txt"};
  ofstream outstream2 {outfile};
  outstream1.precision(17);
  outstream2.precision(17);
  cout << "Writing to " << outfile << "\n";
  cout << "Read " << detvec.size() << " detections and " << pairvec.size() << "pairs.\n";

  // Create time-reversed copies of the planet position
  // and MJD vectors, for use in back-in-time integrations/
  // Reverse the order and the sign of the mjd vector.
  mjd_reverse = mjd;
  for(i=0;i<mjd.size();i++) mjd_reverse[mjd.size()-i-1] = -mjd[i];
  // Reverse the order of the planet vectors
  planetpos_reverse = planetpos;
  for(planetct=0; planetct<planetnum; planetct++) {
    for(i=0;i<mjd.size();i++) {
	temppos[mjd.size()-i-1] = planetpos[planetct*mjd.size() + i];
    }
    for(i=0;i<mjd.size();i++) {
      planetpos_reverse[planetct*mjd.size() + i] = temppos[i];
    }
  }
    // Reverse the direction of the input velocity
    targvel1.x = -targvel1.x;
    targvel1.y = -targvel1.y;
    targvel1.z = -targvel1.z;


  
  for(i=0;i<pairvec.size();i++) {
    // Obtain indices to the detection and heliocentric distance vectors.
    i1=pairvec[i].i1;
    i2=pairvec[i].i2;
    // Project the first point
    RA = detvec[i1].RA;
    Dec = detvec[i1].Dec;
    celestial_to_stateunitLD(RA,Dec,unitbary);
    helioproj01LD(unitbary, observer_barypos[i1], heliodistvec[i1], delta1, barypos1);
    RA = detvec[i2].RA;
    Dec = detvec[i2].Dec;
    celestial_to_stateunitLD(RA,Dec,unitbary);
    helioproj01LD(unitbary, observer_barypos[i2], heliodistvec[i2], delta2, barypos2);
    timediff = (detvec[i2].MJD - detvec[i1].MJD)*SOLARDAY;
    baryvel.x = (barypos2.x - barypos1.x)/timediff;
    baryvel.y = (barypos2.y - barypos1.y)/timediff;
    baryvel.z = (barypos2.z - barypos1.z)/timediff;    
    cout << "pos1, pos2, vel:\n";
    cout << barypos1.x << " " << barypos1.y << " " << barypos1.z << "\n";
    cout << barypos2.x << " " << barypos2.y << " " << barypos2.z << "\n";
    cout << baryvel.x << " " << baryvel.y << " " << baryvel.z << "\n";
    outstream1 << barypos1.x << " " << barypos1.y << " " << barypos1.z << "\n";
    outstream1 << barypos2.x << " " << barypos2.y << " " << barypos2.z << "\n";

    targpos1.x = 0.5L*barypos2.x + 0.5L*barypos1.x;
    targpos1.y = 0.5L*barypos2.y + 0.5L*barypos1.y;
    targpos1.z = 0.5L*barypos2.z + 0.5L*barypos1.z;
    
    // Integrate orbit to the reference time.
    mjdavg = 0.5L*detvec[i1].MJD + 0.5L*detvec[i2].MJD;
    if(mjdref>=mjdavg) {
      // We are integrating forward in time.
      // We use the true velocity.
      targvel1.x = (barypos2.x - barypos1.x)/timediff;
      targvel1.y = (barypos2.y - barypos1.y)/timediff;
      targvel1.z = (barypos2.z - barypos1.z)/timediff;
      //cout << "Integrating from MJD " << mjdavg << " to " << mjdref << "\n";
      //cout << "Starting position " << targpos1.x << " " << targpos1.y << " " << targpos1.z << "\n";
      //cout << "Starting velocity " << targvel1.x << " " << targvel1.y << " " << targvel1.z << "\n";
      integrate_orbit02LD(polyorder,planetnum, mjd, planetmasses, planetpos, mjdavg, targpos1, targvel1, mjdref, targpos2, targvel2);
      //cout << "Output position " << targpos2.x << " " << targpos2.y << " " << targpos2.z << "\n";
      //cout << "Output velocity " << targvel2.x << " " << targvel2.y << " " << targvel2.z << "\n";
      
    } else {
      // We are integrating backward in time.
      // We use the negative of the true velocity.
      targvel1.x = (barypos1.x - barypos2.x)/timediff;
      targvel1.y = (barypos1.y - barypos2.y)/timediff;
      targvel1.z = (barypos1.z - barypos2.z)/timediff;    
      //cout << "Integrating from MJD " << -mjdavg << " to " << -mjdref << "\n";
      //cout << "Starting position " << targpos1.x << " " << targpos1.y << " " << targpos1.z << "\n";
      //cout << "Starting velocity " << targvel1.x << " " << targvel1.y << " " << targvel1.z << "\n";
      integrate_orbit02LD(polyorder,planetnum, mjd_reverse, planetmasses, planetpos_reverse, -mjdavg, targpos1, targvel1, -mjdref, targpos2, targvel2);
      // Reverse the output velocity to align with the true velocity.
      targvel2.x = -targvel2.x;
      targvel2.y = -targvel2.y;
      targvel2.z = -targvel2.z;
      //cout << "Output position " << targpos2.x << " " << targpos2.y << " " << targpos2.z << "\n";
      //cout << "Output velocity " << targvel2.x << " " << targvel2.y << " " << targvel2.z << "\n";
    }
    outstream2 << targpos2.x << " " << targpos2.y << " " << targpos2.z << " ";
    outstream2 << targvel2.x << " " << targvel2.y << " " << targvel2.z << "\n";
  }
  
  return(0);
}
