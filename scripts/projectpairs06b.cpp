// March 11, 2022: projectpairs06b.cpp:
// Like projectpairs06a.cpp, but reads and keeps track of magnitude,
// photometric band, and observatory code for each observation.
//
// March 04, 2022: projectpairs06a.cpp:
// Like projectpairs05a.cpp, but uses a clustering radius that
// depends on geocentric distance. The user-defined clustering
// radius will now be taken to apply to a standard geocentric
// distance of 1.0 AU.
//
// January 24, 2022: projectpairs05a.cpp:
// Like projectpairs04c.cpp, but reads a more sophisticated pair
// file that aggregates tracklets into single effective pairs.
//
// January 07, 2022: projectpairs04c.cpp:
// Like projectpairs04b.cpp, but uses an integerized form of the
// state vectors for the k-d tree and DBSCAN implementations,
// for speed.
//
// January 06, 2022: projectpairs04b.cpp
// Like projectpairs04a.cpp, but streamlined a bunch of kludgy
// aspects of the DBSCAN implementation that were good for debugging
// but not appropriate for production.
//
// projectpairs04a.cpp
// Like projectpairs03a.cpp, but aimed at performing the further
// step of clustering state vectors propagated to the reference
// time in order link asteroid discoveries. In support of doing
// this very fast in cases where we have millions of detections,
// the observers heliocentric coordinates will be supplied in
// the input paired detection file, along with string identifiers
// indicating which detections correspond to the same unique object,
// which are required for testing purposes. NOTE: I first successfully
// implemented the DBSCAN algorithm for projectpairs04a.
//
// December 07, 2021: projectpairs03a.cpp
// Read pair files output by maketrack01a.cpp, and use a
// constant-acceleration model of the heliocentric distance.
// This is the first version that searches a grid
// of heliocentric distance and velocity. The grid is
// defined in a file whose name is supplied in the
// configuration file, which also supplies best-guess
// heliocentric accelerations for each heliocentric distance
// and radial velocity.
// The proceed as in projectpairs01a.cpp:
// Analysis performed:
// 1. Calculate the observer's precise barycentric
//    position at the moment of each observation.
// 2. Use an approximate distance for the target from
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

#define TIMECONVSCALE 4.0L // The characteristic timescale used to convert velocities
                           // to distance units is equal to the full temporal span
                           // divided by TIMECONVSCALE.
#define MINSPAN 1.0 // Temporal span must be at least this large (in days) for a bona fide cluster
#define MINDAYSTEPS 2 // A bona fide cluster must have at least this many intra-point
                      // time intervals greater than INTRANIGHTSTEP days.
#define INTRANIGHTSTEP 0.3 // Minimum interval in days between successive points
                           // in a tracklet, to enable them to be counted as being
                           // on separate nights.
#define INTEGERIZING_SCALE 100.0L // We divide state vectors by this value to integerize
                                   // them. Given a standard integer with a range
                                   // of +/- 2^31 = 2.15e9, this gives the state vectors
                                   // a range of +/- 2.15e11 km = 1400 AU, which is
                                   // comfortably larger than the solar system.
#define TYPECODE_PAIR "P"
#define TYPECODE_TRACKLET "T"

#define MINGEODIST 0.1 // Geocentric distance corresponding to the center of the
                       // smallest logarithmic bin
#define GEOBIN_HALF_WIDTH 1.5 // Logarithmic half-width for bins in geocentric distance.
                              // for example, if set to 2.0, a bin centered on 1.0 AU
                              // will extend from 0.5 to 2.0 AU
#define MAXGEODIST 100.0 // Upper limit on geocentric distance corresponding to
                         // the center of the largest logarithmic bin. The actual
                         // value is guaranteed to be at least half as large.
#define CRAD_REF_GEODIST 1.0 // Value of geocentric distance to which the user-defined
                             // clustering radius is normalized (AU). In general, the
                             // clustering radius is scaled linearly with geocentric distance.

static void show_usage()
{
  cerr << "Usage: projectpairs04a -dets detfile -pairs pairfile -mjd mjdref -obspos observer_position_file -heliodist heliocentric_dist_vel_acc_file -clustrad cluster_radius -out outfile -outrms rmsfile \n";
}
    
int main(int argc, char *argv[])
{
  det_obsmag_indvec o1 = det_obsmag_indvec(0L,0l,0l,0L,0L,0L,"null",0l,"V","I11",0,{});
  vector <det_obsmag_indvec> detvec = {};
  vector <long> pindexvec={};
  vector <vector <long>> pairvec={};
  double mag=0l;
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  string typecode;
  string configfile;
  string indetfile;
  string inpairfile;
  string accelfile;
  string outfile;
  string rmsfile;
  string lnfromfile;
  long double MJD,RA,Dec;
  int reachedeof=0;
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
  long ipt=0;
  int status=0;
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> poscluster;
  vector <point3LD> velcluster;
  vector <point3LD> posclusterB;
  vector <point3LD> velclusterB;
  point3LD Earthrefpos = point3LD(0,0,0);
  point3LD targvel1 = point3LD(0,0,0);
  point3LD targvel2 = point3LD(0,0,0);
  point3LD observerpos = point3LD(0,0,0);
  point3LD unitbary = point3LD(0,0,0);
  point3LD targpos1 = point3LD(0,0,0);
  point3LD targpos2 = point3LD(0,0,0);
  vector <long double> mjd;
  vector <long double> heliodist;
  vector <long double> heliovel;
  vector <long double> helioacc;
  long double mjdstart=0.0;
  long double mjdend=0.0;
  long double mjdavg=0.0;
  long runnum,runct;
  runnum = runct=1;
  vector <long double> heliodistvec;
  vector <long double> rmsvec;
  long double dist=0.0;
  long double pa=0.0;
  long double timediff=0.0;
  int status1=0; int status2=0;
  long double delta1,delta2;
  long double heliodistB=0.0L;
  long double heliovelB=0.0L;
  long double helioaccB=0.0L;
  long double wrms = LARGERR;
  long double wrmsB = LARGERR;
  int accelnum=0;
  int accelct=0;
  int valid_accel_num=0;
  int badpoint=0;
  long double X=0.0L;
  long double Y=0.0L;
  long double Z=0.0L;
  string detid;
  long origind=0;
  long i1=0;
  long i2=0;
  point6LDx2 statevec1 = point6LDx2(0L,0L,0L,0L,0L,0L,0,0);
  point6ix2 stateveci = point6ix2(0,0,0,0,0,0,0,0);
  vector <point6ix2> allstatevecs;
  vector <point6ix2> binstatevecs;
  double cluster_radius = 1e5l;
  KD_point6ix2 kdpoint = KD_point6ix2(stateveci,-1,-1,1,0);
  vector <KD_point6ix2> kdvec;
  long kdroot=0;
  long splitpoint=0;
  long double chartimescale = 1e5L;
  long double minMJD = 1e30L;
  long double maxMJD = -1e30L;
  vector <KD6i_clust> outclusters;
  vector <long> pointind;
  vector <long> pointjunk;
  vector <long double> clustmjd;
  vector <long double> mjdstep;
  long double timespan = 0.0L;
  int numdaysteps=0;
  int clusterct=0;
  int realclusternum=0;
  string rating;
  int trackpointnum=0;
  int trackpointct=0;
  long pairct=0;
  double mingeodist = MINGEODIST;
  double maxgeodist = MAXGEODIST;
  double geologstep = GEOBIN_HALF_WIDTH;
  int geobinct=0;
  double georadcen,georadmin,georadmax,geodist;
  georadcen = georadmin = georadmax = geodist = 0l;
  double clustmetric=0l;

  if(argc!=17)
    {
      show_usage();
      return(1);
    }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-d" || string(argv[i]) == "-dets" || string(argv[i]) == "-det" || string(argv[i]) == "--dets" || string(argv[i]) == "--det" || string(argv[i]) == "--detection" || string(argv[i]) == "--detections") {
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
    } else if(string(argv[i]) == "-obspos" || string(argv[i]) == "-op" || string(argv[i]) == "-obsvec" || string(argv[i]) == "--observer" || string(argv[i]) == "--observer_position" || string(argv[i]) == "--observer_statevec") {
      if(i+1 < argc) {
	//There is still something to read;
	planetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Observer position file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-heliodist" || string(argv[i]) == "-hd" || string(argv[i]) == "-heliodva" || string(argv[i]) == "-hdva" || string(argv[i]) == "--heliodistvelacc" || string(argv[i]) == "--heliodva" || string(argv[i]) == "--observer_statevec") {
      if(i+1 < argc) {
	//There is still something to read;
	accelfile=argv[++i];
	i++;
      }
      else {
	cerr << "Heliocentric distance, velocity, and acceleration\nfile keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clustrad" || string(argv[i]) == "-cr" || string(argv[i]) == "-crad" || string(argv[i]) == "-cluster" || string(argv[i]) == "--cluster_radius" || string(argv[i]) == "--clusterradius" || string(argv[i]) == "--clusterrad") {
      if(i+1 < argc) {
	//There is still something to read;
	cluster_radius=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Heliocentric distance, velocity, and acceleration\nfile keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-rms" || string(argv[i]) == "-outrms" || string(argv[i]) == "-or" || string(argv[i]) == "--outrmsfile" || string(argv[i]) == "--outrms" || string(argv[i]) == "--rmsfile") {
      if(i+1 < argc) {
	//There is still something to read;
	rmsfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output RMS file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
  
  cout.precision(17);  
  cout << "input detection file " << indetfile << "\n";
  cout << "input pair file " << inpairfile << "\n";
  cout << "input observer position file " << planetfile << "\n";
  cout << "input heliocentric distance, radial velocity, and radial acceleration file " << accelfile << "\n";
  cout << "input reference MJD " << mjdref << "\n";
  cout << "output file " << outfile << "\n";

  cout << "Heliocentric ephemeris for Earth is named " << planetfile << "\n";
  mjd={};
  Earthpos={};
  Earthvel={};
  read_horizons_fileLD(planetfile,mjd,Earthpos,Earthvel);
  cout << "Finished reading ephemeris file " << planetfile << "\n";
  // Calculate the position of Earth at the reference time.
  planetpos01LD(mjdref, 5, mjd, Earthpos, Earthrefpos);

  cout << "File with averaged accelerations is named " << accelfile << "\n";
  status=read_accel_fileLD(accelfile,heliodist,heliovel,helioacc);
  if(status!=0) {
    cerr << "ERROR: could not successfully read accleration file " << accelfile << "\n";
    cerr << "read_accel_fileLD returned status =" << status << ".\n";
   return(1);
  } else {
    accelnum = helioacc.size();
    if(heliovel.size() != accelnum || heliodist.size() != accelnum) {
      cerr << "ERROR: size mismatch " << accelnum << " " << heliovel.size() << " " << heliodist.size() << "in vectors read by read_accel_fileLD\n";
      return(1);
    }
    // Convert distances from AU to km, velocities from
    // AU/day to km/day, and accelerations from units of -GMSUN_KM3_SEC2/r^2
    // to units of km/day^2.
    for(accelct=0;accelct<accelnum;accelct++) {
      heliodist[accelct] *= AU_KM;
      heliovel[accelct] *= AU_KM;
      helioacc[accelct] *= -GMSUN_KM3_SEC2*SOLARDAY*SOLARDAY/heliodist[accelct]/heliodist[accelct];
    }
  }
  
  // Read input detection file.
  ifstream instream1 {indetfile};
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << indetfile << "\n";
    return(1);
  }
  while(instream1 >> MJD >> RA >> Dec >> X >> Y >> Z >> detid >> mag >> band >> obscode >> origind) {
    o1=det_obsmag_indvec(MJD,RA,Dec,X,Y,Z,detid,mag,band,obscode,origind,{});
    detvec.push_back(o1);
    if(MJD < minMJD) minMJD = MJD;
    if(MJD > maxMJD) maxMJD = MJD;
  }
  instream1.close();
  chartimescale = (maxMJD - minMJD)*SOLARDAY/TIMECONVSCALE; // Note that the units are seconds.
  
  cout << detvec.size() << " detection records read from " << indetfile << ".\n";
  
  // Read input image pair file
  pairvec={};
  instream1.open(inpairfile);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << inpairfile << "\n";
    return(1);
  }
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    typecode = ""; // Wipe previously read typecode
    // Read the current type code
    instream1 >> typecode;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad() && typecode == TYPECODE_PAIR) {
      // Read a single, isolated pair
      pindexvec={};
      instream1 >> i1 >> i2;
      pindexvec.push_back(i1);
      pindexvec.push_back(i2);
      pairvec.push_back(pindexvec);
    } else if(!instream1.eof() && !instream1.fail() && !instream1.bad() && typecode == TYPECODE_TRACKLET) {
      // Read the representative pair
      pindexvec={};
      instream1 >> i1 >> i2;
      pindexvec.push_back(i1);
      pindexvec.push_back(i2);
      // Read revised RA, Dec for representative detection 1.
      instream1 >> RA >> Dec;
      // Re-assign paired detection RA, Dec
      detvec[i1].RA = RA;
      detvec[i1].Dec = Dec;
      // Read revised RA, Dec for representative detection 2.
      instream1 >> RA >> Dec;
      // Re-assign paired detection RA, Dec
      detvec[i2].RA = RA;
      detvec[i2].Dec = Dec;
      // Read number of points in the tracklet
      instream1 >> trackpointnum;
      for(trackpointct=0; trackpointct<trackpointnum; trackpointct++) {
	instream1 >> ipt;
	if(ipt!=i1 && ipt!=i2) pindexvec.push_back(ipt);
      }
      pairvec.push_back(pindexvec);
    } else if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      cerr << "ERROR: unrecognized pair type code " << typecode << "\n";
      return(1);
    }
  }
  instream1.close();

  ofstream outstream2 {outfile};
  ofstream outstream1 {rmsfile};
  outstream2 << "#ptct,MJD,RA,Dec,idstring,mag,band,obscode,index1,index2,clusternum\n";
  outstream1 << "#clusternum,posRMS,velRMS,totRMS,pairnum,timespan,uniquepoints,daysteps,metric,rating,heliodist,heliovel,helioacc,posX,posY,posZ,velX,velY,velZ\n";
  cout << "Writing to " << outfile << "\n";
  cout << "Read " << detvec.size() << " detections and " << pairvec.size() << " pairs.\n";

  valid_accel_num=0;  
  realclusternum=0;
  for(accelct=0;accelct<accelnum;accelct++) {
    badpoint=0;
    // Calculate approximate heliocentric distances from the
    // input quadratic approximation.
    cout << "Working on grid point " << accelct << ": " << heliodist[accelct] << " " << heliovel[accelct]/SOLARDAY << " " << helioacc[accelct]*1000.0/SOLARDAY/SOLARDAY << "\n";

    heliodistvec={};
    //cout << "Approximate heliocentric distances:\n";
    for(i=0;i<detvec.size();i++)
      {
	delta1 = detvec[i].MJD - mjdref;
	heliodistvec.push_back(heliodist[accelct] + heliovel[accelct]*delta1 + 0.5*helioacc[accelct]*delta1*delta1);
	//cout << detvec[i].MJD << " " << delta1 << " " << heliodistvec[i] << "\n";
	if(heliodistvec[i]<=0.0L) {
	  badpoint=1;
	  break;
	}
      }
    if(badpoint==1) continue;
    if(badpoint==0 && heliodistvec.size()!=detvec.size()) {
      cerr << "ERROR: number of heliocentric distance values does\n";
      cerr << "not match the number of input detections!\n";
      return(1);
    }
    allstatevecs={};
    for(pairct=0; pairct<pairvec.size(); pairct++) {
      badpoint=0;
      //cout << "Working on pair " << i << " of " << pairvec.size() << "\n";
      // Obtain indices to the detection and heliocentric distance vectors.
      i1=pairvec[pairct][0];
      i2=pairvec[pairct][1];
      // Project the first point
      RA = detvec[i1].RA;
      Dec = detvec[i1].Dec;
      celestial_to_stateunitLD(RA,Dec,unitbary);
      observerpos = point3LD(detvec[i1].x,detvec[i1].y,detvec[i1].z);
      status1 = helioproj01LD(unitbary,observerpos, heliodistvec[i1], delta1, targpos1);
      RA = detvec[i2].RA;
      Dec = detvec[i2].Dec;
      celestial_to_stateunitLD(RA,Dec,unitbary);
      observerpos = point3LD(detvec[i2].x,detvec[i2].y,detvec[i2].z);
      status2 = helioproj01LD(unitbary, observerpos, heliodistvec[i2], delta2, targpos2);
      if(status1 == 0 && status2 == 0 && badpoint==0) {
	timediff = (detvec[i2].MJD - detvec[i1].MJD)*SOLARDAY;
	targvel1.x = (targpos2.x - targpos1.x)/timediff;
	targvel1.y = (targpos2.y - targpos1.y)/timediff;
	targvel1.z = (targpos2.z - targpos1.z)/timediff;

	targpos1.x = 0.5L*targpos2.x + 0.5L*targpos1.x;
	targpos1.y = 0.5L*targpos2.y + 0.5L*targpos1.y;
	targpos1.z = 0.5L*targpos2.z + 0.5L*targpos1.z;
      
	// Integrate orbit to the reference time.
	mjdavg = 0.5L*detvec[i1].MJD + 0.5L*detvec[i2].MJD;
	status1 = Keplerint(GMSUN_KM3_SEC2,mjdavg,targpos1,targvel1,mjdref,targpos2,targvel2);
	if(status1 == 0 && badpoint==0) {
	  statevec1 = point6LDx2(targpos2.x,targpos2.y,targpos2.z,chartimescale*targvel2.x,chartimescale*targvel2.y,chartimescale*targvel2.z,pairct,i1);
	  // Note that the multiplication by chartimescale converts velocities in km/sec
	  // to units of km, for apples-to-apples comparison with the positions.
	  stateveci = conv_6LD_to_6i(statevec1,INTEGERIZING_SCALE);
	  allstatevecs.push_back(stateveci);
	} else {
	  //badpoint=1;
	  //cout << "Kepler integration encountered unphysical situation\n";
	  continue;
	}
      } else {
	badpoint=1;
	//cout << "Heliocentric projection found no physical solution\n";
	continue;
      }
    }

    if(allstatevecs.size()<=1) continue; // No clusters possible, skip to the next step.
    cout << "Calculated " << allstatevecs.size() << " state vectors from " << pairvec.size() << " pairs\n";

    // Loop over geocentric bins, selecting the subset of state-vectors
    // in each bin, and running DBSCAN only on those, with clustering radius
    // adjusted accordingly.
    geobinct = 0;
    georadcen = mingeodist*intpowD(geologstep,geobinct);
    while(georadcen<=MAXGEODIST) {
      georadcen = mingeodist*intpowD(geologstep,geobinct);
      georadmin = georadcen/geologstep;
      georadmax = georadcen*geologstep;
      // Load new array of state vectors, limited to those in the current geocentric bin
      binstatevecs={};
      for(i=0;i<allstatevecs.size();i++)
	{
	  // Reverse integerization of the state vector.
	  // This is only possible to a crude approximation, of course.
	  statevec1 = conv_6i_to_6LD(allstatevecs[i],INTEGERIZING_SCALE);
	  // Calculate geocentric distance in AU
	  geodist = sqrt(DSQUARE(statevec1.x-Earthrefpos.x) + DSQUARE(statevec1.y-Earthrefpos.y) + DSQUARE(statevec1.z-Earthrefpos.z))/AU_KM;
	  if(geodist >= georadmin && geodist <= georadmax) {
	    // This state vector is in the geocentric radius bin we are currently considering.
	    // Add it to binstatevecs.
	    binstatevecs.push_back(allstatevecs[i]);
	  }
	}

      cout << "Found " << binstatevecs.size() << " state vectors in geocentric bin from " << georadmin << " to " << georadmax << " AU\n";
      if(binstatevecs.size()<=1) {
	geobinct++;
	continue; // No clusters possible, skip to the next step.
      }
      
      kdvec={};
      kdroot = splitpoint = 0;
      splitpoint=medind_6ix2(binstatevecs,1);
      kdpoint = KD_point6ix2(binstatevecs[splitpoint],-1,-1,1,0);
      kdvec.push_back(kdpoint);
      kdtree_6i01(binstatevecs,1,splitpoint,kdroot,kdvec);
    
      cout << "Created a KD tree with " << kdvec.size() << " branches\n";

      outclusters={};
      int npt=3;
      int clusternum = DBSCAN_6i01(kdvec, cluster_radius*(georadcen/CRAD_REF_GEODIST)/INTEGERIZING_SCALE, npt, INTEGERIZING_SCALE, outclusters);
      cout << "DBSCAN_6i01 finished, with " << clusternum << " = " << outclusters.size() << " clusters found\n";
      for(clusterct=0; clusterct<outclusters.size(); clusterct++) {
	// Scale cluster RMS down to reference geocentric distance
	for(i=0;i<9;i++) outclusters[clusterct].rmsvec[i] *= CRAD_REF_GEODIST/georadcen;
	// Note that RMS is scaled down for more distant clusters, to
	// avoid bias against them in post-processing.
	
	// Map cluster to individual detections.
	// create vector of unique detection indices.
	pointind={};
	for(i=0;i<outclusters[clusterct].numpoints;i++) {
	  pairct=kdvec[outclusters[clusterct].clustind[i]].point.i1;
	  for(j=0; j<pairvec[pairct].size(); j++) {
	    pointind.push_back(pairvec[pairct][j]);
	  }
	}
	// Sort vector of detection indices
	sort(pointind.begin(), pointind.end());
	// Cull out duplicate entries
	pointjunk = pointind;
	pointind={};
	pointind.push_back(pointjunk[0]);
	for(i=1; i<pointjunk.size(); i++) {
	  if(pointjunk[i]!=pointjunk[i-1]) pointind.push_back(pointjunk[i]);
	}

	// Load vector of detection MJD's
	clustmjd = {};
	for(i=0; i<pointind.size(); i++) {
	  clustmjd.push_back(detvec[pointind[i]].MJD);
	}
	// Sort vector of MJD's
	sort(clustmjd.begin(), clustmjd.end());
	timespan = clustmjd[clustmjd.size()-1] - clustmjd[0];
	// Load vector of MJD steps
	mjdstep={};
	for(i=1; i<clustmjd.size(); i++) {
	  mjdstep.push_back(clustmjd[i]-clustmjd[i-1]);
	}
	// Count steps large enough to suggest a daytime period between nights.
	numdaysteps=0;	
	for(i=0; i<mjdstep.size(); i++) {
	  if(mjdstep[i]>INTRANIGHTSTEP) numdaysteps++;
	}
	//cout << "Unique pts: " << pointind.size() << " span: " << timespan << " daysteps: " << numdaysteps << "\n";
	// Does cluster pass the criteria for a linked detection?
	if(timespan >= MINSPAN && numdaysteps >= MINDAYSTEPS) {
	  realclusternum++;
	  //cout << "Cluster passes discovery criteria: will be designated as cluster " << realclusternum << "\n";
	  // Check whether cluster is composed purely of detections from
	  // a single simulated object (i.e., would be a real discovery) or is a mixture
	  // of detections from two or more different simulated objects (i.e., spurious).
	  rating="PURE";
	  for(i=0; i<pointind.size(); i++) {
	    if(i>0 && stringnmatch01(detvec[pointind[i]].idstring,detvec[pointind[i-1]].idstring,SHORTSTRINGLEN)!=0) rating="MIXED";
	  }
	  // Write all individual detections in this cluster to the output cluster file
	  for(i=0; i<pointind.size(); i++) {
	    outstream2  << fixed << setprecision(6) << intzero01i(i,4) << "," << detvec[pointind[i]].MJD << "," << detvec[pointind[i]].RA << "," << detvec[pointind[i]].Dec << "," << detvec[pointind[i]].idstring << ",";
	    outstream2  << fixed << setprecision(3) << detvec[pointind[i]].mag << "," << detvec[pointind[i]].band << "," << detvec[pointind[i]].obscode << "," << pointind[i] << "," << detvec[pointind[i]].index << "," << realclusternum << "\n";
	  }
	  // Write summary line to rms file
	  clustmetric = double(pointind.size())*double(numdaysteps)*timespan/outclusters[clusterct].rmsvec[8];
	  outstream1  << fixed << setprecision(3) << realclusternum << "," << outclusters[clusterct].rmsvec[6] << "," << outclusters[clusterct].rmsvec[7] << "," << outclusters[clusterct].rmsvec[8] << "," << outclusters[clusterct].numpoints << ",";
	  outstream1  << fixed << setprecision(6) << timespan << "," << pointind.size() << "," << numdaysteps  << "," << clustmetric << "," << rating << ",";
	  outstream1  << fixed << setprecision(3) << heliodist[accelct]/AU_KM << "," << heliovel[accelct]/SOLARDAY << ",";
	  outstream1  << fixed << setprecision(6) << helioacc[accelct]*1000.0/SOLARDAY/SOLARDAY << ",";
	  outstream1  << fixed << setprecision(3)  << outclusters[clusterct].meanvec[0] << "," << outclusters[clusterct].meanvec[1] << "," << outclusters[clusterct].meanvec[2] << ",";
	  outstream1  << fixed << setprecision(6) << outclusters[clusterct].meanvec[3]/chartimescale << ","   << outclusters[clusterct].meanvec[4]/chartimescale << "," << outclusters[clusterct].meanvec[5]/chartimescale << "\n";
	}
      }
      // Move on to the next bin in geocentric distance
      geobinct++;
    }
  }
  return(0);
}