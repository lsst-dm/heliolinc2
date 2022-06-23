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

static void show_usage()
{
  cerr << "Usage: projectpairs04a -dets detfile -pairs pairfile -mjd mjdref -obspos observer_position_file -heliodist heliocentric_dist_vel_acc_file -clustrad cluster_radius -out outfile -outrms rmsfile \n";
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
  int status=0;
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> poscluster;
  vector <point3LD> velcluster;
  vector <point3LD> posclusterB;
  vector <point3LD> velclusterB;
  point3LD targvel1 = point3LD(0,0,0);
  point3LD targvel2 = point3LD(0,0,0);
  point3LD outpos = point3LD(0,0,0);
  point3LD unitbary = point3LD(0,0,0);
  point3LD targpos1 = point3LD(0,0,0);
  point3LD targpos2 = point3LD(0,0,0);
  vector <point3LD> observer_heliopos;
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
  vector <string> det_id_vec;
  long origind=0;
  vector <long> orig_index;
  long i1=0;
  long i2=0;
  point6LDx2 statevec1 = point6LDx2(0L,0L,0L,0L,0L,0L,0,0);
  point6ix2 stateveci = point6ix2(0,0,0,0,0,0,0,0);
  vector <point6ix2> allstatevecs;
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
  while(instream1 >> MJD >> RA >> Dec >> X >> Y >> Z >> detid >> origind) {
    o1=det_bsc(MJD,RA,Dec);
    detvec.push_back(o1);
    outpos = point3LD(X,Y,Z);
    observer_heliopos.push_back(outpos);
    det_id_vec.push_back(detid);
    orig_index.push_back(origind);
    if(MJD < minMJD) minMJD = MJD;
    if(MJD > maxMJD) maxMJD = MJD;
  }
  instream1.close();
  chartimescale = (maxMJD - minMJD)*SOLARDAY/TIMECONVSCALE; // Note that the units are seconds.
  
  cout << "Checking size-match of vectors read from " << indetfile << ":\n";
  cout << detvec.size() << " " << observer_heliopos.size() << " " << det_id_vec.size() << " " << orig_index.size() << "\n";
  
  // Read input image pair file
  instream1.open(inpairfile);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << inpairfile << "\n";
    return(1);
  }
  while(instream1 >> i1 >> i2) {    
    onepair = longpair(i1,i2);
    pairvec.push_back(onepair);
  }
  instream1.close();

  // Test: print out time-sorted detection table.
  ofstream outstream2 {outfile};
  ofstream outstream1 {rmsfile};
  cout << "Writing to " << outfile << "\n";
  cout << "Read " << detvec.size() << " detections and " << pairvec.size() << " pairs.\n";

  valid_accel_num=0;  
  for(accelct=0;accelct<accelnum;accelct++) {
    badpoint=0;
    // Calculate approximate heliocentric distances from the
    // input quadratic approximation.
    cout << "Working on grid point " << accelct << ": " << heliodist[accelct] << " " << heliovel[accelct]/SOLARDAY << " " << helioacc[accelct]*1000.0/SOLARDAY/SOLARDAY << "\n";
    outstream2 << fixed << setprecision(6) << "Working on grid point " << accelct << ": " << heliodist[accelct] << " " << heliovel[accelct]/SOLARDAY << " " << helioacc[accelct]*1000.0/SOLARDAY/SOLARDAY << "\n";

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
    for(i=0;i<pairvec.size();i++) {
      badpoint=0;
      //cout << "Working on pair " << i << " of " << pairvec.size() << "\n";
      // Obtain indices to the detection and heliocentric distance vectors.
      i1=pairvec[i].i1;
      i2=pairvec[i].i2;
      // Project the first point
      RA = detvec[i1].RA;
      Dec = detvec[i1].Dec;
      celestial_to_stateunitLD(RA,Dec,unitbary);
      status1 = helioproj01LD(unitbary, observer_heliopos[i1], heliodistvec[i1], delta1, targpos1);
      RA = detvec[i2].RA;
      Dec = detvec[i2].Dec;
      celestial_to_stateunitLD(RA,Dec,unitbary);
      status2 = helioproj01LD(unitbary, observer_heliopos[i2], heliodistvec[i2], delta2, targpos2);
      if(status1 == 0 && status2 == 0 && badpoint==0) {
	timediff = (detvec[i2].MJD - detvec[i1].MJD)*SOLARDAY;
	targvel1.x = (targpos2.x - targpos1.x)/timediff;
	targvel1.y = (targpos2.y - targpos1.y)/timediff;
	targvel1.z = (targpos2.z - targpos1.z)/timediff;
	//cout << "pos1, pos2, vel:\n";
	//cout << targpos1.x << " " << targpos1.y << " " << targpos1.z << "\n";
	//cout << targpos2.x << " " << targpos2.y << " " << targpos2.z << "\n";
	//cout << targvel1.x << " " << targvel1.y << " " << targvel1.z << "\n";
	//outstream1 << targpos1.x << " " << targpos1.y << " " << targpos1.z << "\n";
	//outstream1 << targpos2.x << " " << targpos2.y << " " << targpos2.z << "\n";

	targpos1.x = 0.5L*targpos2.x + 0.5L*targpos1.x;
	targpos1.y = 0.5L*targpos2.y + 0.5L*targpos1.y;
	targpos1.z = 0.5L*targpos2.z + 0.5L*targpos1.z;
    
	// Integrate orbit to the reference time.
	mjdavg = 0.5L*detvec[i1].MJD + 0.5L*detvec[i2].MJD;
	status1 = Keplerint(GMSUN_KM3_SEC2,mjdavg,targpos1,targvel1,mjdref,targpos2,targvel2);
	if(status1 == 0 && badpoint==0) {
	  statevec1 = point6LDx2(targpos2.x,targpos2.y,targpos2.z,chartimescale*targvel2.x,chartimescale*targvel2.y,chartimescale*targvel2.z,i1,i2);
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
    
    kdvec={};
    kdroot = splitpoint = 0;
    splitpoint=medind_6ix2(allstatevecs,1);
    kdpoint = KD_point6ix2(allstatevecs[splitpoint],-1,-1,1,0);
    kdvec.push_back(kdpoint);
    kdtree_6i01(allstatevecs,1,splitpoint,kdroot,kdvec);
    
    cout << "Calculated " << allstatevecs.size() << " state vectors from " << pairvec.size() << " pairs\n";
    cout << "Created a KD tree with " << kdvec.size() << " branches\n";

    outclusters={};
    int npt=3;
    int clusternum = DBSCAN_6i01(kdvec, cluster_radius/INTEGERIZING_SCALE, npt, INTEGERIZING_SCALE, outclusters);
    cout << "DBSCAN_6i01 finished, with " << clusternum << " = " << outclusters.size() << " clusters found\n";
    realclusternum=0;
    for(clusterct=0; clusterct<outclusters.size(); clusterct++) {
      // Calculate some cluster statistics.
      //cout << "Working on cluster " << clusterct << " with " << outclusters[clusterct].numpoints << " = " << outclusters[clusterct].clustind.size() << " points.\n";

      // Map cluster to individual detections.
      // create vector of unique detection indices.
      pointind={};
      for(i=0;i<outclusters[clusterct].numpoints;i++) {
	  pointind.push_back(kdvec[outclusters[clusterct].clustind[i]].point.i1);
	  pointind.push_back(kdvec[outclusters[clusterct].clustind[i]].point.i2);
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
	  if(i>0 && det_id_vec[pointind[i]] != det_id_vec[pointind[i-1]]) rating="MIXED";
	}
	// Write cluster to output file.
	outstream2 << "Accelct " << accelct << " Cluster " << realclusternum << " with " << outclusters[clusterct].numpoints << " = " << outclusters[clusterct].clustind.size() << " points, " << rating << ".\n";
	outstream2 << fixed << setprecision(6) << "Unique pts: " << pointind.size() << " span: " << timespan << " daysteps: " << numdaysteps << "\n";
	//cout << "Cluster pos RMS: " << outclusters[clusterct].rmsvec[0] << " " << outclusters[clusterct].rmsvec[1] << " " << outclusters[clusterct].rmsvec[2] <<  " total pos " << outclusters[clusterct].rmsvec[6] << "\n";
	//cout << "Cluster vel RMS: " << outclusters[clusterct].rmsvec[3] << " " << outclusters[clusterct].rmsvec[4] << " " << outclusters[clusterct].rmsvec[5] <<  " total vel " << outclusters[clusterct].rmsvec[7] << "\n";
	//cout << "Cluster total RMS: " << outclusters[clusterct].rmsvec[8] << "\n";
	outstream2 << fixed << setprecision(3) << "Cluster pos RMS: " << outclusters[clusterct].rmsvec[0] << " " << outclusters[clusterct].rmsvec[1] << " " << outclusters[clusterct].rmsvec[2] <<  " total pos " << outclusters[clusterct].rmsvec[6] << "\n";
	outstream2 << fixed << setprecision(3) << "Cluster vel RMS: " << outclusters[clusterct].rmsvec[3] << " " << outclusters[clusterct].rmsvec[4] << " " << outclusters[clusterct].rmsvec[5] <<  " total vel " << outclusters[clusterct].rmsvec[7] << "\n";
	outstream2 << "Cluster total RMS: " << outclusters[clusterct].rmsvec[8] << "\n";
	// Write individual detections to output file
	for(i=0; i<pointind.size(); i++) {
	  outstream2  << fixed << setprecision(6) << i << " " << detvec[pointind[i]].MJD << " " << detvec[pointind[i]].RA << " " << detvec[pointind[i]].Dec << " " << det_id_vec[pointind[i]] << " " << pointind[i] << " " << orig_index[pointind[i]] << "\n";
	}
	outstream2 << "\n";
	//cout << "\n";
	// Write summary line to rms file
	outstream1  << fixed << setprecision(6) << "Accelct " << accelct << " Cluster " << realclusternum << " : " << outclusters[clusterct].numpoints << " " << pointind.size() << " " << timespan << " " << numdaysteps << " " << outclusters[clusterct].rmsvec[6] << " " << outclusters[clusterct].rmsvec[7] << " " << outclusters[clusterct].rmsvec[8] << " " << rating << "\n";
      }
    }
  }
  return(0);
}
