// September 14, 2022: heliolinc.cpp (created on this date by copying
// an earlier prototype file called projectpairs06d.cpp).
//
// Implements in C++ (with some modifications) the Heliolinc3D
// algorithm developed by Siegfried Eggl, which in turn was based
// on the original HelioLinC algorithm developed by Matthew Holman
// (Holman et al. 2018, Astronomical Journal, 156, 135).
//
// The asteroid linking problem consists in identifying
// subsets of detections in an input catalog of astronomical
// transient detections that are in fact repeated detections of the
// same previously unknown asteroid. For example, by repeatedly
// observing a 1000 square-degree area of sky over a two-week
// period, an astronomical survey might identify millions of
// cases where an apparent source is detected at a position
// where there is no known star or asteroid. Each detection
// is uniquely identified by its celestial coordinates (e.g.,
// right ascension and declination, or RA and Dec) and the time
// of observation (e.g., modified Julian day, or MJD). These
// detections could be spurious (cosmic ray hits, electronic
// crosstalk, fluctuations of Poisson noise); stationary
// transients (e.g., supernovae or previously unknown
// variable stars); or asteroids. Although thousands of unknown
// asteroids may each have been detected many times, none can
// be said to have been 'discovered' until the linking problem
// has been solved and it is known which sets of detections are
// actually repeated measurements of the same unknown asteroid.
// Then the asteroid has really been discovered; its orbit can
// be calculated; future followup observations can be planned; and
// serendipitous past detections ('precoveries') may be identified.
//
// Heliolinc3D solves the asteroid linking problem by probing a
// number of hypotheses for the time-dependent heliocentric distance
// of undiscovered asteroids that may be present in the input data.
// It requires five input files: one is the detection catalog,
// one is a 'pair' file (see below), and one specifies the hypotheses
// to be probed. The remaining two files are reference files that
// will not typically need to be changed from one execution to another:
// one is a JPL Horizons ephemeris for the Earth, and the other is
// an observatory code file based on the one available from the Minor
// Planet Center at https://minorplanetcenter.net/iau/lists/ObsCodesF.html
// For digestion by the current version of heliolinc, the incomplete
// lines that relate to space-based observatories in this file must
// be deleted. Useful defaults for these last two files,
// Earth1day2020s_02a.txt and ObsCodes.txt, are available in the
// same github repository as this source code.
//
// The input detection catalog and preliminary 'pair' file required
// by heliolinc can be produced by using make_tracklets, available in the
// same github repository as the current code. The input to make_tracklets
// is a .csv detection catalog of very flexible formatting, which is
// described in more detail in the documentation for make_tracklets.
// The pair file lists pairs of entries in the detection catalog that
// might be detections of the same object within a relatively short
// span of time (e.g., 1.5 hours). These 'pairs' can also be 'tracklets'
// consisting of more than two detections, if the input catalog
// includes cases where more than two observations obtained close
// together in time lie along a consistent linear trajectory.
//
// The input file providing the hypotheses that should be probed, where
// each hypothesis proposes a particular dependence of unknown asteroids'
// heliocentric distance as a function of time, has one line per
// hypothesis, and four columns, as follows:
//
// 1. Heliocentric distance (AU)
// 2. Heliocentric radial velocity (AU/day)
// 3. Normalization (0: skip this line. Positive: don't skip).
// 4. 2nd time-derivative of heliocentric distance, scaled by -GMsun/r^2
//
// Columns 1, 2, and 4 effectively approximate the time-dependent
// heliocentric distance as a Taylor Series, providing the 0th,
// 1st, and 2nd time-derivatives of the heliocentric distance.
// The forth column can be thought of as an acceleration, but it's
// important to note that it is **not** the vector acceleration
// away from the sun, which always -GMsun/r^2 for any orbit. Instead,
// it's the second derivative of the distance from the sun, which is
// zero for a circular orbit, and reaches -GMsun/r^2 (1.0 in the
// normalized units used by heliolinc), only in the limit of an orbit
// with zero angular momentum (i.e., with eccentricty exactly 1.0)
//
//
//
//
// DEVELOPMENT NOTES IN REVERSE CHRONOLOGICAL ORDER:
//
// June 22, 2022: besides the viewing geometries
// handled by all earlier versions, also handles the previously-ignored
// geometry, possible only at solar elongation less than 90 degrees,
// where the line-of-sight from the observer intersects the
// the heliocentric sphere defined by the asteroid distance hypothesis
// on the near side of the sun -- that is, the asteroid is at the
// point where the line-of-sight vector pierces the heliocentric
// sphere from the **outside**. An alternative specification for
// this geometry is that the Sun-object-observer phase angle is
// **greater** than 90 degrees. It corresponds to a second solution
// of the quadratic equation for the intersection of a line with
// a sphere. This second solution was ignored in previous versions,
// but the current program handles **both** solutions correctly.
// Hence, it is not a supplement but a replacement for previous
// versions.
//
// March 15, 2022: reads the pairdets file in csv format
// with a header. Also replaces the
// daysteps parameter with obsnights = daysteps+1 (on output).
// The parameter daysteps is still used internally exactly as before.
//
// Note well that the format of the pair file (as opposed to the
// pairdets file) has NOT changed. The pair file has an unavoidably
// specialized format because of the need to handle both pairs and
// tracklets. This specialized format does not lend itself to the
// header+csv convention. Also, the pair file is not needed by any other
// program except projectpairs itself, so its specialized format
// is not expected to create compatibility issues for programs external
// to the heliolinc suite.
// 
// March 11, 2022: reads and keeps track of magnitude,
// photometric band, and observatory code for each observation.
//
// March 04, 2022: uses a clustering radius that
// depends on geocentric distance. The user-defined clustering
// radius will now be taken to apply to a standard geocentric
// distance of 1.0 AU.
//
// January 24, 2022: reads a more sophisticated pair
// file that aggregates tracklets into single effective pairs.
//
// January 07, 2022: uses an integerized form of the
// state vectors for the k-d tree and DBSCAN implementations,
// for speed.
//
// January 06, 2022: streamlined a bunch of kludgy
// aspects of the DBSCAN implementation that were good for debugging
// but not appropriate for production.
//
// January 04, 2022: aimed at performing the further
// step of clustering state vectors propagated to the reference
// time in order link asteroid discoveries. In support of doing
// this very fast in cases where we have millions of detections,
// the observers heliocentric coordinates will be supplied in
// the input paired detection file, along with string identifiers
// indicating which detections correspond to the same unique object,
// which are required for testing purposes. NOTE: I first successfully
// implemented the DBSCAN algorithm at this stage of developement
//
// December 07, 2021: Read pair files output by maketrack01a.cpp, and use a
// constant-acceleration model of the heliocentric distance.
// This is the first version that searches a grid
// of heliocentric distance and velocity. The grid is
// defined in a file whose name is supplied in the
// configuration file, which also supplies best-guess
// heliocentric accelerations for each heliocentric distance
// and radial velocity.
//
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
#define DBSCAN_NPT 3  // Default value of npt (min cluster size) for DBSCAN algorithm.

#define DEBUG_A 0

#define MINGEOOBS 0.0l // Geocentric distance (AU) at time of observation that causes objects
                       // to be considered for rejection in order to avoid globs.
#define MINIMPACTPAR 0.0l // Impact parameter with respect to Earth (km) that causes objects
                          // to be rejected, if also within MINGEOOBS, in order to avoid globs.

static void show_usage()
{
  cerr << "Usage: heliolinc -dets detfile -pairs pairfile -mjd mjdref -obspos observer_position_file -heliodist heliocentric_dist_vel_acc_file -clustrad cluster_radius -npt dbscan_npt -minobsnights minobsnights -mintimespan mintimespan -mingeodist minium_geocentric_distance -maxgeodist maximum_geocentric_distance -geologstep logarithmic_step_size_for_geocentric_distance_bins -out outfile -outsum summary_file \n";
  cerr << "\nor, at minimum:\n\n";
  cerr << "heliolinc -dets detfile -pairs pairfile -mjd mjdref -obspos observer_position_file -heliodist heliocentric_dist_vel_acc_file\n";
  cerr << "\nNote that the minimum invocation leaves some things set to defaults\n";
  cerr << "that you may well wish to specify: in particular, the output file names\n";
  
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
  string outfile="hlout.all";
  string rmsfile="hlout.summary";
  string lnfromfile;
  long double MJD,RA,Dec;
  long double mjdref=0L;
  string stest;
  string planetfile;
  long i=0;
  long j=0;
  long ipt=0;
  int status=0;
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> poscluster;
  vector <point3LD> velcluster;
  vector <point3LD> posclusterB;
  vector <point3LD> velclusterB;
  vector <point3LD> targposvec1;
  vector <point3LD> targposvec2;
  point3LD Earthrefpos = point3LD(0,0,0);
  point3LD targvel1 = point3LD(0,0,0);
  point3LD targvel2 = point3LD(0,0,0);
  point3LD observerpos1 = point3LD(0,0,0);
  point3LD observerpos2 = point3LD(0,0,0);
  point3LD unitbary = point3LD(0,0,0);
  point3LD targpos1 = point3LD(0,0,0);
  point3LD targpos2 = point3LD(0,0,0);
  vector <long double> mjd;
  vector <long double> heliodist;
  vector <long double> heliovel;
  vector <long double> helioacc;
  long double mjdavg=0.0;
  vector <long double> heliodistvec;
  vector <long double> rmsvec;
  long double timediff=0.0;
  int status1=0; int status2=0;
  long double delta1;
  vector <long double> deltavec1;
  vector <long double> deltavec2;
  int accelnum=0;
  int num_dist_solutions=0;
  int solnct=0;
  int accelct=0;
  int badpoint=0;
  long double X=0.0L;
  long double Y=0.0L;
  long double Z=0.0L;
  char detid[SHORTSTRINGLEN];
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
  int numobsnights=0;
  int clusterct=0;
  int realclusternum=0;
  int gridpoint_clusternum=0;
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
  int detfilelinect=0;
  int badread=0;
  int startpoint=0;
  int endpoint=0;
  int npt = DBSCAN_NPT;
  int mindaysteps = MINDAYSTEPS;
  int minobsnights = MINDAYSTEPS+1;
  double mintimespan = MINSPAN;
  ifstream instream1;
  vector <long double> mjdvec;
  int default_cluster_radius, default_npt, default_mindaysteps;
  int default_mintimespan, default_mingeodist, default_maxgeodist;
  int default_geologstep,default_outfile,default_rmsfile;
  int default_mingeoobs, default_minimpactpar;
  default_cluster_radius = default_npt = default_mindaysteps = default_mintimespan = 1;
  default_mingeodist = default_maxgeodist = default_geologstep = default_outfile = default_rmsfile = 1;
  default_mingeoobs = default_minimpactpar = 1;
  int verbose=0;
  int glob_warning=0;
  long double impactpar=0.0L;
  long double absvelocity=0.0L;
  double mingeoobs = MINGEOOBS; // Geocentric distance (AU) at time of observation that causes objects
                                // to be considered for rejection in order to avoid globs.
  double minimpactpar =  MINIMPACTPAR; // Impact parameter with respect to Earth (Earth radii) that causes objects
                                       // to be rejected, if also within MINGEOOBS, in order to avoid globs.

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
    } else if(string(argv[i]) == "-m" || string(argv[i]) == "-mjd" || string(argv[i]) == "-mjdref" || string(argv[i]) == "-MJD" || string(argv[i]) == "--mjd" || string(argv[i]) == "--MJD" || string(argv[i]) == "--modifiedjulianday" || string(argv[i]) == "--ModifiedJulianDay") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdref=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Reference MJD keyword supplied with no corresponding argument";
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
    } else if(string(argv[i]) == "-heliodist" || string(argv[i]) == "-hd" || string(argv[i]) == "-heliodva" || string(argv[i]) == "-hdva" || string(argv[i]) == "--heliodistvelacc" || string(argv[i]) == "--heliodva") {
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
	default_cluster_radius = 0;
	i++;
      }
      else {
	cerr << "Clustering radius keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-npt" || string(argv[i]) == "-npoints" || string(argv[i]) == "-minpts" || string(argv[i]) == "-np" || string(argv[i]) == "--npt" || string(argv[i]) == "--dbscan_npt" || string(argv[i]) == "--DBSCANnpt") {
      if(i+1 < argc) {
	//There is still something to read;
	npt=stoi(argv[++i]);
	default_npt = 0;
	i++;
      }
      else {
	cerr << "DBSCAN npt keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minobsnights" || string(argv[i]) == "-mon" || string(argv[i]) == "-minobs" || string(argv[i]) == "-minobsnight" || string(argv[i]) == "--minobsnights" || string(argv[i]) == "--minobsnight" || string(argv[i]) == "-minobsn") {
      if(i+1 < argc) {
	//There is still something to read;
	minobsnights=stoi(argv[++i]);
	mindaysteps = minobsnights-1; // Nights are fenceposts, daysteps are the spaces in between.
	default_mindaysteps = 0;
	i++;
      }
      else {
	cerr << "Min. number of observing nights keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mintimespan" || string(argv[i]) == "-mts" || string(argv[i]) == "-minspan" || string(argv[i]) == "-mintspan" || string(argv[i]) == "--mts" || string(argv[i]) == "--mintimespan" || string(argv[i]) == "--mintspan") {
      if(i+1 < argc) {
	//There is still something to read;
	mintimespan=stod(argv[++i]);
	default_mintimespan = 0;
	i++;
      }
      else {
	cerr << "Minimum time span keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-mingeodist" || string(argv[i]) == "-mingd" || string(argv[i]) == "-mingeo" || string(argv[i]) == "-mingeod" || string(argv[i]) == "--mingeodist" || string(argv[i]) == "--minimum_geocentr_dist") {
      if(i+1 < argc) {
	//There is still something to read;
	mingeodist=stod(argv[++i]);
	default_mingeodist = 0;
	i++;
      }
      else {
	cerr << "Minimum geocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-maxgeodist" || string(argv[i]) == "-maxgd" || string(argv[i]) == "-maxgeo" || string(argv[i]) == "-maxgeod" || string(argv[i]) == "--maxgeodist" || string(argv[i]) == "--maximum_geocentr_dist") {
      if(i+1 < argc) {
	//There is still something to read;
	maxgeodist=stod(argv[++i]);
	default_maxgeodist = 0;
	i++;
      }
      else {
	cerr << "Maximum geocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-geologstep" || string(argv[i]) == "-gls" || string(argv[i]) == "-geostep" || string(argv[i]) == "-glogstep" || string(argv[i]) == "--geologstep" || string(argv[i]) == "--geodistlogstep" || string(argv[i]) == "--geodiststep") {
      if(i+1 < argc) {
	//There is still something to read;
	geologstep=stod(argv[++i]);
	default_geologstep = 0;
	i++;
      }
      else {
	cerr << "Geocentric distance logarithmic step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mingeoobs" || string(argv[i]) == "-mgo" || string(argv[i]) == "-minobsdist" || string(argv[i]) == "-mindistobs" || string(argv[i]) == "--min_geocentric_obsdist" || string(argv[i]) == "--min_observation_distance" || string(argv[i]) == "--mingeoobs") {
      if(i+1 < argc) {
	//There is still something to read;
	mingeoobs=stod(argv[++i]);
	default_mingeoobs = 0;
	i++;
      }
      else {
	cerr << "Minimum geocentric distance at observation keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minimpactpar" || string(argv[i]) == "-mip" || string(argv[i]) == "-minimp" || string(argv[i]) == "-minimppar" || string(argv[i]) == "--minimum_impact_parameter" || string(argv[i]) == "--minimpactpar" || string(argv[i]) == "--min_impact_par") {
      if(i+1 < argc) {
	//There is still something to read;
	minimpactpar=stod(argv[++i]);
	default_minimpactpar = 0;
	i++;
      }
      else {
	cerr << "Minimum impact parameter keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Verbosity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outpair" || string(argv[i]) == "--outpairs") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	default_outfile = 0;
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outsum" || string(argv[i]) == "-sum" || string(argv[i]) == "-rms" || string(argv[i]) == "-outrms" || string(argv[i]) == "-or" || string(argv[i]) == "--outsummaryfile" || string(argv[i]) == "--outsum" || string(argv[i]) == "--sum") {
      if(i+1 < argc) {
	//There is still something to read;
	rmsfile=argv[++i];
	default_rmsfile = 0;
	i++;
      }
      else {
	cerr << "Output summary file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
  if(mindaysteps > npt-1) mindaysteps = npt-1; // Otherwise the low setting of npt is not operative.

  // Do something useful in the specific case that there is
  // an input data file, but no reference MJD

  if(indetfile.size()>0 && mjdref<=0L) {
    // Read input detection file to suggest optimal MJDref 
    instream1.open(indetfile,ios_base::in);
    if(!instream1) {
      cerr << "ERROR: unable to open input file " << indetfile << "\n";
      return(1);
    }
    mjdvec={};
    // Skip header line
    getline(instream1,lnfromfile);
    // Read body of the file
    detfilelinect=0;
    while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      // Read a line from the paired detections file, and store the MJD
      getline(instream1,lnfromfile);
      detfilelinect++;
      badread=0;
      if(lnfromfile.size()>60) {
	// Read MJD
	startpoint=0;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) MJD = stold(stest);
	else badread=1;
	// If there was a file read error, abort.
	if(badread==1) {
	  cerr << "ERROR reading line " << detfilelinect << " of paired detection file" << indetfile << "\n";
	  return(1);
	} else if(MJD<=0L) {
	  cerr << "ERROR reading line " << detfilelinect << " of paired detection file" << indetfile << ":\n";
	  cerr << "non-positive MJD value " << MJD << "\n";
	  return(1);
	}
	// If we reach this point, the line was read OK. Write it to mjdvec
	mjdvec.push_back(MJD);
      } else if(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
	cerr << "WARNING: line " << detfilelinect << " of paired detection file " << indetfile << " was too short\n";
      }
    }
    instream1.close();
    sort(mjdvec.begin(),mjdvec.end());
    cout << "\nERROR: input positive-valued reference MJD is required\n";
    cout << fixed << setprecision(2) << "Suggested value is " << mjdvec[0]*0.5L + mjdvec[mjdvec.size()-1]*0.5L << "\n";
    cout << "based on your input detection catalog " << indetfile << "\n\n\n";
    show_usage();
    return(1);
  }
  
  if(argc<11)
    {
      cerr << "Too few arguments even for minimalist invocation:\n";
      show_usage();
      return(1);
    }
  
  cout.precision(17);  
  cout << "input detection file " << indetfile << "\n";
  cout << "input pair file " << inpairfile << "\n";
  cout << "input observer position file " << planetfile << "\n";
  cout << "input heliocentric hypothesis file " << accelfile << "\n";
  cout << "input reference MJD " << mjdref << "\n";

  // Catch required parameters if missing
  if(indetfile.size()<=0) {
    cout << "\nERROR: input detection file is required\n";
    show_usage();
    return(1);
  } else if(inpairfile.size()<=0) {
    cout << "\nERROR: input pair file is required\n";
    show_usage();
    return(1);
  } else if(planetfile.size()<=0) {
    cout << "\nERROR: input observer position file is required:\n";
    cout << "e.g. Earth1day2020s_02a.txt\n";
    show_usage();
    return(1);
  } else if(accelfile.size()<=0) {
    cout << "\nERROR: input heliocentric hypothesis file is required\n";
    show_usage();
    return(1);
  }


  if(default_cluster_radius==1) cout << "Defaulting to cluster radius = " << cluster_radius << "km\n";
  else cout << "input clustering radius " << cluster_radius << "km\n";
  if(default_npt==1) cout << "Defaulting to DBSCAN npt (min. no. of tracklets in a linkage) = " << npt << "\n";
  else cout << "input DBSCAN npt (min. no. of tracklets in a linkage) is " << npt << "\n";
  if(default_mindaysteps==1) cout << "Defaulting to minimum number of unique nights = " << mindaysteps+1 << "\n";
  else cout << "minimum number of unique nights is " << mindaysteps+1 << "\n";
  if(default_mintimespan==1) cout << "Defaulting to minimum time span for a linkage = " << mintimespan << " days\n";
  else cout << "minimum time span for a linkage is " << mintimespan << " days\n";
  if(default_mingeodist==1) cout << "Defaulting to minimum geocentric distance = " << mingeodist << " AU\n";
  else cout << "minimum geocentric distance is " << mingeodist << " AU\n";
  if(default_maxgeodist==1) cout << "Defaulting to maximum geocentric distance = " << maxgeodist << " AU\n";
  else cout << "maximum geocentric distance is " << maxgeodist << " AU\n";
  if(default_geologstep==1) cout << "Defaulting to logarithmic step size for geocentric distance bins = " << geologstep << "\n";
  else cout << "logarithmic step size for geocentric distance bins is " << geologstep << "\n";
  if(default_mingeoobs==1) cout << "Defaulting to minimum geocentric distance at observation = " << mingeoobs << " AU\n";
  else cout << "Minimum geocentric distance at observation = " << mingeoobs << " AU\n";
  if(default_minimpactpar==1) cout << "Defaulting to minimum impact parameter = " << minimpactpar << " km\n";
  else cout << "Minimum impact parameter is " << minimpactpar << " km\n";
  if(default_outfile==1) cout << "WARNING: using default name " << outfile << " for comprehensive output file\n";
  else cout << "comprehensive output file " << outfile << "\n";
  if(default_rmsfile==1) cout << "WARNING: using default name " << rmsfile << " for summary output file\n";
  else cout << "summary output file " << rmsfile << "\n";

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
    if(long(heliovel.size()) != accelnum || long(heliodist.size()) != accelnum) {
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
  instream1.open(indetfile,ios_base::in);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << indetfile << "\n";
    return(1);
  }
  detvec={};
  // Skip header line
  getline(instream1,lnfromfile);
  cout << "Header line from input paired detection file " << indetfile << ":\n";
  cout << lnfromfile << "\n";
  // Read body of the file
  detfilelinect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the paired detections file, and load an object of class det_obsmag_indvec
    getline(instream1,lnfromfile);
    detfilelinect++;
    badread=0;
    if(lnfromfile.size()>60) {
      // Read MJD, RA, Dec, observer x, y, z
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) MJD = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) RA = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Dec = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) X = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Y = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Z = stold(stest);
      else badread=1;
      // Read the IDstring
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(detid,stest,SHORTSTRINGLEN);
      else badread=1;
      // Read the magnitude
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) mag = stod(stest);
      else badread=1;
      // Read the band and observatory code
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(band,stest,MINSTRINGLEN);
      else badread=1;
       startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(obscode,stest,MINSTRINGLEN);
      else badread=1;
      // Read the original detection index
       startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) origind = stol(stest);
      else badread=1;

      // If there was a file read error, abort.
      if(badread==1) {
	cerr << "ERROR reading line " << detfilelinect << " of paired detection file" << indetfile << "\n";
	return(1);
      }
      // If we reach this point, the line was read OK. Write it to detvec.
      o1=det_obsmag_indvec(MJD,RA,Dec,X,Y,Z,detid,mag,band,obscode,origind,{});
      detvec.push_back(o1);
      if(MJD < minMJD) minMJD = MJD;
      if(MJD > maxMJD) maxMJD = MJD;
    } else if(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      cerr << "WARNING: line " << detfilelinect << " of paired detection file " << indetfile << " was too short\n";
    }
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
  outstream1 << "#clusternum,posRMS,velRMS,totRMS,pairnum,timespan,uniquepoints,obsnights,metric,rating,heliodist,heliovel,helioacc,posX,posY,posZ,velX,velY,velZ\n";
  cout << "Writing to " << outfile << "\n";
  cout << "Read " << detvec.size() << " detections and " << pairvec.size() << " pairs.\n";

  realclusternum=0;
  for(accelct=0;accelct<accelnum;accelct++) {
    gridpoint_clusternum=0;
    badpoint=0;
    // Calculate approximate heliocentric distances from the
    // input quadratic approximation.
    if(verbose>=0) {
      cout << fixed << setprecision(4) << "Working on grid point " << accelct << ": r = " << heliodist[accelct]/AU_KM << " AU, v = " << heliovel[accelct]/SOLARDAY << " km/sec, dv/dt = " << helioacc[accelct]/(GMSUN_KM3_SEC2*SOLARDAY*SOLARDAY/heliodist[accelct]/heliodist[accelct]) << " GMsun/r^2\n";
      cout.precision(17);
    }
    heliodistvec={};
    for(i=0;i<long(detvec.size());i++)
      {
	delta1 = detvec[i].MJD - mjdref;
	heliodistvec.push_back(heliodist[accelct] + heliovel[accelct]*delta1 + 0.5*helioacc[accelct]*delta1*delta1);
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
    for(pairct=0; pairct<long(pairvec.size()); pairct++) {
      badpoint=0;
      // Obtain indices to the detection and heliocentric distance vectors.
      i1=pairvec[pairct][0];
      i2=pairvec[pairct][1];
      // Project the first point
      RA = detvec[i1].RA;
      Dec = detvec[i1].Dec;
      celestial_to_stateunitLD(RA,Dec,unitbary);
      observerpos1 = point3LD(detvec[i1].x,detvec[i1].y,detvec[i1].z);
      targposvec1={};
      deltavec1={};
      status1 = helioproj02LD(unitbary,observerpos1, heliodistvec[i1], deltavec1, targposvec1);
      RA = detvec[i2].RA;
      Dec = detvec[i2].Dec;
      celestial_to_stateunitLD(RA,Dec,unitbary);
      observerpos2 = point3LD(detvec[i2].x,detvec[i2].y,detvec[i2].z);
      targposvec2={};
      deltavec2={};
      status2 = helioproj02LD(unitbary, observerpos2, heliodistvec[i2], deltavec2, targposvec2);
      
      if(status1 > 0 && status2 > 0 && badpoint==0) {
	// Calculate time difference between the observations
	timediff = (detvec[i2].MJD - detvec[i1].MJD)*SOLARDAY;
	// Did helioproj find two solutions in both cases, or only one?
	num_dist_solutions = status1;
	if(num_dist_solutions > status2) num_dist_solutions = status2;
	// Loop over solutions (num_dist_solutions can only be 1 or 2).
	for(solnct=0; solnct<num_dist_solutions; solnct++) {
	  // Begin new stuff added to eliminate 'globs'
	  // These are spurious linkages of unreasonably large numbers (typically tens of thousands)
	  // of detections that arise when the hypothetical heliocentric distance at a time when
	  // many observations are acquired is extremely close to, but slightly greater than,
	  // the heliocentric distance of the observer. Then detections over a large area of sky
	  // wind up with projected 3-D positions in an extremely small volume -- and furthermore,
	  // they all have similar velocities because the very small geocentric distance causes
	  // the inferred velocities to be dominated by the observer's motion and the heliocentric
	  // hypothesis, with only a negligible contribution from the on-sky angular velocity.
	  glob_warning=0;
	  if(deltavec1[solnct]<mingeoobs*AU_KM && deltavec2[solnct]<mingeoobs*AU_KM) {
	    // New-start
	    // Load target positions
	    targpos1 = targposvec1[solnct];
	    targpos2 = targposvec2[solnct];
	    // Calculate positions relative to observer
	    targpos1.x -= observerpos1.x;
	    targpos1.y -= observerpos1.y;
	    targpos1.z -= observerpos1.z;
	    
	    targpos2.x -= observerpos2.x;
	    targpos2.y -= observerpos2.y;
	    targpos2.z -= observerpos2.z;
	    
	    // Calculate velocity relative to observer
	    targvel1.x = (targpos2.x - targpos1.x)/timediff;
	    targvel1.y = (targpos2.y - targpos1.y)/timediff;
	    targvel1.z = (targpos2.z - targpos1.z)/timediff;
   
	    // Calculate impact parameter (past or future).
	    absvelocity = vecabs3LD(targvel1);
	    impactpar = dotprod3LD(targpos1,targvel1)/absvelocity;
	    // Effectively, we've projected targpos1 onto the velocity
	    // vector, and impactpar temporarily holds the magnitude of this projection.
	    // Subtract off the projection of the distance onto the velocity unit vector
	    targpos1.x -= impactpar*targvel1.x/absvelocity;
	    targpos1.y -= impactpar*targvel1.y/absvelocity;
	    targpos1.z -= impactpar*targvel1.z/absvelocity;
	    // Now targpos1 is the impact parameter vector at projected closest approach.
	    impactpar  = vecabs3LD(targpos1); // Now impactpar is really the impact parameter
	    if(impactpar<=minimpactpar) {
	      // The hypothesis implies the object already passed with minimpactpar km of the Earth
	      // in the likely case that minimpactpar has been set to imply an actual impact,
	      // it's not our problem anymore.
	      glob_warning=1;
	    }
	  }
	  if(!glob_warning) {
	    targpos1 = targposvec1[solnct];
	    targpos2 = targposvec2[solnct];
	  
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
	      // Kepler integration encountered unphysical situation.
	    continue;
	    }
	  }
	}
      } else {
	  badpoint=1;
	  // Heliocentric projection found no physical solution.
	  continue;
      }
    }
  
    if(allstatevecs.size()<=1) continue; // No clusters possible, skip to the next step.
    if(verbose>=0) cout << pairvec.size() << " input pairs/tracklets led to " << allstatevecs.size() << " physically reasonable state vectors\n";

    // Loop over geocentric bins, selecting the subset of state-vectors
    // in each bin, and running DBSCAN only on those, with clustering radius
    // adjusted accordingly.
    geobinct = 0;
    georadcen = mingeodist*intpowD(geologstep,geobinct);
    while(georadcen<=maxgeodist) {
      georadcen = mingeodist*intpowD(geologstep,geobinct);
      georadmin = georadcen/geologstep;
      georadmax = georadcen*geologstep;
      // Load new array of state vectors, limited to those in the current geocentric bin
      binstatevecs={};
      for(i=0;i<long(allstatevecs.size());i++)
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

      if(verbose>=1) cout << "Found " << binstatevecs.size() << " state vectors in geocentric bin from " << georadmin << " to " << georadmax << " AU\n";
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
    
      if(verbose>=1) cout << "Created a KD tree with " << kdvec.size() << " branches\n";

      outclusters={};
      int clusternum = DBSCAN_6i01(kdvec, cluster_radius*(georadcen/CRAD_REF_GEODIST)/INTEGERIZING_SCALE, npt, INTEGERIZING_SCALE, outclusters, verbose);
      if(verbose>=1) cout << "DBSCAN_6i01 finished, with " << clusternum << " = " << outclusters.size() << " clusters found\n";
      for(clusterct=0; clusterct<long(outclusters.size()); clusterct++) {
	// Scale cluster RMS down to reference geocentric distance
	if(DEBUG_A >= 1) cout << "scaling outclusters rms for cluster " << clusterct << " out of " << outclusters.size() << "\n";
	fflush(stdout);
	for(i=0;i<9;i++) {
	  if(DEBUG_A >= 1) cout << "scaling rmsvec point " << i << " out of " << outclusters[clusterct].rmsvec.size() << "\n";
	  if(DEBUG_A >= 1) cout << "RMS = " << outclusters[clusterct].rmsvec[i];
	  outclusters[clusterct].rmsvec[i] *= CRAD_REF_GEODIST/georadcen;
	  if(DEBUG_A >= 1) cout << ", scales to " << outclusters[clusterct].rmsvec[i] << "\n";
	}
	// Note that RMS is scaled down for more distant clusters, to
	// avoid bias against them in post-processing.
	
	// Map cluster to individual detections.
	// create vector of unique detection indices.
	if(DEBUG_A >= 1) cout << "Loading pointind for " << outclusters[clusterct].numpoints << " of " << clusterct <<  " out of " << outclusters.size() << "\n";
	fflush(stdout);
	pointind={};
	for(i=0;i<outclusters[clusterct].numpoints;i++) {
	  pairct=kdvec[outclusters[clusterct].clustind[i]].point.i1;
	  for(j=0; j<long(pairvec[pairct].size()); j++) {
	    pointind.push_back(pairvec[pairct][j]);
	  }
	}
	// Sort vector of detection indices
	if(DEBUG_A >= 1) cout << "About to sort pointind\n";
	fflush(stdout);
	sort(pointind.begin(), pointind.end());
	if(DEBUG_A >= 1) cout << "Done sorting pointind\n";
	fflush(stdout);
	// Cull out duplicate entries
	pointjunk = pointind;
	pointind={};
	pointind.push_back(pointjunk[0]);
	for(i=1; i<long(pointjunk.size()); i++) {
	  if(pointjunk[i]!=pointjunk[i-1]) pointind.push_back(pointjunk[i]);
	}
	if(DEBUG_A >= 1) cout << "Done culling pointind\n";
	fflush(stdout);

	// Load vector of detection MJD's
	clustmjd = {};
	for(i=0; i<long(pointind.size()); i++) {
	  clustmjd.push_back(detvec[pointind[i]].MJD);
	}
	if(DEBUG_A >= 1) cout << "Done loading mjds\n";
	fflush(stdout);
	// Sort vector of MJD's
	sort(clustmjd.begin(), clustmjd.end());
	if(DEBUG_A >= 1) cout << "Done sorting mjds\n";
	fflush(stdout);
	timespan = clustmjd[clustmjd.size()-1] - clustmjd[0];
	if(DEBUG_A >= 1) cout << "Timespan = " << timespan << "\n";
	fflush(stdout);
	// Load vector of MJD steps
	mjdstep={};
	for(i=1; i<long(clustmjd.size()); i++) {
	  mjdstep.push_back(clustmjd[i]-clustmjd[i-1]);
	}
	if(DEBUG_A >= 1) cout << "Done loading vector of steps\n";
	fflush(stdout);
	// Count steps large enough to suggest a daytime period between nights.
	numdaysteps=0;	
	for(i=0; i<long(mjdstep.size()); i++) {
	  if(mjdstep[i]>INTRANIGHTSTEP) numdaysteps++;
	}
	if(DEBUG_A >= 1) cout << "Unique pts: " << pointind.size() << " span: " << timespan << " daysteps: " << numdaysteps << "\n";
	// Does cluster pass the criteria for a linked detection?
	if(timespan >= mintimespan && numdaysteps >= mindaysteps) {
	  realclusternum++;
	  gridpoint_clusternum++;
	  if(DEBUG_A >= 1) cout << "Loading good cluster " << realclusternum << " with timespan " << timespan << " and numdaysteps " << numdaysteps << "\n";
	  fflush(stdout);
	  numobsnights = numdaysteps+1;
	  if(DEBUG_A >= 1) cout << "Cluster passes discovery criteria: will be designated as cluster " << realclusternum << "\n";
	  // Check whether cluster is composed purely of detections from
	  // a single simulated object (i.e., would be a real discovery) or is a mixture
	  // of detections from two or more different simulated objects (i.e., spurious).
	  rating="PURE";
	  for(i=0; i<long(pointind.size()); i++) {
	    if(i>0 && stringnmatch01(detvec[pointind[i]].idstring,detvec[pointind[i-1]].idstring,SHORTSTRINGLEN)!=0) rating="MIXED";
	  }
	  if(DEBUG_A >= 1) cout << "Rating is found to be " << rating << "\n";
	  fflush(stdout);
	  // Write all individual detections in this cluster to the output cluster file
	  for(i=0; i<long(pointind.size()); i++) {
	    if(DEBUG_A >= 1) cout << "Writing point " << i << " out of " << pointind.size() << " to cluster file\n";
	    fflush(stdout);
	    outstream2  << fixed << setprecision(6) << intzero01i(i,4) << "," << detvec[pointind[i]].MJD << "," << detvec[pointind[i]].RA << "," << detvec[pointind[i]].Dec << "," << detvec[pointind[i]].idstring << ",";
	    outstream2  << fixed << setprecision(3) << detvec[pointind[i]].mag << "," << detvec[pointind[i]].band << "," << detvec[pointind[i]].obscode << "," << pointind[i] << "," << detvec[pointind[i]].index << "," << realclusternum << "\n";
	  }
	  if(DEBUG_A >= 1) cout << "Finished writing to cluster file\n";
	  fflush(stdout);
	  // Write summary line to rms file
	  clustmetric = double(pointind.size())*double(numobsnights)*timespan/outclusters[clusterct].rmsvec[8];
	  outstream1  << fixed << setprecision(3) << realclusternum << "," << outclusters[clusterct].rmsvec[6] << "," << outclusters[clusterct].rmsvec[7] << "," << outclusters[clusterct].rmsvec[8] << "," << outclusters[clusterct].numpoints << ",";
	  outstream1  << fixed << setprecision(6) << timespan << "," << pointind.size() << "," << numobsnights  << "," << clustmetric << "," << rating << ",";
	  outstream1  << fixed << setprecision(3) << heliodist[accelct]/AU_KM << "," << heliovel[accelct]/SOLARDAY << ",";
	  outstream1  << fixed << setprecision(6) << helioacc[accelct]*1000.0/SOLARDAY/SOLARDAY << ",";
	  outstream1  << fixed << setprecision(3)  << outclusters[clusterct].meanvec[0] << "," << outclusters[clusterct].meanvec[1] << "," << outclusters[clusterct].meanvec[2] << ",";
	  outstream1  << fixed << setprecision(6) << outclusters[clusterct].meanvec[3]/chartimescale << ","   << outclusters[clusterct].meanvec[4]/chartimescale << "," << outclusters[clusterct].meanvec[5]/chartimescale << "\n";
	  if(DEBUG_A >= 1) cout << "Finished writing to RMS file\n";
	  fflush(stdout);
	}
      }
      // Move on to the next bin in geocentric distance
      geobinct++;
    }
    if(verbose>=0) cout << "Identified " << gridpoint_clusternum << " candidate linkages\n";
  }
  return(0);
}
