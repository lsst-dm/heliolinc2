// June 30, 2023: attempts to parallelize the heliocentric
// hypothesis loop.
//
// September 14, 2022: created on this date by copying
// an earlier prototype file called projectpairs06d.cpp.
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

static void show_usage()
{
  cerr << "Usage: heliolinc_ompkd -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection file -mjd mjdref -obspos observer_position_file -heliodist heliocentric_dist_vel_acc_file -clustrad clustrad -clustchangerad min_distance_for_cluster_scaling -npt dbscan_npt -minobsnights minobsnights -mintimespan mintimespan -mingeodist minimum_geocentric_distance -maxgeodist maximum_geocentric_distance -geologstep logarithmic_step_size_for_geocentric_distance_bins -mingeoobs min_geocentric_dist_at_observation(AU) -minimpactpar min_impact_parameter(km) -useunivar 1_for_univar_0_for_fgfunc -vinf max_v_inf  -outsum summary_file -clust2det clust2detfile -verbose verbosity\n";
  cerr << "\nor, at minimum:\n\n";
  cerr << "heliolinc_ompkd -dets detfile -trk2det tracklet-to-detection file -mjd mjdref -obspos observer_position_file -heliodist heliocentric_dist_vel_acc_file\n";
  cerr << "\nNote that the minimum invocation leaves some things set to defaults\n";
  cerr << "that you may well wish to specify: in particular, the output file names\n";
  
}
    
int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  vector <hldet> pairdets = {};
  vector <hlimage> image_log;
  vector <tracklet> tracklets;
  vector <longpair> trk2det;
  vector <hlradhyp> radhyp;
  vector <EarthState> earthpos;
  HeliolincConfig config;
  vector <hlclust> outclust;
  vector <longpair> clust2det;
  string imfile,pairdetfile,trackletfile,trk2detfile,planetfile,accelfile;
  string sumfile = "sumfile_test.csv";
  string clust2detfile = "clust2detfile_test.csv";
  int default_clustrad, default_clustchangerad, default_npt, default_minobsnights;
  int default_mintimespan, default_mingeodist, default_maxgeodist;
  int default_geologstep,default_clust2detfile,default_sumfile;
  int default_mingeoobs, default_minimpactpar;
  int default_use_univar, default_max_v_inf;
  default_clustrad = default_clustchangerad = default_npt = default_minobsnights = 1;
  default_mintimespan = 1;
  default_mingeodist = default_maxgeodist = default_geologstep = 1;
  default_clust2detfile = default_sumfile = 1;
  default_mingeoobs = default_minimpactpar = 1;
  default_use_univar = default_max_v_inf = 1;
  ofstream outstream1;
  long i=0;
  long clustct=0;
  int status=0;  
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" || string(argv[i]) == "-detpairs") {
      if(i+1 < argc) {
	//There is still something to read;
        pairdetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input paired detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
	//There is still something to read;
	imfile=argv[++i];
	i++;
      }
      else {
	cerr << "Image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-tracklets" || string(argv[i]) == "-tracklet" || string(argv[i]) == "-tf" || string(argv[i]) == "-trkfile" || string(argv[i]) == "-trackletfile" || string(argv[i]) == "--trackletfile") {
      if(i+1 < argc) {
	//There is still something to read;
	trackletfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input pair file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trk2det" || string(argv[i]) == "-t2d" || string(argv[i]) == "-tracklet2detection" || string(argv[i]) == "-trk2detfile" || string(argv[i]) == "--tracklet2detection") {
      if(i+1 < argc) {
	//There is still something to read;
	trk2detfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input pair file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-m" || string(argv[i]) == "-mjd" || string(argv[i ]) == "-mjdref" || string(argv[i]) == "-MJD" || string(argv[i]) == "--mjd" || string(argv[i]) == "--MJD" || string(argv[i]) == "--modifiedjulianday" || string(argv[i]) == "--ModifiedJulianDay") {
      if(i+1 < argc) {
	//There is still something to read;
	config.MJDref=stod(argv[++i]);
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
    } else if(string(argv[i]) == "-clustrad" || string(argv[i]) == "-cr" || string(argv[i]) == "-crad" || string(argv[i]) == "-cluster" || string(argv[i]) == "--clustrad" || string(argv[i]) == "--clusterradius" || string(argv[i]) == "--clusterrad") {
      if(i+1 < argc) {
	//There is still something to read;
	config.clustrad=stod(argv[++i]);
	default_clustrad = 0;
	i++;
      }
      else {
	cerr << "Clustering radius keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clustchangerad" || string(argv[i]) == "-clustchangerad" || string(argv[i]) == "-clustchangerad") {
      if(i+1 < argc) {
	//There is still something to read;
	config.clustchangerad = stod(argv[++i]); 
	default_clustchangerad = 0;
	i++;
      }
      else {
	cerr << "Transition distance for cluster scaling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-npt" || string(argv[i]) == "-npoints" || string(argv[i]) == "-minpts" || string(argv[i]) == "-np" || string(argv[i]) == "--npt" || string(argv[i]) == "--dbscan_npt" || string(argv[i]) == "--DBSCANnpt") {
      if(i+1 < argc) {
	//There is still something to read;
	config.dbscan_npt=stoi(argv[++i]);
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
	config.minobsnights=stoi(argv[++i]);
	default_minobsnights = 0;
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
	config.mintimespan=stod(argv[++i]);
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
	config.mingeodist=stod(argv[++i]);
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
	config.maxgeodist=stod(argv[++i]);
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
	config.geologstep=stod(argv[++i]);
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
	config.mingeoobs=stod(argv[++i]);
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
	config.minimpactpar=stod(argv[++i]);
	default_minimpactpar = 0;
	i++;
      }
      else {
	cerr << "Minimum impact parameter keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-useunivar" || string(argv[i]) == "-use_univar" || string(argv[i]) == "-univar" || string(argv[i]) == "-universalvar") {
      if(i+1 < argc) {
	//There is still something to read;
	config.use_univar=stoi(argv[++i]);
	default_use_univar = 0;
	i++;
      }
      else {
	cerr << "Minimum impact parameter keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-vinf" || string(argv[i]) == "-maxvinf" || string(argv[i]) == "-max_v_inf") {
      if(i+1 < argc) {
	//There is still something to read;
	config.max_v_inf=stod(argv[++i]);
	default_max_v_inf = 0;
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
	config.verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Verbosity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outsum" || string(argv[i]) == "-sum" || string(argv[i]) == "-rms" || string(argv[i]) == "-outrms" || string(argv[i]) == "-or" || string(argv[i]) == "--outsummaryfile" || string(argv[i]) == "--outsum" || string(argv[i]) == "--sum") {
      if(i+1 < argc) {
	//There is still something to read;
	sumfile=argv[++i];
	default_sumfile = 0;
	i++;
      }
      else {
	cerr << "Output summary file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clust2det" || string(argv[i]) == "-c2d" || string(argv[i]) == "-clust2detfile" || string(argv[i]) == "--clust2detfile" || string(argv[i]) == "--clust2det" || string(argv[i]) == "--cluster_to_detection") {
      if(i+1 < argc) {
	//There is still something to read;
	clust2detfile=argv[++i];
	default_clust2detfile = 0;
	i++;
      }
      else {
	cerr << "Output cluster-to-detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	config.verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Verbosity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
  if(config.minobsnights > config.dbscan_npt) config.minobsnights = config.dbscan_npt; // Otherwise the low setting of dbscan_npt is not operative.

  // Do something useful in the specific case that there is
  // an input data file, but no reference MJD
  if(pairdetfile.size()>0 && config.MJDref<=0L) {
    // Read input detection file to suggest optimal MJDref
    read_pairdet_file(pairdetfile, detvec, config.verbose);
    sort(detvec.begin(), detvec.end(), early_hldet());
    cout << "\nERROR: input positive-valued reference MJD is required\n";
    cout << "MJD range is " << detvec[0].MJD << " to " << detvec[detvec.size()-1].MJD << "\n";
    cout << fixed << setprecision(2) << "Suggested reference value is " << detvec[0].MJD*0.5L + detvec[detvec.size()-1].MJD*0.5L << "\n";
    cout << "based on your input detection catalog " << pairdetfile << "\n\n\n";
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
  cout << "input image file " << imfile << "\n";
  cout << "input detection file " << pairdetfile << "\n";
  cout << "input tracklet file " << trackletfile << "\n";
  cout << "input trk2det file " << trk2detfile << "\n";
  cout << "input observer position file " << planetfile << "\n";
  cout << "input heliocentric hypothesis file " << accelfile << "\n";
  cout << "input reference MJD " << config.MJDref << "\n";

  // Catch required parameters if missing
  if(pairdetfile.size()<=0) {
    cout << "\nERROR: input detection file is required\n";
    show_usage();
    return(1);
  } else if(trackletfile.size()<=0) {
    cout << "\nERROR: input tracklet file is required\n";
    show_usage();
    return(1);
  } else if(trk2detfile.size()<=0) {
    cout << "\nERROR: input trk2det file is required\n";
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

  // Catch case where max v_inf > 0 but universal variables are not set.
  if(config.max_v_inf>0.0l && config.use_univar<=0) {
    cerr << "ERROR: unbound orbits being probed (max_v_inf = " << config.max_v_inf << "),\n";
    cerr << "but code is not set to use universal variables (use_univar = " << config.use_univar << ".\n";
    show_usage();
    return(1);
  }

  if(default_clustrad==1) cout << "Defaulting to cluster radius = " << config.clustrad << "km\n";
  else cout << "input clustering radius " << config.clustrad << "km\n";
  if(default_clustchangerad==1) cout << "Defaulting to min. geocentric distance for cluster scaling = " << config.clustchangerad << "AU\n";
  else cout << "Min. geocentric distance for cluster scaling is " << config.clustchangerad << "AU\n";
  cout << "Minimum cluster radius, which will apply for all geocentric distances less\n";
  cout << "than " << config.clustchangerad << "AU, is " << config.clustrad*config.clustchangerad/REF_GEODIST << "km\n";
  if(default_npt==1) cout << "Defaulting to DBSCAN npt (min. no. of tracklets in a linkage) = " << config.dbscan_npt << "\n";
  else cout << "input DBSCAN npt (min. no. of tracklets in a linkage) is " << config.dbscan_npt << "\n";
  if(default_minobsnights==1) cout << "Defaulting to minimum number of unique nights = " << config.minobsnights << "\n";
  else cout << "minimum number of unique nights is " << config.minobsnights << "\n";
  if(default_mintimespan==1) cout << "Defaulting to minimum time span for a linkage = " << config.mintimespan << " days\n";
  else cout << "minimum time span for a linkage is " << config.mintimespan << " days\n";
  if(default_mingeodist==1) cout << "Defaulting to minimum geocentric distance = " << config.mingeodist << " AU\n";
  else cout << "minimum geocentric distance is " << config.mingeodist << " AU\n";
  if(default_maxgeodist==1) cout << "Defaulting to maximum geocentric distance = " << config.maxgeodist << " AU\n";
  else cout << "maximum geocentric distance is " << config.maxgeodist << " AU\n";
  if(default_geologstep==1) cout << "Defaulting to logarithmic step size for geocentric distance bins = " << config.geologstep << "\n";
  else cout << "logarithmic step size for geocentric distance bins is " << config.geologstep << "\n";
  if(default_mingeoobs==1) cout << "Defaulting to minimum geocentric distance at observation = " << config.mingeoobs << " AU\n";
  else cout << "Minimum geocentric distance at observation = " << config.mingeoobs << " AU\n";
  if(default_minimpactpar==1) cout << "Defaulting to minimum impact parameter = " << config.minimpactpar << " km\n";
  else cout << "Minimum impact parameter is " << config.minimpactpar << " km\n";
  if(default_use_univar==1) cout << "For Keplerian integration, defaulting to f and g functions\nrather than universal variables\n";
  else if(config.use_univar>0) cout << "Using universal variables for Keplerian integration\n";
  else cout << "Using f and g functions for Keplerian integration\n";
  if(default_max_v_inf==1) cout << "Defaulting to maximum v_inf relative to the sun = " << config.max_v_inf << " km/sec\n";
  else cout << "Maximum v_inf relative to the sun is " << config.max_v_inf << " km\n";
  if(default_sumfile==1) cout << "WARNING: using default name " << sumfile << " for summary output file\n";
  else cout << "summary output file " << sumfile << "\n";
  if(default_clust2detfile==1) cout << "WARNING: using default name " << clust2detfile << " for output clust2det file\n";
  else cout << "output clust2det file " << clust2detfile << "\n";

  cout << "Heliocentric ephemeris for Earth is named " << planetfile << "\n";
  status = read_horizons_csv(planetfile, earthpos);
  if(status!=0) {
    cerr << "ERROR: could not successfully Earth ephemeris file " << planetfile << "\n";
    cerr << "read_horizons_csv returned status = " << status << ".\n";
   return(1);
  } 

  cout << "File with heliocentric motion hypothesis is named " << accelfile << "\n";
  status=read_radhyp_file(accelfile, radhyp, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read accleration file " << accelfile << "\n";
    cerr << "read_radhyp_file returned status = " << status << ".\n";
   return(1);
  }
  
  detvec={};
  status=read_pairdet_file(pairdetfile, detvec, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << pairdetfile << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << detvec.size() << " data lines from paired detection file " << pairdetfile << "\n";
  
  image_log={};
  status=read_image_file2(imfile, image_log);
  if(status!=0) {
    cerr << "ERROR: could not successfully read image file " << imfile << "\n";
    cerr << "read_image_file2 returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << image_log.size() << " data lines from image file " << imfile << "\n";
  
  tracklets={};
  status=read_tracklet_file(trackletfile, tracklets, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read tracklet file " << trackletfile << "\n";
    cerr << "read_tracklet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << tracklets.size() << " data lines from tracklet file " << trackletfile << "\n";
  
  trk2det={};
  status=read_longpair_file(trk2detfile, trk2det, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read trk2det file " << trk2detfile << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << trk2det.size() << " data lines from trk2det file " << trk2detfile << "\n";
  
  status=heliolinc_alg_ompkd4(image_log, detvec, tracklets, trk2det, radhyp, earthpos, config, outclust, clust2det);
  if(status!=0) {
    cerr << "ERROR: heliolinc_alg failed with status " << status << "\n";
    return(status);
  } 
  
  outstream1.open(sumfile);
  cout << "Writing " << outclust.size() << " lines to output cluster-summary file " << sumfile << "\n";
  outstream1 << "#clusternum,posRMS,velRMS,totRMS,astromRMS,pairnum,timespan,uniquepoints,obsnights,metric,rating,heliohyp0,heliohyp1,heliohyp2,posX,posY,posZ,velX,velY,velZ,orbit_a,orbit_e,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count\n";
  for(clustct=0 ; clustct<long(outclust.size()); clustct++) {
    outstream1 << fixed << setprecision(3) << outclust[clustct].clusternum << "," << outclust[clustct].posRMS << "," << outclust[clustct].velRMS << "," << outclust[clustct].totRMS << ",";
    outstream1 << fixed << setprecision(4) << outclust[clustct].astromRMS << ",";
    outstream1 << fixed << setprecision(6) << outclust[clustct].pairnum << "," << outclust[clustct].timespan << "," << outclust[clustct].uniquepoints << "," << outclust[clustct].obsnights << "," << outclust[clustct].metric << "," << outclust[clustct].rating << ",";
    outstream1 << fixed << setprecision(6) << outclust[clustct].heliohyp0 << "," << outclust[clustct].heliohyp1 << "," << outclust[clustct].heliohyp2 << ",";
    outstream1 << fixed << setprecision(1) << outclust[clustct].posX << "," << outclust[clustct].posY << "," << outclust[clustct].posZ << ",";
    outstream1 << fixed << setprecision(4) << outclust[clustct].velX << "," << outclust[clustct].velY << "," << outclust[clustct].velZ << ",";
    outstream1 << fixed << setprecision(6) << outclust[clustct].orbit_a << "," << outclust[clustct].orbit_e << "," << outclust[clustct].orbit_MJD << ",";
    outstream1 << fixed << setprecision(1) << outclust[clustct].orbitX << "," << outclust[clustct].orbitY << "," << outclust[clustct].orbitZ << ",";
    outstream1 << fixed << setprecision(4) << outclust[clustct].orbitVX << "," << outclust[clustct].orbitVY << "," << outclust[clustct].orbitVZ << "," << outclust[clustct].orbit_eval_count << "\n";
  }
  outstream1.close();
  outstream1.open(clust2detfile);
  cout << "Writing " << clust2det.size() << " lines to output clust2det file " << clust2detfile << "\n";
  outstream1 << "#clusternum,detnum\n";
  for(clustct=0 ; clustct<long(clust2det.size()); clustct++) {
    outstream1 << clust2det[clustct].i1 << "," << clust2det[clustct].i2 << "\n";
  }
  outstream1.close();
  
  return(0);
}
