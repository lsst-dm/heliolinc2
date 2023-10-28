// October 18, 2023: heliovane_danby.cpp:
// Like heliolinc_danby and similar  programs to which it is closely
// related, heliovane_danby searches for asteroids using hypotheses
// about their heliocentric motion. However, while the heliolinc
// programs use hypotheses about the heliocentric distance and radial
// velocity (the heliolinc parameter space), and hence look for
// asteroids on heliocentric spheres, heliovane uses hypotheses
// about the rotation of asteroids in heliocentric ecliptic longitude,
// and hence looks for asteroids in planes of constant heliocentric
// ecliptic longitude -- i.e., vanes extending outward from the sun.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: heliovane_danby -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection file -mjd mjdref -obspos observer_position_file -heliolon heliocentric_longitude_hypothesis_file -clustrad clustrad -npt dbscan_npt -minobsnights minobsnights -mintimespan mintimespan -mingeodist minimum_geocentric_distance -maxgeodist maximum_geocentric_distance -geologstep logarithmic_step_size_for_geocentric_distance_bins -minsunelong minimum_solar_elongation(deg) -maxsunelong maximum_solar_elongation(deg) -min_incid_angle min_angle_of_incidence_between_observer_to_target_vector_and_heliocentric_vane(deg) -maxheliodist maximum_heliocentric_radius(AU) -mingeoobs min_geocentric_dist_at_observation(AU) -minimpactpar min_impact_parameter(km) -useunivar 1_for_univar_0_for_fgfunc -vinf max_v_inf -outsum summary_file -clust2det clust2detfile -verbose verbosity\n";
  cerr << "\nor, at minimum:\n\n";
  cerr << "heliovane_danby -dets detfile -trk2det tracklet-to-detection file -mjd mjdref -obspos observer_position_file -heliodist heliocentric_dist_vel_acc_file\n";
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
  vector <hlradhyp> lambdahyp;
  vector <EarthState> earthpos;
  HeliovaneConfig config;
  vector <hlclust> outclust;
  vector <longpair> clust2det;
  string imfile,pairdetfile,trackletfile,trk2detfile,planetfile,lambdafile;
  string sumfile = "sumfile_test.csv";
  string clust2detfile = "clust2detfile_test.csv";
  int default_clustrad, default_npt, default_minobsnights;
  int default_mintimespan, default_mingeodist, default_maxgeodist;
  int default_geologstep,default_clust2detfile,default_sumfile;
  int default_mingeoobs, default_minimpactpar;
  int default_use_univar, default_max_v_inf;
  default_clustrad = default_npt = default_minobsnights = default_mintimespan = 1;
  default_mingeodist = default_maxgeodist = default_geologstep = default_clust2detfile = default_sumfile = 1;
  default_mingeoobs = default_minimpactpar = 1;
  default_use_univar = default_max_v_inf = 1;
  ofstream outstream1;
  long i=0;
  long clustct=0;
  int status=0;  
  int default_minsunelong = 0;
  int default_maxsunelong = 0;
  int default_min_incid_angle = 0;
  int default_maxheliodist = 0;
  
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
	cerr << "Output paired detection file keyword supplied with no corresponding argument\n";
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
	cerr << "Output pair file keyword supplied with no corresponding argument\n";
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
	cerr << "Output pair file keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-heliolon" || string(argv[i]) == "-heliolambda" || string(argv[i]) == "-heliolam" || string(argv[i]) == "-hlam" || string(argv[i]) == "--heliolon" || string(argv[i]) == "--heliolambda") {
      if(i+1 < argc) {
	//There is still something to read;
	lambdafile=argv[++i];
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
    } else if(string(argv[i]) == "-useunivar" || string(argv[i]) == "-use_univar" || string(argv[i]) == "-univar" || string(argv[i]) == "-universalvar") {
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
    } else if(string(argv[i]) == "-minsunelong" || string(argv[i]) == "-minelong" || string(argv[i]) == "-minse") {
      if(i+1 < argc) {
	//There is still something to read;
	config.minsunelong=stod(argv[++i]);
	default_minsunelong = 0;
	i++;
      }
      else {
	cerr << "Minimum solar elongation keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxsunelong" || string(argv[i]) == "-maxelong" || string(argv[i]) == "-maxse") {
      if(i+1 < argc) {
	//There is still something to read;
	config.maxsunelong=stod(argv[++i]);
	default_maxsunelong = 0;
	i++;
      }
      else {
	cerr << "Maximum solar elongation keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-min_incid_angle" || string(argv[i]) == "-min_ang_incid" || string(argv[i]) == "-minincid" || string(argv[i]) == "-minincidang" || string(argv[i]) == "-minincidangle" || string(argv[i]) == "-minangincid") {
      if(i+1 < argc) {
	//There is still something to read;
	config.min_incid_angle=stod(argv[++i]);
	default_min_incid_angle = 0;
	i++;
      }
      else {
	cerr << "Minimum angle of incidence keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxheliodist" || string(argv[i]) == "-maxhelio" || string(argv[i]) == "-maxhdist" || string(argv[i]) == "-maxheliorad" || string(argv[i]) == "-maxhrad") {
      if(i+1 < argc) {
	//There is still something to read;
	config.maxheliodist=stod(argv[++i]);
	default_maxheliodist = 0;
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
  cout << "input heliocentric longitude hypothesis file " << lambdafile << "\n";
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
  } else if(lambdafile.size()<=0) {
    cout << "\nERROR: input heliocentric longitude hypothesis file is required\n";
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
  if(default_minsunelong==1) cout << "Defaulting to minimum solar elongation = " << config.minsunelong << " deg\n";
  else cout << "Minimum solar elongation is " << config.minsunelong << " deg\n";
  if(default_maxsunelong==1) cout << "Defaulting to maximum solar elongation = " << config.maxsunelong << " deg\n";
  else cout << "Maximum solar elongation is " << config.maxsunelong << " deg\n";
  if(default_min_incid_angle==1) cout << "Defaulting to minimum angle of incidence = " << config.min_incid_angle << " deg\n";
  else cout << "Minimum angle of incidence is " << config.min_incid_angle << " deg\n";
  if(default_maxheliodist==1) cout << "Defaulting to maximum heliocentric distance = " << config.maxheliodist << " AU\n";
  else cout << "Maximum heliocentric distance is " << config.maxheliodist << " AU\n";
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

  cout << "File with heliocentric longitude hypotheses is named " << lambdafile << "\n";
  status=read_radhyp_file(lambdafile, lambdahyp, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read accleration file " << lambdafile << "\n";
    cerr << "read_radhyp_file returned status = " << status << ".\n";
   return(1);
  }
  if(config.verbose>=1) {
    for(i=0;i<long(lambdahyp.size());i++) {
      cout << "Lambda matrix check: " << lambdahyp[i].HelioRad << " " << lambdahyp[i].R_dot << " " << lambdahyp[i].R_dubdot << "\n";
    }
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
  
  status=heliovane_alg_danby(image_log, detvec, tracklets, trk2det, lambdahyp, earthpos, config, outclust, clust2det);
  if(status!=0) {
    cerr << "ERROR: heliovane_alg_danby failed with status " << status << "\n";
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
