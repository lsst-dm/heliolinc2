// February 26, 2024: clusterprobe.cpp:
// Based on heliolinc_kd.cpp, this stripped-down, experimental
// program will probe the minimum viable clustering radius as
// a function of heliocentric distance and hypothesis mismatch.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define TIMECONVSCALE 4.0l // The characteristic timescale used to convert velocities
                           // to distance units is equal to the full temporal span
                           // divided by TIMECONVSCALE.

static void show_usage()
{
  cerr << "Usage: clusterprobe -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection file -mjd mjdref -obspos observer_position_file -heliodist heliocentric_dist_vel_acc_file -outsum summary_file -clust2det clust2detfile -innea 1 to use inner Earth geom -verbose verbosity\n";
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
  string imfile,pairdetfile,trackletfile,trk2detfile,planetfile,accelfile;
  string sumfile = "sumfile_test.csv";
  int default_sumfile = 1;
  ofstream outstream1;
  long i=0;
  int status=0;
  int use_innea=0;
  
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
    } else if(string(argv[i]) == "-innea" || string(argv[i]) == "-inner" || string(argv[i]) == "-inner_earth" || string(argv[i]) == "-innEa" || string(argv[i]) == "--inner_Earth" || string(argv[i]) == "--innea" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	use_innea=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Inner near-Earth geometry keyword supplied with no corresponding argument\n";
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


  if(default_sumfile==1) cout << "WARNING: using default name " << sumfile << " for summary output file\n";

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
  
  point3d Earthrefpos = point3d(0l,0l,0l);
  long imnum = image_log.size();
  long pairnum = tracklets.size();
  long trk2detnum = trk2det.size();
  long accelnum = radhyp.size();
  long accelct=0;

  vector <double> heliodist;
  vector <double> heliovel;
  vector <double> helioacc;
  vector <point6dx2> allstatevecs;
   
  if(imnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty image catalog\n";
    return(1);
  } else if(pairnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty tracklet array\n";
    return(1);
  } else if(trk2detnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty trk2det array\n";
    return(1);
  } else if(accelnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty heliocentric hypothesis array\n";
    return(1);
  }
  
  double MJDmin = image_log[0].MJD;
  double MJDmax = image_log[imnum-1].MJD;
  if(config.MJDref<MJDmin || config.MJDref>MJDmax) {
    // Input reference MJD is invalid. Suggest a better value before exiting.
    cerr << "\nERROR: reference MJD, supplied as " << config.MJDref << ",\n";
    cerr << "must fall in the time interval spanned by the data (" << MJDmin << " to " << MJDmax << "\n";
    cerr << fixed << setprecision(2) << "Suggested value is " << MJDmin*0.5l + MJDmax*0.5l << "\n";
    cout << "based on your input image catalog\n";
    return(2);
  }

  double chartimescale = (MJDmax - MJDmin)*SOLARDAY/TIMECONVSCALE; // Note that the units are seconds.
  Earthrefpos = earthpos01(earthpos, config.MJDref);

  // Convert heliocentric radial motion hypothesis matrix
  // from units of AU, AU/day, and GMSun/R^2
  // to units of km, km/day, and km/day^2.
  heliodist = heliovel = helioacc = {};
  for(accelct=0;accelct<accelnum;accelct++) {
    heliodist.push_back(radhyp[accelct].HelioRad * AU_KM);
    heliovel.push_back(radhyp[accelct].R_dot * AU_KM);
    helioacc.push_back(radhyp[accelct].R_dubdot * (-GMSUN_KM3_SEC2*SOLARDAY*SOLARDAY/heliodist[accelct]/heliodist[accelct]));
  }

  // Begin master loop over heliocentric hypotheses
  outstream1.open(sumfile);
  outstream1 << "#heliohyp0(AU) heliohyp1(AU/day) heliohyp2(GM/r^2) heliodist(km) heliovel(km/sec) helioacc(m/sec) minrad vradmin pradmin minvrad minprad numpoints\n";
  for(accelct=0;accelct<accelnum;accelct++) {
    // Covert all tracklets into state vectors at the reference time, under
    // the assumption that the heliocentric distance hypothesis is correct.

    if(use_innea==1) {
      status = trk2statevec_clusterprobe_innea(image_log, tracklets, heliodist[accelct], heliovel[accelct], helioacc[accelct], chartimescale, allstatevecs, config.MJDref);
    } else {
      status = trk2statevec_clusterprobe(image_log, tracklets, heliodist[accelct], heliovel[accelct], helioacc[accelct], chartimescale, allstatevecs, config.MJDref);
    }
    if(status==1) {
      cerr << "WARNING: hypothesis " << accelct << ": " << radhyp[accelct].HelioRad << " " << radhyp[accelct].R_dot << " " << radhyp[accelct].R_dubdot << " led to\nnegative heliocentric distance or other invalid result: SKIPPING\n";
      continue;
    } else if(status==2) {
      // This is a weirder error case and is fatal.
      cerr << "Fatal error case from trk2statevec.\n";
      return(3);
    }
    // If we get here, trk2statevec probably ran OK.
    if(allstatevecs.size()<=1) continue; // No clusters possible, skip to the next step.
    cout << pairnum << " input pairs/tracklets led to " << allstatevecs.size() << " physically reasonable state vectors\n";
    // Determine the minimum radius that will enclose the whole cluster,
    // if centered on the most central point.
    double minrad,vradmin,pradmin,minvrad,minprad;
    minrad = vradmin = pradmin = minvrad = minprad = LARGERR3;
    for(long i=0; i<long(allstatevecs.size()); i++) {
      double dist,vdist,pdist,maxdist,maxvdist,maxpdist;
      dist = vdist = pdist = 0.0;
      maxdist = maxvdist = maxpdist = 0.0l;
      for(long j=0; j<long(allstatevecs.size()); j++) {
	if(j!=i) {
	  dist = sqrt(DSQUARE(allstatevecs[i].x - allstatevecs[j].x) + DSQUARE(allstatevecs[i].y - allstatevecs[j].y) + DSQUARE(allstatevecs[i].z - allstatevecs[j].z) + DSQUARE(allstatevecs[i].vx - allstatevecs[j].vx) + DSQUARE(allstatevecs[i].vy - allstatevecs[j].vy) + DSQUARE(allstatevecs[i].vz - allstatevecs[j].vz));
	  vdist = sqrt(DSQUARE(allstatevecs[i].vx - allstatevecs[j].vx) + DSQUARE(allstatevecs[i].vy - allstatevecs[j].vy) + DSQUARE(allstatevecs[i].vz - allstatevecs[j].vz));
	  pdist = sqrt(DSQUARE(allstatevecs[i].x - allstatevecs[j].x) + DSQUARE(allstatevecs[i].y - allstatevecs[j].y) + DSQUARE(allstatevecs[i].z - allstatevecs[j].z));
	  if(dist>maxdist) {
	    maxdist = dist;
	  }
	  if(pdist>maxpdist) maxpdist = pdist;
	  if(vdist>maxvdist) maxvdist = vdist;
	}
      }
      // Now maxdist is the full 6D distance from point i to the most
      // distant other point. pdistmax and vdistmax are the 3D
      // position and velocity distances from point i to this same point.
      // Meanwhile, maxpdist is the 3D position distance from point i
      // to the other point that is most distance in 3D position space,
      // regardless of velocity, and maxvdist is the analogous value
      // for velocity. The reason pdistmax and vdistmax are distanct
      // (potentially) from maxpdist and maxvdist is that the point that
      // is most distant overall might not be the most distant in position
      // or velocity space individually.
      if(maxdist<minrad) {
	minrad = maxdist;
	pradmin = maxpdist;
	vradmin = maxvdist;
      }
      if(maxpdist<minprad) minprad=maxpdist;
      if(maxvdist<minvrad) minvrad=maxvdist;
    }
    outstream1  << fixed << setprecision(8) << radhyp[accelct].HelioRad << " " << radhyp[accelct].R_dot << " " << radhyp[accelct].R_dubdot << " " << heliodist[accelct] << " " << heliovel[accelct]/SOLARDAY << " " <<  helioacc[accelct]/SOLARDAY/SOLARDAY*1000.0 << " ";
    outstream1 << fixed << setprecision(2) << minrad << " " << vradmin << " " << pradmin << " " << minvrad << " " << minprad << " " << allstatevecs.size() << "\n";
  }

  outstream1.close();
  return(0);
}
