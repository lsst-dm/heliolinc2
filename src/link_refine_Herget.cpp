// November 04, 2022: link_refine_Herget.cpp:
// Like link_refine, but incorporates orbit-fitting using the
// Method of Herget, using the routine Hergertfit01, which does
// a downhill simplex fit in the 2-D space of geocentric distances
// at the instants of the first and last observations.
//
// September 14, 2022: link_refine.cpp (created on this date by copying
// an earlier prototype file called ancluster04a.cpp).
//
// Refine linkage files produced by heliolinc, eliminating cases
// of overlapping linkages by identifying sets of mutually overlapping
// linkages, choosing the best one in each such set using a quality
// metric, and discarding the rest. This enables heliolinc to be
// used aggressively to find all plausible linkages without concerns
// about redundancy, since link_refine is guaranteed to produce
// a non-redundant set with optimized quality metrics.
//
// The choice of the quality metric is obviously of great importance
// for the ultimate level of linkage success, because a badly
// chosen quality metric can cause a poor choice to be made from
// a selection of mutually overlapping linkages, such that a bad
// linkage is retained and good ones (if any) are rejected.
//
// This program uses a quality metric that involves a polynomial
// fit to the sky positions (RA and Dec) as a function of time.
// This is better than alternatives using only inferred state vectors,
// if and only the observations are all from one observatory.
// If observations are from multiple observatories, especially
// ones widely separated across the Earth, the alternative program
// link_refine_multisite should be used. It relies on inferred state
// vectors and hence is less powerful for a single site, but it
// is not confused by parallax with multiple sites. The current
// program is NOT suitable for use in combining data from multiple
// observing sites.
//
// June 28, 2022: ancluster04a.cpp:
// Like ancluster03c.cpp, but additionally fits a polynomial
// of order tracknightcount to the data, and includes the
// inverse residual of this polynomial in the quality metric.
//
// March 11, 2022: ancluster03c.cpp:
// Like ancluster03b, but reads input
// files with a new, standardized csv format.
// This requires both the cluster and the RMS file, where
// previously only the cluster file was required.
//
// February 18, 2022: ancluster03b.cpp: Like ancluster03a.cpp,
// but has some improvements to reduce memory requirements.
// Also, rejects all clusters that have multiple measurements
// at the same time, since this is impossible for a real object.
//
// January 17, 2022: ancluster03a.cpp
// Read cluster files produced by the good version of projectpairs04c.cpp,
// and analyze them. The good version of projectpairs04c.cpp fixes
// the omission of indices in the output cluster files. The old version
// of the current program, ancluster01a.cpp, had some
// time- and memory-intensive cleverness to cope with the absence
// of the indices. This new version requires the indices to be
// present, and uses them to proceed much more efficiently.
// This version differs from ancluster02a.cpp in that it sorts primarily
// on clusters rather than detections: it starts with the very best cluster
// (based on clustermetric), and rejects all overlapping clusters.
// Then it proceeds to the best surviving cluster, etc. This is superior
// to the detection-based analysis used by ancluster02a.cpp because
// the detection-based analysis implicitly assumed that all detections
// are good. The new, cluster-based analysis should work much better
// in cases where some detections are spurious.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MAXCLUSTRMS 1.0e5
#define DEBUG 1
#define FTOL_HERGET_SIMPLEX 1e-5L
#define ESCAPE_SCALE 0.99L // If the input velocity is above escape velocity, we
                           // scale it down by this amount.

static void show_usage()
{
  cerr << "Usage: link_refine_Herget -pairdet pairdet_file -lflist link_file_list -mjd mjdref -simptype simplex_type -ptpow point_num_exponent -nightpow night_num_exponent -timepow timespan_exponent -rmspow astrom_rms_exponent -maxrms maxrms -outfile outfile -outrms rmsfile -verbose verbosity\n\nOR, at minimum:\nlink_refine_Herget -pairdet pairdet_file -lflist link_file_list -mjd mjdref -outfile outfile -outrms rmsfile\n";
}

int main(int argc, char *argv[])
{
  int i=1;
  int j=1;
  int i1=0;
  int i2=0;
  int detct=0;
  string pairdetfile,clusterlist,clusterfile,rmsfile,outfile,outrmsfile,stest;
  int rmslinect=0;
  int clustlinect=0;
  string lnfromfile;
  int clusterfilenum=0;
  int clusterfilect=0;
  vector <string> clusternames;
  vector <string> rmsnames;
  long double mjdref=0L;
  string w1,w2,w3;
  int ptct=0;
  ifstream instream1;
  ifstream instream2;
  ofstream outstream1;
  ofstream outstream2;
  double maxrms = MAXCLUSTRMS;
  vector <det_obsmag_indvec> detvec;
  vector <det_obsmag_indvec> clusterdets;
  det_obsmag_indvec dsv = det_obsmag_indvec(0l,0l,0l,0L,0L,0L,"null",0.0,"V","I11",1,{});
  clusteran04 clustan = clusteran04(0,{},0,0.0,0,0,0.0,"null",{},{},{});
  long double X, Y, Z;
  X = Y = Z = 0L;
  long double MJD = 0L;
  double RA = 0;
  double Dec = 0;
  char detid[SHORTSTRINGLEN];
  double mag = 0;
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  long origind=0;
  vector <clusteran04> clustanvec;
  int clusterct=0;
  int clusterct2=0;
  int clusternum=0;
  vector <float> rmsvec;
  float posrms,velrms,totrms;
  posrms = velrms = totrms = 0;
  int pairnum=0;
  float timespan=0.0;
  int ptnum=0;
  int obsnights=0;
  float clustmetric=0;
  char rating[SHORTSTRINGLEN];
  vector <float> heliopar;
  float heliodist,heliovel,helioacc;
  heliodist = heliovel = helioacc = 0l;
  vector <double> statevecs;
  double x,y,z,vx,vy,vz;
  x=y=z=vx=vy=vz=0.0;
  vector <int> clustind;
  int startpoint=0;
  int endpoint=0;
  int detfilelinect=0;
  int ptpow = 1;      // Set the exponent to be used for the number of unique
                      // points in calculating the cluster quality metric.
  int nightpow = 1;   // Exponent for number of unique nights in quality metric.
  int timepow = 0;    // Exponent for total time span in quality metric.
  int rmspow = 2;     // Exponent for astrometric RMS (sqrt(variance)) in quality metric.
  
  float_index findex = float_index(0L,0);
  vector <float_index> clustanvec2;
  int goodclusternum=0;
  int istimedup=0;
  int badread=0;
  ldouble_index ldi = ldouble_index(0L,0);
  vector <ldouble_index> sortclust;

  // Variables for orbit fitting
  point3LD onepoint = point3LD(0.0L,0.0L,0.0L);
  vector <point3LD> observerpos;
  point3LD startpos = point3LD(0.0L,0.0L,0.0L);
  point3LD startvel = point3LD(0.0L,0.0L,0.0L);
  point3LD endpos = point3LD(0.0L,0.0L,0.0L);
  point3LD endvel = point3LD(0.0L,0.0L,0.0L);
  vector <long double> obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit;
  long double geodist1,geodist2, v_escape, v_helio, astromrms, chisq;
  long double ftol = FTOL_HERGET_SIMPLEX;
  long double simplex_scale = SIMPLEX_SCALEFAC;
  int simptype=0;
  int verbose=0;
  
  if(argc!=11 && argc!=13 && argc!=15 && argc!=17 && argc!=19 && argc!=21 && argc!=23 && argc!=25) {
    show_usage();
    return(1);
  }
  
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pd" || string(argv[i]) == "-pdet" || string(argv[i]) == "--pairdet" || string(argv[i]) == "--paireddetections" || string(argv[i]) == "--pairdetfile" || string(argv[i]) == "--pairdetections") {
      if(i+1 < argc) {
	//There is still something to read;
	pairdetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input pairdet file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-lflist" || string(argv[i]) == "-clist" || string(argv[i]) == "-cl" || string(argv[i]) == "-clustlist" || string(argv[i]) == "--clusterlist" || string(argv[i]) == "--clist" || string(argv[i]) == "--clusterl" || string(argv[i]) == "--clustlist") {
      if(i+1 < argc) {
	//There is still something to read;
	clusterlist=argv[++i];
	i++;
      }
      else {
	cerr << "Input cluster list keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-m" || string(argv[i]) == "-mjd" || string(argv[i]) == "-mjdref" || string(argv[i]) == "-MJD" || string(argv[i]) == "--mjd" || string(argv[i]) == "--MJD" || string(argv[i]) == "--modifiedjulianday" || string(argv[i]) == "--ModifiedJulianDay") {
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
    } else if(string(argv[i]) == "-of" || string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--fileout" || string(argv[i]) == "--fileoutput") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-simptype" || string(argv[i]) == "-stype" || string(argv[i]) == "-simpt" || string(argv[i]) == "-simplextype" || string(argv[i]) == "-simplex_type" || string(argv[i]) == "--simptype" || string(argv[i]) == "--simplex_type") {
      if(i+1 < argc) {
	//There is still something to read;
	simptype=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Simplex type keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-ptpow" || string(argv[i]) == "-ptexp" || string(argv[i]) == "-pointpow" || string(argv[i]) == "-pointnumpow" || string(argv[i]) == "-pointnum_exp" || string(argv[i]) == "--ptpow" || string(argv[i]) == "--pointnum_power") {
      if(i+1 < argc) {
	//There is still something to read;
	ptpow=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "keyword for ptpow, the power to which the number of unique points\nis raised when calcuating the quality metric\nsupplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-nightpow" || string(argv[i]) == "-nightexp" || string(argv[i]) == "-npow" || string(argv[i]) == "-nightnumpow" || string(argv[i]) == "-nightnum_exp" || string(argv[i]) == "--nightpow" || string(argv[i]) == "--nightnum_power") {
      if(i+1 < argc) {
	//There is still something to read;
	nightpow=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "keyword for nightpow, the power to which the number of unique\n observing nights is raised when calcuating the quality metric\nsupplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timepow" || string(argv[i]) == "-timeexp" || string(argv[i]) == "-timexp" || string(argv[i]) == "-tpow" || string(argv[i]) == "-timespan_pow" || string(argv[i]) == "--timepow" || string(argv[i]) == "--timespan_power") {
      if(i+1 < argc) {
	//There is still something to read;
	timepow=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "keyword for timepow, the power to which the total timespan is raised when calcuating\nthe quality metric supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rmspow" || string(argv[i]) == "-rmsexp" || string(argv[i]) == "-rmspower" || string(argv[i]) == "-rpow" || string(argv[i]) == "-astromrms_pow" || string(argv[i]) == "--rmspow" || string(argv[i]) == "--astrom_rms_power") {
      if(i+1 < argc) {
	//There is still something to read;
	rmspow=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "keyword for rmspow, the power to which the astrometric\nRMS residual is raised when calcuatingthe quality metric\nsupplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxrms" || string(argv[i]) == "-mr" || string(argv[i]) == "-mrms" || string(argv[i]) == "-rms" || string(argv[i]) == "-maxr" || string(argv[i]) == "--maxrms" || string(argv[i]) == "--maximumrms") {
      if(i+1 < argc) {
	//There is still something to read;
	maxrms=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Maximum RMS keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rms" || string(argv[i]) == "-sum" || string(argv[i]) == "-outsum" || string(argv[i]) == "-outrms" || string(argv[i]) == "-or" || string(argv[i]) == "--outrmsfile" || string(argv[i]) == "--outrms" || string(argv[i]) == "--rmsfile") {
      if(i+1 < argc) {
	//There is still something to read;
	outrmsfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output RMS file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--VERB" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--verbose") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "keyword for verbosity supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
    
  cout << "input paired detection file " << pairdetfile << "\n";
  cout << "input cluster list file " << clusterlist << "\n";
  cout << "Reference MJD: " << mjdref << "\n";
  cout << "Maximum RMS in km: " << maxrms << "\n";
  cout << "In calculating the cluster quality metric, the number of\nunique points will be raised to the power of " << ptpow << ";\n";
  cout << "the number of unique nights will be raised to the power of " << nightpow << ";\n";
  cout << "the total timespan will be raised to the power of " << timepow << ";\n";
  cout << "and the astrometric RMS will be raised to the power of (negative) " << rmspow << "\n";
  cout << "output cluster file " << outfile << "\n";
  cout << "output rms file " << outrmsfile << "\n";

  // Read paired detection file
  instream1.open(pairdetfile,ios_base::in);
  if(!instream1) {
    cerr << "can't open input file " << pairdetfile << "\n";
    return(1);
  }
  detvec={};
  // Skip header line
  getline(instream1,lnfromfile);
  cout << "Header line from input paired detection file " << pairdetfile << ":\n";
  cout << lnfromfile << "\n";
  // Read body of the file
  detfilelinect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the paired detections file, and load an object of class det_obsmag_indvec
    getline(instream1,lnfromfile);
    detfilelinect++;
    badread=0;
    //while(instream1 >> MJD >> RA >> Dec >> X >> Y >> Z >> detid >> mag >> band >> obscode >> origind) {
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
	cerr << "ERROR reading line " << detfilelinect << " of paired detection file" << pairdetfile << "\n";
	return(1);
      }
      // If we reach this point, the line was read OK. Write it to detvec.
      dsv=det_obsmag_indvec(MJD,RA,Dec,X,Y,Z,detid,mag,band,obscode,origind,{});
      dsv.indvec = {};
      detvec.push_back(dsv);
    } else if(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      cerr << "WARNING: line " << detfilelinect << " of paired detection file " << pairdetfile << " was too short\n";
    }
  }
  instream1.close();
  cout << "Successfully read " << detvec.size() << " detections from " << pairdetfile << "\n";
  
  instream1.open(clusterlist,ios_base::in);
  if(!instream1) {
    cerr << "can't open input file " << clusterlist << "\n";
    return(1);
  }
  while(instream1 >> clusterfile >> rmsfile)
    {
      if(clusterfile.size()>0 && rmsfile.size()>0) {
	clusternames.push_back(clusterfile);
	rmsnames.push_back(rmsfile);
      }	
    }
  instream1.close();
  clusterfilenum = clusternames.size();
  cout << "Read names of " << clusterfilenum << " cluster files from " << clusterlist << "\n";
  
  outstream1.open(outfile,ios_base::out);
  outstream2.open(outrmsfile,ios_base::out);
  outstream1 << "#ptct,MJD,RA,Dec,idstring,mag,band,obscode,index1,index2,clusternum\n";
  outstream2 << "#clusternum,posRMS,astromRMS,totRMS,pairnum,timespan,uniquepoints,obsnights,metric,rating,heliodist,heliovel,helioacc,posX,posY,posZ,velX,velY,velZ,orbit_a,orbit_e,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count\n";

  // Read cluster files, loading clusters.
  for(clusterfilect=0; clusterfilect<clusterfilenum; clusterfilect++) {
    instream1.open(clusternames[clusterfilect],ios_base::in);
    instream2.open(rmsnames[clusterfilect],ios_base::in);
    if(!instream1) {
      cerr << "can't open cluster file " << clusternames[clusterfilect] << "\n";
      return(2);
    } else cout << "Reading cluster file " << clusternames[clusterfilect] << "\n";
    if(!instream2) {
      cerr << "can't open rms file " << rmsnames[clusterfilect] << "\n";
      return(2);
    } else cout << "Reading rms file " << rmsnames[clusterfilect] << "\n";
    // Skip header lines
    getline(instream1,lnfromfile);
    getline(instream2,lnfromfile);
    rmslinect=clustlinect=0;
    while(!instream1.bad() && !instream1.fail() && !instream1.eof() && !instream2.bad() && !instream2.fail() && !instream2.eof()) {
      // Read a line from the rms file, and load an object of type clusteran04
      getline(instream2,lnfromfile);
      rmslinect++;
      badread=0;
      if(lnfromfile.size()>40 && !instream2.bad() && !instream2.fail() && !instream2.eof()) {
	// Read cluster index number;
	startpoint=0;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) clusterct = stoi(stest);
	else badread=1;
 	// Read three RMS values
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) posrms = stof(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) velrms = stof(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) totrms = stof(stest);
	else badread=1;

	// Read the integer pairnum
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) pairnum = stoi(stest);
	else badread=1;
	
	// Read the float timespan
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) timespan = stof(stest);
	else badread=1;
	
	// Read the integers ptnum and obsnights
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) ptnum = stoi(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) obsnights = stoi(stest);
	else badread=1;

	//Read and discard a value for clustermetric
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	
	// Recalculate the float clustermetric
	clustmetric = intpowD(double(ptnum),ptpow)*intpowD(double(obsnights),nightpow)*intpowD(timespan,timepow);
	// Note that the value of clustermetric just calculated
	// will later be divided by the reduced chi-square value of the
	// astrometric fit, before it is ultimately used as a selection criterion.
	
	// read the string rating
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) stringncopy01(rating,stest,SHORTSTRINGLEN);

	// Read three heliocentric hypothesis parameters.
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) heliodist = stof(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) heliovel = stof(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) helioacc = stof(stest);
	else badread=1;
	
	// Read the six elements of the mean state vector
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) x = stof(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) y = stof(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) z = stof(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) vx = stof(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) vy = stof(stest);
	else badread=1;
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) vz = stof(stest);
	else badread=1;

	// Note: the detection index vector is left empty, for now.
     
	// If there was a file read error, abort.
	if(badread==1) {
	  cerr << "ERROR reading line " << rmslinect << " of rms file " << rmsnames[clusterfilect] << "\n";
	  return(1);
	}
      
	// Load vectors
	rmsvec = {};
	rmsvec.push_back(posrms);
	rmsvec.push_back(velrms);
	rmsvec.push_back(totrms); // rmsvec contents: [0] position RMS scatter,
	                          // [1] velocity RMS scatter, [2] total RMS scatter.
	heliopar = {};
	heliopar.push_back(heliodist);
	heliopar.push_back(heliovel);
	heliopar.push_back(helioacc);
	statevecs = {};
	statevecs.push_back(x);
	statevecs.push_back(y);
	statevecs.push_back(z);
	statevecs.push_back(vx);
	statevecs.push_back(vy);
	statevecs.push_back(vz);
	
	clustan = clusteran04(clusterct,rmsvec,pairnum,timespan,ptnum,obsnights,clustmetric,rating,heliopar,statevecs,{});

	// Read the associated lines from the cluster file, but retain only the
	// indices to the det vector.
	for(i=0;i<ptnum;i++) {
	  if (!instream1.bad() && !instream1.fail() && !instream1.eof()) {
	    // Read a line from the cluster file.
	    getline(instream1,lnfromfile);
	    clustlinect++;
	    while(lnfromfile.size()<40 && !instream1.bad() && !instream1.fail() && !instream1.eof()) {
	      cerr << "WARNING: line " << clustlinect << " of cluster file " << clusternames[clusterfilect] << " is too short\n";
	      // Read another line, maybe there's just a blank one.
	      getline(instream1,lnfromfile);
	      clustlinect++;
	    }
	    badread=0;
	    if(lnfromfile.size()>40 && !instream1.bad() && !instream1.fail() && !instream1.eof()) {
	      // Read and discard the first eight quantities: ptct, MJD, RA, Dec, mag, band, and obscode.
	      startpoint=0;
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) startpoint = endpoint+1;
	      else badread=1;
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) startpoint = endpoint+1;
	      else badread=1;
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) startpoint = endpoint+1;
	      else badread=1;
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) startpoint = endpoint+1;
	      else badread=1;
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) startpoint = endpoint+1;
	      else badread=1;
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) startpoint = endpoint+1;
	      else badread=1;
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) startpoint = endpoint+1;
	      else badread=1;
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) startpoint = endpoint+1;
	      else badread=1;
	      // Read the essential quantity: the index to the detection vector
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) {
		startpoint = endpoint+1;
		i1 = stoi(stest);
		clustan.clustind.push_back(i1);
	      } else badread=1;
	      // Read and discard index2, which is currently irrelevant.
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) startpoint = endpoint+1;
	      else badread=1;
	      // Read clusterct, and check it against the index read from the rms file
	      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	      if(endpoint>0) {
		startpoint = endpoint+1;
		i2 = stoi(stest);
		if(i2 != clustan.clusterct) {
		  cerr << "ERROR: cluster count mismatch at rms line " << rmslinect << ", cluster file line " << clustlinect << ": " << i2 << " != " << clustan.clusterct << "\n";
		  return(1);
		}
	      } else badread=1;
	      // If there was a file read error, abort.
	      if(badread==1) {
		cerr << "ERROR reading line " << clustlinect << " of cluster file " << clusternames[clusterfilect] << "\n";
		return(1);
	      }
	    }
	    // Close if-statement checking if we're still reading valid lines from cluster file.
	  } else {
	      cerr << "WARNING: unsuccessful read at line " << clustlinect << " of cluster file " << clusternames[clusterfilect] << "\n";
	  }
	  // Close loop over points in the cluster
	}
	// See if it's a good cluster
	if(clustan.rmsvec[2]<=maxrms) {
	  // Check for duplicate MJDs
	  clusterdets={};
	  for(ptct=0; ptct<ptnum; ptct++) {
	    i1 = clustan.clustind[ptct];
	    clusterdets.push_back(detvec[i1]);
	  }
	  sort(clusterdets.begin(), clusterdets.end(), early_det_obsmag_indvec());

	  istimedup=0;
	  for(ptct=1; ptct<ptnum; ptct++) {
	    if(clusterdets[ptct-1].MJD == clusterdets[ptct].MJD && stringnmatch01(clusterdets[ptct-1].obscode,clusterdets[ptct].obscode,3)==0) istimedup=1;
	  }
	  if(istimedup==0) {
	    // The cluster is good.
	    // Perform orbit fitting using the method of Herget, to get astrometric residuals
	    // Load observational vectors
	    observerpos = {};
	    obsMJD = obsRA = obsDec = sigastrom = fitRA = fitDec = resid = orbit = {};
	    for(ptct=0; ptct<ptnum; ptct++) {
	      obsMJD.push_back(clusterdets[ptct].MJD);
	      obsRA.push_back(clusterdets[ptct].RA);
	      obsDec.push_back(clusterdets[ptct].Dec);
	      sigastrom.push_back(1.0L);
	      onepoint = point3LD(clusterdets[ptct].x,clusterdets[ptct].y,clusterdets[ptct].z);
	      observerpos.push_back(onepoint);
	    }
	    // Use mean state vectors to estimate positions
	    startpos.x = clustan.statevecs[0];
	    startpos.y = clustan.statevecs[1];
	    startpos.z = clustan.statevecs[2];
	    startvel.x = clustan.statevecs[3];
	    startvel.y = clustan.statevecs[4];
	    startvel.z = clustan.statevecs[5];
	    // Check if the velocity is above escape
	    v_escape = sqrt(2.0L*GMSUN_KM3_SEC2/vecabs3LD(startpos));
	    v_helio = vecabs3LD(startvel);
	    if(v_helio>=v_escape) {
	      cerr << "WARNING: mean state vector velocity was " << v_helio/v_escape << " times higher than solar escape\n";
	      startvel.x *= ESCAPE_SCALE*v_escape/v_helio;
	      startvel.y *= ESCAPE_SCALE*v_escape/v_helio;
	      startvel.z *= ESCAPE_SCALE*v_escape/v_helio;
	    }
	    // Calculate position at first observation
	    Keplerint(GMSUN_KM3_SEC2, mjdref, startpos, startvel, obsMJD[0], endpos, endvel);
	    // Find vector relative to the observer by subtracting off the observer's position.
	    endpos.x -= observerpos[0].x;
	    endpos.y -= observerpos[0].y;
	    endpos.z -= observerpos[0].z;
	    geodist1 = vecabs3LD(endpos)/AU_KM;
	    // Calculate position at last observation
	    Keplerint(GMSUN_KM3_SEC2, mjdref, startpos, startvel, obsMJD[ptnum-1], endpos, endvel);
	    endpos.x -= observerpos[ptnum-1].x;
	    endpos.y -= observerpos[ptnum-1].y;
	    endpos.z -= observerpos[ptnum-1].z;
	    geodist2 = vecabs3LD(endpos)/AU_KM;
	    simplex_scale = SIMPLEX_SCALEFAC;
	    chisq = Hergetfit01(geodist1, geodist2, simplex_scale, simptype, ftol, 1, ptnum, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
	    // orbit vector contains: semimajor axis [0], eccentricity [1],
	    // mjd at epoch [2], the state vectors [3-8], and the number of
	    // orbit evaluations (~iterations) required to reach convergence [9].

	    chisq /= (long double)ptnum; // Now it's the reduced chi square value
	    astromrms = sqrt(chisq); // This gives the actual astrometric RMS in arcseconds if all the
	                             // entries in sigastrom are 1.0. Otherwise it's a measure of the
	                             // RMS in units of the typical uncertainty.
	    // Include this astrometric RMS value in the cluster metric and the RMS vector
	    clustan.rmsvec.push_back(astromrms); // rmsvec[3]: astrometric rms in arcsec.
	    clustan.clustermetric /= intpowD(astromrms,rmspow); // Under the default value rmspow=2, this is equivalent
	                                                        // to dividing by the chi-square value rather than just
	                                                        // the astrometric RMS, which has the desireable effect of
	                                                        // prioritizing low astrometric error even more.
	    for(i=0;i<10;i++) clustan.statevecs.push_back(orbit[i]);
	    // statevecs now contains original heliolinc state vecs [0-5],
	    // Herget fit semimajor axis [6], eccentricity [7], MJD at epoch [8],
	    // state vectors [9-14], and number of orbit evaluations (~iterations)
	    // required for fit [15].
	    // Add cluster indices to detvec.
	    for(ptct=0; ptct<ptnum; ptct++) {
	      i1 = clustan.clustind[ptct];
	      detvec[i1].indvec.push_back(clustanvec.size());
	    }
	    // Push new cluster on to clustanvec.
	    clustanvec.push_back(clustan);
	  }
	  // Close if-statement checking the RMS was low enough.
	}	
	// Close if-statement checking if RMS line was long enough
      } else if(!instream2.bad() && !instream2.fail() && !instream2.eof()) {
	cerr << "WARNING: line " << rmslinect << " of rms file " << rmsnames[clusterfilect] << " is too short\n";
      }
      // Close loop checking if we're still successfully reading both files
    }
    instream1.close();
    instream2.close();
    // Close loop on input cluster files.
  }

  clusternum = clustanvec.size();
  cout << "Found " << clusternum << " clusters\n";
  // Load just clustermetric values and indices from
  // clustanvec into clustanvec2
  clustanvec2 = {};
  // Record indices so information won't be lost on sort
  for(clusterct=0; clusterct<clusternum; clusterct++) {
    findex = float_index(clustanvec[clusterct].clustermetric,clusterct);
    clustanvec2.push_back(findex);
  }
  // Sort clustanvec2
  sort(clustanvec2.begin(), clustanvec2.end(), lower_float_index());
  
  // Loop on clusters, starting with the best (highest metric),
  // and eliminating duplicates
  goodclusternum=0;
  for(clusterct2=clusternum-1; clusterct2>=0; clusterct2--) {
    clusterct = clustanvec2[clusterct2].index;
    if(clustanvec[clusterct].uniquepoints>=6 && clustanvec[clusterct].rmsvec[2]<=maxrms) {
      // This is a good cluster not already marked as used.
      goodclusternum++;
      cout << "Accepted good cluster " << goodclusternum << " with metric " << clustanvec[clusterct].clustermetric << "\n";
      // See whether cluster is pure or mixed.
      stringncopy01(rating,"PURE",SHORTSTRINGLEN);
      for(i=0; i<long(clustanvec[clusterct].clustind.size()); i++) {
	if(i>0 && stringnmatch01(detvec[clustanvec[clusterct].clustind[i]].idstring,detvec[clustanvec[clusterct].clustind[i-1]].idstring,SHORTSTRINGLEN) != 0) {
	  stringncopy01(rating,"MIXED",SHORTSTRINGLEN);
	}
      }
      // Figure out the time order of cluster points, so we can write them out in order.
      sortclust = {};
      for(i=0; i<long(clustanvec[clusterct].clustind.size()); i++) {
	i1 = clustanvec[clusterct].clustind[i];
	ldi = ldouble_index(detvec[i1].MJD,i1);
	sortclust.push_back(ldi);
      }
      sort(sortclust.begin(), sortclust.end(), lower_ldouble_index());
   
      // Write all individual detections in this cluster to the output cluster file
      for(i=0; i<long(clustanvec[clusterct].clustind.size()); i++) {
	i1 = sortclust[i].index;
	outstream1  << fixed << setprecision(6) << intzero01i(i,4) << "," << detvec[i1].MJD << "," << detvec[i1].RA << "," << detvec[i1].Dec << "," << detvec[i1].idstring << ",";
	outstream1  << fixed << setprecision(3) << detvec[i1].mag << "," << detvec[i1].band << "," << detvec[i1].obscode << "," << i1 << "," << detvec[i1].index << "," << goodclusternum << "\n";
      }
      // Write summary line to rms file
      outstream2  << fixed << setprecision(3) << goodclusternum << "," << clustanvec[clusterct].rmsvec[0] << "," << clustanvec[clusterct].rmsvec[3] << "," << clustanvec[clusterct].rmsvec[2] << "," << clustanvec[clusterct].pairnum << ",";
      outstream2  << fixed << setprecision(6) << clustanvec[clusterct].timespan << "," << clustanvec[clusterct].clustind.size() << "," << clustanvec[clusterct].obsnights  << "," << clustanvec[clusterct].clustermetric << "," << clustanvec[clusterct].rating << ",";
      outstream2  << fixed << setprecision(3) << clustanvec[clusterct].heliopar[0] << "," << clustanvec[clusterct].heliopar[1] << ",";
      outstream2  << fixed << setprecision(6) << clustanvec[clusterct].heliopar[2] << ",";
      // Original heliolinc position state vector
      outstream2  << fixed << setprecision(3)  << clustanvec[clusterct].statevecs[0] << "," << clustanvec[clusterct].statevecs[1] << "," << clustanvec[clusterct].statevecs[2] << ",";
      // Original heliolinc velocity state vector
      outstream2  << fixed << setprecision(6) << clustanvec[clusterct].statevecs[3] << ","   << clustanvec[clusterct].statevecs[4] << "," << clustanvec[clusterct].statevecs[5] << ",";
      // Orbit fit semimajor axis
      outstream2  << fixed << setprecision(6)  << clustanvec[clusterct].statevecs[6]/AU_KM << ",";
      // Orbit fit eccentricity
      outstream2 << fixed << setprecision(6) << clustanvec[clusterct].statevecs[7] << ",";
      // Orbit fit MJD at epoch
      outstream2 << fixed << setprecision(9) << clustanvec[clusterct].statevecs[8] << ",";
      // Orbit fit position state vector
      outstream2 << fixed << setprecision(3) << clustanvec[clusterct].statevecs[9] << "," << clustanvec[clusterct].statevecs[10] << "," << clustanvec[clusterct].statevecs[11] << ",";
      // Orbit fit velocity state vector
      outstream2 << fixed << setprecision(6) << clustanvec[clusterct].statevecs[12] << "," << clustanvec[clusterct].statevecs[13] << "," << clustanvec[clusterct].statevecs[14] << ",";
      // Total number of orbit evaluations (~iterations) used in orbit fit
      outstream2 << fixed << setprecision(0) << clustanvec[clusterct].statevecs[15] << "\n";
      // FINISHED WRITING SUMMARY LINE TO RMS/SUMMARY FILE
      
      for(i=0;i<long(clustanvec[clusterct].clustind.size());i++) {
	// This point is in the chosen cluster, and cannot be in any other
	detct = clustanvec[clusterct].clustind[i];
	//cout << "Deduplicating detection " << detct << "\n";
	// Mark all the clusters containing this point as already considered
	for(j=0;j<long(detvec[detct].indvec.size());j++) {
	  //cout << "         wiping cluster " << detvec[detct].indvec[j] << "\n";
	  clustanvec[detvec[detct].indvec[j]].uniquepoints = 0;
	}
      }
    } 
  }
  return(0);
}
