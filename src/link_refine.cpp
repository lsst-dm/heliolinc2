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
#define ONE_POINT_PER_IMAGE 1
#define DEBUG 1


static void show_usage()
{
  cerr << "Usage: link_refine -pairdet pairdet_file -lflist link_file_list -maxrms maxrms -outfile outfile -outrms rmsfile\n";
}
    
int main(int argc, char *argv[])
{
  int i=1;
  int j=1;
  int i1=0;
  int i2=0;
  int detnum=0;
  int detct=0;
  char c=0;
  string pairdetfile,clusterlist,clusterfile,rmsfile,outfile,outrmsfile,stest;
  int rmslinect=0;
  int clustlinect=0;
  string lnfromfile;
  int clusterfilenum=0;
  int clusterfilect=0;
  vector <string> clusternames;
  vector <string> rmsnames;
  string w1,w2,w3;
  int ptct=0;
  ifstream instream1;
  ifstream instream2;
  ofstream outstream1;
  ofstream outstream2;
  double tmjd;
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
  int pairnum;
  float timespan;
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
  point3d unitvec = point3d(0l,0l,0l);
  point3d testvec = point3d(0l,0l,0l);
  vector <double> timevec;
  vector <double> sigvec;
  vector <double> xvec;
  vector <double> yvec;
  vector <double> zvec;
  vector <double> avecx;
  vector <double> avecy;
  vector <double> avecz;
  vector <double> xmod;
  vector <double> ymod;
  vector <double> zmod;
  int polyorder = 2;
  double astromrms,norm;
  
  float_index findex = float_index(0L,0);
  vector <float_index> clustanvec2;
  int goodclusternum=0;
  int istimedup=0;
  int badread=0;
  ldouble_index ldi = ldouble_index(0L,0);
  vector <ldouble_index> sortclust;
  

  if(argc!=9 && argc!=11) {
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
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
    
  cout << "input paired detection file " << pairdetfile << "\n";
  cout << "input cluster list file " << clusterlist << "\n";
  cout << "Maximum RMS in km: " << maxrms << "\n";
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
  outstream2 << "#clusternum,posRMS,velRMS,totRMS,pairnum,timespan,uniquepoints,obsnights,metric,rating,heliodist,heliovel,helioacc,posX,posY,posZ,velX,velY,velZ\n";

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
	clustmetric = double(ptnum)*double(obsnights)*timespan/totrms;
	
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
	rmsvec.push_back(totrms);
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
		  cerr << "ERROR: cluster count mismatch at rms line " << rmslinect << ", cluster file line " << clustlinect << "\n";
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
	    // Perform polynomial fitting, to get astrometric residuals
	    timevec = {};
	    sigvec = {};
	    xvec = {};
	    yvec = {};
	    zvec = {};
	    polyorder = clustan.obsnights;
	    for(ptct=0; ptct<ptnum; ptct++) {
	      timevec.push_back(clusterdets[ptct].MJD - clusterdets[0].MJD);
	      sigvec.push_back(1.0);
	      celestial_to_stateunit(clusterdets[ptct].RA, clusterdets[ptct].Dec, unitvec);
	      xvec.push_back(unitvec.x);
	      yvec.push_back(unitvec.y);
	      zvec.push_back(unitvec.z);
	    }
	    // Perform fits to x, y, and z as a function of time
	    polyfit01(xvec, sigvec, timevec, ptnum, polyorder, avecx);
	    polyfit01(yvec, sigvec, timevec, ptnum, polyorder, avecy);
	    polyfit01(zvec, sigvec, timevec, ptnum, polyorder, avecz);
	    // Calculate the astrometric RMS relative to the polynomial fit
	    astromrms = norm = 0l;
	    for(ptct=0; ptct<ptnum; ptct++) {
	      // Use fit coefficients to generate a model vector.
	      testvec.x = avecx[0];
	      testvec.y = avecy[0];
	      testvec.z = avecz[0];
	      for(i=1;i<=polyorder;i++) {
		testvec.x += avecx[i]*intpowD(timevec[ptct], i);
		testvec.y += avecy[i]*intpowD(timevec[ptct], i);
		testvec.z += avecz[i]*intpowD(timevec[ptct], i);
	      }
	      // Normalize the model vector
	      vecnorm3d(testvec);
	      // Reconstruct projected data vector
	      unitvec.x = xvec[ptct];
	      unitvec.y = yvec[ptct];
	      unitvec.z = zvec[ptct];	      
	      // Calculate the difference between the model vector
	      // and the data
	      testvec.x -= unitvec.x;
	      testvec.y -= unitvec.y;
	      testvec.z -= unitvec.z;
	      astromrms += DSQUARE(vecabs3d(testvec)*ASECPRAD);
	      norm += 1.0l;
	    }
	    astromrms = sqrt(astromrms/norm);
	    // Include this astrometric RMS value in the cluster metric and the RMS vector
	    clustan.rmsvec.push_back(astromrms);
	    clustan.clustermetric /= astromrms;
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
      for(i=0; i<clustanvec[clusterct].clustind.size(); i++) {
	if(i>0 && stringnmatch01(detvec[clustanvec[clusterct].clustind[i]].idstring,detvec[clustanvec[clusterct].clustind[i-1]].idstring,SHORTSTRINGLEN) != 0) {
	  stringncopy01(rating,"MIXED",SHORTSTRINGLEN);
	}
      }
      // Figure out the time order of cluster points, so we can write them out in order.
      sortclust = {};
      for(i=0; i<clustanvec[clusterct].clustind.size(); i++) {
	i1 = clustanvec[clusterct].clustind[i];
	ldi = ldouble_index(detvec[i1].MJD,i1);
	sortclust.push_back(ldi);
      }
      sort(sortclust.begin(), sortclust.end(), lower_ldouble_index());
   
      // Write all individual detections in this cluster to the output cluster file
      for(i=0; i<clustanvec[clusterct].clustind.size(); i++) {
	i1 = sortclust[i].index;
	outstream1  << fixed << setprecision(6) << intzero01i(i,4) << "," << detvec[i1].MJD << "," << detvec[i1].RA << "," << detvec[i1].Dec << "," << detvec[i1].idstring << ",";
	outstream1  << fixed << setprecision(3) << detvec[i1].mag << "," << detvec[i1].band << "," << detvec[i1].obscode << "," << i1 << "," << detvec[i1].index << "," << goodclusternum << "\n";
      }
      // Write summary line to rms file
      outstream2  << fixed << setprecision(3) << goodclusternum << "," << clustanvec[clusterct].rmsvec[0] << "," << clustanvec[clusterct].rmsvec[3] << "," << clustanvec[clusterct].rmsvec[2] << "," << clustanvec[clusterct].pairnum << ",";
      outstream2  << fixed << setprecision(6) << clustanvec[clusterct].timespan << "," << clustanvec[clusterct].clustind.size() << "," << clustanvec[clusterct].obsnights  << "," << clustanvec[clusterct].clustermetric << "," << clustanvec[clusterct].rating << ",";
      outstream2  << fixed << setprecision(3) << clustanvec[clusterct].heliopar[0] << "," << clustanvec[clusterct].heliopar[1] << ",";
      outstream2  << fixed << setprecision(6) << clustanvec[clusterct].heliopar[2] << ",";
      outstream2  << fixed << setprecision(3)  << clustanvec[clusterct].statevecs[0] << "," << clustanvec[clusterct].statevecs[1] << "," << clustanvec[clusterct].statevecs[2] << ",";
      outstream2  << fixed << setprecision(6) << clustanvec[clusterct].statevecs[3] << ","   << clustanvec[clusterct].statevecs[4] << "," << clustanvec[clusterct].statevecs[5] << "\n";
      for(i=0;i<clustanvec[clusterct].clustind.size();i++) {
	// This point is in the chosen cluster, and cannot be in any other
	detct = clustanvec[clusterct].clustind[i];
	//cout << "Deduplicating detection " << detct << "\n";
	// Mark all the clusters containing this point as already considered
	for(j=0;j<detvec[detct].indvec.size();j++) {
	  //cout << "         wiping cluster " << detvec[detct].indvec[j] << "\n";
	  clustanvec[detvec[detct].indvec[j]].uniquepoints = 0;
	}
      }
    } 
  }
  return(0);
}
