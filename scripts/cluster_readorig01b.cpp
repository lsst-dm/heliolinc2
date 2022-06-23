// March 15, 2022: cluster_readorig01b.cpp:
// Like cluster_readorig01a, but reads files with the new, header+csv formatting,
// such as those output by projectpairs06c.cpp, ancluster03c.cpp, or related
// programs. Note that this means that an rms file is required as well as
// a cluster file.
// 
// March 04, 2022: cluster_readorig01a.cpp
// Read cluster files produced by the good version of projectpairs05a.cpp,
// ancluster03b.cpp, or related programs using the same format, and use
// the recorded line numbers to pull input lines from the original input
// file, and create a new file giving all the clusters in terms of the
// orginal input lines.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MAXCLUSTRMS 1.0e5
#define ONE_POINT_PER_IMAGE 1


static void show_usage()
{
  cerr << "Usage: cluster_readorig01b -orig original_input -clust clusterfile -rms rmsfile -outfile outfile\n";
}
    
int main(int argc, char *argv[])
{
  int i=1;
  int linect=0;
  int linenum=0;
  string origfile,clusterfile,rmsfile,outfile,rating,stest;
  string lnfromfile;
  string w1,w2;
  double posrms,velrms,totrms;
  posrms = velrms = totrms = 0l;
  double heliodist, heliovel, helioacc;
  heliodist = heliovel = helioacc = 0l;
  double x,y,z,vx,vy,vz;
  x = y = z = vz = vy = vz = 0l;
  vector <string> inputlines;
  ifstream instream1;
  ifstream instream2;
  ofstream outstream1;
  vector <double> mjdvec;
  float clustmetric=0;
  int ptnum=0;
  int pairnum=0;
  int ptct=0;
  int nightobs=0;
  double timespan=0l;
  int obsnights=0;
  int globalct = 0;
  int ijunk1 = 0;
  int ijunk2 = 0;
  long i1 = 0;
  long i2 = 0;
  double mjd,RA,Dec;
  mjd = RA = Dec = 0l;
  string detid;
  int rmslinect,clustlinect;
  int badread=0;
  int startpoint=0;
  int endpoint=0;

  if(argc<9) {
    show_usage();
    return(1);
  }
  
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-orig" || string(argv[i]) == "-o" || string(argv[i]) == "-ori" || string(argv[i]) == "--original" || string(argv[i]) == "-org") {
      if(i+1 < argc) {
	//There is still something to read;
        origfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input original file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-cf" || string(argv[i]) == "-c" || string(argv[i]) == "-clust" || string(argv[i]) == "--clusterfile" || string(argv[i]) == "--cfile" || string(argv[i]) == "-cfile" || string(argv[i]) == "--clustfile") {
      if(i+1 < argc) {
	//There is still something to read;
	clusterfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input cluster file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-rf" || string(argv[i]) == "-r" || string(argv[i]) == "-rms" || string(argv[i]) == "-rmsfile" || string(argv[i]) == "--rmsfile" || string(argv[i]) == "-rfile" || string(argv[i]) == "--rfile") {
      if(i+1 < argc) {
	//There is still something to read;
	rmsfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input cluster file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-of" || string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--fileout" || string(argv[i]) == "--fileoutput") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
      }
    }
  }
    
  cout << "input original data file " << origfile << "\n";
  cout << "input cluster file " << clusterfile << "\n";
  cout << "input rms file " << rmsfile << "\n";
  cout << "output cluster file " << outfile << "\n";

  // Read original input file into a string vector
  instream1.open(origfile,ios_base::in);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << origfile << "\n";
    return(1);
  }
  linect=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    inputlines.push_back(lnfromfile);
    linect++;
  }
  linenum = linect;
  cout << "File " << origfile << " appears to have " << linenum << " lines\n";
  instream1.close();

  // Open cluster file
  instream1.open(clusterfile,ios_base::in);
  instream2.open(rmsfile,ios_base::in);
  if(!instream1) {
    cerr << "can't open cluster file " << clusterfile << "\n";
    return(2);
  } else cout << "Reading cluster file " << clusterfile << "\n";
  if(!instream2) {
    cerr << "can't open rms file " << rmsfile << "\n";
    return(2);
  } else cout << "Reading rms file " << rmsfile << "\n";
  // Skip header lines
  getline(instream1,lnfromfile);
  getline(instream2,lnfromfile);
  rmslinect=clustlinect=0;

  // Open output file
  outstream1.open(outfile,ios_base::out);

  rmslinect=clustlinect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof() && !instream2.bad() && !instream2.fail() && !instream2.eof()) {
    // Read a line from the rms file
    getline(instream2,lnfromfile);
    rmslinect++;
    badread=0;
    if(lnfromfile.size()>40) {
      // Read cluster index number;
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) globalct = stoi(stest);
      else badread=1;
      cout << globalct << " = globalct\n";
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
	
      // Read the double timespan
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) timespan = stod(stest);
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
	
      // Read the float clustmetric
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) clustmetric = stof(stest);
      else badread=1;
	
      // read the string rating
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) rating=stest;

      // Read three heliocentric hypothesis parameters.
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) heliodist = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) heliovel = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) helioacc = stod(stest);
      else badread=1;
	
      // Read the six elements of the mean state vector
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) x = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) y = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) z = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) vx = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) vy = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) vz = stod(stest);
      else badread=1;
     
      // If there was a file read error, abort.
      if(badread==1) {
	cerr << "ERROR reading line " << rmslinect << " of rms file " << rmsfile << "\n";
	return(1);
      }
      
      // Read the associated lines from the cluster file, but retain only the
      // indices to the original input file
      for(i=0;i<ptnum;i++) {
	if (!instream1.bad() && !instream1.fail() && !instream1.eof()) {
	  // Read a line from the cluster file.
	  getline(instream1,lnfromfile);
	  clustlinect++;
	  while(lnfromfile.size()<40 && !instream1.bad() && !instream1.fail() && !instream1.eof()) {
	    cerr << "WARNING: line " << clustlinect << " of cluster file " << clusterfile << " is too short\n";
	    // Read another line, maybe there's just a blank one.
	    getline(instream1,lnfromfile);
	    clustlinect++;
	  }
	  badread=0;
	  if(lnfromfile.size()>40) {
	    // Read and discard the first nine quantities: ptct, MJD, RA, Dec, mag, band, obscode, and i1.
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
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    // Read the essential quantity: the index to the original file
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) {
	      startpoint = endpoint+1;
	      i1 = stoi(stest);
	    } else badread=1;
	    cout << i1 << " = origindex\n";
	    // Read clusterct, and check it against the index read from the rms file
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) {
	      startpoint = endpoint+1;
	      i2 = stoi(stest);
	      cout << i2 << " = cluster index\n";
	      if(i2 != globalct) {
		cerr << "ERROR: cluster count mismatch at rms line " << rmslinect << ", cluster file line " << clustlinect << "\n";
		return(1);
	      } else {
		// All good. Right this detection to the output file.
		outstream1 << inputlines[i1-1] << ",";
		// Additional output: globalct, heliodist, heliovel, helioacc, ptnum, timespan, nightobs, totrms, clustermetric
		outstream1 << fixed << setprecision(6) << globalct << "," << heliodist << "," << heliovel << "," << helioacc << "," << ptnum << ",";
		outstream1 << fixed << setprecision(6) << timespan << "," << nightobs << "," << totrms << "," << clustmetric << "\n";
	      }
	    } else badread=1;
	    // If there was a file read error, abort.
	    if(badread==1) {
	      cerr << "ERROR reading line " << clustlinect << " of cluster file " << clusterfile << "\n";
	      return(1);
	    }
	  }
	  // Close if-statement checking if we're still reading valid lines from cluster file.
	} else {
	  cerr << "WARNING: unsuccessful read at line " << clustlinect << " of cluster file " << clusterfile << "\n";
	}
	// Close loop over points in the cluster
      }
    }
  }
  instream1.close();
  instream2.close();
  outstream1.close();

  return(0);
}
