// February 02, 2023: mpcobs_cull01.cpp:
// Cull down a file of observations of minor planets so that it contains
// only observations within a specific range of dates (MJDs).
// Observations corresponding to space observatories such as C51 are
// also automatically eliminated.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: mpcobs_match -mpcfile mpc_file -mjdrange MJDstart MJDend -mpcout mpc_output\n";
}

int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  char c='0';
  int reachedeof=0;
  int lct=0;
  string mpcfile,mpcout;
  string stest;
  double mjdstart=0.0;
  double mjdend=1.0e7;
  ifstream instream1;
  ofstream outstream1;
  string lnfromfile;
  string obscode;
  double MJD;
  vector <string> badcodes;
  int isbadcode=0;
  
  badcodes.push_back("C51");
  
  if(argc<0 && argc<1) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-mpcfile" || string(argv[i]) == "-mpcin" || string(argv[i]) == "-mpcobs" || string(argv[i]) == "--mpcobs" || string(argv[i]) == "-mpcinfile" || string(argv[i]) == "--mpcfile" || string(argv[i]) == "--mpcin" || string(argv[i]) == "--mpcinfile") {
      if(i+1 < argc) {
	//There is still something to read;
	mpcfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input mpc file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdrange" || string(argv[i]) == "-MJDrange" || string(argv[i]) == "-mjd" || string(argv[i]) == "-MJD" || string(argv[i]) == "--mjdrange" || string(argv[i]) == "--MJDrange" || string(argv[i]) == "--MJDRANGE" || string(argv[i]) == "--MJD") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdstart = stod(argv[++i]);
      }
      else {
	cerr << "Input MJD range keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      } if(i+1 < argc) {
	//There is still something to read;
	mjdend = stod(argv[++i]);
	i++;
      } else {
	cerr << "Input MJD range keyword supplied only one of two required arguments\n";
	show_usage();
	return(1);
      } 
    } else if(string(argv[i]) == "-mpcout" || string(argv[i]) == "-mpcoutfile" || string(argv[i]) == "--mpcoutfile" || string(argv[i]) == "--mpcout" || string(argv[i]) == "-mpco") {
      if(i+1 < argc) {
	//There is still something to read;
	mpcout=argv[++i];
	i++;
      }
      else {
	cerr << "Output mpc file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  // Check for required arguments
  if(mpcfile.size()<=0) {
    cerr << "ERROR: you have not specified an input mpc obs file\n";
    show_usage();
    return(1);
  }
  if(mpcout.size()<=0) {
    cerr << "ERROR: you have not specified an output matched mpc file\n";
    show_usage();
    return(1);
  }
  
  // Read file of mpc observations
  instream1.open(mpcfile,ios_base::in);
  outstream1.open(mpcout);
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the mpc file
    getline(instream1,lnfromfile);
    if(!instream1.bad() && !instream1.fail() && !instream1.eof() && lnfromfile.size()>=80) {
      // Extract the MJD
      MJD = mpc80_mjd(lnfromfile);
      // Read the obscode
      obscode={};
      for(i=77;i<80;i++) {
	obscode.push_back(lnfromfile[i]);
      }
      if(MJD>=mjdstart && MJD<=mjdend) {
	// See if the code is OK
	isbadcode=0;
	for(i=0; i<badcodes.size(); i++) {
	  if(obscode.compare(badcodes[i])==0) isbadcode=1;
	}
	if(!isbadcode) {
	  // This line is good: write it to the output file
	  outstream1 << lnfromfile << "\n";
	}
      }
    }
  }
  instream1.close();
  outstream1.close();
  return(0);
}

