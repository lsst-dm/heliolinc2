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
  cerr << "Usage: cluster_readorig01a.cpp -orig original_input -clust clusterfile -outfile outfile\n";
}
    
int main(int argc, char *argv[])
{
  int i=1;
  int linect=0;
  int linenum=0;
  string origfile,clusterfile,outfile;
  string lnfromfile;
  string w1,w2;
  double xrms,yrms,zrms,posrms;
  xrms = yrms = zrms = posrms = 0l;
  double vxrms,vyrms,vzrms,velrms,totrms;
  vxrms = vyrms = vzrms = velrms = totrms = 0l;
  double heliodist, heliovel, helioacc;
  vector <string> inputlines;
  ifstream instream1;
  ofstream outstream1;
  vector <float> heliopar;
  vector <float> rmsvec;
  vector <double> mjdvec;
  float clustmetric=0;
  int ptnum=0;
  int ptct=0;
  double timespan=0l;
  int daysteps=0;
  int globalct = 0;
  int ijunk1 = 0;
  int ijunk2 = 0;
  long i1 = 0;
  long i2 = 0;
  double mjd,RA,Dec;
  mjd = RA = Dec = 0l;
  string detid;

  if(argc<7) {
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
  if(!instream1) {
    cerr << "can't open cluster file " << clusterfile << "\n";
    return(2);
  } else cout << "Reading cluster file " << clusterfile << "\n";

  // Open output file
  outstream1.open(outfile,ios_base::out);

  // Begin reading cluster file
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Find out if we are reading a line introducing a grid point, like this:
    // Working on grid point 0: 216916911.500000 -2468364.855000 -18185.943176
    // or only one introducing a cluster, like this:
    // outclust_01bmb22.txt 965  Cluster 6794 	with 12 = 12 points, Globalct 6077 metric 0.329127 PURE
    instream1 >> w1;
    if(w1 == "Working") {
      // introducing a grid point. Read heliocentric dist, vel, acc.
      instream1 >> w1 >> w2; // read, "on grid"
      instream1 >> w1 >> w2; // readl, "point nnnn:"
      instream1 >> heliodist >> heliovel >> helioacc;
      continue;
    } else {
      // Introducing a cluster. Read and discard a bunch of stuff,
      // because all we care about right now is the global cluster count
      // and the cluster metric.
      instream1 >> ijunk1 >> w1 >> ijunk2 >> w2; // expected format: "965  Cluster 6794 with" 
      instream1 >> ijunk1 >> w1 >> ijunk2 >> w2; // expected format: "12 = 12 points,"
      // Now for the stuff we actually care about
      instream1 >> w1 >> globalct >> w2 >> clustmetric; // expected format: "Globalct 6077 metric 0.329127"
      getline(instream1,lnfromfile);
    }
    // Read third line to get number of unique points, time span, and daysteps
    // Format of line to be read:
    // Unique pts: 136 span: 3.018265 daysteps: 2
    instream1 >> w1 >> w2 >> ptnum;
    instream1 >> w1 >> timespan >> w2 >> daysteps;
    //cout << "Unique: " << ptnum << " " << timespan << " " << daysteps << "\n";
    // Read position RMS line.
    // Format of line to be read:
    // Cluster pos RMS: 69955.132 20488.722 22187.227 total pos 76195.678
    instream1 >> w1 >> w2;
    instream1 >> w1 >> xrms >> yrms >> zrms;
    instream1 >> w1 >>  w2 >> posrms;
    // Read velocity RMS line.
    // Format of line to be read:
    // Cluster vel RMS: 30801.886 19634.520 17591.633 total vel 40543.015
    instream1 >> w1 >> w2;
    instream1 >> w1 >> vxrms >> vyrms >> vzrms;
    instream1 >> w1 >>  w2 >> velrms;
    // Read total RMS line.
    // Format of line to be read:
    // Cluster total RMS: 86310.587
    instream1 >> w1 >> w2;
    instream1 >> w1 >> totrms;
    rmsvec={};
    rmsvec.push_back(xrms);
    rmsvec.push_back(yrms);
    rmsvec.push_back(zrms);
    rmsvec.push_back(posrms);
    rmsvec.push_back(vxrms);
    rmsvec.push_back(vyrms);
    rmsvec.push_back(vzrms);
    rmsvec.push_back(velrms);
    rmsvec.push_back(totrms);
    //cout << fixed << setprecision(6) << "Cluster total RMS: " << rmsvec[8] << "\n";
    // Read detection lines
    if(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      for(ptct=0; ptct<ptnum; ptct++) {
	// Format of line to be read
	// 0 59929.300257 144.917102 1.510424 S100nkXka 3644298 4169268
	instream1 >> w1 >> mjd;
	instream1 >> RA >> Dec >> detid;
	instream1 >> i1 >> i2;

	// Write it to the output file.
	outstream1 << inputlines[i2-1] << ",";
	// Additional output: globalct, heliodist, heliovel, helioacc, ptnum, timespan, daysteps, totrms, clustermetric
	outstream1 << fixed << setprecision(6) << globalct << "," << heliodist << "," << heliovel << "," << helioacc << "," << ptnum << ",";
	outstream1 << fixed << setprecision(6) << timespan << "," << daysteps << "," << totrms << "," << clustmetric << "\n";
      }
    }
  }
  instream1.close();
  outstream1.close();

  return(0);
}
