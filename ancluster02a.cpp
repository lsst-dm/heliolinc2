// January 12, 2022: ancluster02a.cpp
// Read cluster files produced by the good version of projectpairs04c.cpp,
// and analyze them. The good version of projectpairs04c.cpp fixes
// the omission of indices in the output cluster files. The old version
// of the current program, ancluster01a.cpp, had some
// time- and memory-intensive cleverness to cope with the absence
// of the indices. This new version requires the indices to be
// present, and uses them to proceed much more efficiently.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MAXBINSEARCH 50
#define MAXCLUSTRMS 1.0e5

class det_svec { // Detection with string ID, integer idex, and vector of
                 // indices. The vector is intended to contain indices for
                 // all the clusters of which this detection is a memeber.
public:
  double mjd;
  double RA;
  double Dec;
  string detid;
  int index;
  vector <int> indvec;
  det_svec(double mjd, double RA, double Dec, string detid, int index, vector <int> indvec) :mjd(mjd), RA(RA), Dec(Dec), detid(detid), index(index), indvec(indvec) { }
};
  
class clusteran01{ // Analysis of a heliolinc cluster (candidate discovery).
public:
  string recfile;
  int accelct;
  int clusternum;
  int numpoints;
  double timespan;
  int daysteps;
  vector <double> heliopar;
  vector <double> rmsvec; // pos(3), vel(3), postotal, veltotal, total 
  vector <int> clustind;
  clusteran01( string recfile, int accelct, int clusternum, int numpoints, double timespan, int daysteps, vector <double> rmsvec, vector <double> heliopar, vector <int> clustind) :recfile(recfile), accelct(accelct), clusternum(clusternum), numpoints(numpoints), timespan(timespan), daysteps(daysteps), rmsvec(rmsvec), heliopar(heliopar), clustind(clustind) {}
};


static void show_usage()
{
  cerr << "Usage: ancluster02a -pairdet pairdet_file -clist clusterlist -outfile outfile -outrms rmsfile\n";
}
    
int main(int argc, char *argv[])
{
  int i=1;
  int j=1;
  int i1=0;
  int i2=0;
  int detnum=0;
  int detct=0;
  string pairdetfile,clusterlist,clusterfile,outfile,rmsfile,sjunk;
  string lnfromfile;
  vector <string> clusternames;
  int clusternum = 0;
  int clusterct = 0;
  int clusternum2 = 0;
  int clusternum3 = 0;
  string w1,w2,w3;
  int ptnum=0;
  int ptct=0;
  double timespan=0l;
  int daysteps=0;
  int accelct=0;
  ifstream instream1;
  ofstream outstream1;
  ofstream outstream2;
  string detidstring;
  double tmjd;
  vector <det_svec> detvec;
  det_svec dsv = det_svec(0l,0l,0l,"",1,{});
  clusteran01 clustan = clusteran01("",0,0,0,0l,0,{},{},{});
  vector <clusteran01> clustanvec;
  double xrms,yrms,zrms,posrms;
  xrms = yrms = zrms = posrms = 0l;
  double vxrms,vyrms,vzrms,velrms,totrms;
  vxrms = vyrms = vzrms = velrms = totrms = 0l;
  double heliodist,heliovel,helioacc;
  heliodist = heliovel = helioacc = 0l;
  vector <double> heliopar;
  vector <double> rmsvec;
  int bestclust=0;
  double clustmetric=0;
  double bestclustmetric=0;
  int goodclusternum=0;
  string rating;

  if(argc!=9) {
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
    } else if(string(argv[i]) == "-clist" || string(argv[i]) == "-cl" || string(argv[i]) == "-clustlist" || string(argv[i]) == "--clusterlist" || string(argv[i]) == "--clist" || string(argv[i]) == "--clusterl" || string(argv[i]) == "--clustlist") {
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
    
  cout << "input paired detection file " << pairdetfile << "\n";
  cout << "input cluster list file " << clusterlist << "\n";
  cout << "output cluster file " << outfile << "\n";
  cout << "output rms file " << rmsfile << "\n";

  // Read paired detection file
  instream1.open(pairdetfile,ios_base::in);
  if(!instream1) {
    cerr << "can't open input file " << pairdetfile << "\n";
    return(1);
  }
  detvec={};
  while(instream1 >> dsv.mjd >> dsv.RA >> dsv.Dec >> w1 >> w2 >> w3 >> dsv.detid >> dsv.index) {
    dsv.indvec = {};
    detvec.push_back(dsv);
  }
  instream1.close();
  cout << "Successfully read " << detvec.size() << " detections from " << pairdetfile << "\n";
  
  instream1.open(clusterlist,ios_base::in);
  if(!instream1) {
    cerr << "can't open input file " << clusterlist << "\n";
    return(1);
  }
  while(instream1 >> sjunk)
    {
      if(sjunk.size()>0) clusternames.push_back(sjunk);
    }
  instream1.close();
  clusternum = clusternames.size();
  cout << "Read names of " << clusternum << " cluster files from " << clusterlist << "\n";
  
  outstream1.open(outfile,ios_base::out);
  outstream2.open(rmsfile,ios_base::out);

  // Read cluster files, loading clusters.
  clusternum3=0;
  for(clusterct=0; clusterct<clusternum; clusterct++) {
    instream1.open(clusternames[clusterct],ios_base::in);
    if(!instream1) {
      cerr << "can't open cluster file " << clusternames[clusterct] << "\n";
      return(2);
    } else cout << "Reading cluster file " << clusternames[clusterct] << "\n";
    while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      // Find out if we are reading a line introducing a grid point, like this:
      // Working on grid point 0: 216916911.500000 -2468364.855000 -18185.943176
      // or only one introducing a cluster, like this:
      // Accelct 0 Cluster 1 with 7257 = 7257 points.
      instream1 >> w1;
      if(w1 == "Working") {
	// introducing a grid point. Read heliocentric dist, vel, acc.
	instream1 >> w1 >> w2; // read, "on grid"
	instream1 >> w1 >> w2; // readl, "point nnnn:"
	instream1 >> heliodist >> heliovel >> helioacc;
	heliopar = {};
	heliopar.push_back(heliodist);
	heliopar.push_back(heliovel);
	heliopar.push_back(helioacc);
	//cout << "Heliopar: " << heliopar[0] << " " << heliopar[1] << " " << heliopar[2] << "\n";
	continue;
      } else {
	// introducing a cluster. Read accelct and clusternum2.
	instream1 >> accelct >> w1 >> clusternum2;
	//cout << "Accelct " << accelct << " Cluster " << clusternum2 << "\n";
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
      // Define clusteran01 object, except for final clustind vector.
      clustan = clusteran01(clusternames[clusterct],accelct,clusternum2,ptnum,timespan,daysteps,rmsvec,heliopar,{});
      // Read detection lines
      if(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
	for(ptct=0; ptct<ptnum; ptct++) {
	  // Format of line to be read
	  // 0 59929.300257 144.917102 1.510424 S100nkXka 3644298 4169268
	  instream1 >> w1 >> dsv.mjd;
	  instream1 >> dsv.RA >> dsv.Dec >> dsv.detid;
	  instream1 >> i1 >> i2;
	  // Add detection to cluster index vector.
	  clustan.clustind.push_back(i1);
	  // Add cluster index to detvec
	  detvec[i1].indvec.push_back(clustanvec.size());
	}
	// Push new cluster on to clustanvec.
	clustanvec.push_back(clustan);
      }
    }
    instream1.close();
  }

  cout << "Found " << clustanvec.size() << " clusters\n";

  // Loop on detections, finding the best cluster
  // for each detection.
  goodclusternum=0;
  for(detct=0; detct<detvec.size(); detct++)
    {
      if(detvec[detct].index>0) { // Detection not already marked as used.
	cout  << fixed << setprecision(6) << "Detection " << detct << ": " << detvec[detct].mjd << " " << detvec[detct].detid << " " << detvec[detct].indvec.size() << "\n";
	// Loop on clusters for this detection
	bestclust=-1;
	clustmetric=bestclustmetric=0.0;
	for(i=0;i<detvec[detct].indvec.size();i++) {
	  //cout  << fixed << setprecision(6) << "Considering cluster " << detvec[detct].indvec[i] << " pts, span, days, rms, metric: " << clustanvec[detvec[detct].indvec[i]].numpoints << " " << clustanvec[detvec[detct].indvec[i]].timespan << " " << clustanvec[detvec[detct].indvec[i]].daysteps;
	  //cout   << fixed << setprecision(3) << " " << clustanvec[detvec[detct].indvec[i]].rmsvec[8];
	  
	  // Best cluster is the one that maximizes clustmetric,
	  // provided it meets the basic conditions.
	  if(clustanvec[detvec[detct].indvec[i]].numpoints>=6 && clustanvec[detvec[detct].indvec[i]].rmsvec[8]<=MAXCLUSTRMS) {
	    clustmetric = double(clustanvec[detvec[detct].indvec[i]].numpoints)*double(clustanvec[detvec[detct].indvec[i]].daysteps)*clustanvec[detvec[detct].indvec[i]].timespan/clustanvec[detvec[detct].indvec[i]].rmsvec[8];
	    //cout   << fixed << setprecision(8) << " " << clustmetric << "\n";
	    if(clustmetric>bestclustmetric) {
	      bestclustmetric=clustmetric;
	      bestclust=detvec[detct].indvec[i];
	    }
	  }
	}
	if(bestclust>=0) {
	  goodclusternum++;
	  // A good cluster has been found.
	  cout  << fixed << setprecision(6) << "Best cluster is " << bestclust << " pts, span, days, rms, metric: " << clustanvec[bestclust].numpoints << " " << clustanvec[bestclust].timespan << " " << clustanvec[bestclust].daysteps;
	  cout  << fixed << setprecision(3) << " " << clustanvec[bestclust].rmsvec[8];
	  cout  << fixed << setprecision(8) << " " << bestclustmetric << "\n";
	  // See whether cluster is pure or mixed.
	  rating="PURE";
	  for(i=0; i<clustanvec[bestclust].clustind.size(); i++) {
	    if(i>0 && detvec[clustanvec[bestclust].clustind[i]].detid != detvec[clustanvec[bestclust].clustind[i-1]].detid) rating="MIXED";
	  }
	  // Write it to the output file.
	  outstream1 << "Working from heliolinc distvelacc point ";
	  outstream1 << fixed << setprecision(6) << clustanvec[bestclust].heliopar[0] << " " << clustanvec[bestclust].heliopar[1] << " " << clustanvec[bestclust].heliopar[2] << "\n";
	  outstream1  << fixed << setprecision(6) << clustanvec[bestclust].recfile << " " << clustanvec[bestclust].accelct << " Cluster " << clustanvec[bestclust].clusternum << " with " << clustanvec[bestclust].numpoints << " = " << clustanvec[bestclust].clustind.size() << " points, Globalct " << goodclusternum << " metric " << bestclustmetric << " " << rating << "\n";
	  outstream1 << fixed << setprecision(6) << "Unique pts: " << clustanvec[bestclust].numpoints << " span: " << clustanvec[bestclust].timespan << " daysteps: " << clustanvec[bestclust].daysteps << "\n";
	  outstream1 << fixed << setprecision(3) << "Cluster pos RMS: " << clustanvec[bestclust].rmsvec[0] << " "  << clustanvec[bestclust].rmsvec[1] << " "  << clustanvec[bestclust].rmsvec[2] << " total pos " << clustanvec[bestclust].rmsvec[6] << "\n";  
	  outstream1 << fixed << setprecision(3) << "Cluster vel RMS: " << clustanvec[bestclust].rmsvec[3] << " "  << clustanvec[bestclust].rmsvec[4] << " "  << clustanvec[bestclust].rmsvec[5] << " total vel " << clustanvec[bestclust].rmsvec[7] << "\n";  
	  outstream1 << fixed << setprecision(3) << "Cluster total RMS: " << clustanvec[bestclust].rmsvec[8] << "\n";
	  for(i=0; i<clustanvec[bestclust].clustind.size(); i++)
	    {
	      j = clustanvec[bestclust].clustind[i];
	      outstream1 << fixed << setprecision(6) << i << " " << detvec[j].mjd << " " << detvec[j].RA << " " << detvec[j].Dec << " " << detvec[j].detid << " " << j << " " << detvec[j].index << "\n";
	    }
	  outstream1 << "\n";
	  // Write summary line to the RMS file
	  outstream2 << "Globalct " << goodclusternum << " : ";
	  outstream2 << fixed << setprecision(6) << clustanvec[bestclust].heliopar[0] << " " << clustanvec[bestclust].heliopar[1] << " " << clustanvec[bestclust].heliopar[2] << " : ";
	  outstream2  << fixed << setprecision(6) << clustanvec[bestclust].numpoints << " " << clustanvec[bestclust].timespan << " " << clustanvec[bestclust].daysteps << " ";
	  outstream2 << fixed << setprecision(3) << clustanvec[bestclust].rmsvec[6] << " "  << clustanvec[bestclust].rmsvec[7] << " "  << clustanvec[bestclust].rmsvec[8] << " " << rating << " " << bestclustmetric << "\n";  

	  // Mark all the clusters containing this point as already considered
	  for(i=0;i<detvec[detct].indvec.size();i++) {
	    clustanvec[detvec[detct].indvec[i]].numpoints = 0;
	  }
	  // Mark all of the points in the chosen cluster as already considered
	  for(i=0;i<clustanvec[bestclust].clustind.size();i++) {
	    detvec[clustanvec[bestclust].clustind[i]].index = 0;
	  }
	} else {
	  cout << "No good cluster was found containing this detection\n";
	}
      }
    }

  return(0);
}
