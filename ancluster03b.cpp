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


static void show_usage()
{
  cerr << "Usage: ancluster03a -pairdet pairdet_file -clist clusterlist -maxrms maxrms -outfile outfile -outrms rmsfile\n";
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
  int clusterct2 = 0;
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
  double maxrms = MAXCLUSTRMS;
  vector <det_svec> detvec;
  vector <det_svec> clusterdets;
  det_svec dsv = det_svec(0l,0l,0l,"",1,{});
  clusteran03 clustan = clusteran03("",0,0,0,0l,0,{},{},{},0l);
  vector <clusteran03> clustanvec;
  float_index findex = float_index(0L,0);
  vector <float_index> clustanvec2;
  double xrms,yrms,zrms,posrms;
  xrms = yrms = zrms = posrms = 0l;
  double vxrms,vyrms,vzrms,velrms,totrms;
  vxrms = vyrms = vzrms = velrms = totrms = 0l;
  double heliodist,heliovel,helioacc;
  heliodist = heliovel = helioacc = 0l;
  vector <float> heliopar;
  vector <float> rmsvec;
  vector <double> mjdvec;
  float clustmetric=0;
  int goodclusternum=0;
  int istimedup=0;
  string rating;

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
  cout << "Maximum RMS in km: " << maxrms << "\n";
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
      // Calculate cluster metric.
      clustmetric = double(ptnum)*double(daysteps)*timespan/totrms;
      // Define clusteran03 object, except for final clustind vector.
      clustan = clusteran03(clusternames[clusterct],accelct,clusternum2,ptnum,timespan,daysteps,rmsvec,heliopar,{},clustmetric);
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
	}
	// See if it's a good cluster
	if(clustan.rmsvec[8]<=maxrms) {
	  // Check for duplicate MJDs
	  mjdvec={};
	  for(ptct=0; ptct<ptnum; ptct++) {
	    i1 = clustan.clustind[ptct];
	    mjdvec.push_back(detvec[i1].mjd);
	  }
	  sort(mjdvec.begin(),mjdvec.end());
	  istimedup=0;
	  for(ptct=1; ptct<ptnum; ptct++) {
	    if(mjdvec[ptct-1] == mjdvec[ptct]) istimedup=1;
	  }
	  if(istimedup==0) {
	    // The cluster is good.
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
	// Close if-statement checking that we're still reading the file.
      }
      // Close while-loop reading the file.
    }
    instream1.close();
    // Close loop on input cluster files.
  }

  clusternum3 = clustanvec.size();
  cout << "Found " << clusternum3 << " clusters\n";
  // Load just clustermetric values and indices from
  // clustanvec into clustanvec2
  clustanvec2 = {};
  // Record indices so information won't be lost on sort
  for(clusterct=0; clusterct<clusternum3; clusterct++) {
    findex = float_index(clustanvec[clusterct].clustermetric,clusterct);
    clustanvec2.push_back(findex);
  }
  // Sort clustanvec2
  sort(clustanvec2.begin(), clustanvec2.end(), lower_float_index());
  
  // Loop on clusters, starting with the best (highest metric),
  // and eliminating duplicates
  goodclusternum=0;
  for(clusterct2=clusternum3-1; clusterct2>=0; clusterct2--) {
    clusterct = clustanvec2[clusterct2].index;
    cout  << fixed << setprecision(6) << "Cluster " << clusterct << ": of " << clusternum3 << ": " << clustanvec[clusterct].numpoints << " " << clustanvec[clusterct].timespan << " " << clustanvec[clusterct].daysteps << " " << clustanvec[clusterct].rmsvec[8] << " " << clustanvec[clusterct].clustermetric << "\n";
    if(clustanvec[clusterct].numpoints>=6 && clustanvec[clusterct].rmsvec[8]<=maxrms) {
      // This is a good cluster not already marked as used.
      goodclusternum++;
      cout << "Accepted as good cluster " << goodclusternum << "\n";
      // See whether cluster is pure or mixed.
      rating="PURE";
      for(i=0; i<clustanvec[clusterct].clustind.size(); i++) {
	if(i>0 && detvec[clustanvec[clusterct].clustind[i]].detid != detvec[clustanvec[clusterct].clustind[i-1]].detid) rating="MIXED";
      }
      // Write it to the output file.
      outstream1 << "Working from heliolinc distvelacc point ";
      outstream1 << fixed << setprecision(6) << clustanvec[clusterct].heliopar[0] << " " << clustanvec[clusterct].heliopar[1] << " " << clustanvec[clusterct].heliopar[2] << "\n";
      outstream1  << fixed << setprecision(6) << clustanvec[clusterct].recfile << " " << clustanvec[clusterct].accelct << " Cluster " << clustanvec[clusterct].clusternum << " with " << clustanvec[clusterct].numpoints << " = " << clustanvec[clusterct].clustind.size() << " points, Globalct " << goodclusternum << " metric " << clustanvec[clusterct].clustermetric << " " << rating << "\n";
      outstream1 << fixed << setprecision(6) << "Unique pts: " << clustanvec[clusterct].numpoints << " span: " << clustanvec[clusterct].timespan << " daysteps: " << clustanvec[clusterct].daysteps << "\n";
      outstream1 << fixed << setprecision(3) << "Cluster pos RMS: " << clustanvec[clusterct].rmsvec[0] << " "  << clustanvec[clusterct].rmsvec[1] << " "  << clustanvec[clusterct].rmsvec[2] << " total pos " << clustanvec[clusterct].rmsvec[6] << "\n";
      outstream1 << fixed << setprecision(3) << "Cluster vel RMS: " << clustanvec[clusterct].rmsvec[3] << " "  << clustanvec[clusterct].rmsvec[4] << " "  << clustanvec[clusterct].rmsvec[5] << " total vel " << clustanvec[clusterct].rmsvec[7] << "\n";  
      outstream1 << fixed << setprecision(3) << "Cluster total RMS: " << clustanvec[clusterct].rmsvec[8] << "\n";
      // Load vector of detections in this cluster (needed only for time-sorting).
      clusterdets={};
      for(i=0; i<clustanvec[clusterct].clustind.size(); i++) {
	j = clustanvec[clusterct].clustind[i];
	clusterdets.push_back(detvec[j]);
	clusterdets[i].index=j;
      }
      sort(clusterdets.begin(), clusterdets.end(), early_det_svec());
      // Write detections to output file in time-sorted order.
      for(i=0; i<clustanvec[clusterct].clustind.size(); i++) {
	j = clusterdets[i].index;
	outstream1 << fixed << setprecision(6) << i << " " << detvec[j].mjd << " " << detvec[j].RA << " " << detvec[j].Dec << " " << detvec[j].detid << " " << j << " " << detvec[j].index << "\n";
      }
      outstream1 << "\n";
      // Write summary line to the RMS file
      outstream2 << "Globalct " << goodclusternum << " : ";
      outstream2 << fixed << setprecision(6) << clustanvec[clusterct].heliopar[0] << " " << clustanvec[clusterct].heliopar[1] << " " << clustanvec[clusterct].heliopar[2] << " : ";
      outstream2  << fixed << setprecision(6) << clustanvec[clusterct].numpoints << " " << clustanvec[clusterct].timespan << " " << clustanvec[clusterct].daysteps << " ";
      outstream2 << fixed << setprecision(3) << clustanvec[clusterct].rmsvec[6] << " "  << clustanvec[clusterct].rmsvec[7] << " "  << clustanvec[clusterct].rmsvec[8] << " " << rating << " ";
      outstream2 << fixed << setprecision(9) << clustanvec[clusterct].clustermetric << "\n";
      for(i=0;i<clustanvec[clusterct].clustind.size();i++) {
	// This point is in the chosen cluster, and cannot be in any other
	detct = clustanvec[clusterct].clustind[i];
	//cout << "Deduplicating detection " << detct << "\n";
	// Mark all the clusters containing this point as already considered
	for(j=0;j<detvec[detct].indvec.size();j++) {
	  //cout << "         wiping cluster " << detvec[detct].indvec[j] << "\n";
	  clustanvec[detvec[detct].indvec[j]].numpoints = 0;
	}
      }
    } else cout << "Rejected as bad or else redundant with a better cluster\n";
  }

  return(0);
}
