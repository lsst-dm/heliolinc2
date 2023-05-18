// May 15, 2023: crosslink_clust2det
// Given a set of pairdets, cluster summary, and clust2det files produced
// by running link_refine_Herget_new on ITF data, link objects across different
// files by detecting cluster overlap through the tracklet identifiers

#include "solarsyst_dyn_geo01.h"
#include "cmath"

string char2string(char *inarray, int size)
{
  int i=0;
  string stest="";
  while(i<size && inarray[i]!=0) {
    stest = stest + inarray[i];
    i++;
  }
  //for(i=0; i<size; i++) {
  //  stest = stest + inarray[i];
  //}
  return(stest);
}

static void show_usage()
{
  cerr << "Usage: crosslink_clustdet -setlist list of file sets -trackdiv tracklet_division_time -maxrms max astromeric RMS (arcsec) -out output file -outmpc output MPC80 file -outid output ID files\n";
}

int main(int argc, char *argv[])
{
  string sumfile,clust2detfile,setlist;
  string outfile = "crossout_test01";
  string outMPC = "crossMPC_test01";
  string outid = "crossID_test01";
  int stringsize=0;
  vector <hlclust> inclustvec;
  vector  <longpair> inclust2det;
  vector <hldet> detvec;
  vector <hldet> clustvec;
  vector <hldet> clustvec2;
  vector <hldet> trackvec;
  vector <hldet> cluster_detvec_main = {};
  vector <hlclust> inclust_main = {};
  long detnum = 0;
  long c2dnum = 0;
  long clustnum = 0;
  long detct=0;
  long c2dct=0;
  long clustct=0;
  long clustct2=0;
  vector <hldet> clustdet = {};
  hldet o1 = hldet(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, "null", "V", "500", 0, 0, 0);
 
  double angvel,crosstrack,alongtrack,PA,poleRA,poleDec;
  angvel=crosstrack=alongtrack=PA=poleRA=poleDec=0.0;
  double arc,timespan;
  arc = timespan = 0;
  vector <double> angvelvec;
  vector <double> GCRvec;
  vector <double> PAvec;
  vector <double> arcvec;
  vector <double> timespanvec;
  long tracknum=0;
  double nightstep = 3.0l/24.0l;
  string pairdetfile,stest;
  ifstream instream1;
  ofstream outstream1;
  ofstream outstream2;
  ofstream outstream3;
  int status=0;
  long i=0;
  long j=0;
  int bandlen=0;
  int verbose=0;
  long max_known_obj=0;
  double avg_det_qual=0.0l;
  string stest1,stest2;
  int year,month;
  double day;
  int rahr,ramin;
  int decdeg,decmin;
  double Dec,rasec,decsec;
  string signstring;
  double maxrms = 1.0l;
  int filect=0;
  vector <long> loadindex;
  vector <string> loadstring;
  vector <vector <long>> cluster_index;
  vector <vector <string>> cluster_IDvec;
  int nomatch=1;
  long BAmatch=0;
  long mergecount=0;
  long totalct=0;

  if(argc<3) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-setlist" || string(argv[i]) == "-list" || string(argv[i]) == "-set" || string(argv[i]) == "-sl" || string(argv[i]) == "--setlist" || string(argv[i]) == "--set" || string(argv[i]) == "--list" || string(argv[i]) == "--filesetlist") {
      if(i+1 < argc) {
	//There is still something to read;
	setlist=argv[++i];
	i++;
      }
      else {
	cerr << "Input setlist file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-trackdiv" || string(argv[i]) == "-nightstep") {
      if(i+1 < argc) {
	//There is still something to read;
	nightstep = stod(argv[++i]);
	nightstep /= 24.0l;
	i++;
      }
      else {
	cerr << "Input cluster list keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-maxrms" || string(argv[i]) == "-rms") {
      if(i+1 < argc) {
	//There is still something to read;
	maxrms = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input cluster list keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outclust" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outclust" || string(argv[i]) == "--output_cluster_file") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outMPC" || string(argv[i]) == "-MPC" || string(argv[i]) == "-outmpc" || string(argv[i]) == "-mpc" || string(argv[i]) == "--outMPC" ) {
      if(i+1 < argc) {
	//There is still something to read;
	outMPC=argv[++i];
	i++;
      }
      else {
	cerr << "Output MPC file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outid" || string(argv[i]) == "-outID" || string(argv[i]) == "-ID" || string(argv[i]) == "-id") {
      if(i+1 < argc) {
	//There is still something to read;
	outid=argv[++i];
	i++;
      }
      else {
	cerr << "Output ID file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose=stoi(argv[++i]);
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

  // Catch required parameters if missing
  if(setlist.size()<=0) {
    cout << "\nERROR: input list of file sets is required\n";
    show_usage();
    return(1);
  }

  // Open list file
  instream1.open(setlist);
  while(instream1 >> pairdetfile >> sumfile >> clust2detfile) {
    // Echo input files
    cout << "Input file set: " << pairdetfile << " " << sumfile << " " << clust2detfile << "\n";
    // Read input files
    // Read paired detection file
    detvec={};
    status=read_pairdet_file(pairdetfile, detvec, verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read paired detection file " << pairdetfile << "\n";
      cerr << "read_pairdet_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << detvec.size() << " data lines from paired detection file " << pairdetfile << "\n";
  
    // Read cluster summary file
    inclustvec={};
    status=read_clustersum_file(sumfile, inclustvec, verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read input cluster summary file " << sumfile << "\n";
      cerr << "read_clustersum_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << inclustvec.size() << " data lines from cluster summary file " << sumfile << "\n";
    inclust2det={};
    // Read cluster-to-detection file
    status=read_longpair_file(clust2detfile, inclust2det, verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read cluster-to-detection file " << clust2detfile << "\n";
      cerr << "read_longpair_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << inclust2det.size() << " data lines from cluster-to-detection file " << clust2detfile << "\n";

    detnum = detvec.size();
    c2dnum = inclust2det.size();
    clustnum = inclustvec.size();
    detct=0;
    c2dct=0;
    clustct=0;
    for(clustct=0 ; clustct<clustnum; clustct++) {
      if(inclustvec[clustct].astromRMS < maxrms) {
	// This cluster is good, add it to the master vector
	inclust_main.push_back(inclustvec[clustct]);
	// Add all of its points to the clust2det vector
	while(inclust2det[c2dct].i1==clustct && c2dct<c2dnum) {
	  // This point belongs to the cluster under consideration.
	  // Add it to the master vector
	  detct = inclust2det[c2dct].i2;
	  o1 = detvec[detct];
	  o1.index = inclust_main.size()-1;
	  cluster_detvec_main.push_back(o1);
	  c2dct++;
	}
      } else {
	// This cluster is no good. Skip over all the corresponding points
	while(inclust2det[c2dct].i1==clustct && c2dct<c2dnum) c2dct++;
      }
    }
    filect++;
    cout << "Finished reading file set number " << filect << ": vector sizes now " << inclust_main.size() << " and " << cluster_detvec_main.size() << "\n";
  }
  instream1.close();
  cout << "Finished reading all of the file sets\n";
  // Load vectors
  cluster_index={};
  cluster_IDvec={};
  detnum = cluster_detvec_main.size();
  clustnum = inclust_main.size();
  detct=0;
  clustct=0;
  for(clustct=0 ; clustct<clustnum; clustct++) {
    loadindex={};
    loadstring={};
    cout << "Working on cluster " << clustct << ", detct = " << detct << ", index = " << cluster_detvec_main[detct].index << "\n";
    while(cluster_detvec_main[detct].index==clustct && detct<detnum) {
      // This point belongs to the cluster under consideration.
      stringsize = sizeof(cluster_detvec_main[detct].idstring)/sizeof(char);      
      stest = char2string(cluster_detvec_main[detct].idstring, stringsize);
      cout << "copying " << cluster_detvec_main[detct].idstring << " to " << stest << "\n";
      if(loadstring.size()<=0) loadstring.push_back(stest);
      else {
	nomatch=1;
	for(i=0;i<long(loadstring.size());i++) {
	  if(loadstring[i] == stest) nomatch=0; // It's a match
	}
	if(nomatch==1) {
	  // This string tracklet identifier hasn't previously been seen
	  loadstring.push_back(stest);
	}
      }
      detct++;
    }
    loadindex.push_back(clustct);
    cluster_index.push_back(loadindex);
    cluster_IDvec.push_back(loadstring);
  }
  cout << "Finished loading cluster_index and cluster_IDvec\n";
  // Search the master vectors for duplicates
  for(clustct=0 ; clustct<clustnum; clustct++) {
    if(cluster_index[clustct].size()>0) {
      // Primary cluster is worth considering, hasn't been subsumed
      // to another cluster
      for(clustct2=clustct+1; clustct2<clustnum; clustct2++) {
	if(cluster_index[clustct2].size()>0) {
	  // clustct2 hasn't yet been subsumed either
	  BAmatch=0;
	  // How many of the IDs from clustct2 match one from clustct?
	  for(i=0; i<long(cluster_IDvec[clustct2].size()); i++) {
	    nomatch=1;
	    for(j=0; j<long(cluster_IDvec[clustct].size()); j++) {
	      if(cluster_IDvec[clustct2][i]==cluster_IDvec[clustct][j]) {
		nomatch = 0; // It's a match
	      }
	    }
	    if(nomatch==0) {
	      BAmatch+=1;
	    }
	  }
	  if(BAmatch == long(cluster_IDvec[clustct2].size())) {
	    // All of the points in clustct2 were matched:
	    // the clusters are entirely redundant
	    cluster_index[clustct2] = {}; // Wipe the cluster index to mark redundancy
	  } else if(BAmatch>0 && BAmatch<long(cluster_IDvec[clustct2].size())) {
	    // Some of the IDstrings were matched, but not all of them
	    // Add the non-redundant IDstrings to the primary cluster
	    for(i=0; i<long(cluster_IDvec[clustct2].size()); i++) {
	      nomatch=1;
	      for(j=0; j<long(cluster_IDvec[clustct].size()); j++) {
		if(cluster_IDvec[clustct2][i]==cluster_IDvec[clustct][j]) {
		  nomatch = 0; // It's a match
		}
	      }
	      if(nomatch==1) {
		// We have a non-redundant one here. Add it to the primary cluster
		cluster_IDvec[clustct].push_back(cluster_IDvec[clustct2][i]);
		if(cluster_index[clustct][cluster_index[clustct].size()-1] != clustct2) {
		  cluster_index[clustct].push_back(clustct2);
		}
	      }
	    }
	    // Now clustct2 has been subsumed into clustct
	    cluster_index[clustct2] = {}; // Wipe the cluster index to mark redundancy
	  } else if(BAmatch > long(cluster_IDvec[clustct2].size())) {
	    cerr << "ERROR: logically impossible case between clusters " << clustct << " and " << clustct2 << "\n";
	    return(1);
	  }
	}
      }
    }
  }
  
  outstream1.open(outfile);
  outstream2.open(outMPC);
  outstream3.open(outid);
  // Write data for cross-linked clusters first
  for(clustct=0 ; clustct<clustnum; clustct++) {
    if(cluster_index[clustct].size()>1) {
      // Write out all the tracklet strings
      for(i=0; i<long(cluster_IDvec[clustct].size()); i++) {
	outstream3 << cluster_IDvec[clustct][i] << " ";
      }
      outstream3 << "\n";
      mergecount+=1;
      totalct+=1;

      // Load cluster-specific detection vector for this cluster
      clustvec={};
      detct=0;
      for(i=0; i<long(cluster_index[clustct].size()); i++) {
	clustct2 = cluster_index[clustct][i];
	while(detct<detnum && cluster_detvec_main[detct].index < clustct2) detct++;
	while(detct<detnum && cluster_detvec_main[detct].index == clustct2) {
	  clustvec.push_back(cluster_detvec_main[detct]);
	  detct++;
	}
      }
      // Time-sort and de-duplicate clustvec
      sort(clustvec.begin(), clustvec.end(), early_hldet());
      clustvec2={};
      clustvec2.push_back(clustvec[0]);
      for(i=1; i<long(clustvec.size()); i++) {
	if(clustvec[i].MJD != clustvec2[clustvec2.size()-1].MJD && clustvec[i].RA != clustvec2[clustvec2.size()-1].RA) {
	  // The new point is not redundant: add it.
	  clustvec2.push_back(clustvec[i]);
	}
      }
      clustvec = clustvec2;
      // Run some analytics on this detection cluster
      avg_det_qual = 0.0l;
      max_known_obj=0;
      for(i=0; i<long(clustvec.size()); i++) {
	avg_det_qual += double(clustvec[i].det_qual);
	if(clustvec[i].known_obj > max_known_obj) max_known_obj = clustvec[i].known_obj;
      }
      avg_det_qual/=double(clustvec.size());
      // Loop over clustvec to extract individual tracklets
      trackvec={};
      trackvec.push_back(clustvec[0]);
      angvelvec = GCRvec = PAvec = timespanvec = arcvec = {};
      for(i=1; i<long(clustvec.size()); i++) {
	if((clustvec[i].MJD - clustvec[i-1].MJD) < nightstep) {
	  // Add a new point to this tracklet
	  trackvec.push_back(clustvec[i]);
	} else {
	  // A tracklet is finished. Analyze it.
	  tracknum = trackvec.size();
	  if(tracknum>1) {
	    greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
	    angvelvec.push_back(angvel);
	    GCRvec.push_back(sqrt(DSQUARE(crosstrack)+DSQUARE(alongtrack)));
	    PAvec.push_back(PA);
	    timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
	    arc = timespan*angvel;
	    arcvec.push_back(arc*3600.0l);
	    timespanvec.push_back(timespan*24.0l);
	    // Wipe trackvec, and load the next point of the next tracklet.
	    trackvec = {};
	    trackvec.push_back(clustvec[i]);
	  }
	}
      }
      if(tracknum>1) {
	// Handle a final tracklet.
	greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
	angvelvec.push_back(angvel);
	GCRvec.push_back(sqrt(DSQUARE(crosstrack)+DSQUARE(alongtrack)));
	PAvec.push_back(PA);
	timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
	arc = timespan*angvel;
	arcvec.push_back(arc*3600.0l);
	timespanvec.push_back(timespan*24.0l);
	// Wipe trackvec
	trackvec = {};
      }
      // Sort all of the tracklet statistics vectors
      tracknum = angvelvec.size();
      sort(angvelvec.begin(), angvelvec.end());
      sort(GCRvec.begin(), GCRvec.end());
      sort(PAvec.begin(), PAvec.end());
      sort(timespanvec.begin(), timespanvec.end());
      sort(arcvec.begin(), arcvec.end());
    
      outstream1 << "\n#clusternum,posRMS,velRMS,totRMS,astromRMS,timespan,uniquepoints,obsnights,metric,orbit_a,orbit_e,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count,avg_det_qual,max_known_obj,minvel,maxvel,minGCR,maxGCR,minpa,maxpa,mintimespan,maxtimespan,minarc,maxarc,contributors\n";
      outstream1 << fixed << setprecision(3) << inclustvec[clustct].clusternum << "," << inclustvec[clustct].posRMS << "," << inclustvec[clustct].velRMS << "," << inclustvec[clustct].totRMS << ",";
      outstream1 << fixed << setprecision(4) << inclustvec[clustct].astromRMS << ",";
      outstream1 << fixed << setprecision(6) << inclustvec[clustct].timespan << "," << inclustvec[clustct].uniquepoints << "," << inclustvec[clustct].obsnights << "," << inclustvec[clustct].metric << ",";
      outstream1 << fixed << setprecision(6) << inclustvec[clustct].orbit_a << "," << inclustvec[clustct].orbit_e << "," << inclustvec[clustct].orbit_MJD << ",";
      outstream1 << fixed << setprecision(1) << inclustvec[clustct].orbitX << "," << inclustvec[clustct].orbitY << "," << inclustvec[clustct].orbitZ << ",";
      outstream1 << fixed << setprecision(4) << inclustvec[clustct].orbitVX << "," << inclustvec[clustct].orbitVY << "," << inclustvec[clustct].orbitVZ << "," << inclustvec[clustct].orbit_eval_count << ",";
      outstream1 << fixed << setprecision(1) << avg_det_qual << "," << max_known_obj << ",";
      outstream1 << fixed << setprecision(6) << angvelvec[0] << "," << angvelvec[tracknum-1] << "," << GCRvec[0] << "," << GCRvec[tracknum-1] << "," << PAvec[0] << "," << PAvec[tracknum-1] << "," << timespanvec[0] << "," << timespanvec[tracknum-1] << "," << arcvec[0] << "," << arcvec[tracknum-1] << "," << cluster_index[clustct].size() << "\n";
    
      // Write this data to the output file
      outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
      for(i=0; i<long(clustvec.size()); i++) {
	outstream1 << fixed << setprecision(7) << clustvec[i].MJD << "," << clustvec[i].RA << "," << clustvec[i].Dec << ",";
	outstream1 << fixed << setprecision(4) << clustvec[i].mag << ",";
	outstream1 << fixed << setprecision(2) << clustvec[i].trail_len << "," << clustvec[i].trail_PA << ",";
	outstream1 << fixed << setprecision(4) << clustvec[i].sigmag << ",";
	outstream1 << fixed << setprecision(3) << clustvec[i].sig_across << "," << clustvec[i].sig_along << ",";
	outstream1 << clustvec[i].image << "," << clustvec[i].idstring << "," << clustvec[i].band << ",";
	outstream1 << clustvec[i].obscode << "," << clustvec[i].known_obj << ",";
	outstream1 << clustvec[i].det_qual << "," << clustvec[i].index << "\n";
      }
    	
      outstream2 << "\n#clusternum,posRMS,velRMS,totRMS,astromRMS,timespan,uniquepoints,obsnights,metric,orbit_a,orbit_e,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count,avg_det_qual,max_known_obj,minvel,maxvel,minGCR,maxGCR,minpa,maxpa,mintimespan,maxtimespan,minarc,maxarc,contributors\n";
      outstream2 << fixed << setprecision(3) << inclustvec[clustct].clusternum << "," << inclustvec[clustct].posRMS << "," << inclustvec[clustct].velRMS << "," << inclustvec[clustct].totRMS << ",";
      outstream2 << fixed << setprecision(4) << inclustvec[clustct].astromRMS << ",";
      outstream2 << fixed << setprecision(6) << inclustvec[clustct].timespan << "," << inclustvec[clustct].uniquepoints << "," << inclustvec[clustct].obsnights << "," << inclustvec[clustct].metric << ",";
      outstream2 << fixed << setprecision(6) << inclustvec[clustct].orbit_a << "," << inclustvec[clustct].orbit_e << "," << inclustvec[clustct].orbit_MJD << ",";
      outstream2 << fixed << setprecision(1) << inclustvec[clustct].orbitX << "," << inclustvec[clustct].orbitY << "," << inclustvec[clustct].orbitZ << ",";
      outstream2 << fixed << setprecision(4) << inclustvec[clustct].orbitVX << "," << inclustvec[clustct].orbitVY << "," << inclustvec[clustct].orbitVZ << "," << inclustvec[clustct].orbit_eval_count << ",";
      outstream2 << fixed << setprecision(1) << avg_det_qual << "," << max_known_obj << ",";
      outstream2 << fixed << setprecision(6) << angvelvec[0] << "," << angvelvec[tracknum-1] << "," << GCRvec[0] << "," << GCRvec[tracknum-1] << "," << PAvec[0] << "," << PAvec[tracknum-1] << "," << timespanvec[0] << "," << timespanvec[tracknum-1] << "," << arcvec[0] << "," << arcvec[tracknum-1] << "," << cluster_index[clustct].size() << "\n";
    
      // Write this data to the output file
      outstream2 << "#MPC 80-column formatted observations:\n";
      for(i=0; i<long(clustvec.size()); i++) {
	// The temporary ID will just be 'a' followed by a 6-digit number
	if(totalct<=999999) stest1 = to_string(totalct);
	else {
	  stest1 = to_string(totalct%1000000);
	}
	while(stest1.size()<6) {
	  stest1 = "0" + stest1;
	}
	outstream2 << "     a" << stest1 << "  C";
	mjd2mpcdate(clustvec[i].MJD,year,month,day);
	outstream2 << year << " ";
	if(month<10) outstream2 << "0";
	outstream2 << month << " ";
	day = round(day*1000000.0l);
	day /= 1000000.0l;
	if(day<10.0l) outstream2  << fixed << setprecision(6) << "0";
	outstream2  << fixed << setprecision(6) << day;
	// Convert RA, Dec from decimal degrees to sexagesimal format.
	rahr = int(clustvec[i].RA/15.0l);
	ramin = int(clustvec[i].RA*4.0l - double(rahr)*60.0l);
	rasec = clustvec[i].RA*240.0l - double(rahr)*3600.0l - double(ramin)*60.0l;
	rasec = round(rasec*1000.0l);
	rasec /= 1000.0l;
	if(clustvec[i].Dec>=0) {
	  signstring="+";
	  Dec = clustvec[i].Dec;
	} else {
	  signstring="-";
	  Dec = -clustvec[i].Dec;
	}
	decdeg = int(Dec);
	decmin = int(Dec*60.0l - double(decdeg)*60.0l);
	decsec = Dec*3600.0l - double(decdeg)*3600.0l - double(decmin)*60.0l;
	decsec = round(decsec*100.0l);
	decsec /= 100.0l;
	// Write out RA and Dec with appropriate zero-padding.
	if(rahr<10) outstream2 << "0";
	outstream2 << rahr << " ";
	if(ramin<10) outstream2 << "0";
	outstream2 << ramin << " ";
	if(rasec<10.0l) outstream2 << "0";
	outstream2 << fixed << setprecision(3) << rasec << signstring;
	if(decdeg<10) outstream2 << "0";
	outstream2 << decdeg << " ";
	if(decmin<10) outstream2 << "0";
	outstream2 << decmin << " ";
	if(decsec<10.0l) outstream2 << "0";
	outstream2 << fixed << setprecision(2) << decsec << "         ";
	// Write out magnitude and band.
	outstream2 << fixed << setprecision(1) << clustvec[i].mag << " " << clustvec[i].band;
	// Add correct number of spaces after band.
	bandlen = j = 0;
	while(j<MINSTRINGLEN && clustvec[i].band[j]!='\0') {
	  bandlen++;
	  j++;
	}	    
	for(j=0;j<7-bandlen;j++) outstream2 << " ";
	// Write out obscode
	outstream2 << clustvec[i].obscode << "\n";
      }
    }
  }
  outstream3.close();
  // Now write the data for de-duplicated clusters that were not cross-linked
  for(clustct=0 ; clustct<clustnum; clustct++) {
    if(cluster_index[clustct].size()==1) {
      totalct++;
      // Load cluster-specific detection vector for this cluster
      clustvec={};
      detct=0;
      clustct2 = cluster_index[clustct][0];
      while(detct<detnum && cluster_detvec_main[detct].index < clustct2) detct++;
      while(detct<detnum && cluster_detvec_main[detct].index == clustct2) {
	clustvec.push_back(cluster_detvec_main[detct]);
	detct++;
      }
      // Time-sort clustvec
      sort(clustvec.begin(), clustvec.end(), early_hldet());
      // Run some analytics on this detection cluster
      avg_det_qual = 0.0l;
      max_known_obj=0;
      for(i=0; i<long(clustvec.size()); i++) {
	avg_det_qual += double(clustvec[i].det_qual);
	if(clustvec[i].known_obj > max_known_obj) max_known_obj = clustvec[i].known_obj;
      }
      avg_det_qual/=double(clustvec.size());
      // Loop over clustvec to extract individual tracklets
      trackvec={};
      trackvec.push_back(clustvec[0]);
      angvelvec = GCRvec = PAvec = timespanvec = arcvec = {};
      for(i=1; i<long(clustvec.size()); i++) {
	if((clustvec[i].MJD - clustvec[i-1].MJD) < nightstep) {
	  // Add a new point to this tracklet
	  trackvec.push_back(clustvec[i]);
	} else {
	  // A tracklet is finished. Analyze it.
	  tracknum = trackvec.size();
	  if(tracknum>1) {
	    greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
	    angvelvec.push_back(angvel);
	    GCRvec.push_back(sqrt(DSQUARE(crosstrack)+DSQUARE(alongtrack)));
	    PAvec.push_back(PA);
	    timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
	    arc = timespan*angvel;
	    arcvec.push_back(arc*3600.0l);
	    timespanvec.push_back(timespan*24.0l);
	    // Wipe trackvec, and load the next point of the next tracklet.
	    trackvec = {};
	    trackvec.push_back(clustvec[i]);
	  }
	}
      }
      if(tracknum>1) {
	// Handle a final tracklet.
	greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
	angvelvec.push_back(angvel);
	GCRvec.push_back(sqrt(DSQUARE(crosstrack)+DSQUARE(alongtrack)));
	PAvec.push_back(PA);
	timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
	arc = timespan*angvel;
	arcvec.push_back(arc*3600.0l);
	timespanvec.push_back(timespan*24.0l);
	// Wipe trackvec
	trackvec = {};
      }
      // Sort all of the tracklet statistics vectors
      tracknum = angvelvec.size();
      sort(angvelvec.begin(), angvelvec.end());
      sort(GCRvec.begin(), GCRvec.end());
      sort(PAvec.begin(), PAvec.end());
      sort(timespanvec.begin(), timespanvec.end());
      sort(arcvec.begin(), arcvec.end());
    
      outstream1 << "\n#clusternum,posRMS,velRMS,totRMS,astromRMS,timespan,uniquepoints,obsnights,metric,orbit_a,orbit_e,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count,avg_det_qual,max_known_obj,minvel,maxvel,minGCR,maxGCR,minpa,maxpa,mintimespan,maxtimespan,minarc,maxarc,contributors\n";
      outstream1 << fixed << setprecision(3) << inclustvec[clustct].clusternum << "," << inclustvec[clustct].posRMS << "," << inclustvec[clustct].velRMS << "," << inclustvec[clustct].totRMS << ",";
      outstream1 << fixed << setprecision(4) << inclustvec[clustct].astromRMS << ",";
      outstream1 << fixed << setprecision(6) << inclustvec[clustct].timespan << "," << inclustvec[clustct].uniquepoints << "," << inclustvec[clustct].obsnights << "," << inclustvec[clustct].metric << ",";
      outstream1 << fixed << setprecision(6) << inclustvec[clustct].orbit_a << "," << inclustvec[clustct].orbit_e << "," << inclustvec[clustct].orbit_MJD << ",";
      outstream1 << fixed << setprecision(1) << inclustvec[clustct].orbitX << "," << inclustvec[clustct].orbitY << "," << inclustvec[clustct].orbitZ << ",";
      outstream1 << fixed << setprecision(4) << inclustvec[clustct].orbitVX << "," << inclustvec[clustct].orbitVY << "," << inclustvec[clustct].orbitVZ << "," << inclustvec[clustct].orbit_eval_count << ",";
      outstream1 << fixed << setprecision(1) << avg_det_qual << "," << max_known_obj << ",";
      outstream1 << fixed << setprecision(6) << angvelvec[0] << "," << angvelvec[tracknum-1] << "," << GCRvec[0] << "," << GCRvec[tracknum-1] << "," << PAvec[0] << "," << PAvec[tracknum-1] << "," << timespanvec[0] << "," << timespanvec[tracknum-1] << "," << arcvec[0] << "," << arcvec[tracknum-1] << "," << cluster_index[clustct].size() << "\n";
    
      // Write this data to the output file
      outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
      for(i=0; i<long(clustvec.size()); i++) {
	outstream1 << fixed << setprecision(7) << clustvec[i].MJD << "," << clustvec[i].RA << "," << clustvec[i].Dec << ",";
	outstream1 << fixed << setprecision(4) << clustvec[i].mag << ",";
	outstream1 << fixed << setprecision(2) << clustvec[i].trail_len << "," << clustvec[i].trail_PA << ",";
	outstream1 << fixed << setprecision(4) << clustvec[i].sigmag << ",";
	outstream1 << fixed << setprecision(3) << clustvec[i].sig_across << "," << clustvec[i].sig_along << ",";
	outstream1 << clustvec[i].image << "," << clustvec[i].idstring << "," << clustvec[i].band << ",";
	outstream1 << clustvec[i].obscode << "," << clustvec[i].known_obj << ",";
	outstream1 << clustvec[i].det_qual << "," << clustvec[i].index << "\n";
      }
    	
      outstream2 << "\n#clusternum,posRMS,velRMS,totRMS,astromRMS,timespan,uniquepoints,obsnights,metric,orbit_a,orbit_e,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count,avg_det_qual,max_known_obj,minvel,maxvel,minGCR,maxGCR,minpa,maxpa,mintimespan,maxtimespan,minarc,maxarc,contributors\n";
      outstream2 << fixed << setprecision(3) << inclustvec[clustct].clusternum << "," << inclustvec[clustct].posRMS << "," << inclustvec[clustct].velRMS << "," << inclustvec[clustct].totRMS << ",";
      outstream2 << fixed << setprecision(4) << inclustvec[clustct].astromRMS << ",";
      outstream2 << fixed << setprecision(6) << inclustvec[clustct].timespan << "," << inclustvec[clustct].uniquepoints << "," << inclustvec[clustct].obsnights << "," << inclustvec[clustct].metric << ",";
      outstream2 << fixed << setprecision(6) << inclustvec[clustct].orbit_a << "," << inclustvec[clustct].orbit_e << "," << inclustvec[clustct].orbit_MJD << ",";
      outstream2 << fixed << setprecision(1) << inclustvec[clustct].orbitX << "," << inclustvec[clustct].orbitY << "," << inclustvec[clustct].orbitZ << ",";
      outstream2 << fixed << setprecision(4) << inclustvec[clustct].orbitVX << "," << inclustvec[clustct].orbitVY << "," << inclustvec[clustct].orbitVZ << "," << inclustvec[clustct].orbit_eval_count << ",";
      outstream2 << fixed << setprecision(1) << avg_det_qual << "," << max_known_obj << ",";
      outstream2 << fixed << setprecision(6) << angvelvec[0] << "," << angvelvec[tracknum-1] << "," << GCRvec[0] << "," << GCRvec[tracknum-1] << "," << PAvec[0] << "," << PAvec[tracknum-1] << "," << timespanvec[0] << "," << timespanvec[tracknum-1] << "," << arcvec[0] << "," << arcvec[tracknum-1] << "," << cluster_index[clustct].size() << "\n";
    
      // Write this data to the output file
      outstream2 << "#MPC 80-column formatted observations:\n";
      for(i=0; i<long(clustvec.size()); i++) {
	// The temporary ID will just be 'a' followed by a 6-digit number
	if(totalct<=999999) stest1 = to_string(totalct);
	else {
	  stest1 = to_string(totalct%1000000);
	}
	while(stest1.size()<6) {
	  stest1 = "0" + stest1;
	}
	outstream2 << "     a" << stest1 << "  C";
	mjd2mpcdate(clustvec[i].MJD,year,month,day);
	outstream2 << year << " ";
	if(month<10) outstream2 << "0";
	outstream2 << month << " ";
	day = round(day*1000000.0l);
	day /= 1000000.0l;
	if(day<10.0l) outstream2  << fixed << setprecision(6) << "0";
	outstream2  << fixed << setprecision(6) << day;
	// Convert RA, Dec from decimal degrees to sexagesimal format.
	rahr = int(clustvec[i].RA/15.0l);
	ramin = int(clustvec[i].RA*4.0l - double(rahr)*60.0l);
	rasec = clustvec[i].RA*240.0l - double(rahr)*3600.0l - double(ramin)*60.0l;
	rasec = round(rasec*1000.0l);
	rasec /= 1000.0l;
	if(clustvec[i].Dec>=0) {
	  signstring="+";
	  Dec = clustvec[i].Dec;
	} else {
	  signstring="-";
	  Dec = -clustvec[i].Dec;
	}
	decdeg = int(Dec);
	decmin = int(Dec*60.0l - double(decdeg)*60.0l);
	decsec = Dec*3600.0l - double(decdeg)*3600.0l - double(decmin)*60.0l;
	decsec = round(decsec*100.0l);
	decsec /= 100.0l;
	// Write out RA and Dec with appropriate zero-padding.
	if(rahr<10) outstream2 << "0";
	outstream2 << rahr << " ";
	if(ramin<10) outstream2 << "0";
	outstream2 << ramin << " ";
	if(rasec<10.0l) outstream2 << "0";
	outstream2 << fixed << setprecision(3) << rasec << signstring;
	if(decdeg<10) outstream2 << "0";
	outstream2 << decdeg << " ";
	if(decmin<10) outstream2 << "0";
	outstream2 << decmin << " ";
	if(decsec<10.0l) outstream2 << "0";
	outstream2 << fixed << setprecision(2) << decsec << "         ";
	// Write out magnitude and band.
	outstream2 << fixed << setprecision(1) << clustvec[i].mag << " " << clustvec[i].band;
	// Add correct number of spaces after band.
	bandlen = j = 0;
	while(j<MINSTRINGLEN && clustvec[i].band[j]!='\0') {
	  bandlen++;
	  j++;
	}	    
	for(j=0;j<7-bandlen;j++) outstream2 << " ";
	// Write out obscode
	outstream2 << clustvec[i].obscode << "\n";
      }
    }
  }
  cout << "Wrote out " << mergecount << " merged linkages and " << totalct << " total linkages\n";


  
  outstream1.close();
  outstream2.close();
  return(0);
}
