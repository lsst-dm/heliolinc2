// April 26, 2023: parse_clust2det_new

#include "solarsyst_dyn_geo01.h"
#include "cmath"


static void show_usage()
{
  cerr << "Usage: parse_clust2det_new -pairdet pairdet_file -insum input cluster summary file -clust2det input cluster-to-detection file -trackdiv tracklet_division_time -out output file\n";
}

int main(int argc, char *argv[])
{
  string sumfile,clust2detfile,outfile;
  vector <hlclust> inclustvec;
  vector  <longpair> inclust2det;
  vector <hldet> detvec;
  vector <hldet> cluster_detvec;
  vector <hldet> clustvec;
  vector <hldet> trackvec;
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
  int status=0;
  long detnum=0;
  long detct=0;
  long clustct=0;
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

  if(argc<9) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-pairdets" || string(argv[i]) == "-pairdet" || string(argv[i]) == "-pd" || string(argv[i]) == "-pdet" || string(argv[i]) == "--pairdet" || string(argv[i]) == "--paireddetections" || string(argv[i]) == "--pairdetfile" || string(argv[i]) == "--pairdetections") {
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
    } else if(string(argv[i]) == "-insum" || string(argv[i]) == "-inclust" || string(argv[i]) == "-clust" || string(argv[i]) == "-sum" || string(argv[i]) == "--input_summary" || string(argv[i]) == "--input_cluster" || string(argv[i]) == "--input_cluster_file" || string(argv[i]) == "--input_summary_file") {
      if(i+1 < argc) {
	//There is still something to read;
	sumfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input cluster summary file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clust2det" || string(argv[i]) == "-c2d" || string(argv[i]) == "-inc2d" || string(argv[i]) == "-input_c2d" || string(argv[i]) == "--input_clust2det" ) {
      if(i+1 < argc) {
	//There is still something to read;
	clust2detfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input cluster list keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trackdiv" || string(argv[i]) == "-nightstep") {
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
    }  else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
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
  if(pairdetfile.size()<=0) {
    cout << "\nERROR: input paired detection file is required\n";
    show_usage();
    return(1);
  } else if(sumfile.size()<=0) {
    cout << "\nERROR: input cluster summary file is required\n";
    show_usage();
    return(1);
  } else if(clust2detfile.size()<=0) {
    cout << "\nERROR: input cluster-to-detection file is required\n";
    show_usage();
    return(1);
  } else if(outfile.size()<=0) {
    cout << "\nERROR: output filename is required\n";
    show_usage();
    return(1);
  }

  cout << "input paired detection file " << pairdetfile << "\n";
  cout << "input cluster summary file " << sumfile << "\n";
  cout << "input cluster-to-detection file " << clust2detfile << "\n";
  cout << "output file " << outfile << "\n";


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

  status = parse_clust2det(detvec, inclust2det, cluster_detvec);
  if(status!=0) {
    cerr << "ERROR: parse_clust2det failed with error status " << status << "\n";
    return(status);
  }
  detnum = cluster_detvec.size();
  cout << "Wrote cluster detection vector with " << detnum << "entries\n";
  
  outstream1.open(outfile);
  cout << "Writing " << inclustvec.size() << " clusters to output file " << outfile << "\n";
  detct=0;
  for(clustct=0 ; clustct<long(inclustvec.size()); clustct++) {
    // Load cluster-specific detection vector for this cluster
    clustvec={};
    if(detct<detnum) {
      if(cluster_detvec[detct].index!=clustct) {
	cerr << "ERROR: cluster counting mismatch, clustct = " << clustct << ", detct = " << detct << ", index = " << cluster_detvec[detct].index << "\n";
	return(2);
      }
      while(detct<detnum && cluster_detvec[detct].index==clustct) {
	clustvec.push_back(cluster_detvec[detct]);
	detct++;
      }
      cout << clustvec.size() << " points found for cluster " << clustct << "\n";
    }
    if(clustvec.size()>0) {
      // Time-sort clustvec, just to be sure
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

      outstream1 << "\n#clusternum,posRMS,velRMS,totRMS,astromRMS,timespan,uniquepoints,obsnights,metric,orbit_a,orbit_e,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count,avg_det_qual,max_known_obj,minvel,maxvel,minGCR,maxGCR,minpa,maxpa,mintimespan,maxtimespan,minarc,maxarc\n";
      outstream1 << fixed << setprecision(3) << inclustvec[clustct].clusternum << "," << inclustvec[clustct].posRMS << "," << inclustvec[clustct].velRMS << "," << inclustvec[clustct].totRMS << ",";
      outstream1 << fixed << setprecision(4) << inclustvec[clustct].astromRMS << ",";
      outstream1 << fixed << setprecision(6) << inclustvec[clustct].timespan << "," << inclustvec[clustct].uniquepoints << "," << inclustvec[clustct].obsnights << "," << inclustvec[clustct].metric << ",";
      outstream1 << fixed << setprecision(6) << inclustvec[clustct].orbit_a << "," << inclustvec[clustct].orbit_e << "," << inclustvec[clustct].orbit_MJD << ",";
      outstream1 << fixed << setprecision(1) << inclustvec[clustct].orbitX << "," << inclustvec[clustct].orbitY << "," << inclustvec[clustct].orbitZ << ",";
      outstream1 << fixed << setprecision(4) << inclustvec[clustct].orbitVX << "," << inclustvec[clustct].orbitVY << "," << inclustvec[clustct].orbitVZ << "," << inclustvec[clustct].orbit_eval_count << ",";
      outstream1 << fixed << setprecision(1) << avg_det_qual << "," << max_known_obj << ",";
      outstream1 << fixed << setprecision(6) << angvelvec[0] << "," << angvelvec[tracknum-1] << "," << GCRvec[0] << "," << GCRvec[tracknum-1] << "," << PAvec[0] << "," << PAvec[tracknum-1] << "," << timespanvec[0] << "," << timespanvec[tracknum-1] << "," << arcvec[0] << "," << arcvec[tracknum-1] << "\n";
	
      // Write this data to the output file
      outstream1 << "#MPC 80-column formatted observations:\n";
      for(i=0; i<long(clustvec.size()); i++) {
	// The temporary ID will just be 'a' followed by a 6-digit number
	if(clustvec[i].index<=999999) stest1 = to_string(clustvec[i].index);
	else {
	  stest1 = to_string(clustvec[i].index%1000000);
	}
	while(stest1.size()<6) {
	  stest1 = "0" + stest1;
	}
	outstream1 << "     a" << stest1 << "  C";
	mjd2mpcdate(clustvec[i].MJD,year,month,day);
	outstream1 << year << " ";
	if(month<10) outstream1 << "0";
	outstream1 << month << " ";
	day = round(day*1000000.0l);
	day /= 1000000.0l;
	if(day<10.0l) outstream1  << fixed << setprecision(6) << "0";
	outstream1  << fixed << setprecision(6) << day;
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
	if(rahr<10) outstream1 << "0";
	outstream1 << rahr << " ";
	if(ramin<10) outstream1 << "0";
	outstream1 << ramin << " ";
	if(rasec<10.0l) outstream1 << "0";
	outstream1 << fixed << setprecision(3) << rasec << signstring;
	if(decdeg<10) outstream1 << "0";
	outstream1 << decdeg << " ";
	if(decmin<10) outstream1 << "0";
	outstream1 << decmin << " ";
	if(decsec<10.0l) outstream1 << "0";
	outstream1 << fixed << setprecision(2) << decsec << "         ";
	// Write out magnitude and band.
	outstream1 << fixed << setprecision(1) << clustvec[i].mag << " " << clustvec[i].band;
	// Add correct number of spaces after band.
	bandlen = j = 0;
	while(j<MINSTRINGLEN && clustvec[i].band[j]!='\0') {
	  bandlen++;
	  j++;
	}	    
	for(j=0;j<7-bandlen;j++) outstream1 << " ";
	// Write out obscode
	outstream1 << clustvec[i].obscode << "\n";
      }
    }
  }
  outstream1.close();
  
  return(0);
}
