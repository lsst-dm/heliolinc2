// April 26, 2023: parse_clust2det_new

#include "solarsyst_dyn_geo01.h"
#include "cmath"


static void show_usage()
{
  cerr << "Usage: parse_trk2det_MPC80 -pairdet pairdet_file -trk2det input tracklet-to-detection file -out output file\n";
}

int main(int argc, char *argv[])
{
  string trk2detfile,outfile;
  vector  <longpair> intrk2det;
  vector <hldet> detvec;
  vector <hldet> trackvec;
  double angvel,crosstrack,alongtrack,PA,poleRA,poleDec;
  angvel=crosstrack=alongtrack=PA=poleRA=poleDec=0.0l;
  double arc,timespan;
  arc = timespan = 0;
  double GCR = 0.0l;
  vector <double> magvec;
  long tracknum=0;
  string pairdetfile,stest;
  ifstream instream1;
  ofstream outstream1;
  int status=0;
  long detnum=0;
  long detct=0;
  long trkct=0;
  long i=0;
  int verbose=0;
  long max_known_obj=0;
  double avg_det_qual=0.0l;
  double magrange,magmean,magrms;
  magrange = magmean = magrms = 0.0l;
  
  
  if(argc!=7) {
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
    } else if(string(argv[i]) == "-trk2det" || string(argv[i]) == "-t2d" || string(argv[i]) == "-tracklet2detection" || string(argv[i]) == "-trk2detfile" || string(argv[i]) == "--tracklet2detection") {
      if(i+1 < argc) {
	//There is still something to read;
	trk2detfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input trk2det file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "--outfile") {
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
  } else if(trk2detfile.size()<=0) {
    cout << "\nERROR: input tracklet-to-detection file is required\n";
    show_usage();
    return(1);
  } else if(outfile.size()<=0) {
    cout << "\nERROR: output filename is required\n";
    show_usage();
    return(1);
  }

  cout << "input paired detection file " << pairdetfile << "\n";
  cout << "input tracklet-to-detection file " << trk2detfile << "\n";
  cout << "output file " << outfile << "\n";


  // Read paired detection file
  detvec={};
  status=read_pairdet_file(pairdetfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << pairdetfile << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  detnum = detvec.size();
  cout << "Read " << detnum << " data lines from paired detection file " << pairdetfile << "\n";
  
  intrk2det={};
  // Read tracklet-to-detection file
  status=read_longpair_file(trk2detfile, intrk2det, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read tracklet-to-detection file " << trk2detfile << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
    return(1);
  }
  cout << "Read " << intrk2det.size() << " data lines from tracklet-to-detection file " << trk2detfile << "\n";

  outstream1.open(outfile);
  detct=0;
  trkct=0;
  // Load tracklet-specific detection vector for this cluster
  trackvec={};
  for(i=0 ; i<long(intrk2det.size()); i++) {
    detct = intrk2det[i].i2;
    cout << "i=" << i << ", trk2det: " << intrk2det[i].i1 << "," << intrk2det[i].i2 << " trkct=" << trkct << ", detct=" << detct << ", detnum=" << detnum << "\n";
    if(detct>=0 && detct<detnum && trkct==intrk2det[i].i1) {
      trackvec.push_back(detvec[detct]);
    } else if(trkct!=intrk2det[i].i1 || i>=long(intrk2det.size()-1)) {
      // We just finished a tracklet. Write it out
      cout << trackvec.size() << " points found for tracklet " << trkct << "\n";
      if(trackvec.size()>0) {
	// Time-sort trackvec, just to be sure
	sort(trackvec.begin(), trackvec.end(), early_hldet());
	// Run some analytics on this tracklet
	avg_det_qual = 0.0l;
	max_known_obj=0;
	magvec={};
	for(long j=0; j<long(trackvec.size()); j++) {
	  avg_det_qual += double(trackvec[j].det_qual);
	  if(trackvec[j].known_obj > max_known_obj) max_known_obj = trackvec[j].known_obj;
	  if(trackvec[j].mag>0.0l) magvec.push_back(trackvec[j].mag);
	}
	avg_det_qual/=double(trackvec.size());
	angvel = GCR = PA = timespan = arc = 0.0l;
	tracknum = trackvec.size();
	if(tracknum>1) {
	  greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
	  if(tracknum>2) GCR = sqrt(DSQUARE(crosstrack)+DSQUARE(alongtrack));
	  else GCR = 0.0l;
	  timespan = (trackvec[trackvec.size()-1].MJD - trackvec[0].MJD)*24.0l;
	  arc = timespan*angvel*3600.0l;
	  // Sort magvec
	  sort(magvec.begin(), magvec.end());
	  dmeanrms01(magvec, &magmean, &magrms);
	  // Magrange will be the full max-min
	  magrange = magvec[magvec.size()-1] - magvec[0];
	} else if(tracknum==1) {
	  // The 'tracklet' is a singleton. Set all tracket vectors to error codes or zero.
	  angvel = -1.0l;
	  GCR = -1.0l;
	  PA = -999.0l;
	  arc = 0.0l;
	  timespan = 0.0l;
	  magmean = trackvec[0].mag;
	  magrms = magrange = 99.9;
	}
	// Write the summary data to the output file
	outstream1 << "\n#trknum,timespan,pointnum,angvel,PA,crosstrack,alongtrack,GCR,arc,avg_det_qual,max_known_obj,stringID,magmean,magrms,magrange\n";
	outstream1 << fixed << setprecision(6) << trkct << "," << timespan << "," << tracknum << ",";
	outstream1 << fixed << setprecision(3) << angvel << "," << PA << "," << crosstrack << "," << alongtrack << "," << GCR << "," << arc << ",";
	outstream1 << fixed << setprecision(1) << avg_det_qual << "," << max_known_obj << ",";	
	outstream1 << trackvec[0].idstring << "," << magmean << "," << magrms << "," << magrange << "\n";

	// Write the individual points to the output file in MPC80 format.
	string stest1;
	int year,month;
	double day;
	int rahr,ramin;
	int decdeg,decmin;
	double Dec,rasec,decsec;;
	string signstring;
	
	outstream1 << "#MPC 80-column formatted observations:\n";
	for(long j=0; j<long(trackvec.size()); j++) {
	  // The temporary ID will just be 'a' followed by a 6-digit number
	  if(trkct<=999999) stest1 = to_string(trkct);
	  else {
	    stest1 = to_string(trkct%1000000);
	  }
	  while(stest1.size()<6) {
	    stest1 = "0" + stest1;
	  }
	  outstream1 << "     a" << stest1 << "  C";
	  mjd2mpcdate(trackvec[j].MJD,year,month,day);
	  outstream1 << year << " ";
	  if(month<10) outstream1 << "0";
	  outstream1 << month << " ";
	  day = round(day*1000000.0l);
	  day /= 1000000.0l;
	  if(day<10.0l) outstream1  << fixed << setprecision(6) << "0";
	  outstream1  << fixed << setprecision(6) << day;
	  // Convert RA, Dec from decimal degrees to sexagesimal format.
	  rahr = int(trackvec[j].RA/15.0l);
	  ramin = int(trackvec[j].RA*4.0l - double(rahr)*60.0l);
	  rasec = trackvec[j].RA*240.0l - double(rahr)*3600.0l - double(ramin)*60.0l;
	  rasec = round(rasec*1000.0l);
	  rasec /= 1000.0l;
	  if(trackvec[j].Dec>=0) {
	    signstring="+";
	    Dec = trackvec[j].Dec;
	  } else {
	    signstring="-";
	    Dec = -trackvec[j].Dec;
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
	  outstream1 << fixed << setprecision(1) << trackvec[j].mag << " " << trackvec[j].band;
	  // Add correct number of spaces after band.
	  int bandlen,k;
	  bandlen = k = 0;
	  while(k<MINSTRINGLEN && trackvec[j].band[k]!='\0') {
	    bandlen++;
	    k++;
	  }    
	  for(k=0;k<7-bandlen;k++) outstream1 << " ";
	  // Write out obscode
	  outstream1 << trackvec[j].obscode << "\n";
	}
      }
      trackvec = {}; // Wipe trackvec
      if(detct>=0 && detct<detnum) {
	// Load the next point of the next tracklet.
	trkct=intrk2det[i].i1;
	trackvec.push_back(detvec[detct]);
      }
    }
  }
  outstream1.close();
  
  return(0);
}
