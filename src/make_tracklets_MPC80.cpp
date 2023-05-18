// May 10, 2023: make_tracklets_MPC80.cpp:
// Based on make_tracklets_new, but designed to ingest the Minor Planet
// Center's Isolated Tracklet File (itf.txt). The ITF data is already
// trackletized in MPC 80-column formate.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define NUMPOS 3

#define DEBUG 1
      
static void show_usage()
{
  cerr << "Usage: make_tracklets_MPC80 -dets detfile -outimgs output image file/ \n";
  cerr << "-pairdets paired detection file -tracklets tracklet file -trk2det output tracklet-to-detection file/ \n";
  cerr << "-earth earthfile -obscode obscodefile\n";
}


int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  vector <hldet> detvec_fixed = {};
  vector <hldet> pairdets = {};
  vector <observatory> observatory_list = {};
  vector <hlimage> img_log = {};
  vector <longpair> trk2det ={};
  vector <tracklet> tracklets={};
  vector <point3d> Earthpos;
  vector <point3d> Earthvel;
  vector <double> EarthMJD;
  string lnfromfile;
  int status = 0;
  long i = 0;
  int imct=0;
  string indetfile;
  string outimfile;
  string earthfile;
  string obscodefile;
  string pairdetfile="pairdetfile01.csv";
  string trackletfile="trackletfile01.csv";
  string trk2detfile="trk2detfile01.csv";
  ifstream instream1;
  ofstream outstream1;
  string stest;
  int outimfile_set;
  outimfile_set = 0;
  int pairdetfile_default,trackletfile_default,trk2detfile_default;  
  pairdetfile_default = trackletfile_default = trk2detfile_default = 1;
  int verbose=0;
  
  if(argc<7)
    {
      show_usage();
      return(1);
    }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-d" || string(argv[i]) == "-dets" || string(argv[i]) == "-det" || string(argv[i]) == "--dets" || string(argv[i]) == "--det" || string(argv[i]) == "--detection" || string(argv[i]) == "--detections") {
      if(i+1 < argc) {
	//There is still something to read;
	indetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outim" || string(argv[i]) == "-outimg" || string(argv[i]) == "-outimgs" || string(argv[i]) == "--outimages" || string(argv[i]) == "--outimage" || string(argv[i]) == "--outimgs" || string(argv[i]) == "--outimg") {
      if(i+1 < argc) {
	//There is still something to read;
	outimfile=argv[++i];
	outimfile_set = 1;
	i++;
      }
      else {
	cerr << "Output image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" || string(argv[i]) == "-detpairs") {
      if(i+1 < argc) {
	//There is still something to read;
        pairdetfile=argv[++i];
	pairdetfile_default = 0;
	i++;
      }
      else {
	cerr << "Output paired detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-tracklets" || string(argv[i]) == "-tracklet" || string(argv[i]) == "-pair" || string(argv[i]) == "-pairfile" || string(argv[i]) == "-pairs") {
      if(i+1 < argc) {
	//There is still something to read;
	trackletfile=argv[++i];
	trackletfile_default = 0;
	i++;
      }
      else {
	cerr << "Output pair file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trk2det" || string(argv[i]) == "-t2d" || string(argv[i]) == "-tracklet2detection" || string(argv[i]) == "-trk2detfile" || string(argv[i]) == "--tracklet2detection") {
      if(i+1 < argc) {
	//There is still something to read;
	trk2detfile=argv[++i];
	trk2detfile_default = 0;
	i++;
      }
      else {
	cerr << "Output pair file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-earth" || string(argv[i]) == "-e" || string(argv[i]) == "-Earth" || string(argv[i]) == "--earthfile" || string(argv[i]) == "--Earthfile" || string(argv[i]) == "--earth" || string(argv[i]) == "--Earth") {
      if(i+1 < argc) {
	//There is still something to read;
	earthfile=argv[++i];
	i++;
      }
      else {
	cerr << "Earth file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obscode" || string(argv[i]) == "-obs" || string(argv[i]) == "-oc" || string(argv[i]) == "-obscodes" || string(argv[i]) == "--obscode" || string(argv[i]) == "--obscodes" || string(argv[i]) == "--observatorycodes") {
      if(i+1 < argc) {
	//There is still something to read;
	obscodefile=argv[++i];
	i++;
      }
      else {
	cerr << "Observatory code file keyword supplied with no corresponding argument\n";
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
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }

  if(indetfile.size()<=0) {
    cerr << "Please supply an input detection file:\n\n";
    show_usage();
    return(1);
  }

  if(earthfile.size()<=0) {
    cerr << "Please supply a heliocentric ephemeris file for the Earth:\n\n";
    show_usage();
    return(1);
  }
  
  if(obscodefile.size()<=0) {
    cerr << "Please supply an observatory code file:\n\n";
    show_usage();
    return(1);
  }
  
  cout << "\nInput detection file is called " << indetfile << "\n\n";
  
  if(outimfile_set ==1) cout << "output image file = " << outimfile << "\n";
  else cout << "No output image file specified: required image information will only be used internally.\n";
  
  cout << "Observatory code file " << obscodefile << "\n";

  if(pairdetfile_default == 0) cout << "Output paired detection file will be called " << pairdetfile << "\n";
  else cout << "Defaulting to output paired detection file name = " << pairdetfile << "\n";

  if(trackletfile_default == 0) cout << "Output tracklet file will be called " << trackletfile << "\n";
  else cout << "Defaulting to output tracklet file name = " << trackletfile << "\n";
 
  if(trk2detfile_default == 0) cout << "Output tracklet-to-detection mapping file will be called " << trk2detfile << "\n";
  else cout << "Defaulting to output tracklet-to-detection file name = " << trk2detfile << "\n";
  
  cout << "Heliocentric ephemeris file for the Earth: " << earthfile << "\n";

  // Read observatory code file
  status = read_obscode_file2(obscodefile, observatory_list, verbose);
  if(status!=0) {
    cerr << "ERROR reading observatory code file " << obscodefile << "\n";
    return(1);
  }
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile << "\n";
  if(DEBUG>=2) {
    for(i=0;i<long(observatory_list.size());i++) {
      cout << observatory_list[i].obscode << " " << observatory_list[i].obslon << " " << observatory_list[i].plxcos << " " << observatory_list[i].plxsin << "\n";
    }
  }

  // Read input detection file.
  status = read_detection_file_MPC80(indetfile, detvec);
  
  if(status==0) { 
    cout << "Input file " << indetfile << " read successfully to the end.\n";
  }
  else if(status==1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(status==2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else {
    cerr << "Warning: unknown file read problem\n";
    return(3);
  }

  // Record indices to preserve the original order after a sort
  for(i=0; i<long(detvec.size()); i++) {
    detvec[i].index=i;
    if(detvec[i].index<0 || detvec[i].index>=long(detvec.size())) {
      cerr << "Error loading detvec.index, i, index = " << i << " " << detvec[i].index << "\n";
      return(2);
    }
  }
  
  // Save a native-order version of the detection vector
  detvec_fixed = detvec;
  
  for(i=0; i<long(detvec.size()); i++) {
    if(detvec[i].index<0 || detvec[i].index>=long(detvec.size())) {
      cerr << "Pre-sort sanity check fail: detvec.index, i, index = " << i << " " << detvec[i].index << "\n";
      return(3);
    }
  }

  // time-sort the detection vector
  if(DEBUG==1) cout << "Preparing to time-sort the detection vector\n";
  sort(detvec.begin(), detvec.end(), early_hldet());

  for(i=0; i<long(detvec.size()); i++) {
    if(detvec[i].index<0 || detvec[i].index>=long(detvec.size())) {
      cerr << "Post-sort sanity check fail: detvec.index, i, index = " << i << " " << detvec[i].index << "\n";
      return(4);
    }
  }
 
  EarthMJD={};
  Earthpos={};
  Earthvel={};
  read_horizons_csv(earthfile,EarthMJD,Earthpos,Earthvel);
  cout << "Finished reading heliocentric ephemeris file " << earthfile << " for the Earth.\n";

  if(DEBUG==1) cout << "Preparing to load the image table\n";
  status = load_image_table(img_log, detvec, observatory_list, EarthMJD, Earthpos, Earthvel);
  if(DEBUG==1) cout << "Loaded the image table\n";

  for(i=0; i<long(detvec.size()); i++) {
    if(detvec[i].index<0 || detvec[i].index>=long(detvec.size())) {
      cerr << "Post-image sanity check fail: detvec.index, i, index = " << i << " " << detvec[i].index << "\n";
      return(5);
    }
  }
  
  if(DEBUG>=2) {
    // Test: print out time-sorted detection table.
    outstream1.open("testjunk01.txt");
    for(i=0;i<long(detvec.size());i++) {
      outstream1 << fixed << setprecision(6) << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec <<  " " << detvec[i].mag << " " << detvec[i].band <<  " " << detvec[i].obscode << "\n";
    }
    outstream1.close();
  }
  
  if(outimfile.size()>0)
    {
      // Write and print image log table
      outstream1.open(outimfile);
      for(imct=0;imct<long(img_log.size());imct++)
	{
	  outstream1 << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
	  outstream1 << fixed << setprecision(6) << " " << img_log[imct].Dec << " " << img_log[imct].obscode << " ";
	  outstream1 << fixed << setprecision(1) << img_log[imct].X << " " << img_log[imct].Y << " " << img_log[imct].Z << " ";
	  outstream1 << fixed << setprecision(4) << img_log[imct].VX << " " << img_log[imct].VY << " " << img_log[imct].VZ << " ";
	  outstream1 << img_log[imct].startind << " " << img_log[imct].endind << "\n";
	}
      outstream1.close();
    }

  cout << "Running remake_tracklets\n";
  remake_tracklets(detvec, detvec_fixed, img_log, tracklets, trk2det, verbose);
  cout << "Finished remake_tracklets\n";

  // Write paired detection file
  cout << "Writing paired detection file with " << detvec.size() << " lines\n";
  outstream1.open(pairdetfile);
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(i=0;i<long(detvec.size());i++) {
    outstream1 << fixed << setprecision(7) << detvec[i].MJD << "," << detvec[i].RA << "," << detvec[i].Dec << ",";
    outstream1 << fixed << setprecision(4) << detvec[i].mag << ",";
    outstream1 << fixed << setprecision(2) << detvec[i].trail_len << "," << detvec[i].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << detvec[i].sigmag << ",";
    outstream1 << fixed << setprecision(3) << detvec[i].sig_across << "," << detvec[i].sig_along << ",";
    outstream1 << detvec[i].image << "," << detvec[i].idstring << "," << detvec[i].band << ",";
    outstream1 << detvec[i].obscode << "," << detvec[i].known_obj << ","; 
    outstream1 << detvec[i].det_qual << "," << detvec[i].index << "\n"; 
  }
  outstream1.close();

  // Write tracklet file
  cout << "Writing tracklet file with " << tracklets.size() << " lines\n";
  outstream1.open(trackletfile);
  outstream1 << "#Image1,RA1,Dec1,Image2,RA2,Dec2,npts,trk_ID\n";
  for(i=0;i<long(tracklets.size());i++) {
    outstream1 << fixed << setprecision(7) << tracklets[i].Img1 << "," << tracklets[i].RA1 << "," << tracklets[i].Dec1 << ",";
    outstream1 << fixed << setprecision(7) << tracklets[i].Img2 << "," << tracklets[i].RA2 << "," << tracklets[i].Dec2 << ",";
    outstream1 << tracklets[i].npts << "," << tracklets[i].trk_ID << "\n"; 
  }
  outstream1.close();

   // Write trk2det file
  cout << "Writing trk2det file with " << trk2det.size() << " lines\n";
  outstream1.open(trk2detfile);
  outstream1 << "#trk_ID,detnum\n";
  for(i=0;i<long(trk2det.size());i++) {
    outstream1 << trk2det[i].i1 << "," << trk2det[i].i2 << "\n"; 
  }
  outstream1.close();

  return(0);
}

