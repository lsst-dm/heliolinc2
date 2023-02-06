// February 01, 2023: mpcobs_match01.cpp

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define MAGCOL 6
#define OBSCODECOL 8
#define COLS_TO_READ 5

#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds

static void show_usage()
{
  cerr << "Usage: mpcobs_match -mpcfile mpc_file -detfile detection_file -colformat colformat_file -mjdrange MJDstart MJDend -qrange query_range -mpcout mpc_output -detout detection_output -mpcNM mpc_nomatch -detNM det_nomatch\n";
}

int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  char c='0';
  int reachedeof=0;
  int lct=0;
  string mpcfile,detfile,mpcout,detout,mpc_nomatch,det_nomatch;
  string colformatfile,stest;
  double mjdstart=0.0;
  double mjdend=1.0e7;
  ifstream instream1;
  ofstream outstream1;
  ofstream outstream2;
  string lnfromfile;
  string object,band;
  string obscode;
  string stringobscode;
  double MJD, RA, Dec, mag;
  vector <double> mpc_MJD, mpc_RA, mpc_Dec, mpc_mag;
  vector <string> mpc_obscode;
  vector <string> mpc_object;
  vector <double> det_MJD, det_RA, det_Dec, det_mag;
  vector <string> det_obscode;
  vector <string> mpclines;
  vector <string> detlines;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int magcol = MAGCOL;
  int obscodecol = OBSCODECOL;
  int mjdread,raread,decread,magread,obscoderead;
  mjdread = raread = decread = magread = obscoderead = 0;
  int colreadct=0;
  point3d_index point0 = point3d_index(0l,0l,0l,0);
  vector <point3d_index> mpcvec;
  vector <point3d_index> detvec;
  KD_point3d_index kdroot = KD_point3d_index(point0,0,0,1,0);
  vector <KD_point3d_index> det_kdtree;
  KD_point3d_index lp1 = KD_point3d_index(point0,0,0,1,0);
  KD_point3d_index rp1 = KD_point3d_index(point0,0,0,1,0);
  int dim=1;
  long medpt;
  vector <long> indexvec;
  double qrange = 0.0005l;
  int matchpoint,matchnum;
  double timetol = IMAGETIMETOL/SOLARDAY;
  
  if(argc<0 && argc<1) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-mpcfile" || string(argv[i]) == "-mpcin" || string(argv[i]) == "-mpcobs" || string(argv[i]) == "-mpcinfile" || string(argv[i]) == "--mpcobs" || string(argv[i]) == "--mpcfile" || string(argv[i]) == "--mpcin" || string(argv[i]) == "--mpcinfile") {
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
    } else if(string(argv[i]) == "-detfile" || string(argv[i]) == "-detin" || string(argv[i]) == "-di" || string(argv[i]) == "-detinfile" || string(argv[i]) == "--detfile" || string(argv[i]) == "--detectionfile" || string(argv[i]) == "--detection_input" || string(argv[i]) == "--detection_input_file") {
      if(i+1 < argc) {
	//There is still something to read;
	detfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format"  || string(argv[i]) == "-col" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" || string(argv[i]) == "--cformat") {
      if(i+1 < argc) {
	//There is still something to read;
	colformatfile=argv[++i];
	i++;
      }
      else {
	cerr << "Column format file keyword supplied with no corresponding argument\n";
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
    }  else if(string(argv[i]) == "-qrange" || string(argv[i]) == "-range" || string(argv[i]) == "-queryrange" || string(argv[i]) == "-QR" || string(argv[i]) == "-qr" || string(argv[i]) == "--queryrange") {
      if(i+1 < argc) {
	//There is still something to read;
	qrange = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input MJD range keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-detout" || string(argv[i]) == "-detfileout" || string(argv[i]) == "-do" || string(argv[i]) == "-detoutfile" || string(argv[i]) == "--detout" || string(argv[i]) == "--detection_outfile" || string(argv[i]) == "--detection_output" || string(argv[i]) == "--detection_output_file") {
      if(i+1 < argc) {
	//There is still something to read;
	detout=argv[++i];
	i++;
      }
      else {
	cerr << "Output detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-mpcNM" || string(argv[i]) == "-mpcNM" || string(argv[i]) == "-mpc_nomatch") {
      if(i+1 < argc) {
	//There is still something to read;
	mpc_nomatch=argv[++i];
	i++;
      }
      else {
	cerr << "Output mpc unmatched file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-detNM" || string(argv[i]) == "-det_nomatch") {
      if(i+1 < argc) {
	//There is still something to read;
	det_nomatch=argv[++i];
	i++;
      }
      else {
	cerr << "Output detection unmatched file keyword supplied with no corresponding argument\n";
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
  if(detfile.size()<=0) {
    cerr << "ERROR: you have not specified an input detection file\n";
    show_usage();
    return(1);
  }
  if(mpcout.size()<=0) {
    cerr << "ERROR: you have not specified an output matched mpc file\n";
    show_usage();
    return(1);
  }
   if(detout.size()<=0) {
    cerr << "ERROR: you have not specified an output matched detection file\n";
    show_usage();
    return(1);
  }
   if(mpc_nomatch.size()<=0) {
    cerr << "ERROR: you have not specified an output unmatched mpc file\n";
    show_usage();
    return(1);
  }
  if(det_nomatch.size()<=0) {
    cerr << "ERROR: you have not specified an output unmatched detection file\n";
    show_usage();
    return(1);
  }

  // Read the column formatting file, if any
  if(colformatfile.size()>0)
    {
      instream1.open(colformatfile);
      if(!instream1)  {
	cerr << "ERROR: unable to open input file " << colformatfile << "\n";
	return(1);
      }
      colreadct=0;
      while(!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct<COLS_TO_READ) {
	instream1 >> stest;
	if(stest == "MJDCOL") {
	  instream1 >> mjdcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "RACOL") {
	  instream1 >> racol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "DECCOL") {
	  instream1 >> deccol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "MAGCOL") {
	  instream1 >> magcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "OBSCODECOL") {
	  instream1 >> obscodecol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else {
	  cout << "WARNING: unrecognized string " << stest << " read from column formatting file\n";
	}
      }
      instream1.close();
      if(colreadct<COLS_TO_READ) {
	cout << "WARNING: only " << colreadct << " column specifications, of " << COLS_TO_READ << " expected, were read from column format file " << colformatfile << ".\n";
      }
    }

  cout << "Column specifications:\n";
  cout << "MJDCOL " << mjdcol << "\n";
  cout << "RACOL " << racol << "\n";
  cout << "DECCOL " << deccol << "\n";
  cout << "MAGCOL " << magcol << "\n";
  cout << "OBSCODECOL " << obscodecol << "\n";

  // Read input detection file.
  instream1.open(detfile);
  if(!instream1) {
    cerr << "can't open input file " << detfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  lct++;
  //cout << lnfromfile << "\n";
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    lct++;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    i=0;
    j = 0;
    c='0';
    mjdread = raread = decread = magread = obscoderead = 0;
    while(i<lnfromfile.size() && lnfromfile.size()>=30 && reachedeof == 0) {
      // Note check on line length: it is completely impossible for a
      // line containing all the required quantities at minimum plausible
      // precision to be less than 30 characters long.
      c='0';
      stest="";
      while(i<lnfromfile.size() && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==mjdcol) {
	MJD=stold(stest);
	mjdread=1;
      } else if(j==racol) {
	RA=stold(stest);
	raread=1;
      } else if(j==deccol) {
	Dec=stold(stest);
	decread=1;
      } else if(j==magcol) {
	mag=stod(stest);
	magread=1;
      } else if(j==obscodecol) {
	obscode=stest;
	obscoderead=1;
      }
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0 && lnfromfile.size()>=30) {
      if(!mjdread) {
	cerr << "ERROR: MJD not read from line " << detlines.size()+1 << " of input detection file " << detfile << "!\n";
	return(2);
      }
      if(!raread) {
	cerr << "ERROR: RA not read from line " << detlines.size()+1 << " of input detection file " << detfile << "!\n";
	return(2);
      }
      if(!decread) {
	cerr << "ERROR: Dec not read from line " << detlines.size()+1 << " of input detection file " << detfile << "!\n";
	return(2);
      }
      if(!magread) {
	mag = 99.999;
	cout << "WARNING: magnitude not read from line " << detlines.size()+1 << " of input detection file " << detfile << ".\n";
	cout << "magnitude will be set to 99.999\n";
      }
      if(!obscoderead) {
	cerr << "ERROR: observatory code not read from line " << detlines.size()+1 << " of input detection file " << detfile << "!\n";
	return(2);
      }
      if(MJD>mjdstart && MJD<mjdend) {
	detlines.push_back(lnfromfile);
	det_MJD.push_back(MJD);
	det_RA.push_back(RA);
	det_Dec.push_back(Dec);
	det_mag.push_back(mag);
	det_obscode.push_back(obscode);
	//cout << "Detection read: " << MJD << " " << Dec << " " << RA << " " << mag << " " << obscode << "\n";
      }
    }
  }
  instream1.close();
  if(reachedeof==1) { 
    cout << "Input file " << detfile << " read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  
  // Read file of mpc observations
  instream1.open(mpcfile,ios_base::in);
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the mpc file
    getline(instream1,lnfromfile);
    if(!instream1.bad() && !instream1.fail() && !instream1.eof() && lnfromfile.size()>=80) {
      // Parse out the interesting infromation.
      mpc80_parseline(lnfromfile, object, &MJD, &RA, &Dec, &mag, band, obscode);
      if(MJD>mjdstart && MJD<mjdend) {
	mpclines.push_back(lnfromfile);
	mpc_MJD.push_back(MJD);
	mpc_RA.push_back(RA);
	mpc_Dec.push_back(Dec);
	mpc_mag.push_back(mag);
	mpc_obscode.push_back(obscode);
	mpc_object.push_back(object);
      }
    }
  }

  // Load vectors for the KD trees
  // Detection
  for(i=0; i<detlines.size(); i++) {
    point0=point3d_index(det_MJD[i],det_RA[i],det_Dec[i],i);
    detvec.push_back(point0);
  }
  // mpc
  for(i=0; i<mpclines.size(); i++) {
    point0=point3d_index(mpc_MJD[i],mpc_RA[i],mpc_Dec[i],i);
    mpcvec.push_back(point0);
  }
  // Create KD-tree of detections
  dim=1;
  medpt = medind_3d_index(detvec,dim);
  kdroot = KD_point3d_index(detvec[medpt],-1,-1,1,0);
  det_kdtree={};
  det_kdtree.push_back(kdroot);
  kdtree_3d_index(detvec, dim, medpt, 0, det_kdtree);
  cout << "From input vector with " << detvec.size() << " points, made kdtree of size " << det_kdtree.size() << "\n";
  
  // Loop over the mpc vector, writing matched output
  outstream1.open(mpcout);
  outstream2.open(detout);
  for(i=0; i<mpclines.size(); i++) {
    indexvec = {};
    kdrange_3d_index(det_kdtree, mpcvec[i], qrange, indexvec);
    if(indexvec.size()>0) {
      cout << "For point " << i << " at " << mpcvec[i].x << " " << mpcvec[i].y << " " << mpcvec[i].z << " " << mpcvec[i].index << ", found " << indexvec.size() << " possible matches\n";
    }
    matchnum=0;
    for(j=0; j<indexvec.size(); j++) {
      // Map back to input detection vectors
      matchpoint = det_kdtree[indexvec[j]].point.index;
      // See if it's really a match
      cout << "Possible match: " << mpc_MJD[i] << " " << mpc_RA[i] << " " << mpc_Dec[i] << " " << mpc_mag[i] << " " << mpc_obscode[i] << "\n";
      cout << "And: " << det_MJD[matchpoint] << " " << det_RA[matchpoint] << " " << det_Dec[matchpoint] << " " << det_mag[matchpoint] << " " << det_obscode[matchpoint] << "\n";
      cout << "Comparison: " << det_obscode[matchpoint].compare(mpc_obscode[i]) << " " << SOLARDAY*fabs(det_MJD[matchpoint]-mpc_MJD[i]) << " < " << timetol*SOLARDAY << "?\n";
	
      if(det_obscode[matchpoint].compare(mpc_obscode[i])==0 && fabs(det_MJD[matchpoint]-mpc_MJD[i])<timetol) {
	// Yes, really a match
	matchnum++;
	// Write stuff to the mpc match file
	outstream1 << mpclines[i] << " " << det_MJD[matchpoint] << " " << det_RA[matchpoint] << " " << det_Dec[matchpoint] << " " << det_mag[matchpoint] << " " << det_obscode[matchpoint] << "\n";
	// Write stuff to the detection match
	outstream2 << detlines[matchpoint] << "," << mpc_object[i] << "," << mpc_MJD[i] << "," << mpc_RA[i] << "," << mpc_Dec[i] << "," << mpc_mag[i] << "," << mpc_obscode[i] << "," << mpclines[i] << "\n";
	// Wipe out the detection line
	detlines[matchpoint] = {};
      }
      if(matchnum>0) {
	// Wipe out the mpc line
	mpclines[i] = {};
      }
    }
  }
  outstream1.close();
  outstream2.close();
  
  // Write files of unmatched objects
  outstream1.open(mpc_nomatch);
  outstream2.open(det_nomatch);
  for(i=0; i<mpclines.size(); i++) {
    if(mpclines[i].size()>0) {
      outstream1 << mpclines[i] << " " << mpc_object[i] << " " << mpc_MJD[i] << " " << mpc_RA[i] << " " << mpc_Dec[i] << " " << mpc_mag[i] << " " << mpc_obscode[i] << "\n";
    }
  }
  for(i=0; i<detlines.size(); i++) {
    if(detlines[i].size()>0) {
      outstream2 << detlines[i] << "\n";
    }
  }
  outstream1.close();
  outstream2.close();

 return(0);
}

