// February 14, 2024: make_trailed_tracklets.cpp:
// Like make_tracklets_new.cpp, but rather than implicitly assuming
// that all detections are point sources, it uses trail length and
// orientation to construct only tracklets composed of trailed detections
// whose length and orientation are consistent with the motion implied
// by the tracklet itself. To do this, requires image exposure time
// and also two parameters specifying the matching tolerance for
// trail length and orientation, 
//
// Description of related program make_tracklets_new:
// Differs from the original C++ make_tracklets.cpp in that it uses
// the same core algorithmic routines as the new, python-wrapped codes.
// This requires abandoning some of the elegant efficiency of the
// original C++ codes' I/O and file structure, but it is necessary
// to enable consistent maintenance and developement of both python-wrapped
// and pure C++ into the future.


#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define NUMPOS 3

#define IDCOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define MAGCOL 5
#define BANDCOL 6
#define OBSCODECOL 7
#define TRAILPACOL 8
#define TRAILLENCOL 9
#define COLS_TO_READ 14
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME (1.5/24.0) // Default max inter-image time interval
                           // for tracklets, in days.
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)
#define MAX_GCR 0.5 // Default maximum Great Circle Residual allowed for a valid tracklet
#define DEBUG 0
#define DEBUGA 0
#define DEBUGB 0

      
static void show_usage()
{
  cerr << "Usage: make_trailed_tracklets -dets detfile -imgs imfile -outimgs output image file/ \n";
  cerr << "-pairdets paired detection file -tracklets tracklet file -trk2det output tracklet-to-detection file/ \n";
  cerr << "-colformat column format file -imrad image radius(deg)/ \n";
  cerr << "-maxtime max inter-image time interval (hr) -mintime min inter-image time interval (hr)/ \n";
  cerr << "-maxGCR maximum GRC -mintrkpts min. num. of tracklet points/ \n";
  cerr << "-minvel minimum angular velocity (deg/day) -maxvel maximum angular velocity (deg/day)/ \n";
  cerr << "-minarc minimum total angular arc (arcsec) -siglenscale fractional trail length matching tolerance/ \n";
  cerr << "-sigpascale trail orientation matching tolerance (arcsec; will be divided by trail length/ \n";
  cerr << "to get matching tolerance in radians) -exptime exposure time (seconds) -earth earthfile/ \n";
  cerr << "-obscode obscodefile -forcerun\n";
  cerr << "\nor, at minimum\n\n";
  cerr << "make_tracklets -dets detfile -earth earthfile -obscode obscodefile\n";
  cerr << "Note well that the minimum invocation will leave a bunch of things\n";
  cerr << "set to defaults that may not be what you want.\n";
}


int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
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
  string inimfile;
  string outimfile;
  string earthfile;
  string obscodefile;
  string colformatfile;
  string pairdetfile="pairdetfile01.csv";
  string trackletfile="trackletfile01.csv";
  string trk2detfile="trk2detfile01.csv";
  int idcol=IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int magcol = MAGCOL;
  int bandcol = BANDCOL;
  int obscodecol = OBSCODECOL;
  int trail_PA_col = TRAILPACOL;
  int trail_len_col = TRAILLENCOL;
  int sigmag_col, sig_across_col, sig_along_col, known_obj_col, det_qual_col;
  trail_len_col = trail_PA_col = sigmag_col = sig_across_col = -1;
  sig_along_col = known_obj_col = det_qual_col = -1;
  int colreadct=0;
  ifstream instream1;
  ofstream outstream1;
  string stest;
  int inimfile_set,outimfile_set,colformatfile_set;
  inimfile_set = outimfile_set = colformatfile_set = 0;
  int pairdetfile_default,trackletfile_default,trk2detfile_default,imagerad_default;
  int maxtime_default,mintime_default,minvel_default,maxvel_default;
  int maxgcr_default,minarc_default,mintrkpts_default;
  int exptime_default,siglenscale_default,sigpascale_default;
  MakeTrackletsConfig config;
  
  pairdetfile_default = trackletfile_default = trk2detfile_default = imagerad_default = 1;
  maxtime_default = mintime_default = minvel_default = maxvel_default = 1;
  maxgcr_default = minarc_default = mintrkpts_default = 1;
  exptime_default = siglenscale_default = sigpascale_default = 1;
  
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
    } else if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
	//There is still something to read;
	inimfile=argv[++i];
	inimfile_set = 1;
	i++;
      }
      else {
	cerr << "Image file keyword supplied with no corresponding argument\n";
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
    }  else if(string(argv[i]) == "-trk2det" || string(argv[i]) == "-t2d" || string(argv[i]) == "-tracklet2detection" || string(argv[i]) == "-trk2detfile" || string(argv[i]) == "--tracklet2detection") {
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
    } else if(string(argv[i]) == "-imrad" || string(argv[i]) == "-imagerad" ) {
      if(i+1 < argc) {
	//There is still something to read;
        config.imagerad=stod(argv[++i]);
	imagerad_default = 0;
	i++;
	if(!isnormal(config.imagerad) || config.imagerad<=0.0) {
	  cerr << "Error: invalid image radius (" << config.imagerad << " deg) supplied.\n";
	  cerr << "Image radius must be strictly positive!\n";
	  return(2);
	}
      }
      else {
	cerr << "Output image radius keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxtime") {
      if(i+1 < argc) {
	//There is still something to read;
        config.maxtime=stod(argv[++i]);
	maxtime_default = 0;
	i++;
	if(isnormal(config.maxtime) && config.maxtime>0.0) {
	  config.maxtime/=24.0; // Convert from hours to days.
	} else {
	  cerr << "Error: invalid maximum inter-image time interval\n";
	  cerr << "(" << config.maxtime << " hr) supplied: must be strictly positive.\n";
	  return(2);
	}      
      } else {
	cerr << "Maximum inter-image time interval\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-mintime") {
      if(i+1 < argc) {
	//There is still something to read;
        config.mintime=stod(argv[++i]);
	mintime_default = 0;
	i++;	
	if((isnormal(config.mintime) || config.mintime==0.0) && config.mintime>=0.0) {
	  config.mintime/=24.0; // Convert from hours to days
	} else {
	  cerr << "Error: invalid minimum inter-image time interval\n";
	  cerr << "(" << config.mintime << " hr) supplied: must be non-negative.\n";
	  return(2);
	}      
      } else {
	cerr << "Minimum inter-image time interval\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minvel") {
      if(i+1 < argc) {
	//There is still something to read;
        config.minvel=stod(argv[++i]);
	minvel_default = 0;
	i++;
	if(!isnormal(config.minvel) && config.minvel!=0.0l) {
	  cerr << "Error: invalid minimum angular velocity\n";
	  cerr << "(" << config.minvel << "deg/day) supplied.\n";
	  return(2);
	}
      }
      else {
	cerr << "Minimum angular velocity\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxvel") {
      if(i+1 < argc) {
	//There is still something to read;
        config.maxvel=stod(argv[++i]);
	maxvel_default = 0;
	i++;
	if(!isnormal(config.maxvel) || config.maxvel<=0.0) {
	  cerr << "Error: invalid maximum angular velocity\n";
	  cerr << "(" << config.maxvel << "deg/day) supplied: must be strictly positive.\n";
	  return(2);
	}
      }
      else {
	cerr << "Maximum angular velocity\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxGCR" || string(argv[i]) == "-maxgcr" ) {
      if(i+1 < argc) {
	//There is still something to read;
        config.maxgcr=stod(argv[++i]);
	maxgcr_default = 0;
	i++;
	if(!isnormal(config.maxgcr) || config.maxgcr<=0.0) {
	  cerr << "Error: invalid maximum Great Circle residual\n";
	  cerr << "(" << config.maxgcr << " arcsec) supplied: must be strictly positive.\n";
	  return(2);
	}
      }
      else {
	cerr << "Output maximum Great Circle Residual\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minarc") {
      if(i+1 < argc) {
	//There is still something to read;
        config.minarc=stod(argv[++i]);
	minarc_default = 0;
	i++;
	if(!isnormal(config.minarc) && config.minarc!=0.0l) {
	  cerr << "Error: invalid minimum angular arc\n";
	  cerr << "(" << config.minarc << " arcsec) supplied.\n";
	  return(2);
	}
      }
      else {
	cerr << "Minimum angular arc\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-siglenscale") {
      if(i+1 < argc) {
	//There is still something to read;
        config.siglenscale=stod(argv[++i]);
	siglenscale_default = 0;
	i++;
	if(!isnormal(config.siglenscale) || config.siglenscale<=0.0l) {
	  cerr << "Error: invalid fractional trail length uncertainty parameter\n";
	  cerr << "(" << config.siglenscale << ") supplied.\n";
	  return(2);
	}
      }
      else {
	cerr << "Fractional trail length uncertainty parameter\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-sigpascale") {
      if(i+1 < argc) {
	//There is still something to read;
        config.sigpascale=stod(argv[++i]);
	sigpascale_default = 0;
	i++;
	if(!isnormal(config.sigpascale) || config.sigpascale<=0.0l) {
	  cerr << "Error: invalid trail position angle uncertainty parameter\n";
	  cerr << "(" << config.sigpascale << " arcsec) supplied.\n";
	  return(2);
	}
      }
      else {
	cerr << "Trail position angle uncertainty parameter\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-exptime") {
      if(i+1 < argc) {
	//There is still something to read;
        config.exptime=stod(argv[++i]);
	exptime_default = 0;
	i++;
	if(!isnormal(config.exptime) || config.exptime<=0.0l) {
	  cerr << "Error: invalid exposure time\n";
	  cerr << "(" << config.exptime << " seconds) supplied.\n";
	  return(2);
	}
      }
      else {
	cerr << "Exposure time keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-mintrkpts" || string(argv[i]) == "-mpt" || string(argv[i]) == "-mintrackpts" || string(argv[i]) == "-minpts" || string(argv[i]) == "--minimumtrackletpoints" || string(argv[i]) == "--mintrackpoints" || string(argv[i]) == "--mintrackletpoints") {
      if(i+1 < argc) {
	//There is still something to read;
	config.mintrkpts=stoi(argv[++i]);
	mintrkpts_default = 0;
	i++;
      }
      else {
	cerr << "Min. tracklet points keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-obscode" || string(argv[i]) == "-obs" || string(argv[i]) == "-oc" || string(argv[i]) == "-obscodes" || string(argv[i]) == "--obscode" || string(argv[i]) == "--obscodes" || string(argv[i]) == "--observatorycodes") {
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
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format"  || string(argv[i]) == "-col" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" || string(argv[i]) == "--cformat") {
      if(i+1 < argc) {
	//There is still something to read;
	colformatfile=argv[++i];
	colformatfile_set = 1;
	i++;
      }
      else {
	cerr << "Column format file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	config.verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Verbosity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-forcerun" || string(argv[i]) == "-force"  || string(argv[i]) == "-fr" || string(argv[i]) == "-f" || string(argv[i]) == "--force" || string(argv[i]) == "--forcerun") {
      config.forcerun=1;
      i++;
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
    cerr << "Please supply a observatory code file:\n\n";
    show_usage();
    return(1);
  }
  
  if(config.mintrkpts<2) config.mintrkpts=2;
  
  cout << "\nInput detection file is called " << indetfile << "\n\n";
  
  if(inimfile_set == 1) cout << "Input image file = " << inimfile << "\n";
  else cout << "No input image file specified: image catalog will be generated internally.\n";	

  if(outimfile_set ==1) cout << "output image file = " << outimfile << "\n";
  else cout << "No output image file specified: required image information will only be used internally.\n";

  if(colformatfile_set == 1) cout << "column formatting file = " << colformatfile << "\n";
  else {
    cout << "No column formatting file specified for input file " << indetfile << "\n";
    cout << "The following default format will be assumed:\n";
    cout << "String identifer in column " << idcol << "\n";
    cout << "Modified Julian Day in column " << mjdcol << "\n";
    cout << "Right Ascension (RA) in column " << racol << "\n";
    cout << "Declination (Dec) in column " << deccol << "\n";
    cout << "Magnitude in column " << magcol << "\n";
    cout << "Photometric band in column " << bandcol << "\n";
    cout << "Observatory code in column " << obscodecol << "\n";
    cout << "Trail celestial position angle (PA, degrees) in column " << trail_PA_col << "\n";
    cout << "Trail length (arcsec) in column " << trail_len_col << "\n";
  }
  
  cout << "Observatory code file " << obscodefile << "\n";

  if(pairdetfile_default == 0) cout << "Output paired detection file will be called " << pairdetfile << "\n";
  else cout << "Defaulting to output paired detection file name = " << pairdetfile << "\n";

  if(trackletfile_default == 0) cout << "Output tracklet file will be called " << trackletfile << "\n";
  else cout << "Defaulting to output tracklet file name = " << trackletfile << "\n";
 
  if(trk2detfile_default == 0) cout << "Output tracklet-to-detection mapping file will be called " << trk2detfile << "\n";
  else cout << "Defaulting to output tracklet-to-detection file name = " << trk2detfile << "\n";
  
  cout << "Heliocentric ephemeris file for the Earth: " << earthfile << "\n";
  
  if(imagerad_default == 0) cout << "Image radius = " << config.imagerad << " degrees.\n";
  else cout << "Defaulting to image radius = " <<  config.imagerad << " degrees.\n";
  
  if(maxtime_default == 0) cout << "Max time interval = " << config.maxtime*24.0 << " hours.\n";
  else cout << "Defaulting to max time interval = " << config.maxtime*24.0 << " hours.\n";
  if(mintime_default ==0) cout << "Min time interval = " << config.mintime*24.0 << " hours.\n";
  else cout << "Defaulting to min time interval = " << config.mintime*24.0 << " hours.\n";

  if(minvel_default == 0) cout << "Min angular velocity = " << config.minvel << " deg/day.\n";
  else cout << "Defaulting to min angular velocity = " << config.minvel << " deg/day.\n";
  if(maxvel_default == 0) cout << "Maximum angular velocity = " << config.maxvel << " deg/day.\n";
  else cout << "Defaulting to maximum angular velocity = " << config.maxvel << " deg/day.\n";
  if(mintrkpts_default == 0) cout << "Minimum number of points per tracklet = " << config.mintrkpts << "\n";
  else cout << "Defaulting to minimum number of points per tracklet = " << config.mintrkpts << "\n";
  if(minarc_default == 0) cout << "Minimum tracklet length = " << config.minarc << " arcsec.\n";
  else cout << "Defaulting to minimum tracklet length = " << config.minarc << " arcsec.\n";
  if(maxgcr_default == 0) cout << "Maximum tracklet Great Circle residual = " << config.maxgcr << " arcsec.\n";
  else cout << "Defulting to maximum tracklet Great Circle residual = " << config.maxgcr << " arcsec.\n";
  if(siglenscale_default == 0) cout << "Fractional uncertainty on trail length = " << config.siglenscale << "\n";
  else cout << "Defulting to fractional uncertainty on trail length = " << config.siglenscale << "\n";
  if(sigpascale_default == 0) cout << "Trail position angle uncertainty parameter = " << config.sigpascale << " arcsec.\n";
  else cout << "Defulting to trail position angle uncertainty parameter = " << config.sigpascale << " arcsec.\n";
  if(exptime_default == 0) cout << "Exposure time = " << config.exptime << " seconds.\n";
  else cout << "Defulting to exposure time = " << config.exptime << " seconds.\n";
  
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
	} else if(stest == "TRAILLENCOL") {
	  instream1 >> trail_len_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "TRAILPACOL") {
	  instream1 >> trail_PA_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "SIGMAGCOL") {
	  instream1 >> sigmag_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "SIGACROSSCOL") {
	  instream1 >> sig_across_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "SIGALONGCOL") {
	  instream1 >> sig_along_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if (stest == "IDCOL") {
	  instream1 >> idcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "BANDCOL") {
	  instream1 >> bandcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "OBSCODECOL") {
	  instream1 >> obscodecol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "KNOWNOBJCOL") {
	  instream1 >> known_obj_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "DETQUALCOL") {
	  instream1 >> det_qual_col;
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
  cout << "IDCOL " << idcol << "\n";
  cout << "MJDCOL " << mjdcol << "\n";
  cout << "RACOL " << racol << "\n";
  cout << "DECCOL " << deccol << "\n";
  cout << "MAGCOL " << magcol << "\n";
  cout << "BANDCOL " << bandcol << "\n";
  cout << "OBSCODECOL " << obscodecol << "\n";
  cout << "TRAILLENCOL " << trail_len_col << "\n";
  cout << "TRAILPACOL " << trail_PA_col << "\n";

  // Read observatory code file
  status = read_obscode_file2(obscodefile, observatory_list, config.verbose);
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
  status = read_detection_filemt2(indetfile, mjdcol, racol, deccol, magcol, idcol, bandcol, obscodecol, trail_len_col, trail_PA_col, sigmag_col, sig_across_col, sig_along_col, known_obj_col, det_qual_col, detvec, config.verbose, config.forcerun);
  
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

  if(DEBUGB==1) cout << "Preparing to time-sort the detection vector\n";
  
  // time-sort the detection vector
  sort(detvec.begin(), detvec.end(), early_hldet());
  
  if(DEBUGB==1) cout << "Time-sorted the detection vector\n";

  // Get image information, if there is an image file
  if(inimfile.size()>0) {
    cout << "About to read image file " << inimfile << "\n";
    status = read_image_file(inimfile, img_log);
    cout << "Read image file with " << img_log.size() << " lines\n";
    if(status==0) {
      // time-sort the image file
      sort(img_log.begin(), img_log.end(), early_hlimage());
    } else {
      cerr << "Warning: failed to read supplied image file " << inimfile << "\n";
      cerr << "Constructing image table by inference from input detections instead\n";
      img_log={};
    }
  }
  EarthMJD={};
  Earthpos={};
  Earthvel={};
  read_horizons_csv(earthfile,EarthMJD,Earthpos,Earthvel);
  cout << "Finished reading heliocentric ephemeris file " << earthfile << " for the Earth.\n";

  if(DEBUGB==1) cout << "Preparing to load the image table\n";
  status = load_image_table(img_log, detvec, observatory_list, EarthMJD, Earthpos, Earthvel);
  if(DEBUGB==1) cout << "Loaded the image table\n";
  // Replace any invalid exposure times with the default (or user-supplied constant) value.
  long exp_resetnum=0;
  for(i=0;i<long(img_log.size());i++) {
    if(img_log[i].exptime<=0.0l) {
      cout << "Correcting exposure time on image " << i << ": " << img_log[i].MJD << " " << img_log[i].RA << " " << img_log[i].Dec << " " << img_log[i].obscode << ", exptime was " << img_log[i].exptime << "\n";
      img_log[i].exptime = config.exptime;
      exp_resetnum++;
    }
  }
  cout << "Exposure time was corrected for " << exp_resetnum << " out of " << img_log.size() << " images\n";

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
	  outstream1 << fixed << setprecision(8) << img_log[imct].MJD << " " << img_log[imct].RA;
	  outstream1 << fixed << setprecision(8) << " " << img_log[imct].Dec << " " << img_log[imct].obscode << " ";
	  outstream1 << fixed << setprecision(1) << img_log[imct].X << " " << img_log[imct].Y << " " << img_log[imct].Z << " ";
	  outstream1 << fixed << setprecision(4) << img_log[imct].VX << " " << img_log[imct].VY << " " << img_log[imct].VZ << " ";
	  outstream1 << img_log[imct].startind << " " << img_log[imct].endind << "\n";
	}
      outstream1.close();
    }

  make_trailed_tracklets(detvec, img_log, config, pairdets, tracklets, trk2det);

  // Write paired detection file
  cout << "Writing paired detection file with " << pairdets.size() << " lines\n";
  outstream1.open(pairdetfile);
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(i=0;i<long(pairdets.size());i++) {
    outstream1 << fixed << setprecision(7) << pairdets[i].MJD << "," << pairdets[i].RA << "," << pairdets[i].Dec << ",";
    outstream1 << fixed << setprecision(4) << pairdets[i].mag << ",";
    outstream1 << fixed << setprecision(2) << pairdets[i].trail_len << "," << pairdets[i].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << pairdets[i].sigmag << ",";
    outstream1 << fixed << setprecision(3) << pairdets[i].sig_across << "," << pairdets[i].sig_along << ",";
    outstream1 << pairdets[i].image << "," << pairdets[i].idstring << "," << pairdets[i].band << ",";
    outstream1 << pairdets[i].obscode << "," << pairdets[i].known_obj << ","; 
    outstream1 << pairdets[i].det_qual << "," << pairdets[i].index << "\n"; 
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

