// March 15, 2022: maketrack04c.cpp:
// Like maketrack04b.cpp, but writes output pairdets file in csv format
// with a header.
//
// Note that the format of the pair file (as opposed
// to the pairdets file) has NOT changed. The pair file has an unavoidably
// specialized format (because of the need to handle both pairs and
// tracklets), and this specialized format does not lend itself to the
// header+csv convention. Also, the pair file is not needed by any other
// program except projectpairs (AKA heliolinc), so its specialized format
// is not expected to create compatibility issues for programs external
// to the heliolinc suite.
//
// March 08, 2022: maketrack04b.cpp:
// Like maketrack04a.cpp, but uses MPC observatory codes to figure
// out the latitude and longitudes of observing sites, and allows
// input data to come from multiple observatories. Besides the
// observatory code, also reads in magnitude and photometric band
// (i.e., filter) information. Accepts an input column formatting file
// specifying in which columns of the input observation files these
// various quantities are to be found.
// Example column formatting file:
//
// IDCOL 1
// MJDCOL 3
// RACOL 6
// DECCOL 8
// MAGCOL 32
// BANDCOL 26
// OBSCODECOL 38
// 
// January 18, 2022: maketrack04a.cpp:
// Description of ancestor program maketrack03a.cpp:
// Like maketrack02b.cpp, but saves more information in the paired
// detection file. Specifically, it saves the observers heliocentric
// (or barycentric) X, Y, Z, position, and a string identifier for each
// detection. Whether the observer positions in the output paired
// detection file are heliocentric or barycentric depends on the input
// Earth position file. It is expected that in general, heliocentric
// observer positions will be more useful because we will be using a
// simple, 2-body Keplerian model for orbit propagation.
//
// November 15, 2021: maketrack02a.cpp
// Like maketrack01b.cpp, but uses a k-d tree to speed up pairing.
// 
// November 10, 2021: maketrack01b.cpp
// Like maketrack01a.cpp, but projects celestial positions onto an
// x-y plane centered on the primary image prior to matching.
// This means that the distance calculation can use fast 2-D Cartesian
// distance rather than spherical geometry, which speeds up the
// program by a factor of a few.
// 
// Description of ancestor program maketrack01a.cpp:
// Read a csv file of detections of astronomical objects, test
// crude n^2 pairing based on image matching.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define NUMPOS 3

#define IDCOL 1
#define MJDCOL 3
#define RACOL 6
#define DECCOL 8
#define MAGCOL 32
#define BANDCOL 26
#define OBSCODECOL 38
#define COLS_TO_READ 7
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME (1.5/24.0) // Default max inter-image time interval
                           // for tracklets, in days.
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)
#define MAX_GCR 0.5 // Default maximum Great Circle Residual allowed for a valid tracklet
#define DEBUG 0

      
static void show_usage()
{
  cerr << "Usage: maketrack04b -dets detfile -imgs imfile -outimgs output image file/ \n";
  cerr << "-pairs pairfile -pairdets paired detection file -colformat column format file/ \n";
  cerr << "-imrad image radius(deg) -maxtime max inter-image time interval (hr)/ \n";
  cerr << "-mintime min inter-image time interval (hr) -maxGCR maximum GRC -mintrkpts min. num. of tracklet points/\n";
  cerr << "-minvel minimum angular velocity (deg/day) -maxvel maximum angular velocity (deg/day) \n";
  cerr << "-minarc minimum total angular arc (arcsec) -earth earthfile -obscode obscodefile\n";
  cerr << "\nor, at minimum\n\n";
  cerr << "maketrack04b -dets detfile -earth earthfile -obscode obscodefile\n";
  cerr << "Note well that the minimum invocation will leave a bunch of things\n";
  cerr << "set to defaults that may not be what you want.\n";
}
    
int main(int argc, char *argv[])
{
  det_obsmag_indvec o1 = det_obsmag_indvec(0L,0l,0l,0L,0L,0L,"null",0l,"V","I11",0,{});
  vector <det_obsmag_indvec> detvec = {};
  vector <det_obsmag_indvec> pairdets = {};
  vector <det_obsmag_indvec> ppset = {};
  observatory obs1 = observatory("Ill",0l,0l,0l);
  vector <observatory> observatory_list = {};
  img_log03 imlog = img_log03(0.0,0.0,0.0,"I11",0,0);
  vector <img_log03> img_log_tmp = {};
  vector <img_log03> img_log = {};
  longpair onepair = longpair(0,0);
  vector <longpair> pairvec ={};
  point3d p3 = point3d(0,0,0);
  point3d p3avg = point3d(0,0,0);
  vector <point3LD> Earthpos;
  vector <point3LD> Earthvel;
  vector <point3LD> observer_heliopos;
  vector <long double> EarthMJD;
  point3LD outpos = point3LD(0,0,0);
  double tdelt = 0;
  double mjdmean = 0;
  double mjdnorm = 0;
  char idstring[SHORTSTRINGLEN];
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  string lnfromfile;
  int status = 0;
  int i = 0;
  int j = 0;
  int k = 0;
  int imct=0;
  int imctp=0;
  int imnum=0;
  long detnum=0;
  long num_dets=0;
  long detct=0;
  int startind=0;
  int endind=0;
  int reachedeof = 0;
  char c='0';
  long double MJD,RA,Dec;
  MJD = RA = Dec = 0.0L;
  double mag = 0l;
  double maxvel = MAXVEL; // Max angular velocity in deg/day
  double minvel = 0.0l; // Min angular velocity in deg/day
  double minarc = 0.0l; // Min total angular arc in arcseconds
  double angvel = 0.0l;
  double maxtime = MAXTIME; // Max time interval a tracklet could span,
                            // in days.
  double maxdist = MAXVEL*MAXTIME; // Max angular distance a tracklet
                                   // could span, in degrees.
  double imrad = IMAGERAD; // radius from image center to most distant corner (deg).
  string indetfile;
  string inimfile;
  string outimfile;
  string earthfile;
  string obscodefile;
  string colformatfile;
  string outpairfile="outpairfile01.txt";
  string pairdetfile="pairdetfile01.txt";
  double obslon = 289.26345L;
  double plxcos = 0.865020L;
  double plxsin = -0.500901L;
  long lct=0;
  long_index ppn = long_index(0,0);
  vector <long_index> pair_partner_num={};
  vector <long_index> tracklet_check={};
  double dt,dtref,dx,dy;
  dt = dtref = dx = dy = 0.0l;
  xy_index xyind=xy_index(0.0, 0.0, 0);
  vector <xy_index> axyvec = {};
  double dist,pa;
  dist=pa=0.0;
  int dettarg=0;
  vector <double> timevec;
  vector <double> xvec;
  vector <double> yvec;
  vector <long> detindexvec;
  int biggest_tracklet=-1;
  int tracklet_size=0;
  double slopex,slopey,interceptx,intercepty,worsterr;
  vector <double> fiterr = {};
  vector <double> fiterr2 = {};
  int worstpoint=-1;
  int istracklet=0;
  int rp1,rp2,instep;
  rp1=rp2=instep=0;
  double outra1,outra2,outdec1,outdec2;
  outra1 = outra2 = outdec1 = outdec2 = 0.0l;
  point3d_index p3di =  point3d_index(0l,0l,0l,0);
  vector <point3d_index> track_mrdi_vec;
  double mintime = IMAGETIMETOL/SOLARDAY;
  int trkptnum,istimedup=1;
  double maxgcr = MAX_GCR;
  int idcol=IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int magcol = MAGCOL;
  int bandcol = BANDCOL;
  int obscodecol = OBSCODECOL;
  int colreadct=0;
  ifstream instream1;
  ofstream outstream1;
  string stest;
  int mintrkpts=2;
  
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
	i++;
      }
      else {
	cerr << "Output image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-pairs" || string(argv[i]) == "-pairfile") {
      if(i+1 < argc) {
	//There is still something to read;
	outpairfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output pair file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" || string(argv[i]) == "-detpairs") {
      if(i+1 < argc) {
	//There is still something to read;
        pairdetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output paired detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-imrad") {
      if(i+1 < argc) {
	//There is still something to read;
        imrad=stod(argv[++i]);
	i++;
	if(!isnormal(imrad) || imrad<=0.0) {
	  cerr << "Error: invalid image radius (" << imrad << " deg) supplied.\n";
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
        maxtime=stod(argv[++i]);
	i++;
	if(isnormal(maxtime) && maxtime>0.0) {
	  maxtime/=24.0; // Convert from hours to days.
	} else {
	  cerr << "Error: invalid maximum inter-image time interval\n";
	  cerr << "(" << maxtime << " hr) supplied: must be strictly positive.\n";
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
        mintime=stod(argv[++i]);
	i++;	
	if((isnormal(mintime) || mintime==0.0) && mintime>=0.0) {
	  mintime/=24.0; // Convert from hours to days
	} else {
	  cerr << "Error: invalid minimum inter-image time interval\n";
	  cerr << "(" << mintime << " hr) supplied: must be non-negative.\n";
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
        minvel=stod(argv[++i]);
	i++;
	if(!isnormal(minvel) && minvel!=0.0l) {
	  cerr << "Error: invalid minimum angular velocity\n";
	  cerr << "(" << minvel << "deg/day) supplied.\n";
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
        maxvel=stod(argv[++i]);
	i++;
	if(!isnormal(maxvel) || maxvel<=0.0) {
	  cerr << "Error: invalid maximum angular velocity\n";
	  cerr << "(" << maxvel << "deg/day) supplied: must be strictly positive.\n";
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
        maxgcr=stod(argv[++i]);
	i++;
	if(!isnormal(maxgcr) || maxgcr<=0.0) {
	  cerr << "Error: invalid maximum Great Circle residual\n";
	  cerr << "(" << maxgcr << " arcsec) supplied: must be strictly positive.\n";
	  return(2);
	}
      }
      else {
	cerr << "Output maximum Great Circle Residual\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-minarc") {
      if(i+1 < argc) {
	//There is still something to read;
        minarc=stod(argv[++i]);
	i++;
	if(!isnormal(minarc) && minarc!=0.0l) {
	  cerr << "Error: invalid minimum angular arc\n";
	  cerr << "(" << minarc << " arcsec) supplied.\n";
	  return(2);
	}
      }
      else {
	cerr << "Minimum angular arc\nkeyword supplied with no corresponding argument\n";
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
	mintrkpts=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Earth file keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" || string(argv[i]) == "--cformat") {
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
  
  if(mintrkpts<2) mintrkpts=2;
  
  cout << "indet file " << indetfile << "\n";
  cout << "inimage file " << inimfile << "\n";
  cout << "column formatting file " << colformatfile << "\n";
  cout << "observatory code file " << obscodefile << "\n";
  cout << "output image file " << outimfile << "\n";
  cout << "pairfile file " << outpairfile << "\n";
  cout << "paired detection file " << pairdetfile << "\n";
  cout << "Heliocentric ephemeris file for the Earth: " << earthfile << "\n";
  cout << "image radius " << imrad << " degrees\n";
  cout << "max time interval " << maxtime*24.0 << " hours\n";
  cout << "min time interval " << mintime*24.0 << " hours\n";
  cout << "minvel " << minvel << " deg/day\n";
  cout << "maxvel " << maxvel << " deg/day\n";
  cout << "minimum number of points per tracklet " << mintrkpts << "\n";
  cout << "minarc " << minarc << " arcsec\n";
  cout << "maxGCR " << maxgcr << " arcsec\n";
  maxdist = maxtime*maxvel;
  
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
	if (stest == "IDCOL") {
	  instream1 >> idcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "MJDCOL") {
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
	} else if(stest == "BANDCOL") {
	  instream1 >> bandcol;
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
  cout << "IDCOL " << idcol << "\n";
  cout << "MJDCOL " << mjdcol << "\n";
  cout << "RACOL " << racol << "\n";
  cout << "DECCOL " << deccol << "\n";
  cout << "MAGCOL " << magcol << "\n";
  cout << "BANDCOL " << bandcol << "\n";
  cout << "OBSCODECOL " << obscodecol << "\n";

  // Read observatory code file
  instream1.open(obscodefile);
  if(!instream1) {
    cerr << "can't open input file " << obscodefile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  while (!instream1.eof() && !instream1.fail() && !instream1.bad())
    {
      instream1 >> stest;
      stringncopy01(obscode,stest,MINSTRINGLEN);
      instream1 >> obslon;
      instream1 >> plxcos;
      instream1 >> plxsin;
      obs1 = observatory(obscode,obslon,plxcos,plxsin);
      observatory_list.push_back(obs1);
      // Skip the rest of the line
      getline(instream1,lnfromfile);
    }
  instream1.close();
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile << ":\n";
  
  if(DEBUG>=2) {
    for(i=0;i<observatory_list.size();i++) {
      cout << observatory_list[i].obscode << " " << observatory_list[i].obslon << " " << observatory_list[i].plxcos << " " << observatory_list[i].plxsin << "\n";
    }
  }
  
  // Read input detection file.
  instream1.open(indetfile);
  if(!instream1) {
    cerr << "can't open input file " << indetfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  lct++;
  cout << lnfromfile << "\n";
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
      if(j==idcol) stringncopy01(idstring,stest,SHORTSTRINGLEN);
      if(j==mjdcol) MJD=stold(stest);
      else if(j==racol) RA=stold(stest);
      else if(j==deccol) Dec=stold(stest);
      else if(j==magcol) mag=stod(stest);
      else if(j==bandcol) stringncopy01(band,stest,MINSTRINGLEN);
      else if(j==obscodecol) stringncopy01(obscode,stest,MINSTRINGLEN);
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0 && lnfromfile.size()>=30) {
      // cout<<"MJD, RA, Dec: " << MJD-floor(MJD) << " " << RA << " " << Dec << "\n";
      o1=det_obsmag_indvec(MJD,RA,Dec,0L,0L,0L,idstring,mag,band,obscode,-lct,{});
      detvec.push_back(o1);
    }
  }
  instream1.close();
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";

  cout << "Last two obscodes: " << detvec[detvec.size()-2].obscode << " and " << detvec[detvec.size()-1].obscode << "\n"; 
  
  // time-sort the detection vector
  sort(detvec.begin(), detvec.end(), early_det_obsmag_indvec());
  
  cout << "Last two obscodes: " << detvec[detvec.size()-2].obscode << " and " << detvec[detvec.size()-1].obscode << "\n"; 

  //Print out sorted detvec for debugging.
  //for(i=0;i<detvec.size();i++) {
  //   cout << "det# " << i << ": " << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << " " << detvec[i].idstring << " " << detvec[i].mag << " " << detvec[i].band << " " << detvec[i].obscode << " " << detvec[i].index << "\n";
  // }
  
  // Get image information.
  if(inimfile.size()>0) {
    // Read input image file: MJD, RA, Dec, obscode:
    instream1.open(inimfile);
    if(!instream1) {
      cerr << "can't open input file " << inimfile << "\n";
      return(1);
    }
    reachedeof=0;
    while(reachedeof==0) {
      getline(instream1,lnfromfile);
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
      else if(instream1.eof()) reachedeof=1; //End of file, fine.
      else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
      else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
      i=0;
      j = 0;
      c='0';
      MJD=0.0;
      while(i<lnfromfile.size() && reachedeof == 0) {
	stest="";
	c='0';
	while(i<lnfromfile.size() && c!=',' && c!=' ' && c!='\n' && c!=EOF) {
	  // We allow the file to be delimited by comma or space.
	  c=lnfromfile[i];
	  if(c!=',' && c!=' ' && c!='\n' && c!=EOF) stest.push_back(c);
	  i++;
	}
	// We just finished reading something
	j++;
	if(j==1) MJD=stod(stest); // We assume we have MJD, RA, Dec, obscode
	else if(j==2) RA=stod(stest);
	else if(j==3) Dec=stod(stest);
	else if(j==4) stringncopy01(obscode,stest,MINSTRINGLEN);
      }
      if((reachedeof == 0 || reachedeof == 1) && MJD>0.0) {
	// Requirement of MJD>0.0 tests that we read a plausibly
	// valid line.
	imlog=img_log03(MJD,RA,Dec,obscode,0,0);
	img_log_tmp.push_back(imlog);
      }
    }
    if(reachedeof==1) { 
      cout << "File read successfully to the end.\n";
    }
    else if(reachedeof==-1) cout << "Warning: file read failed\n";
    else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
    else cout << "Warning: unknown file read problem\n";
    // time-sort the image file
    sort(img_log_tmp.begin(), img_log_tmp.end(), early_imlg3());
    // find the indices in the time-sorted detection file
    // that correspond to the earliest and latest detections
    // on each image, and load these values into imglog02.
    detct=0;
    for(imct=0;imct<img_log_tmp.size();imct++) {
      while(detct<detvec.size() && detvec[detct].MJD < img_log_tmp[imct].MJD-IMAGETIMETOL/SOLARDAY) detct++; //Not on any image
      if(detct<detvec.size() && fabs(detvec[detct].MJD-img_log_tmp[imct].MJD)<=IMAGETIMETOL/SOLARDAY) {
	// This should be the first detection on image imct.
	img_log_tmp[imct].startind = detct;
	while(detct<detvec.size() && fabs(detvec[detct].MJD-img_log_tmp[imct].MJD)<=IMAGETIMETOL/SOLARDAY) detct++; //Still on this same image
	// This should be the first detection on the next image
	img_log_tmp[imct].endind = detct;
      }
      if(img_log_tmp[imct].startind >= 0 && img_log_tmp[imct].endind > 0) {
	img_log.push_back(img_log_tmp[imct]);
      }
    }
    instream1.close();
  } else {
    // No input image file was supplied: we have to create one from
    // the sorted detection file.
    mjdnorm = 1.0;
    mjdmean = detvec[0].MJD;
    startind=0;
    for(i=1;i<detvec.size();i++) {
      tdelt = detvec[i].MJD - detvec[i-1].MJD;
      if(tdelt < IMAGETIMETOL/SOLARDAY && stringnmatch01(detvec[i].obscode,detvec[i-1].obscode,3)==0) {
	//This point corresponds to the same image as the previous one.
	mjdmean += detvec[i].MJD;
	mjdnorm += 1.0;
      }
      else {
	//Now we are considering a new image.
	//Calculate the meanmjd of the previous image, for which
	// we have now seen all points.
	//Record the current detct i as the detection index just
	//after the end of the previous image
	endind=i;
	if(isnormal(mjdnorm)) mjdmean /= mjdnorm;
	else mjdmean = 0.0;
	//Load it into the vector with mean MJD for all images,
	// and increment image count.
	imlog = img_log03(mjdmean,0.0,0.0,detvec[endind-1].obscode,startind,endind);
	img_log.push_back(imlog);
	// Set up for the next image, starting with detvec[i].MJD;
	mjdmean = detvec[i].MJD;
	mjdnorm = 1.0;
	startind=i;
      }
    }
    cout << "OK here, now for final image\n";
    fflush(stdout);
    // Account for the final image.
    if(isnormal(mjdnorm)) {
      endind=i;
      mjdmean /= mjdnorm;
      //Load it into the vector with mean MJD for all images,
      // and increment image count.
      imlog = img_log03(mjdmean,0.0,0.0,detvec[endind-1].obscode,startind,endind);
      img_log.push_back(imlog);
    }
    cout << "Done with final image\n";
    fflush(stdout);

    //We've now loaded the mean MJDs and the starting and ending
    //detection table indices for each image; it still remains to
    //get the mean RA and Dec.
   
    detnum = detvec.size();
    imnum = img_log.size();
    cout << img_log.size() << " unique images were identified.\n";
    cout << "Given our total of " << detvec.size() << " detections,\n";
    cout << "we have " << double(detvec.size())/double(img_log.size()) << " detections per image, on average\n";

    // Find the number of detections and the average RA, Dec on each image.
    // We perform the average after projection onto the unit circle, to
    // avoid wrapping issues.
    detct=imct=0;
    while( imct<imnum && detct<detnum ) {
      vector <det_obsmag_indvec> imobs = {};
      num_dets=0;
      p3avg = point3d(0,0,0);
      while(detct<detnum && detvec[detct].MJD < img_log[imct].MJD + IMAGETIMETOL/SOLARDAY) {
	num_dets++; //Keep count of detections on this image
	imobs.push_back(detvec[detct]); //Load vector of observations for this image
	p3 =  celeproj01(detvec[detct].RA,detvec[detct].Dec); // Project current detection
	p3avg.x += p3.x;
	p3avg.y += p3.y;
	p3avg.z += p3.z; // Average projected coords
	detct++;
      }
      // If we got here, we must just have finished with an image.
      // Finish the average
      if(num_dets>0) {
	p3avg.x /= double(num_dets);
	p3avg.y /= double(num_dets);
	p3avg.z /= double(num_dets);
	i=celedeproj01(p3avg, &img_log[imct].RA, &img_log[imct].Dec);
	if(i==0) ; // All is well.
	else if(i==1) {
	  cout << "Warning: vector of zeros fed to celedeproj01\n";
	  img_log[imct].RA = img_log[imct].Dec = 0.0;
	}
	else if(i==2) {
	  cout << "Warning: impossible z value " << p3avg.z << " fed to celedeproj01\n";
	  img_log[imct].RA = img_log[imct].Dec = 0.0;
	}
	else {
	  cout << "Warning: unspecified failure from celedeproj01 with\n";
	  cout << "input " << p3avg.x << " " << p3avg.y << " " << p3avg.z << "\n";
	  img_log[imct].RA = img_log[imct].Dec = 0.0;
	}
      }
      if(DEBUG>=1) {
	cout << "Image " << imct << " of " << img_log.size() << ": " << num_dets << " = " << img_log[imct].endind-img_log[imct].startind ;
	cout << " from " << img_log[imct].startind << " to " << img_log[imct].endind << " of " << detvec.size() << ".\n";
	fflush(stdout);
      }	
      imct++;
    }
  }

  EarthMJD={};
  Earthpos={};
  Earthvel={};
  read_horizons_fileLD(earthfile,EarthMJD,Earthpos,Earthvel);
  cout << "Finished reading heliocentric ephemeris file " << earthfile << " for the Earth.\n";

  if(DEBUG>=2) {
    // Test: print out time-sorted detection table.
    outstream1.open("testjunk01.txt");
    for(i=0;i<detvec.size();i++) {
      outstream1 << fixed << setprecision(6) << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec <<  " " << detvec[i].mag << " " << detvec[i].band <<  " " << detvec[i].obscode << "\n";
    }
    outstream1.close();
  }
  
  if(outimfile.size()>0)
    {
      // Write and print image log table
      outstream1.open(outimfile);
      for(imct=0;imct<img_log.size();imct++)
	{
	  //	  cout << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
	  //	  cout << " " << img_log[imct].Dec << " " << img_log[imct].startind << " " << img_log[imct].endind << " ";
	  //	  cout << img_log[imct].endind-img_log[imct].startind << "\n";
	  outstream1 << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
	  outstream1 << fixed << setprecision(6) << " " << img_log[imct].Dec << " " << img_log[imct].obscode << " ";
	  outstream1 << img_log[imct].startind << " " << img_log[imct].endind << "\n";
	}
      outstream1.close();
    }
  // Calculate observer's heliocentric position at the time of each image.
  observer_heliopos={};
  for(imct=0;imct<img_log.size();imct++)
    {
      if(imct==0 || (imct>0 && stringnmatch01(img_log[imct].obscode,img_log[imct-1].obscode,3)==0)) {
	// Observatory has changed: get observatory coordinates for this image.
	status = obscode_lookup(observatory_list,img_log[imct].obscode,obslon,plxcos,plxsin);
	if(status>0) {
	  cerr << "ERROR: obscode_lookup failed for observatory code " << img_log[imct].obscode << "\n";
	  return(3);
	}
      }
      observer_barycoords01LD(img_log[imct].MJD, 5, obslon, plxcos, plxsin, EarthMJD, Earthpos, outpos);
      observer_heliopos.push_back(outpos);
    }
  // assign observer heliocentric position to each detection.
  for(imct=0;imct<img_log.size();imct++) {
    for(i=img_log[imct].startind;i<img_log[imct].endind;i++) {
      detvec[i].x = observer_heliopos[imct].x;
      detvec[i].y = observer_heliopos[imct].y;
      detvec[i].z = observer_heliopos[imct].z;
      // Check that the obscodes match between image and detection
      if(stringnmatch01(detvec[i].obscode,img_log[imct].obscode,3)!=0) {
	cout << "ERROR: obscode mismatch (" << detvec[i].obscode << " vs. " << img_log[imct].obscode << " between image " << imct << " and detection " << i << "\n";
	return(4);
      }
      detvec[i].indvec = {}; // Just making sure the index vectors are empty at this point.
    }
  }

  // PERFORM PAIRING
  long pdct=0; // count of detections that have been paired
  long pairct=0; // count of actual pairs
  // Loop over images for image A
  for(imct=0;imct<img_log.size();imct++) {
    int apct=0;
    int adetct=0;
    // See if there are any images that might match
    vector <int> imagematches = {};
    int imatchcount = 0;
    int imtarg=imct+1;
    while(imtarg<img_log.size() && img_log[imtarg].MJD < img_log[imct].MJD + maxtime) {
      // See if the images are close enough on the sky.
      double timediff = img_log[imtarg].MJD-img_log[imct].MJD;
      if(!isnormal(timediff) || timediff<0.0) {
	cerr << "WARNING: Negative time difference " << timediff << " encountered between images " << imct << " and " << imtarg << "\n";
      }
      double imcendist = distradec01(img_log[imct].RA, img_log[imct].Dec, img_log[imtarg].RA, img_log[imtarg].Dec);
      if(imcendist<2.0*imrad+maxvel*timediff && timediff>=mintime) {
	if(DEBUG>=1) cout << "  pairs may exist between images " << imct << " and " << imtarg << ": dist = " << imcendist << ", timediff = " << timediff << "\n";
	imagematches.push_back(imtarg);
      }
      imtarg++;
    }
    if(DEBUG>=1) cout << "Looking for pairs for image " << imct << ": " << imagematches.size() << " later images are worth searching\n";
    if(imagematches.size()>0) {
      // Search is worth doing. Project all the detections
      // on image A.
      xyind=xy_index(0.0, 0.0, 0);
      axyvec = {};
      dist=pa=0.0;
      dettarg=0;
      for(detct=img_log[imct].startind ; detct<img_log[imct].endind ; detct++) {
	distradec02(img_log[imct].RA, img_log[imct].Dec,detvec[detct].RA,detvec[detct].Dec,&dist,&pa);
	xyind = xy_index(dist*sin(pa/DEGPRAD),dist*cos(pa/DEGPRAD),detct);
	axyvec.push_back(xyind);
	if((!isnormal(xyind.x) && xyind.x!=0) || (!isnormal(xyind.y) && xyind.y!=0)) {
	  cerr << "nan-producing input: ra1, dec1, ra2, dec2, dist, pa:\n";
	  cerr << img_log[imct].RA << " " << img_log[imct].Dec << " " << detvec[detct].RA << " " << detvec[detct].Dec << " " << dist << " " << pa << " " << xyind.x << " " << xyind.x << "\n";
	  //return(6);
	}
      }
      // Loop over images with potential matches (image B's)
      for(imatchcount=0;imatchcount<imagematches.size();imatchcount++)
      {
	imtarg = imagematches[imatchcount];
	double range = (img_log[imtarg].MJD-img_log[imct].MJD)*maxvel;
	vector <xy_index> bxyvec = {};
	// Project all detections on image B
	for(dettarg=img_log[imtarg].startind ; dettarg<img_log[imtarg].endind ; dettarg++) {
	  distradec02(img_log[imct].RA, img_log[imct].Dec,detvec[dettarg].RA,detvec[dettarg].Dec,&dist,&pa);
	  xyind = xy_index(dist*sin(pa/DEGPRAD),dist*cos(pa/DEGPRAD),dettarg);
	  bxyvec.push_back(xyind);
	}
	// Create k-d tree of detections on image B (imtarg).
	int dim=1;
	xy_index xyi = bxyvec[0];
	kdpoint root = kdpoint(xyi,-1,-1,dim);
	kdpoint lp1 = kdpoint(xyi,-1,-1,dim);
	kdpoint rp1 = kdpoint(xyi,-1,-1,dim);
	kdpoint kdtest = kdpoint(xyi,-1,-1,dim);
	vector <kdpoint> kdvec ={};
	long medpt;
	medpt = medindex(bxyvec,dim);
	root = kdpoint(bxyvec[medpt],-1,-1,1);
	kdvec.push_back(root);
	kdtest=kdvec[0];
	kdtree01(bxyvec,dim,medpt,0,kdvec);
	// Loop over detections on image A
	if(DEBUG>=1) cout << "Looking for pairs between " << axyvec.size() << " detections on image " << imct << " and " << kdvec.size() << " on image " << imtarg << "\n";
	for(detct=0 ; detct<axyvec.size() ; detct++) {
	  vector <long> indexvec = {};
	  if((isnormal(axyvec[detct].x) || axyvec[detct].x==0) && (isnormal(axyvec[detct].y) || axyvec[detct].y==0)) {
	     kdrange01(kdvec,axyvec[detct].x,axyvec[detct].y,range,indexvec);
	  } else {
	    cerr << "WARNING: detection " << detct << " on image " << imct << " not normal: " << axyvec[detct].x << " " << axyvec[detct].y << "\n";
	  }
	  int matchnum=indexvec.size();
	  long matchpt=0;
	  int matchct=0;
	  if(matchnum>0) {
	    // Record image A detection as paired, if not already recorded.
	    if(detvec[axyvec[detct].index].index<0) {
	      //This detection has not yet been paired with any other.
	      detvec[axyvec[detct].index].index *= -1; // Mark as paired by changing to positive sign.
	      pairdets.push_back(detvec[axyvec[detct].index]); // Load into paired detection vector
	      pairdets[pdct].indvec={}; // Make sure index vector is currently empty.
	      detvec[axyvec[detct].index].index = pdct; // Re-assign index to apply to paired detection vector
	      pdct++; // Increment count of paired detections
	      adetct++;
	      if(pdct!=pairdets.size()) cerr << "\nWARNING: PAIRED DETECTION MISMATCH: " << pdct << " vs " << pairdets.size() << "\n\n";
	    }
	    // Record image B detections
	    for(matchct=0;matchct<matchnum;matchct++) {
	      matchpt = indexvec[matchct];
	      if(detvec[kdvec[matchpt].point.index].index<0) {
		//This detection has not yet been paired with any other.
		detvec[kdvec[matchpt].point.index].index *= -1; // Mark as paired by changing to positive sign
		pairdets.push_back(detvec[kdvec[matchpt].point.index]); // Load into paired detection vector
		pairdets[pdct].indvec={}; // Make sure index vector is currently empty.
		detvec[kdvec[matchpt].point.index].index = pdct; // Re-assign index to apply to paired detection vector
		pdct++; // Increment count of paired detections
		if(pdct!=pairdets.size()) cerr << "\nWARNING: PAIRED DETECTION MISMATCH: " << pdct << " vs " << pairdets.size() << "\n\n";
	      }
	      // Write index values for both components of the
	      // new pair to the pair vector, regardless of whether
	      // the index values are pre-existing or newly assigned.
	      onepair = longpair(detvec[axyvec[detct].index].index,detvec[kdvec[matchpt].point.index].index);
	      if(DEBUG>=1) cout << "Writing pair " << detvec[axyvec[detct].index].index << ", " << detvec[kdvec[matchpt].point.index].index << "\n";
	      pairvec.push_back(onepair);
	      pairct++;
	      apct++;
	      // Load index of each detection into the paired index vector of the other
	      pairdets[detvec[axyvec[detct].index].index].indvec.push_back(detvec[kdvec[matchpt].point.index].index);
	      pairdets[detvec[kdvec[matchpt].point.index].index].indvec.push_back(detvec[axyvec[detct].index].index);
	    }
	    // Close if-statement checking if image A detection was matched to anything.
	  }
	  // Close loop over detections on source image (image A)
	}
	// Close loop over image B candidates
      }
      // Close if-statement checking if any images could match image A      
    }
    if(DEBUG>=1) cout << "Image " << imct << ": found " << adetct << " newly paired detections and a total of " << apct << " pairs.\n";
    // Close loop over images for image A
  }
  if(DEBUG>=1) cout << "Test count of paired detections: " << pdct << " " << pairdets.size() << "\n";
  if(DEBUG>=1) cout << "Test count of pairs: " << pairct << " " << pairvec.size() << "\n";
  
  // Write paired detections vector to file
  cout << "Writing paired detections file\n";
  outstream1.open(pairdetfile);
  outstream1 << "#MJD,RA,Dec,observerX,observerY,observerZ,stringID,mag,band,obscode,origindex\n";
  for(i=0;i<pairdets.size();i++) {
    outstream1 << fixed << setprecision(7) << pairdets[i].MJD << "," << pairdets[i].RA << "," << pairdets[i].Dec << ",";
    outstream1 << fixed << setprecision(3) << pairdets[i].x << "," << pairdets[i].y << "," << pairdets[i].z << ",";
    outstream1 << fixed << setprecision(3) << pairdets[i].idstring << "," << pairdets[i].mag << "," << pairdets[i].band << ",";
    outstream1 << pairdets[i].obscode << "," << pairdets[i].index << "\n"; 
  }
  outstream1.close();

  // Load a vector storing the number of pair-partners found for each detection.
  for(i=0;i<pairdets.size();i++) {
    ppn = long_index(pairdets[i].indvec.size(),i);
    pair_partner_num.push_back(ppn);
  }
  // Sort the new vector by number of pair-partners
  sort(pair_partner_num.begin(), pair_partner_num.end(), lower_long_index());
  
  // Analyze paired detections in order of decreasing number of partners.
  // At the same time, write the output pair file, distinguishing
  // between real pairs comprising just two detections (indicated
  // with the letter P), and effective pairs that are really tracklets
  // created by fitting and averaging many points that lie along a
  // consistent trajectory.
  cout << "Constructing tracklets, and writing pairs to output file\n";
  outstream1.open(outpairfile);
  
  for(i=pairdets.size()-1; i>=0 ;i--) {
    pdct=pair_partner_num[i].index;
    istracklet=0; // Assume there is no tracklet unless one is confirmed to exist.
    if(pairdets[pdct].indvec.size() > mintrkpts-1) {
      if(DEBUG>=2) {
	cout << "Working on detection " << i << " = " << pdct << " with " << pair_partner_num[i].lelem << " = " << pairdets[pdct].indvec.size() << " pair partners:\n";
	for(j=0; j<pairdets[pdct].indvec.size(); j++) {
	  cout << pairdets[pdct].indvec[j] << ", ";
	}
	cout << "\n";
      }
      // The corresponding detection is paired with more than one
      // other detection.
      // Project all of these pairs relative to detection pdct,
      // storing x,y projected coordinates in axyvec.
      axyvec={};
      ppset={};
      for(j=0; j<pairdets[pdct].indvec.size(); j++) {
	detct = pairdets[pdct].indvec[j];
	if(pairdets[detct].indvec.size()>0) {
	  // Detection detct hasn't already been allocated to a tracklet.
	  distradec02(pairdets[pdct].RA, pairdets[pdct].Dec, pairdets[detct].RA, pairdets[detct].Dec, &dist, &pa);
	  dist *= 3600.0L; // Convert distance from degrees to arcsec.
	  xyind = xy_index(dist*sin(pa/DEGPRAD),dist*cos(pa/DEGPRAD),detct);
	  axyvec.push_back(xyind);
	  ppset.push_back(pairdets[detct]); // We need this vector (of type det_obsmag_indvec)
	                                    // mainly just to have some way to store the
	                                    // indices of mutually consistent pair partners
	                                    // on the next step
	}
      }
      if(DEBUG>=2) cout << "Loaded axyvec and ppset vectors OK, with sizes " << axyvec.size() << " and " << ppset.size() << "\n";
      if(axyvec.size() != ppset.size()) {
	cerr << "ERROR: vectors of projected and original\n";
	cerr << "pair partner candidates do not have the same length!\n";
	cerr << axyvec.size() << " != " << ppset.size() << "\n";
	return(3);
      }
      // Perform n^2 search to find the largest cluster of consistent points.
      // Load all the pair-partners for detection pdct into a vector
      // (this is mainly just to be able to use the index vectors)
      for(j=0; j<axyvec.size(); j++) {
	dtref = ppset[j].MJD - pairdets[pdct].MJD;
	if(dtref == 0) {
	  cerr << "ERROR: paired detections with no time separation!\n";
	  return(4);
	}
	// Make sure corresponding index vector in ppset is empty
	ppset[j].indvec = {};
	// Count consistent pair partners
	if(DEBUG>=2) cout << "Counting consistent pair partners\n";
	for(k=0; k<axyvec.size(); k++) {
	  if(j!=k) {
	    dt = ppset[k].MJD - pairdets[pdct].MJD;
	    dx = axyvec[k].x - axyvec[j].x*(dt/dtref);
	    dy = axyvec[k].y - axyvec[j].y*(dt/dtref);
	    dist = sqrt(dx*dx + dy*dy);
	    if(DEBUG>=2) cout << "Detection " << axyvec[j].index << ":" << axyvec[k].index << " dist = " << dist << "\n";
	    if(dist < 2.0*maxgcr) {
	      ppset[j].indvec.push_back(k);
	    }
	  }
	}
      }
      // Find the largest set of pair-partners lying along a line.
      if(DEBUG>=2) cout << "Find the largest set of pair-partners lying along a line.\n";
      biggest_tracklet=-1;
      tracklet_size=0;
      if(DEBUG>=2) cout << "size = " << ppset.size() << "\n";
      for(j=0; j<ppset.size(); j++) {
	if(DEBUG>=2) cout << j << ":" << ppset.size()-1 << " size = " << ppset[j].indvec.size() << " ";
	if(ppset[j].indvec.size()+2 > tracklet_size) {
	  tracklet_size = ppset[j].indvec.size()+2; //We add one for pdct, one for j, to get actual tracklet size
	  biggest_tracklet = j;
	  if(DEBUG>=2) cout << "bt = " << biggest_tracklet << ", size = " << tracklet_size << "\n";
	} else if(DEBUG>=2) cout << "not the biggest\n";
      }
      if(DEBUG>=2) cout << "Biggest tracklet is " << biggest_tracklet << ", which corresponds to " << axyvec[biggest_tracklet].index << ", with size " << tracklet_size << "\n";
      istracklet=0; // Assume there is no tracklet until one is confirmed to exist.
      if(tracklet_size <= mintrkpts) {
	istracklet=0;
      } else {
	// Perform linear fits to x and y vs time.
	// Load all the points from the biggest potential tracklet.
	track_mrdi_vec={}; // We need this vector purely so we can do a time-sort.
	                   // mrdi stands for MJD, RA, Dec, index
	// Load the reference point
	p3di = point3d_index(0.0l,0.0l,0.0l,pdct);
	track_mrdi_vec.push_back(p3di);
	// Load anchor point corresponding to biggest_tracklet
	p3di = point3d_index(ppset[biggest_tracklet].MJD - pairdets[pdct].MJD,axyvec[biggest_tracklet].x,axyvec[biggest_tracklet].y,axyvec[biggest_tracklet].index);
	track_mrdi_vec.push_back(p3di);
	// Load the other points
 	for(j=0; j<ppset[biggest_tracklet].indvec.size(); j++) {
	  p3di = point3d_index(ppset[ppset[biggest_tracklet].indvec[j]].MJD - pairdets[pdct].MJD,axyvec[ppset[biggest_tracklet].indvec[j]].x,axyvec[ppset[biggest_tracklet].indvec[j]].y,axyvec[ppset[biggest_tracklet].indvec[j]].index);
	  track_mrdi_vec.push_back(p3di);
	  //timevec.push_back(ppset[ppset[biggest_tracklet].indvec[j]].MJD - pairdets[pdct].MJD);
	  //xvec.push_back(axyvec[ppset[biggest_tracklet].indvec[j]].x);
	  //yvec.push_back(axyvec[ppset[biggest_tracklet].indvec[j]].y);
	  //detindexvec.push_back(axyvec[ppset[biggest_tracklet].indvec[j]].index);
	}
	// Sort track_mrdi_vec by time.
	sort(track_mrdi_vec.begin(), track_mrdi_vec.end(), lower_point3d_index_x());
	// Load time, x, y, and index vectors from sorted track_mrdi_vec.
	timevec=xvec=yvec={};
	detindexvec={};
	for(j=0;j<track_mrdi_vec.size();j++)
	  {
	    timevec.push_back(track_mrdi_vec[j].x);
	    xvec.push_back(track_mrdi_vec[j].y);
	    yvec.push_back(track_mrdi_vec[j].z);
	    detindexvec.push_back(track_mrdi_vec[j].index);
	  }
 	if(DEBUG>=2) {
	  cout << "First iteration linear fit vectors:\n";
	  for(j=0; j<timevec.size(); j++) {
	    cout << detindexvec[j] << " " << timevec[j] << " " << xvec[j] << " " << yvec[j] << "\n";
	  }
	}

	// Perform fit to projected x coordinate as a function of time
	linfituw01(timevec, xvec, slopex, interceptx);
 	// Perform fit to projected y coordinate as a function of time
	linfituw01(timevec, yvec, slopey, intercepty);
	// Load vector of residuals
	fiterr = {};
	for(j=0; j<timevec.size(); j++) {
	  fiterr.push_back(sqrt(DSQUARE(timevec[j]*slopex+interceptx-xvec[j]) + DSQUARE(timevec[j]*slopey+intercepty-yvec[j])));
	}
	// Ditch duplicate times, if there are any
	istimedup=1; // Guilty until proven innocent
	while(istimedup==1 && timevec.size()>=mintrkpts+1) {
	  istimedup=0;
	  j=1;
	  while(j<timevec.size() && istimedup==0) {
	    if(fabs(timevec[j] - timevec[j-1]) < IMAGETIMETOL/SOLARDAY) {
	      istimedup=1; // Point j and j-1 are time-duplicates.
	      // Mark for rejection whichever one has the largest fitting error
	      if(fiterr[j]>=fiterr[j-1]) worstpoint = j; 
	      else worstpoint = j-1;
	    }
	    j++;
	  }
	  if(istimedup==1) {
	    // Reject the bad point
	    trkptnum=timevec.size();
	    for(j=worstpoint; j<trkptnum-1; j++) {
	      timevec[j] = timevec[j+1];
	      xvec[j] = xvec[j+1];
	      yvec[j] = yvec[j+1];
	      detindexvec[j] = detindexvec[j+1];
	    }
	    trkptnum--;
	    timevec.resize(trkptnum);
	    xvec.resize(trkptnum);
	    yvec.resize(trkptnum);
	    detindexvec.resize(trkptnum);
	    // Re-do linear fit
	    // Perform fit to projected x coordinate as a function of time
	    linfituw01(timevec, xvec, slopex, interceptx);
	    // Perform fit to projected y coordinate as a function of time
	    linfituw01(timevec, yvec, slopey, intercepty);
	    // Load vector of residuals
	    fiterr = {};
	    for(j=0; j<timevec.size(); j++) {
	      fiterr.push_back(sqrt(DSQUARE(timevec[j]*slopex+interceptx-xvec[j]) + DSQUARE(timevec[j]*slopey+intercepty-yvec[j])));
	    }
	  }
	}
	// Find worst error.  
	worsterr = 0.0l;
	for(j=0; j<timevec.size(); j++) {
	  if(fiterr[j]>worsterr) {
	    worsterr = fiterr[j];
	    worstpoint = j;
	  }
	}
	// Reject successive points until either there are only three left
	// or the worst error drops below maxgcr.
	while(worsterr>maxgcr && timevec.size()>3 && timevec.size()>=mintrkpts) {
	  // Reject the worst point
	  trkptnum=timevec.size();
	  for(j=worstpoint; j<trkptnum-1; j++) {
	    timevec[j] = timevec[j+1];
	    xvec[j] = xvec[j+1];
	    yvec[j] = yvec[j+1];
	    detindexvec[j] = detindexvec[j+1];
	  }
	  trkptnum--;
	  timevec.resize(trkptnum);
	  xvec.resize(trkptnum);
	  yvec.resize(trkptnum);
	  detindexvec.resize(trkptnum);	  
	  // Perform fit to projected x coordinate as a function of time
	  linfituw01(timevec, xvec, slopex, interceptx);
	  // Perform fit to projected y coordinate as a function of time
	  linfituw01(timevec, yvec, slopey, intercepty);
	  // Load vector of residuals
	  fiterr = {};
	  for(j=0; j<timevec.size(); j++) {
	    fiterr.push_back(sqrt(DSQUARE(timevec[j]*slopex+interceptx-xvec[j]) + DSQUARE(timevec[j]*slopey+intercepty-yvec[j])));
	  }
	  // Find worst error.  
	  worsterr = 0.0l;
	  for(j=0; j<timevec.size(); j++) {
	    if(fiterr[j]>worsterr) {
	      worsterr = fiterr[j];
	      worstpoint = j;
	    }
	  }
	}
	if(worsterr<=maxgcr && timevec.size()>=3 && timevec.size()>=mintrkpts) {
	  // We succeeded in finding a tracklet with no time-duplicates, and
	  // no outliers beyond maxgcr. Prepare to write it to the pair file.
	  // Select points that will represent this tracklet.
	  instep = (timevec.size()-1)/4;
	  rp1 = instep;
	  rp2 = timevec.size()-1-instep;
	  if(rp1==rp2) {
	    cerr << "ERROR: both representative points for a tracklet are the same!\n";
	    cerr << "size, instep, rp1, rp2: " << timevec.size() << " " << instep << " " << rp1 << " " << rp2 << "\n";
	    return(4);
	  }
	  // Calculate angular velocity in deg/day. The slope values
	  // correspond to velocities in arcsec/day.
	  angvel = sqrt(slopex*slopex + slopey*slopey)/3600.0l;
	  
	  // Determine improved RA, Dec based on tracklet fit for the representative points
	  // Calculated projected x, y at rp1
	  dx = timevec[rp1]*slopex + interceptx;
	  dy = timevec[rp1]*slopey + intercepty;
	  // Calculate equivalent celestial position angle.
	  if(dx==0l && dy>=0l) pa = 0.0l;
	  else if(dx==0l && dy<0l) pa = M_PI;
	  else if(dx>0l) pa = M_PI/2.0l - atan(dy/dx);
	  else if(dx<0l) pa = 3.0l*M_PI/2.0l - atan(dy/dx);
	  else {
	    cerr << "ERROR: logical impossibility while trying to solve for PA\n";
	    cerr << "dx = " << dx << " dy = " << dy << "\n";
	  }
	  dist = sqrt(dx*dx + dy*dy)/3600.0l; // renders distance in degrees, not arcsec.
	  pa*=DEGPRAD; // position angle in degrees, not radians.
	  arc2cel01(pairdets[pdct].RA, pairdets[pdct].Dec, dist, pa, outra1, outdec1);
	  if(!isnormal(outra1)) {
	    cerr << "NAN WARNING: dx, dy, dist, pa: " << dx << " " << dy << " " << dist << " " << pa << "\n";
	  }
	  // Calculated projected x, y at rp2
	  dx = timevec[rp2]*slopex + interceptx;
	  dy = timevec[rp2]*slopey + intercepty;
	  // Calculate equivalent celestial position angle.
	  if(dx==0l && dy>=0l) pa = 0.0l;
	  else if(dx==0l && dy<0l) pa = M_PI;
	  else if(dx>0l) pa = M_PI/2.0l - atan(dy/dx);
	  else if(dx<0l) pa = 3.0l*M_PI/2.0l - atan(dy/dx);
	  else {
	    cerr << "ERROR: logical impossibility while trying to solve for PA\n";
	    cerr << "dx = " << dx << " dy = " << dy << "\n";
	  }
	  dist = sqrt(dx*dx + dy*dy)/3600.0l; // renders distance in degrees, not arcsec.
	  pa*=DEGPRAD; // position angle in degrees, not radians.
	  arc2cel01(pairdets[pdct].RA, pairdets[pdct].Dec, dist, pa, outra2, outdec2);
	  // Calculate total angular arc
	  distradec02(outra1, outdec1, outra2, outdec2, &dist, &pa);
	  dist *= 3600.0l;
	  if(dist>=minarc && angvel>=minvel && angvel<=maxvel) {
	    // Write out representative pair, followed by RA, Dec and the total number of constituent points
	    // representative pair
	    outstream1 << "T " << detindexvec[rp1] << " " <<  detindexvec[rp2] << " ";
	    // RA1, Dec1, RA2, Dec2
	    outstream1 << fixed << setprecision(6) << outra1 << " " << outdec1 << " " <<  outra2 << " " << outdec2 << " ";
	    // Number of points in final, refined tracklet.
	    outstream1 << detindexvec.size() << "\n";
	    // Now write out the detection indices for this full number of points,
	    // and wipe all the associated index vectors.
	    for(j=0; j<detindexvec.size(); j++) {
	      outstream1 << detindexvec[j] << "\n";
	      pairdets[detindexvec[j]].indvec = {};
	    }
	    istracklet=1;
	    // Close if-statement confirming that a bona fide,
	    // aligned tracklet was found and written to the output file.
	  } else {
	    istracklet=0;
	    cout << "A tracklet was rejected: arc = " << dist << " < " << minarc << " or angvel = " << angvel << " < " << minvel << "\n";
	  }
	} else istracklet=0;
	// Close else-statement confirming there was a candidate for
	// being an aligned tracklet.
      }
      // Close if-statement checking that detection i has more than
      // one pair-partner, and hence COULD be part of a tracklet
    } else istracklet=0;
    if(istracklet==0 && mintrkpts==2) {
      // Write out all the pairs as normal
      for(j=0; j<pairdets[pdct].indvec.size(); j++) {
	k=pairdets[pdct].indvec[j];
	// Calculate angular arc and angular velocity
	distradec02(pairdets[pdct].RA,pairdets[pdct].Dec,pairdets[k].RA,pairdets[k].Dec, &dist, &pa);
	angvel = dist/fabs(pairdets[pdct].MJD-pairdets[k].MJD); // Degrees per day
	dist *= 3600.0l; // Arcseconds
	if(pairdets[k].indvec.size()>0 && k>pdct && angvel>=minvel && dist>=minarc && angvel<=maxvel) {
	  outstream1 << "P " <<  pdct << " " <<  k << "\n";
	} else if(angvel<minvel || dist<minarc) {
	  cout << "A pair was rejected: arc = " << dist << " < " << minarc << " or angvel = " << angvel << " < " << minvel << "\n";
	}
      }
    }
    // Close loop over all detections
  }
  outstream1.close();
  return(0);
}
