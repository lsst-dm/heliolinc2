// December 17, 2021: maketrack03a.cpp
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
#define MINOBSINTERVAL 1.0 // Minimum time-between-images in seconds
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME 1.5 // Default max inter-image time interval
                    // for tracklets, in hours (will be converted
                    // to days before use).
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)


      
static void show_usage()
{
  cerr << "Usage: maketrack03a -dets detfile -imgs imfile -outimgs output image file/ \n";
  cerr << "-pairs pairfile -pairdets paired detection file/ \n";
  cerr << "-imrad image radius(deg) -maxtime max inter-image time interval (hr)/ \n";
  cerr << "-maxvel maximum angular velocity (deg/day) -earth earthfile -obscoords lon plxcos plxsin\n";
  cerr << "\nor, at minimum\n\n";
  cerr << "maketrack03a -dets detfile -earth earthfile\n";
  cerr << "Note well that the minimum invocation will leave a bunch of things\n";
  cerr << "set to defaults that may not be what you want.\n";
}
    
int main(int argc, char *argv[])
{
  det_OC_index o1 = det_OC_index(0L,0L,0L,0L,0L,0L,"null",0);
  vector <det_OC_index> detvec = {};
  vector <det_OC_index> pairdets = {};
  img_log02 imlog = img_log02(0.0,0.0,0.0,0,0);
  vector <img_log02> img_log = {};
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
  string idstring;
  string lnfromfile;
  int i = 0;
  int j = 0;
  int imct=0;
  int imctp=0;
  int imnum=0;
  long detnum=0;
  long num_dets=0;
  long detct=0;
  long detctp=0;
  int startind=0;
  int endind=0;
  int reachedeof = 0;
  char c='0';
  long double MJD,RA,Dec;
  MJD = RA = Dec = 0.0L;
  double maxvel = MAXVEL; // Max angular velocity in deg/day
  double maxtime = MAXTIME; // Max time interval a tracklet could span,
                            // in days.
  double maxdist = MAXVEL*MAXTIME/24.0; // Max angular distance a tracklet
                                   // could span, in degrees.
  double imrad = IMAGERAD; // radius from image center to most distant corner (deg).
  string indetfile;
  string inimfile;
  string outimfile;
  string earthfile;
  string outpairfile="outpairfile01.txt";
  string pairdetfile="pairdetfile01.txt";
  long double obslon = 289.26345L;
  long double plxcos = 0.865020L;
  long double plxsin = -0.500901L;
  long lct=0;

  if(argc<5)
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
    } else if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-img" || string(argv[i]) == "--imgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
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
	if(!isnormal(maxtime) || maxtime<=0.0) {
	  cerr << "Error: invalid maximum inter-image time interval\n";
	  cerr << "(" << maxtime << " hr) supplied: must be strictly positive.\n";
	  return(2);
	}      
      } else {
	cerr << "Output maximum inter-image time interval\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxvel") {
      if(i+1 < argc) {
	//There is still something to read;
        maxtime=stod(argv[++i]);
	i++;
	if(!isnormal(maxvel) || maxvel<=0.0) {
	  cerr << "Error: invalid maximum angular velocity\n";
	  cerr << "(" << maxvel << "deg/day) supplied: must be strictly positive.\n";
	  return(2);
	}
      }
      else {
	cerr << "Output maximum inter-image time interval\nkeyword supplied with no corresponding argument\n";
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
    }  else if(string(argv[i]) == "-o" || string(argv[i]) == "-obs" || string(argv[i]) == "-obscoords" || string(argv[i]) == "-obs" || string(argv[i]) == "--obscoords" || string(argv[i]) == "--observatory" || string(argv[i]) == "--observatorycoords") {
      if(i+1 < argc) {
	//There is still something to read;
	obslon=stod(argv[++i]);
      }
      else {
	cerr << "Observatory coordinates keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	plxcos=stod(argv[++i]);
      }
      else {
	cerr << "Observatory coordinates keyword supplied with only one of the three required arguments\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	plxsin=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Observatory coordinates keyword supplied with only two of the three required arguments\n";
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
  
  cout << "indet file " << indetfile << "\n";
  cout << "inimage file " << inimfile << "\n";
  cout << "output image file " << outimfile << "\n";
  cout << "pairfile file " << outpairfile << "\n";
  cout << "paired detection file " << pairdetfile << "\n";
  cout << "Heliocentric ephemeris file for the Earth: " << earthfile << "\n";
  cout << "image radius " << imrad << "\n";
  cout << "max time interval " << maxtime << "\n";
  cout << "maxvel " << maxvel << "\n";
  maxtime/=24.0; /*Unit conversion from hours to days*/
  maxdist = maxtime*maxvel;
  
  ifstream instream1 {indetfile};
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
    while(i<lnfromfile.size() && reachedeof == 0) {
      string stest;
      c='0';
      while(i<lnfromfile.size() && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==IDCOL) idstring=stest;
      if(j==MJDCOL) MJD=stold(stest);
      else if(j==RACOL) RA=stold(stest);
      else if(j==DECCOL) Dec=stold(stest);
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0) {
      // cout<<"MJD, RA, Dec: " << MJD-floor(MJD) << " " << RA << " " << Dec << "\n";
      o1=det_OC_index(MJD,RA,Dec,0L,0L,0L,idstring,-lct);
      detvec.push_back(o1);
    }
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  // time-sort the detection vector
  sort(detvec.begin(), detvec.end(), early_det_OC_index());

  // Get image information.
  if(inimfile.size()>0)
    {
      // Read input image file: MJD, RA, Dec:
        ifstream instream1 {inimfile};
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
	    string stest;
	    c='0';
	    while(i<lnfromfile.size() && c!=',' && c!=' ' && c!='\n' && c!=EOF) {
	      // We allow the file to be delimited by comma or space.
	      c=lnfromfile[i];
	      if(c!=',' && c!=' ' && c!='\n' && c!=EOF) stest.push_back(c);
	      i++;
	    }
	    // We just finished reading something
	    j++;
	    if(j==1) MJD=stod(stest); // We assume we have MJD, RA, Dec.
	    else if(j==2) RA=stod(stest);
	    else if(j==3) Dec=stod(stest);
	  }
	  if((reachedeof == 0 || reachedeof == 1) && MJD>0.0) {
	    // Requirement of MJD>0.0 tests that we read a plausibly
	    // valid line.
	    imlog=img_log02(MJD,RA,Dec,0,0);
	    img_log.push_back(imlog);
	  }
	}
	if(reachedeof==1) { 
	  cout << "File read successfully to the end.\n";
	}
	else if(reachedeof==-1) cout << "Warning: file read failed\n";
	else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
	else cout << "Warning: unknown file read problem\n";
	// time-sort the image file
	sort(img_log.begin(), img_log.end(), early_imlg2());
	// find the indices in the time-sorted detection file
	// that correspond to the earliest and latest detections
	// on each image, and load these values into imglog02.
	detct=0;
	for(imct=0;imct<img_log.size();imct++) {
	  while(detct<detvec.size() && detvec[detct].MJD < img_log[imct].MJD-IMAGETIMETOL/SOLARDAY) detct++; //Not on any image
	  if(detct<detvec.size() && fabs(detvec[detct].MJD-img_log[imct].MJD)<=IMAGETIMETOL/SOLARDAY) {
	    // This should be the first detection on image imct.
	    img_log[imct].startind = detct;
	    while(detct<detvec.size() && fabs(detvec[detct].MJD-img_log[imct].MJD)<=IMAGETIMETOL/SOLARDAY) detct++; //Still on this same image
	    // This should be the first detection on the next image
	    img_log[imct].endind = detct;
	  }
	}
    } else {
    // No input image file was supplied: we have to create one from
    // the sorted detection file.
    mjdnorm = 1.0;
    mjdmean = detvec[0].MJD;
    startind=0;
    for(i=1;i<detvec.size();i++) {
      tdelt = detvec[i].MJD - detvec[i-1].MJD;
      if(tdelt < MINOBSINTERVAL/SOLARDAY) {
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
	imlog = img_log02(mjdmean,0.0,0.0,startind,endind);
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
      imlog = img_log02(mjdmean,0.0,0.0,startind,endind);
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
      vector <det_OC_index> imobs = {};
      num_dets=0;
      p3avg = point3d(0,0,0);
      while(detct<detnum && detvec[detct].MJD < img_log[imct].MJD + MINOBSINTERVAL/SOLARDAY) {
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
      cout << "Image " << imct << " of " << img_log.size() << ": " << num_dets << " = " << img_log[imct].endind-img_log[imct].startind ;
      cout << " from " << img_log[imct].startind << " to " << img_log[imct].endind << " of " << detvec.size() << ".\n";
      fflush(stdout);
      imct++;
    }
  }

  EarthMJD={};
  Earthpos={};
  Earthvel={};
  read_horizons_fileLD(earthfile,EarthMJD,Earthpos,Earthvel);
  cout << "Finished reading heliocentric ephemeris file " << earthfile << " for the Earth.\n";

  // Test: print out time-sorted detection table.
  ofstream outstream1 {"testjunk01.txt"};
  for(i=0;i<detvec.size();i++) {
    outstream1 << fixed << setprecision(6) << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << "\n";
  }
  
  if(outimfile.size()>0)
    {
      // Write and print image log table
      ofstream outstream2 {outimfile};
      for(imct=0;imct<img_log.size();imct++)
	{
	  //	  cout << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
	  //	  cout << " " << img_log[imct].Dec << " " << img_log[imct].startind << " " << img_log[imct].endind << " ";
	  //	  cout << img_log[imct].endind-img_log[imct].startind << "\n";
	  outstream2 << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
	  outstream2 << " " << img_log[imct].Dec << " " << img_log[imct].startind << " ";
	  outstream2 << img_log[imct].endind << "\n";
	}
    }
  // Calculate observer's heliocentric position at the time of each image.
  observer_heliopos={};
  for(imct=0;imct<img_log.size();imct++)
    {
      observer_barycoords01LD(img_log[imct].MJD, 5, obslon, plxcos, plxsin, EarthMJD, Earthpos, outpos);
      observer_heliopos.push_back(outpos);
    }
  // assign observer heliocentric position to each detection.
  for(imct=0;imct<img_log.size();imct++) {
    for(i=img_log[imct].startind;i<=img_log[imct].endind;i++) {
      detvec[i].x = observer_heliopos[imct].x;
      detvec[i].y = observer_heliopos[imct].y;
      detvec[i].z = observer_heliopos[imct].z;
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
	cerr << "ERROR: Negative time difference " << timediff << " encountered between images " << imct << " and " << imtarg << "\n";
	cerr << "testprog02a ABORTING!\n" ;
      }
      double imcendist = distradec01(img_log[imct].RA, img_log[imct].Dec, img_log[imtarg].RA, img_log[imtarg].Dec);
      if(imcendist<2.0*imrad+maxvel*timediff) {
	cout << "  pairs may exist between images " << imct << " and " << imtarg << ": dist = " << imcendist << ", timediff = " << timediff << "\n";
	imagematches.push_back(imtarg);
      }
      imtarg++;
    }
    cout << "Looking for pairs for image " << imct << ": " << imagematches.size() << " later images are worth searching\n";
    if(imagematches.size()>0) {
      // Search is worth doing. Project all the detections
      // on image A.
      xy_index xyind=xy_index(0.0, 0.0, 0);
      vector <xy_index> axyvec = {};
      double dist,pa;
      int dettarg=0;
      for(detct=img_log[imct].startind ; detct<img_log[imct].endind ; detct++) {
	distradec02(img_log[imct].RA, img_log[imct].Dec,detvec[detct].RA,detvec[detct].Dec,&dist,&pa);
	xyind = xy_index(dist*sin(pa/DEGPRAD),dist*cos(pa/DEGPRAD),detct);
	axyvec.push_back(xyind);
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
	cout << "Looking for pairs between " << axyvec.size() << " detections on image " << imct << " and " << kdvec.size() << " on image " << imtarg << "\n";
	for(detct=0 ; detct<axyvec.size() ; detct++) {
	  vector <long> indexvec = {};
	  if(isnormal(axyvec[detct].x) &&isnormal(axyvec[detct].y)) {
	     kdrange01(kdvec,axyvec[detct].x,axyvec[detct].y,range,indexvec);
	  } else {
	    cout << "WARNING: detection " << detct << " on image " << imct << " not normal: " << axyvec[detct].x << " " << axyvec[detct].y << "\n";
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
		detvec[kdvec[matchpt].point.index].index = pdct; // Re-assign index to apply to paired detection vector
		pdct++; // Increment count of paired detections
		if(pdct!=pairdets.size()) cerr << "\nWARNING: PAIRED DETECTION MISMATCH: " << pdct << " vs " << pairdets.size() << "\n\n";
	      }
	      // Write index values for both components of the
	      // new pair to the pair vector, regardless of whether
	      // the index values are pre-existing or newly assigned.
	      onepair = longpair(detvec[axyvec[detct].index].index,detvec[kdvec[matchpt].point.index].index);
	      pairvec.push_back(onepair);
	      pairct++;
	      apct++;
	    }
	    // Close if-statement checking if image A detection was matched to anything.
	  }
	  // Close loop over detections on source image (image A)
	}
	// Close loop over image B candidates
      }
      // Close if-statement checking if any images could match image A      
    }
    cout << "Image " << imct << ": found " << adetct << " newly paired detections and a total of " << apct << " pairs.\n";
    // Close loop over images for image A
  }
  cout << "Test count of paired detections: " << pdct << " " << pairdets.size() << "\n";
  cout << "Test count of pairs: " << pairct << " " << pairvec.size() << "\n";

  // Write paired detections vector to file
  cout << "Writing paired detections file\n";
  ofstream outstream3 {pairdetfile};
  outstream3.precision(12);
  for(i=0;i<pairdets.size();i++) {
    outstream3 << pairdets[i].MJD << " " << pairdets[i].RA << " " << pairdets[i].Dec << " ";
    outstream3 << pairdets[i].x << " " << pairdets[i].y << " " << pairdets[i].z << " ";
    outstream3 << pairdets[i].idstring << " " << pairdets[i].index << "\n"; 
  }
  // Write pair vector to file
  cout << "Writing pairs to output file\n";
  ofstream outstream4 {outpairfile};
  for(i=0;i<pairvec.size();i++) {
    outstream4 <<  pairvec[i].i1 << " " <<  pairvec[i].i2 << "\n";
  }
  
  // Write index'ed input detection vector to file as a test.
  cout << "Writing indexed input detection vector as a test\n";
  ofstream outstream5 {"indextest01.txt"};
  for(i=0;i<detvec.size();i++) {
    outstream5 << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << " " << detvec[i].index << "\n"; 
  }

  return(0);
}
