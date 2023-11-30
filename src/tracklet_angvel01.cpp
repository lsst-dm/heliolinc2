// November 28, 2023: tracklet_angvel01.cpp:
// Use the same file-reading machinery as heliolinc_danby to
// read make_tracklets output, translate each tracklet from
// ordinary celestial coordinates into geocentric ecliptic
// coordinates, calculate the angular velocities in these
// coordinates, and then store them in bins of constant
// solar elongation and geocentric ecliptic latitude.
// This converts the tracklets into the same parameter space
// used by mpcorb_anvgel01.cpp, enabling apples-to-apples
// comparison to ensure that tracklet-culling on the basis of
// mpcorb_anvgel01 output will be effective.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: tracklet_angvel01 -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection file -obspos observer_position_file -latrange range +/- for ecliptic latitude -latstep step size for ecliptic latitude -elongstep step size for solar elongation -outfile output file -verbose verbosity\n";
}
    
int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  vector <hldet> pairdets = {};
  vector <hlimage> image_log;
  vector <tracklet> tracklets;
  vector <longpair> trk2det;
  vector <EarthState> earthpos;
  point3d Earth2sun = point3d(0,0,0);
  point2d angvelpt = point2d(0l,0l);
  string imfile,pairdetfile,trackletfile,trk2detfile,planetfile;
  int elongstep=10;
  int latstep=5;
  int latmax = 30;
  ofstream outstream1;
  string outfile;
  string outname;
  long i=0;
  long j=0;
  int status=0;
  long pairnum,pairct,imnum,i1,i2;
  pairnum = pairct = imnum = i1 = i2 = 0;
  double b1,b2,l1,l2; // b is ecliptic latitude, l is longitude.
  b1 = b2 = l1 = l2 = 0.0l;
  double sunlambda, sunbeta, sunelong, lambdavel, betavel;
  sunlambda = sunbeta = sunelong = lambdavel = betavel = 0.0l;
  char nsfx[64];
  string numsuffix;
  vector <int> elongcen;
  vector <int> latcen;
  vector <point2d> velvec;
  vector <vector <point2d>> velmat;
  int lat,elong,latnum,elongnum,latct,elongct;
  double dlatct,delongct;
  int verbose=0;
  double timediff=0.0l;

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" || string(argv[i]) == "-detpairs") {
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
    }  else if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
	//There is still something to read;
	imfile=argv[++i];
	i++;
      }
      else {
	cerr << "Image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-tracklets" || string(argv[i]) == "-tracklet" || string(argv[i]) == "-tf" || string(argv[i]) == "-trkfile" || string(argv[i]) == "-trackletfile" || string(argv[i]) == "--trackletfile") {
      if(i+1 < argc) {
	//There is still something to read;
	trackletfile=argv[++i];
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
	i++;
      }
      else {
	cerr << "Output pair file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obspos" || string(argv[i]) == "-op" || string(argv[i]) == "-obsvec" || string(argv[i]) == "--observer" || string(argv[i]) == "--observer_position" || string(argv[i]) == "--observer_statevec") {
      if(i+1 < argc) {
	//There is still something to read;
	planetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Observer position file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-latrange" || string(argv[i]) == "-latmax" || string(argv[i]) == "-lrange" || string(argv[i]) == "--latrange" || string(argv[i]) == "--latmax" || string(argv[i]) == "-lmax") {
      if(i+1 < argc) {
	//There is still something to read;
	latmax=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Latitude range keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-latstep" || string(argv[i]) == "-lstep" || string(argv[i]) == "-ls" || string(argv[i]) == "--latstep" || string(argv[i]) == "--lstep" ) {
      if(i+1 < argc) {
	//There is still something to read;
	latstep=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Latitude step size keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-elongstep" || string(argv[i]) == "-estep" || string(argv[i]) == "-es" || string(argv[i]) == "--elongstep" || string(argv[i]) == "--estep" ) {
      if(i+1 < argc) {
	//There is still something to read;
	elongstep=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Solar elongation step size keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile") {
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
  
  if(argc<11)
    {
      cerr << "Too few arguments even for minimalist invocation:\n";
      show_usage();
      return(1);
    }
  
  cout.precision(17);  
  cout << "input image file " << imfile << "\n";
  cout << "input detection file " << pairdetfile << "\n";
  cout << "input tracklet file " << trackletfile << "\n";
  cout << "input trk2det file " << trk2detfile << "\n";
  cout << "input observer position file " << planetfile << "\n";

  // Catch required parameters if missing
  if(pairdetfile.size()<=0) {
    cout << "\nERROR: input detection file is required\n";
    show_usage();
    return(1);
  } else if(trackletfile.size()<=0) {
    cout << "\nERROR: input tracklet file is required\n";
    show_usage();
    return(1);
  } else if(trk2detfile.size()<=0) {
    cout << "\nERROR: input trk2det file is required\n";
    show_usage();
    return(1);
  } else if(planetfile.size()<=0) {
    cout << "\nERROR: input observer position file is required:\n";
    cout << "e.g. Earth1day2020s_02a.txt\n";
    show_usage();
    return(1);
  }
  
  // Load vectors holding central ecliptic latitudes and solar elongations
  latnum = 2*latmax/latstep + 1;
  if((latnum-1)*latstep != 2*latmax) {
    cerr << "ERROR: latitude step size " << latstep << " is not evenly divisible\n";
    cerr << "into the range +/- " << latmax << "\n";
    return(1);
  }
  elongnum = 360/elongstep;
  if(elongnum*elongstep != 360) {
    cerr << "ERROR: solar elongation step size " << elongstep << " is not evenly divisible into 360\n";
    return(1);
  }
  cout << "There will be " << latnum << " steps in ecliptic latitude, and " << elongnum << " steps in solar elongation,\n";
  cout << "for a total of " << latnum*elongnum << " steps\n";
  lat=-latmax;
  while(lat<=latmax) {
    elong = elongstep/2;
    while(elong<=(360-elongstep/2)) {
      latcen.push_back(lat);
      elongcen.push_back(elong);
      velmat.push_back(velvec);
      elong += elongstep;
    }
    lat += latstep;
  }
 
  cout << "Heliocentric ephemeris for Earth is named " << planetfile << "\n";
  status = read_horizons_csv(planetfile, earthpos);
  if(status!=0) {
    cerr << "ERROR: could not successfully Earth ephemeris file " << planetfile << "\n";
    cerr << "read_horizons_csv returned status = " << status << ".\n";
   return(1);
  } 

  detvec={};
  status=read_pairdet_file(pairdetfile, detvec, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << pairdetfile << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << detvec.size() << " data lines from paired detection file " << pairdetfile << "\n";
  
  image_log={};
  status=read_image_file2(imfile, image_log);
  if(status!=0) {
    cerr << "ERROR: could not successfully read image file " << imfile << "\n";
    cerr << "read_image_file2 returned status = " << status << ".\n";
   return(1);
  }
  imnum = image_log.size();
  cout << "Read " << imnum << " data lines from image file " << imfile << "\n";
  
  tracklets={};
  status=read_tracklet_file(trackletfile, tracklets, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read tracklet file " << trackletfile << "\n";
    cerr << "read_tracklet_file returned status = " << status << ".\n";
   return(1);
  }
  pairnum = tracklets.size();
  cout << "Read " << pairnum << " data lines from tracklet file " << trackletfile << "\n";
  
  trk2det={};
  status=read_longpair_file(trk2detfile, trk2det, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read trk2det file " << trk2detfile << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << trk2det.size() << " data lines from trk2det file " << trk2detfile << "\n";

  // Main loop over all tracklets
  for(pairct=0; pairct<pairnum; pairct++) {
    // Obtain indices to the image_log and heliocentric distance vectors.
    i1=tracklets[pairct].Img1;
    i2=tracklets[pairct].Img2;
    if(i1<0 || i1>=imnum || i2<0 || i2>=imnum) {
      cerr << "ERROR: image indices " << i1 << ", " << i2 << " not in range 0 to " << imnum << " spanned by image log\n";
      return(2);
    }
    // Convert the representative points to ecliptic coordinates.
    poleswitch02(tracklets[pairct].RA1,tracklets[pairct].Dec1, NEPRA, NEPDEC, 90.0l, l1, b1);
    poleswitch02(tracklets[pairct].RA2,tracklets[pairct].Dec2, NEPRA, NEPDEC, 90.0l, l2, b2);
    // Calculate time difference between the observations
    timediff = (image_log[i2].MJD - image_log[i1].MJD);
    // Find the geocentric ecliptic longitude of the sun.
    Earth2sun.x = -image_log[i1].X;
    Earth2sun.y = -image_log[i1].Y;
    Earth2sun.z = -image_log[i1].Z;
    celedeproj01(Earth2sun, &sunlambda, &sunbeta);
    lambdavel = (l2-l1)/timediff;
    betavel = (b2-b1)/timediff;
    angvelpt = point2d(lambdavel,betavel);
    sunelong = l1 - sunlambda;
    while(sunelong>=360.0l) sunelong-=360.0l;
    while(sunelong<0.0l) sunelong+=360.0l;
    dlatct = (b1+double(latmax))/double(latstep) + 0.5l;
    latct = dlatct;
    delongct = double(sunelong)/double(elongstep);
    elongct = delongct;
    if(verbose>=1) cout << "angvel " << lambdavel << " " << betavel << " intcoords: " << latct << " " << elongct << " " << latct*elongnum + elongct << "\n";
    if(latct>=0 && latct<latnum && elongct>=0 && elongct<elongnum) {
      velmat[latct*elongnum + elongct].push_back(angvelpt);
    }
  }
  for(i=0;i<latnum*elongnum;i++) {
    // Construct the name of the output file
    if(latcen[i]>=0) {
      int dnum=sprintf(nsfx,"%d",latcen[i]);
      numsuffix = nsfx;
      if(dnum<=1) numsuffix = "n0" + numsuffix;
      else numsuffix = "n" + numsuffix;
    } else if(latcen[i]<0) {
      int dnum=sprintf(nsfx,"%d",-latcen[i]);
      numsuffix = nsfx;
      if(dnum<=1) numsuffix = "s0" + numsuffix;
      else numsuffix = "s" + numsuffix;
    }
    outname = outfile + numsuffix + "_";
    int dnum=sprintf(nsfx,"%d",elongcen[i]);
    numsuffix = nsfx;
    if(dnum<=1) numsuffix = "00" + numsuffix;
    else if(dnum==2) numsuffix = "0" + numsuffix;
    outname = outname + numsuffix + ".txt";
    cout << "Writing output file " << i << " of " << velmat.size() << ", named " << outname << "\n";
    outstream1.open(outname);
    for(j=0;j<long(velmat[i].size());j++) {
      outstream1 << velmat[i][j].x << " " << velmat[i][j].y << "\n";
    }
    outstream1.close();
  }

  return(0);
}
