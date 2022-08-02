/// July 20, 2022: sunrelative_eclipcoords01a.cpp
// Given an input file providing MJD, RA, and Dec for some set
// of detections (or image boresights, etc), calculate the RA and Dec
// of the sun (as observed from the geocenter) at each MJD. Transform
// the RA and Dec of the input points and the sun into ecliptic coordinates.
// Shift to sun-relative ecliptic coordinates (still geocentric) by
// subtracting the solar ecliptic longitude from the ecliptic longitude
// of each input point. Output the resulting sun-relative geocentric
// ecliptic coordinates for every point. Also, hammer-project these
// coordinates centered both on the sun and on opposition, and include
// the x,y coordinates of both hammer projections in the output file.
//
// Hence, the output file will have the following columns:
// 1. Input MJD
// 2. Input RA
// 3. Input Dec
// 4. Output sun-relative geocentric ecliptic longitude
// 5. Output geocentric ecliptic latitude.
// 6. Output hammer-projection x with the sun at plot center
// 7. Output hammer-projection y with the sun at plot center
// 8. Output hammer-projection x with opposition at plot center
// 9. Output hammer-projection y with opposition at plot center

#define INTERP_POLYORDER 3 // Order of polynomials for interpolation in planetpos files.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: sunrelative_eclipcoords01a -sunfile sunfile -earthfile earthfile -incoords infile -outfile outfile\n";
}
    
int main(int argc, char *argv[])
{
  string sunfile, earthfile, infile, outfile;
  vector <det_bsc> indets;
  det_bsc det1 = det_bsc(0l,0l,0l);
  int i=0;
  int detnum=0;
  int detct=0;
  double MJD,RA,Dec,detmjd,esdist;
  MJD = RA = Dec = detmjd = esdist = 0l;
  double sunRA,sunDec;
  sunRA = sunDec = 0l;
  vector <double> mjd_ephem;
  vector <double> mjd_test;
  vector <point3d> Sunpos; 
  vector <point3d> Sunvel; 
  vector <point3d> Earthpos; 
  vector <point3d> Earthvel;
  point3d Sunposnow = point3d(0l,0l,0l);
  point3d Earthposnow = point3d(0l,0l,0l);
  double SunEclipLon, SunEclipLat,EclipLon,EclipLat;
  SunEclipLon = SunEclipLat = EclipLon = EclipLat = 0l;
  double lambda,phi,x,y;
  lambda = phi = x = y = 0l;
  
  if(argc<9) {
    show_usage();
    return(1);
  }


  // Parse input values
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-s" || string(argv[i]) == "-sun" || string(argv[i]) == "-sf" || string(argv[i]) == "-sunfile" || string(argv[i]) == "--sunfile" || string(argv[i]) == "--sunephem" || string(argv[i]) == "--SunEphemeris") {
      if(i+1 < argc) {
	//There is still something to read;
	sunfile=argv[++i];
	i++;
      } else {
	cerr << "sunfile keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-e" || string(argv[i]) == "-earth" || string(argv[i]) == "-ef" || string(argv[i]) == "-earthfile" || string(argv[i]) == "--earthfile" || string(argv[i]) == "-earthephem" || string(argv[i]) == "--EarthEphemeris") {
      if(i+1 < argc) {
	//There is still something to read;
	earthfile=argv[++i];
	i++;
      } else {
	cerr << "earthfile keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-indets" || string(argv[i]) == "-incoords" || string(argv[i]) == "-infile" || string(argv[i]) == "-i" || string(argv[i]) == "-ic"  || string(argv[i]) == "-in" || string(argv[i]) == "-if" || string(argv[i]) == "-id" || string(argv[i]) == "--InputDetections") {
      if(i+1 < argc) {
	//There is still something to read;
	infile = argv[++i];
	i++;
      }
      else {
	cerr << "Input detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-of" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "-out" || string(argv[i]) == "-ofile" || string(argv[i]) == "-outf" || string(argv[i]) == "--OutputFile") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile = argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Unrecognized input keyword " << argv[i] << "\n";
      show_usage();
      return(1);
    }
  }

  // Echo input
  cout << "Input ephemeris files: " << sunfile << " and " << earthfile << "\n";
  cout << "Input detection file: " << infile << "\n";
  cout << "Output file name: " << outfile << "\n";
  
  // Read input detection file
  indets = {};
  ifstream instream1 {infile};
  while(instream1 >> MJD >> RA >> Dec) {
    det1 = det_bsc(MJD,RA,Dec);
    indets.push_back(det1);
  }
  instream1.close();
  detnum = indets.size();
  cout << detnum << " lines read from " << infile << "\n";

  // Read input ephemeris files for Earth and the sun
  Sunpos={};
  Sunvel={};
  mjd_ephem={};
  read_horizons_file(sunfile,mjd_ephem,Sunpos,Sunvel);
  Earthpos={};
  Earthvel={};
  mjd_test={};
  read_horizons_file(earthfile,mjd_test,Earthpos,Earthvel);
  if(mjd_test != mjd_ephem) {
    cerr << "ERROR: input ephemeris files for the sun and Earth don't match!\n";
    return(2);
  }

  // Open output file for writing
  ofstream outstream1 {outfile};

  // Calculate RA, Dec for the sun at every point
  for(detct=0 ; detct<detnum ; detct++) {
    // Copy input data to output file
    outstream1 << fixed << setprecision(6) << indets[detct].MJD << " " << indets[detct].RA << " " << indets[detct].Dec << " ";
    
    detmjd = indets[detct].MJD;
    // Barycentric coordinates for the sun
    planetpos01(detmjd, INTERP_POLYORDER, mjd_ephem, Sunpos, Sunposnow);
    // Barycentric coordinates for Earth
    planetpos01(detmjd, INTERP_POLYORDER, mjd_ephem, Earthpos, Earthposnow);
    // Vector from Earth to the sun
    Sunposnow.x -= Earthposnow.x;
    Sunposnow.y -= Earthposnow.y;
    Sunposnow.z -= Earthposnow.z;
    esdist = sqrt(Sunposnow.x*Sunposnow.x + Sunposnow.y*Sunposnow.y + Sunposnow.z*Sunposnow.z);
    Sunposnow.x /= esdist;
    Sunposnow.y /= esdist;
    Sunposnow.z /= esdist;
    stateunit_to_celestial(Sunposnow, sunRA, sunDec);
    // cout << "At MJD " << detmjd << ", the sun is at celestial coordinates " << sunRA << ", " << sunDec << "\n";
  
    // Transform input RA, Dec to ecliptic coordinates
    poleswitch02(indets[detct].RA, indets[detct].Dec, NEPRA, NEPDEC, 90.0l, EclipLon, EclipLat);
    // Transform solar RA, Dec to ecliptic coordinates
    poleswitch02(sunRA, sunDec, NEPRA, NEPDEC, 90.0l, SunEclipLon, SunEclipLat);

    // Difference the input ecliptic longitudes with that of the sun
    lambda = EclipLon - SunEclipLon;
    while(lambda<0l) lambda+=360.0l;
    while(lambda>=360.0l) lambda-=360.0l;
    phi = EclipLat;
    // Write sun-relative ecliptic coords to ouput file
    outstream1 << fixed << setprecision(6) << lambda << " " << phi << " ";
    
    // Carry out hammer projection centered on the sun
    if(lambda>180.0) lambda-=360.0;
    lambda /= -DEGPRAD; /*Negative sign -> East is left*/
    phi /= DEGPRAD;
    x = 2.0*pow(2.0,0.5)*cos(phi)*sin(lambda/2.0)/pow(1.0+cos(phi)*cos(lambda/2.0),0.5);
    y = pow(2.0,0.5)*sin(phi)/pow(1.0+cos(phi)*cos(lambda/2.0),0.5);
    // Write hammer-projected sun-centered coordinates to output file
    outstream1 << fixed << setprecision(6) << x << " " << y << " ";
 
    // Carry out hammer projection centered on opposition
    lambda = EclipLon - SunEclipLon - 180.0l;
    while(lambda<0l) lambda+=360.0l;
    while(lambda>=360.0l) lambda-=360.0l;
    phi = EclipLat;
    if(lambda>180.0) lambda-=360.0;
    lambda /= -DEGPRAD; /*Negative sign -> East is left*/
    phi /= DEGPRAD;
    x = 2.0*pow(2.0,0.5)*cos(phi)*sin(lambda/2.0)/pow(1.0+cos(phi)*cos(lambda/2.0),0.5);
    y = pow(2.0,0.5)*sin(phi)/pow(1.0+cos(phi)*cos(lambda/2.0),0.5);
    // Write hammer-projected opposition-centered coordinates to output file
    outstream1 << fixed << setprecision(6) << x << " " << y << "\n";
  }
  outstream1.close();
  return(0);
}
