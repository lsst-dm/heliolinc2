// November 09, 2021: readpairs01a.cpp
// Read pair files output by maketrack01a.cpp, and calculate
// the angular velocities, etc.

#include "std_lib_facilities.h"
#include "cmath"
#define NUMPOS 3
#define DEGPRAD ((double)180.0/M_PI) /*Degrees per radian*/
#define DSQUARE(x) double(x)*double(x)

#define MJDCOL 3
#define RACOL 6
#define DECCOL 8
#define SOLARDAY double(86400)
#define MINOBSINTERVAL 1.0 // Minimum time-between-images in seconds
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME 1.5 // Default max inter-image time interval
                    // for tracklets, in hours (will be converted
                    // to days before use).
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)

class det_bsc{ // Basic, minimal detection of some astronomical source 
public:
  double MJD;
  double RA;
  double Dec;
  det_bsc(double mjd, double ra, double dec) :MJD(mjd), RA(ra), Dec(dec) { }
};

class det_xindex{ // Detection of some astronomical source with
                  // cross-reference index.
public:
  double MJD;
  double RA;
  double Dec;
  long xindex;
  det_xindex(double mjd, double ra, double dec, long xindex) :MJD(mjd), RA(ra), Dec(dec), xindex(xindex) { }
};

class img_log01{ // Minimalist image log: MJD, RA, Dec, num_dets, imradius
public:
  double MJD;
  double RA;
  double Dec;
  int num_dets;
  double imradius;
  img_log01(double mjd, double ra, double dec, int ndet, double imrad) :MJD(mjd), RA(ra), Dec(dec), num_dets(ndet), imradius(imrad) { }
};

class img_log02{ // Image log that indexes a time-sorted detection vector:
                 // MJD, RA, Dec, starting index, ending index
public:
  double MJD;
  double RA;
  double Dec;
  long startind;
  long endind;
  img_log02(double mjd, double ra, double dec, long startind, long endind) :MJD(mjd), RA(ra), Dec(dec), startind(startind), endind(endind) { }
};

class early_det{
public:
  inline bool operator() (const det_bsc& o1, const det_bsc& o2) {
    return(o1.MJD < o2.MJD);
  }
};

class early_detindx{
public:
  inline bool operator() (const det_xindex& o1, const det_xindex& o2) {
    return(o1.MJD < o2.MJD);
  }
};

class early_imlg2{
public:
  inline bool operator() (const img_log02& i1, const img_log02& i2) {
    return(i1.MJD < i2.MJD);
  }
};

class point3d{ // Double-precision 3-D point
public:
  double x;
  double y;
  double z;
  point3d(double x, double y, double z) :x(x), y(y), z(z) { }
};

class longpair{ // Pair of long integers
public:
  long i1;
  long i2;
  longpair(long i1, long i2) :i1(i1), i2(i2) { }
};

// celeproj01: November 05, 2021
// Given double precision RA, Dec, project
// onto the unit sphere and return an object of class point3d.
// Input coordinates are in degrees, input RA=0, Dec=0
// projects to x=1,y=0,z=0; then y increases for positive
// RA.
point3d celeproj01(double RA, double Dec) {
  return( point3d( cos(RA/DEGPRAD)*cos(Dec/DEGPRAD) , sin(RA/DEGPRAD)*cos(Dec/DEGPRAD), sin(Dec/DEGPRAD)));
};

// celedeproj01: November 05, 2021
// Given a 3-d point (class point3d), de-project it back to
// celestial coordinates: i.e., reverse the process carried out
// by celeproj01.
int celedeproj01(point3d p3, double *RA, double *Dec)
{
  //Normalize the point
  double norm = sqrt(p3.x*p3.x + p3.y*p3.y + p3.z*p3.z);
  if(norm<=0.0) {
    *RA=0.0;
    *Dec=0.0;
    return(1);
  }
  double x = p3.x/norm;
  double y = p3.y/norm;
  double z = p3.z/norm;
  if(fabs(z)<=1.0) *Dec = asin(z)*DEGPRAD;
  else return(2);
  if(y==0 && x<0.0) {
    // y is zero and x is negative
    *RA = 180.0;
    return(0);
  }
  else if(y==0.0) {
    // y is zero and x is zero or positive
    *RA=0.0;
    return(0);
  }
  else if(y>0.0) {
    // y is strictly positive
    *RA = 90.0 - atan(x/y)*DEGPRAD;
    return(0);
  }
  else if(y<0.0) {
    // y is strictly negative
    *RA = 270.0 - atan(x/y)*DEGPRAD;
    return(0);
  }
  else {
    // Weird case, should be impossible
    return(3);
  }
};
    
  
// angspeed01: November 03, 2021
// given two observations of class Obs, which contains
// just the RA, Dec, and MJD of a detection of some astronomical
// object, calculate the angular velocity required to move from
// the first position to the second. Output units are deg/day.
double angspeed01(det_bsc o1, det_bsc o2)
{
  double x1,y1,z1,x2,y2,z2,h,d;
  double deltat = o2.MJD-o1.MJD;
  if(!isnormal(deltat)) throw(runtime_error("Bad time specifiers in angspeed01"));
  x1=cos(o1.Dec/DEGPRAD)*cos(o1.RA/DEGPRAD);
  y1=cos(o1.Dec/DEGPRAD)*sin(o1.RA/DEGPRAD);
  z1=sin(o1.Dec/DEGPRAD);
  x2=cos(o2.Dec/DEGPRAD)*cos(o2.RA/DEGPRAD);
  y2=cos(o2.Dec/DEGPRAD)*sin(o2.RA/DEGPRAD);
  z2=sin(o2.Dec/DEGPRAD);
  h=sqrt(DSQUARE(x1-x2)+DSQUARE(y1-y2)+DSQUARE(z1-z2));
  d=(double)2.0*asin(h/((double)2.0));
  return(d*DEGPRAD/deltat);
}

// distradec01: November 05, 2021
// Given two pairs of RA, Dec coordinates, calculate
// their angular separation on the sky in degrees.
double distradec01(double RA1, double Dec1, double RA2, double Dec2)
{
  double x1,y1,z1,x2,y2,z2,h;
  x1=cos(Dec1/DEGPRAD)*cos(RA1/DEGPRAD);
  y1=cos(Dec1/DEGPRAD)*sin(RA1/DEGPRAD);
  z1=sin(Dec1/DEGPRAD);
  x2=cos(Dec2/DEGPRAD)*cos(RA2/DEGPRAD);
  y2=cos(Dec2/DEGPRAD)*sin(RA2/DEGPRAD);
  z2=sin(Dec2/DEGPRAD);
  h=sqrt(DSQUARE(x1-x2)+DSQUARE(y1-y2)+DSQUARE(z1-z2));
  return(DEGPRAD*2.0*asin(h/2.0));
}


// distradec02: November 08, 2021: Like distradec01, but also
// calculates the celestial position angle from point 1 to point 2.
int distradec02(double ra1,double dec1,double ra2,double dec2,double *dist,double *pa)
{
  double x1,y1,z1,x2,y2,z2,h,d,poleangle,celpa,sinepa,cosinepa,colat1,colat2;

  // Calculate the distance d, in radians, between
  // the two points.
  x1=cos(dec1/DEGPRAD)*cos(ra1/DEGPRAD);
  y1=cos(dec1/DEGPRAD)*sin(ra1/DEGPRAD);
  z1=sin(dec1/DEGPRAD);
  x2=cos(dec2/DEGPRAD)*cos(ra2/DEGPRAD);
  y2=cos(dec2/DEGPRAD)*sin(ra2/DEGPRAD);
  z2=sin(dec2/DEGPRAD);
  h=sqrt(DSQUARE(x1-x2)+DSQUARE(y1-y2)+DSQUARE(z1-z2));
  d=2.0*asin(h/2.0);
  *dist = d*DEGPRAD;

  colat1 = M_PI/2.0 - dec1/DEGPRAD;
  colat2 = M_PI/2.0 - dec2/DEGPRAD;

  // Catch trivial cases
  if(d<=0.0 || d>=M_PI) {
    *pa = 0.0;
    return(0);
  }
  else if(ra1==ra2 && dec1>=dec2) {
    *pa = 180.0;
    return(0);
  }
  else if(ra1==ra2 && dec1<dec2) {
    *pa = 0.0;
    return(0);
  }
  else if(sin(colat1)<=0.0) {
    *pa = ra2;
    return(0);
  }

  // Calculate the difference in RA, paying careful
  // attention to wrapping.
  cosinepa = (cos(colat2) - cos(d)*cos(colat1))/(sin(d)*sin(colat1));
  if(ra1<ra2 && (ra2-ra1)<=180.0) {
    // Simple case, point 2 is east of point 1,
    // so PA should be less than 180 degrees.
    poleangle = (ra2-ra1)/DEGPRAD;
    sinepa = sin(colat2)*sin(poleangle)/sin(d);
    if(cosinepa>=0.0) celpa = asin(sinepa);
    else celpa = M_PI - asin(sinepa);
    *pa = celpa*DEGPRAD;
  }
  else if(ra1<ra2) {
    // Wrapped case with point 2 west of point 1
    // across zero RA: the PA should be greater
    // than 180 degrees.
    poleangle = (ra1+(double)360.0-ra2)/DEGPRAD;
    sinepa = sin(colat2)*sin(poleangle)/sin(d);
    if(cosinepa>=0.0) celpa = asin(sinepa);
    else celpa = M_PI - asin(sinepa);
    *pa = (double)360.0 - celpa*DEGPRAD;
  }
  else if(ra1>ra2 && (ra1-ra2)<=180.0) {
    // Simple case, point 2 is west of point 1,
    // so PA should be greater than 180 degrees.
    poleangle = (ra1-ra2)/DEGPRAD;
    sinepa = sin(colat2)*sin(poleangle)/sin(d);
    if(cosinepa>=0.0) celpa = asin(sinepa);
    else celpa = M_PI - asin(sinepa);
    *pa = (double)360.0 - celpa*DEGPRAD;
  }
  else if(ra1>ra2) {
    // Wrapped case with point 2 east of point 1
    // across zero RA: the PA should be less
    // than 180.0 degrees.
    poleangle = (ra2+(double)360.0-ra1)/DEGPRAD;
    sinepa = sin(colat2)*sin(poleangle)/sin(d);
    if(cosinepa>=0.0) celpa = asin(sinepa);
    else celpa = M_PI - asin(sinepa);
    *pa = celpa*DEGPRAD;
  }
  return(0);
}

static void show_usage()
{
  cerr << "Usage: readpairs01a -dets detfile -pairs pairfile -out outfile image file\n";
}
    
int main(int argc, char *argv[])
{
  det_bsc o1 = det_bsc(0,0,0);
  vector <det_bsc> detvec = {};
  longpair onepair = longpair(0,0);
  vector <longpair> pairvec ={};
  string indetfile;
  string inpairfile;
  string outfile;
  string lnfromfile;
  double MJD,RA, Dec;
  int reachedeof=0;

  if(argc<3)
    {
      show_usage();
      return(1);
    }
  
  int i=1;
  int j=0;
  int c='0';
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
    } else if(string(argv[i]) == "-p" || string(argv[i]) == "-pair" || string(argv[i]) == "-pairs" || string(argv[i]) == "--pairs" || string(argv[i]) == "--pair" || string(argv[i]) == "--pairfile" || string(argv[i]) == "--pairsfile") {
      if(i+1 < argc) {
	//There is still something to read;
	inpairfile=argv[++i];
	i++;
      }
      else {
	cerr << "Pair file keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outpair" || string(argv[i]) == "--outpairs") {
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
    }
  }
    
  cout << "input detection file " << indetfile << "\n";
  cout << "input pair file " << inpairfile << "\n";
  cout << "output file " << outfile << "\n";
  
  ifstream instream1 {indetfile};
  if(!instream1) error("can't open input file ",indetfile);
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
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
      while(i<lnfromfile.size() && c!=' ' && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=' ' && c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==1) MJD=stod(stest);
      else if(j==2) RA=stod(stest);
      else if(j==3) Dec=stod(stest);
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0) {
      // cout<<"MJD, RA, Dec: " << MJD-floor(MJD) << " " << RA << " " << Dec << "\n";
      o1=det_bsc(MJD,RA,Dec);
      detvec.push_back(o1);
    }
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";

  // Read input image pair file
  ifstream instream2 {inpairfile};
  if(!instream2) error("can't open input file ",inpairfile);
  reachedeof=0;
  while(reachedeof==0) {
    getline(instream2,lnfromfile);
    if(!instream2.eof() && !instream2.fail() && !instream2.bad()) ; // Read on.
    else if(instream2.eof()) reachedeof=1; //End of file, fine.
    else if(instream2.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream2.bad()) reachedeof=-2; //Worse problem, warn
    i=0;
    j = 0;
    c='0';
    long i1=0;
    long i2=0;
    while(i<lnfromfile.size() && reachedeof == 0) {
      string stest;
      c='0';
      while(i<lnfromfile.size() && c!=' ' && c!=',' && c!='\n' && c!=EOF) {
	// We allow the file to be delimited by comma or space.
	c=lnfromfile[i];
	if(c!=',' && c!=' ' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==1) i1=stol(stest);
      else if(j==2) i2=stol(stest);
    }
    if((reachedeof == 0 || reachedeof == 1) && i2>0) {
      onepair = longpair(i1,i2);
      pairvec.push_back(onepair);
    }
  }
  if(reachedeof==1) {
      cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  // Test: print out time-sorted detection table.
  ofstream outstream1 {outfile};
  double dist=0.0;
  double pa=0.0;
  double timediff=0.0;
  int i1=0; int i2=0;
  cout << "Read " << detvec.size() << " detections and " << pairvec.size() << "pairs.\n";

  for(i=0;i<pairvec.size();i++) {
    i1=pairvec[i].i1;
    i2=pairvec[i].i2;
    distradec02(detvec[i1].RA,detvec[i1].Dec,detvec[i2].RA,detvec[i2].Dec,&dist,&pa);
    timediff = detvec[i2].MJD - detvec[i1].MJD;
    outstream1 << fixed << setprecision(6) << detvec[i1].MJD << " " << detvec[i1].RA << " " << detvec[i1].Dec << " " << timediff << " " << dist << " " << dist/timediff*sin(pa/DEGPRAD) << " " << dist/timediff*cos(pa/DEGPRAD) << "\n";
    cout << fixed << setprecision(6) << detvec[i1].MJD << " " << detvec[i1].RA << " " << detvec[i1].Dec << " " << timediff << " " << dist << " " << dist/timediff*sin(pa/DEGPRAD) << " " << dist/timediff*cos(pa/DEGPRAD) << "\n";
  }
  
  return(0);
}
