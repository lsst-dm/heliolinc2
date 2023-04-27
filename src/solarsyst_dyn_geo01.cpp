// December 07, 2021: solarsyst_dyn_geo01.cpp
// Library of functions useful for Solar System dynamics and
// geometry. Includes functions originally developed in orbint02a.cpp,
// maketrack02b.cpp, projtest01b.cpp, and other places.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

// stringncopy01: March 09, 2022:
// like library function strncpy, but works the way I want it to.
// Copies first n characters of a string into a character array.
void stringncopy01(char *dest, const string &source, int n)
{
  int i=0;
  int nchar = source.size();
  if(nchar<n) {
    // We have enough space to copy the whole thing.
    for(i=0;i<nchar;i++) {
      dest[i] = source[i];
    }
    dest[nchar]='\0';
  } else {
    // Not enough space: copy first n-1 characters.
    for(i=0;i<n;i++) {
      dest[i] = source[i];
    }
    dest[n-1]='\0';
  }
}
 
// stringnmatch01: March 09, 2022:
// like library function strncmp, but works the way I want it to.
// Compares first n characters of two character arrays
int stringnmatch01(const char *string1, const char *string2, int n)
{
  int i=0;
  while(i<n && string1[i]!='\0' && string2[i]!='\0') {
    if(string1[i]<string2[i]) return(-1);
    else if(string1[i]>string2[i]) return(1);
    i++;
  }
  // If we get here without returning, the strings must have been equal
  return(0);
}

// char1000_getstring: February 20, 2023: given an object of class char1000_index,
// extract just character array and return it as a string
string char1000_getstring(const char1000_index &ci) {
  string s1;
  int i=0;
  while(i<1000 && ci.selem[i]!=0) {
    s1.push_back(ci.selem[i]);
    i++;
  }
  return(s1);
}

// char500_getstring: February 20, 2023: given an object of class char500_index,
// extract just character array and return it as a string
string char500_getstring(const char500_index &ci) {
  string s1;
  int i=0;
  while(i<500 && ci.selem[i]!=0) {
    s1.push_back(ci.selem[i]);
    i++;
  }
  return(s1);
}

// char256_getstring: February 20, 2023: given an object of class char256_index,
// extract just character array and return it as a string
string char256_getstring(const char256_index &ci) {
  string s1;
  int i=0;
  while(i<256 && ci.selem[i]!=0) {
    s1.push_back(ci.selem[i]);
    i++;
  }
  return(s1);
}

// char128_getstring: February 20, 2023: given an object of class char128_index,
// extract just character array and return it as a string
string char128_getstring(const char128_index &ci) {
  string s1;
  int i=0;
  while(i<128 && ci.selem[i]!=0) {
    s1.push_back(ci.selem[i]);
    i++;
  }
  return(s1);
}

// char64_getstring: February 20, 2023: given an object of class char64_index,
// extract just character array and return it as a string
string char64_getstring(const char64_index &ci) {
  string s1;
  int i=0;
  while(i<64 && ci.selem[i]!=0) {
    s1.push_back(ci.selem[i]);
    i++;
  }
  return(s1);
}

// January 20, 2023: changed all the vector and matrix
// allocations to use <nx, <ny rather than <=nx, <=ny.
// The old way made them safe against using 1-indexed
// rather than zero-indexed loops, but it caused unexpected
// blank values that caused trouble in sorting. From now
// on, they can only be used in zero-indexed loops.
void make_ivec(long nx, vector <int> &ivec)
{
  long i=0;
  ivec={};
  for(i=0;i<nx;i++) ivec.push_back(0);
}

void make_imat(int nx, int ny, vector <vector <int>> &imat)
{
  int i=0;
  int j=0;
  vector <int> tvec;
  imat = {};
  
  for(i=0;i<nx;i++) {
    tvec={};
    for(j=0;j<ny;j++) tvec.push_back(0);
    imat.push_back(tvec);
  }
}

void make_lvec(int nx, vector <long> &lvec)
{
  int i=0;
  lvec={};
  for(i=0;i<nx;i++) lvec.push_back(0);
}

void make_lmat(int nx, int ny, vector <vector <long>> &lmat)
{
  int i=0;
  int j=0;
  vector <long> tvec;
  lmat = {};
  
  for(i=0;i<nx;i++) {
    tvec={};
    for(j=0;j<ny;j++) tvec.push_back(0);
    lmat.push_back(tvec);
  }
}

void make_cvec(int nx, vector <char> &cvec)
{
  int i=0;
  cvec={};
  for(i=0;i<nx;i++) cvec.push_back('\0');
}

void make_cmat(int nx, int ny, vector <vector <char>> &cmat)
{
  int i=0;
  int j=0;
  vector <char> tvec;
  cmat = {};
  
  for(i=0;i<nx;i++) {
    tvec={};
    for(j=0;j<ny;j++) tvec.push_back('\0');
    cmat.push_back(tvec);
  }
}

void make_dvec(int nx, vector <double> &dvec)
{
  int i=0;
  dvec={};
  for(i=0;i<nx;i++) dvec.push_back(0.0);
}

void make_dmat(int nx, int ny, vector <vector <double>> &dmat)
{
  int i=0;
  int j=0;
  vector <double> tvec;
  dmat = {};
  
  for(i=0;i<nx;i++) {
    tvec={};
    for(j=0;j<ny;j++) tvec.push_back(0.0);
    dmat.push_back(tvec);
  }
}

void make_LDvec(int nx, vector <long double> &ldvec)
{
  int i=0;
  ldvec={};
  for(i=0;i<nx;i++) ldvec.push_back(0.0);
}

void make_LDmat(int nx, int ny, vector <vector <long double>> &ldmat)
{
  int i=0;
  int j=0;
  vector <long double> tvec;
  ldmat = {};
  
  for(i=0;i<nx;i++) {
    tvec={};
    for(j=0;j<ny;j++) tvec.push_back(0.0);
    ldmat.push_back(tvec);
  }
}

double dotprod3d(point3d p1, point3d p2)
{
  return(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z);
}

long double dotprod3LD(point3LD p1, point3LD p2)
{
  return(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z);
}

double vecabs3d(point3d p1)
{
  return(sqrt(dotprod3d(p1,p1)));
}

long double vecabs3LD(point3LD p1)
{
  return(sqrt(dotprod3LD(p1,p1)));
}

int vecnorm3d(point3d &p1)
{
  double norm = vecabs3d(p1);
  if(isnormal(norm)) {
    p1.x /= norm;
    p1.y /= norm;
    p1.z /= norm;
    return(0);
  } else return(1);
}

int vecnorm3LD(point3LD &p1)
{
  long double norm = vecabs3LD(p1);
  if(isnormal(norm)) {
    p1.x /= norm;
    p1.y /= norm;
    p1.z /= norm;
    return(0);
  } else return(1);
}

point3d crossprod3d(point3d p1, point3d p2)
{
  point3d pout = point3d(0,0,0);
  pout.x = p1.y*p2.z - p1.z*p2.y;
  pout.y = p1.z*p2.x - p1.x*p2.z;
  pout.z = p1.x*p2.y - p1.y*p2.x;
  return(pout);
}

point3LD crossprod3LD(point3LD p1, point3LD p2)
{
  point3LD pout = point3LD(0L,0L,0L);
  pout.x = p1.y*p2.z - p1.z*p2.y;
  pout.y = p1.z*p2.x - p1.x*p2.z;
  pout.z = p1.x*p2.y - p1.y*p2.x;
  return(pout);
}

long double intpowLD(long double x, int p)
{
  int i=0;
  long double y=1.0L;
  for(i=0;i<p;i++) y*=x;
  return(y);
}

long double factorialLD(int p)
{
  int i=0;
  long double y=1.0L;
  for(i=1;i<=p;i++) y*=(long double)i;
  return(y);
}

double intpowD(double x, int p)
{
  int i=0;
  double y=1.0l;
  for(i=0;i<p;i++) y*=x;
  return(y);
}

double factorialD(int p)
{
  int i=0;
  double y=1.0l;
  for(i=1;i<=p;i++) y*=(double)i;
  return(y);
}

// dmean01: January 16, 2023
// Calculate and return the mean of a double-precision vector
double dmean01(const vector <double> &invec)
{
  int pnum=invec.size();
  if(pnum<=0) return(0.0l);
  double mean=0.0l;
  for(int pct=0; pct<pnum; pct++) {
    mean+=invec[pct];
  }
  return(mean/double(pnum));
}

// drms01: January 16, 2023
// Calculate and return the RMS deviation of a
// double-precision vector from its mean
double drms01(const vector <double> &invec)
{
  int pnum=invec.size();
  int pct=0;
  double mean,rms;
  if(pnum<=1) return(-1.0l);

  mean=rms=0.0l;
  // Calculate the mean
  for(pct=0; pct<pnum; pct++) {
    mean+=invec[pct];
  }
  mean/=double(pnum);
  // Calculate the RMS
  for(pct=0; pct<pnum; pct++) {
    rms += DSQUARE(invec[pct]-mean);
  }
  rms /= double(pnum-1);
  return(sqrt(rms));
}

// dmeanrms01: January 16, 2023
// Calculate and return the mean and RMS of
// a double-precision vector.
int dmeanrms01(const vector <double> &invec, double *mean, double *rms)
{
  int pnum=invec.size();
  int pct=0;
  
  if(pnum<=0) return(2);
  *mean = *rms = 0.0l;
  // Calculate the mean
  for(pct=0; pct<pnum; pct++) {
    *mean+=invec[pct];
  }
  *mean /= double(pnum);
  if(pnum==1) return(1);
  
  // Calculate the RMS
  for(pct=0; pct<pnum; pct++) {
    (*rms) += DSQUARE(invec[pct] - *mean);
  }
  *rms = sqrt((*rms)/double(pnum-1));
  return(0);
}


// celeproj01: November 05, 2021
// Given double precision RA, Dec in DEGREES, project
// onto the unit sphere and return an object of class point3d.
// Input coordinates are in degrees, input RA=0, Dec=0
// projects to x=1,y=0,z=0; then y increases for positive
// RA.
point3d celeproj01(double RA, double Dec) {
  return( point3d( cos(RA/DEGPRAD)*cos(Dec/DEGPRAD) , sin(RA/DEGPRAD)*cos(Dec/DEGPRAD), sin(Dec/DEGPRAD)));
};

// celedeproj01: November 05, 2021
// Given a 3-d point (class point3d), de-project it back to
// celestial coordinates IN DEGREES: i.e., reverse the process
// carried out by celeproj01.
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
    
// celeproj01LD: November 05, 2021
// Given double precision RA, Dec in DEGREES, project
// onto the unit sphere and return an object of class point3d.
// Input coordinates are in degrees, input RA=0, Dec=0
// projects to x=1,y=0,z=0; then y increases for positive
// RA.
point3LD celeproj01LD(long double RA, long double Dec) {
  return( point3LD( cos(RA/DEGPRAD)*cos(Dec/DEGPRAD) , sin(RA/DEGPRAD)*cos(Dec/DEGPRAD), sin(Dec/DEGPRAD)));
};

// celedeproj01LD: November 05, 2021
// Given a 3-d point (class point3LD), de-project it back to
// celestial coordinates IN DEGREES: i.e., reverse the process
// carried out by celeproj01.
int celedeproj01LD(point3LD p3, long double *RA, long double *Dec)
{
  //Normalize the point
  long double norm = sqrt(p3.x*p3.x + p3.y*p3.y + p3.z*p3.z);
  if(norm<=0.0) {
    *RA=0.0;
    *Dec=0.0;
    return(1);
  }
  long double x = p3.x/norm;
  long double y = p3.y/norm;
  long double z = p3.z/norm;
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
  double arcsinepa=0.0l;

  // Handle trivial cases.
  if(fabs(ra1-ra2)/DEGPRAD < VSMALLANG) {
    // the two RA values are functionally identical
    if(fabs(dec1-dec2)/DEGPRAD < VSMALLANG) {
      // the two Dec values are functionally identical
      *dist=0.0l;
      *pa=0.0l; // Dummy value for zero distance.
      return(0);
    } else if(dec2<dec1) {
      *dist = dec1-dec2;
      *pa = 180.0l; // Due South
      return(0);
    } else {
      // dec2>dec1 by logical elimination
      *dist = dec2-dec1;
      *pa = 0.0l; // Due North
      return(0);
    }
  } else if(fabs(dec1-dec2)/DEGPRAD < VSMALLANG) {
    // the two Dec values are functionally identical,
    // although the two RA values are not.
    if(ra2<ra1) {
      *dist = ra1-ra2;
      *pa = 270.0l; // Due West.
      return(0);
    } else {
      // ra2>ra1 by logical elimination
      *dist = ra2-ra1;
      *pa = 90.0l; // Due East.
      return(0);
    }
  } else {
    // We have a non-trivial case
    // Calculate the distance d, in radians, between
    // the two points.
    x1=cos(dec1/DEGPRAD)*cos(ra1/DEGPRAD);
    y1=cos(dec1/DEGPRAD)*sin(ra1/DEGPRAD);
    z1=sin(dec1/DEGPRAD);
    x2=cos(dec2/DEGPRAD)*cos(ra2/DEGPRAD);
    y2=cos(dec2/DEGPRAD)*sin(ra2/DEGPRAD);
    z2=sin(dec2/DEGPRAD);
    h=sqrt(DSQUARE(x1-x2)+DSQUARE(y1-y2)+DSQUARE(z1-z2));
    if(h/2.0l <= 1.0l) {
      d=2.0*asin(h/2.0l);
    }
    else {
      if(WARN_INVERSE_TRIG>0) cerr << "WARNING: distradec02 attempting to take arcsine of 1 + " << h/2.0l - 1.0l << "\n";
      d = M_PI/2.0l;
    }
    *dist = d*DEGPRAD;

    colat1 = M_PI/2.0 - dec1/DEGPRAD;
    colat2 = M_PI/2.0 - dec2/DEGPRAD;

    // Calculate the difference in RA, paying careful
    // attention to wrapping.
    cosinepa = (cos(colat2) - cos(d)*cos(colat1))/(sin(d)*sin(colat1));
    if(ra1<ra2 && (ra2-ra1)<=180.0) {
      // Simple case, point 2 is east of point 1,
      // so PA should be less than 180 degrees.
      poleangle = (ra2-ra1)/DEGPRAD;
      sinepa = sin(colat2)*sin(poleangle)/sin(d);
      if(sinepa<=1.0l) {
	arcsinepa = asin(sinepa);
      } else {
	if(WARN_INVERSE_TRIG>0) cerr << "WARNING: distradec02 attempting to take the arcsine of 1 + " << sinepa-1.0l << "\n";
	arcsinepa = M_PI/2.0l;
      }
      if(cosinepa>=0.0) celpa = arcsinepa;
      else celpa = M_PI - arcsinepa;
      *pa = celpa*DEGPRAD;
    }
    else if(ra1<ra2) {
      // Wrapped case with point 2 west of point 1
      // across zero RA: the PA should be greater
      // than 180 degrees.
      poleangle = (ra1+(double)360.0-ra2)/DEGPRAD;
      sinepa = sin(colat2)*sin(poleangle)/sin(d);
      if(sinepa<=1.0l) {
	arcsinepa = asin(sinepa);
      } else {
	if(WARN_INVERSE_TRIG>0) cerr << "WARNING: distradec02 attempting to take the arcsine of 1 + " << sinepa-1.0l << "\n";
	arcsinepa = M_PI/2.0l;
      }
      if(cosinepa>=0.0) celpa = arcsinepa;
      else celpa = M_PI - arcsinepa;
      *pa = (double)360.0 - celpa*DEGPRAD;
    }
    else if(ra1>ra2 && (ra1-ra2)<=180.0) {
      // Simple case, point 2 is west of point 1,
      // so PA should be greater than 180 degrees.
      poleangle = (ra1-ra2)/DEGPRAD;
      sinepa = sin(colat2)*sin(poleangle)/sin(d);
      if(sinepa<=1.0l) {
	arcsinepa = asin(sinepa);
      } else {
	if(WARN_INVERSE_TRIG>0) cerr << "WARNING: distradec02 attempting to take the arcsine of 1 + " << sinepa-1.0l << "\n";
	arcsinepa = M_PI/2.0l;
      }
      if(cosinepa>=0.0) celpa = arcsinepa;
      else celpa = M_PI - arcsinepa;
      *pa = (double)360.0 - celpa*DEGPRAD;
    }
    else if(ra1>ra2) {
      // Wrapped case with point 2 east of point 1
      // across zero RA: the PA should be less
      // than 180.0 degrees.
      poleangle = (ra2+(double)360.0-ra1)/DEGPRAD;
      sinepa = sin(colat2)*sin(poleangle)/sin(d);
      if(sinepa<=1.0l) {
	arcsinepa = asin(sinepa);
      } else {
	if(WARN_INVERSE_TRIG>0) cerr << "WARNING: distradec02 attempting to take the arcsine of 1 + " << sinepa-1.0l << "\n";
	arcsinepa = M_PI/2.0l;
      }
      if(cosinepa>=0.0) celpa = arcsinepa;
      else celpa = M_PI - arcsinepa;
      *pa = celpa*DEGPRAD;
    }
    return(0);
  }
  return(0);
}

long medindex(const vector <xy_index> &xyvec, int dim)
{
  vector <xy_index> xyv = xyvec; //Mutable copy of immutable input vector
  for(long unsigned int i=0; i<xyv.size(); i++) xyv[i].index=i; //Redefine indices
  long medpt = xyv.size()/2;
  if(dim%2==1) sort(xyv.begin(), xyv.end(), xyind_lower_x());
  else sort(xyv.begin(), xyv.end(), xyind_lower_y());
  return(xyv[medpt].index);
}

int splitxy(const vector <xy_index> &xyvec, int dim, long unsigned int splitpoint, vector <xy_index> &left, vector <xy_index> &right)
{
  long unsigned int i=0;
  double xval = xyvec[splitpoint].x;
  double yval = xyvec[splitpoint].y;
  
  if(dim%2==1) {
    // Split on x
    for(i=0 ; i<xyvec.size() ; i++) {
      if(i!=splitpoint && xyvec[i].x<=xval) {
	left.push_back(xyvec[i]);
      } else if(i!=splitpoint) {
	right.push_back(xyvec[i]);
      }
    }
  } else {
    // Split on y
    for(i=0 ; i<xyvec.size() ; i++) {
      if(i!=splitpoint && xyvec[i].y<=yval) {
	left.push_back(xyvec[i]);
      } else if(i!=splitpoint) right.push_back(xyvec[i]);
    }
  }
  return(0);
}

// kdtree01: November 11, 2021:
// Given an input root point, presumed to have null
// right and left branches, load the branches and then
// call kdtree01 on them recursively.
// NOTE THAT THIS IS FOR 2-D KD trees.
int kdtree01(const vector <xy_index> &xyvec, int dim, long unsigned int rootptxy, long unsigned int rootptkd, vector <kdpoint> &kdvec)
{
  int lmed=0;
  int rmed=0;
  int kdct = kdvec.size()-1;
  long leftrootkd=-1;
  long rightrootkd=-1;
  xy_index xyi = xy_index(0.0,0.0,0);
  kdpoint lp = kdpoint(xyi,-1,-1,0);
  kdpoint rp = kdpoint(xyi,-1,-1,0);
  kdpoint kdtest = kdpoint(xyi,-1,-1,0);
  vector <xy_index> leftvec = {};
  vector <xy_index> rightvec = {};

  splitxy(xyvec,dim,rootptxy,leftvec,rightvec);
  
  if(dim==1) dim=2;
  else dim=1;

  if(leftvec.size()==1) {
    // Left branch is just a single leaf
    lp = kdpoint(leftvec[0],-1,-1,dim);
    kdvec.push_back(lp);
    kdct++;
    kdvec[rootptkd].left = kdct;
    kdtest = kdvec[kdct];
  } else if(leftvec.size()<=0) {
    // There is no left branch
    kdvec[rootptkd].left = -1;
  }
  if(rightvec.size()==1) {
    // Right branch is just a single leaf
    rp = kdpoint(rightvec[0],-1,-1,dim);
    kdvec.push_back(rp);
    kdct++;
    kdvec[rootptkd].right = kdct;
    kdtest = kdvec[kdct];
  } else if(rightvec.size()<=0) {
    // There is no right branch
    kdvec[rootptkd].right = -1;
  }
   
 if(leftvec.size()>1) {
    lmed = medindex(leftvec,dim);
    lp = kdpoint(leftvec[lmed],-1,-1,dim);
    kdvec.push_back(lp);
    kdct++;
    kdvec[rootptkd].left = kdct;
    leftrootkd = kdct;
    kdtest = kdvec[kdct];
 }
 
  if(rightvec.size()>1) {
    rmed = medindex(rightvec,dim);
    rp = kdpoint(rightvec[rmed],-1,-1,dim);
    kdvec.push_back(rp);
    kdct++;
    kdvec[rootptkd].right = kdct;
    rightrootkd = kdct;
    kdtest = kdvec[kdct];
  }
  // I moved these down out of the above loops, because I thought
  // that otherwise, a bunch of stuff might get pushed down by the
  // left loop that the right loop didn't know about.
  if(leftvec.size()>1 && leftrootkd>=0) kdtree01(leftvec,dim,lmed,leftrootkd,kdvec);
  else if(leftvec.size()>1 && leftrootkd<0)
    {
      cerr << "Error, kdtree01 finds leftroot less than zero with leftvec.size() = " << leftvec.size() << "\n";
    }
  if(rightvec.size()>1 && rightrootkd>=0) kdtree01(rightvec,dim,rmed,rightrootkd,kdvec);
  else if(rightvec.size()>1 && rightrootkd<0)
    {
      cerr << "Error, kdtree01 finds rightroot less than zero with rightvec.size() = " << rightvec.size() << "\n";
    }

  return(0);
}

// kdrange01: November 15, 2021:
// Given a k-d tree vector kdvec created by kdtree01,
// perform a range-query about the point x,y. Returns
// a vector indexing all of the points in the input k-d tree
// that lie within the specified range of the input coordinates.
// Assumes that kdvec[0] is the root of the k-d tree.
// NOTE THAT THIS IS FOR 2-D KD trees.
int kdrange01(const vector <kdpoint> &kdvec,double x,double y,double range,vector <long> &indexvec)
{
  double rng2 = range*range;
  int notdone=1;
  int dim=1;
  int currentpoint=0;
  int leftpoint=0;
  int rightpoint=0;
  int goleft=0;
  int goright=0;
  double xdiff=0.0;
  double ydiff=0.0;
  vector <long> checkit={};
  long unsigned int checknum=0;

  while(notdone>0) {
    // Climb to the top of the k-d tree, keeping track
    // of potentially interesting unexplored branches
    // in the vector checkit.
    while(leftpoint>=0 || rightpoint>=0) {
      // Previous step did not end on a leaf.
      leftpoint = kdvec[currentpoint].left;
      rightpoint = kdvec[currentpoint].right;
      dim = kdvec[currentpoint].dim;
      if(dim%2==1) {
	xdiff = kdvec[currentpoint].point.x - x;
	// cout << "kxdr: " << kdvec[currentpoint].point.x << " " << x << " " << xdiff << " " << range << "\n";
	goright = (xdiff <= range); // possible hits lie to the left;
	goleft = (xdiff >= -range); // possible hits lie to the right;
	if(goleft && goright) {
	  // Current point might be within range.
	    ydiff = kdvec[currentpoint].point.y - y;
	    if(fabs(ydiff)<=range && (xdiff*xdiff + ydiff*ydiff)<=rng2) {
	      // Current point is within range. Add it to the output vector
	      indexvec.push_back(currentpoint);
	    }
	    if(leftpoint>=0) {
	       //Explore leftward first.
	      currentpoint = leftpoint;
	      if(rightpoint>=0) {
		// Rightward branch will also be explored later
		checknum++;
		if(checknum>checkit.size()) {
		  checkit.push_back(rightpoint);
		}
		else {
		  checkit[checknum-1] = rightpoint;
		}		
	      }
	    }
	    else if(rightpoint>=0) {
	      // Leftward branch is a dead end: explore rightward branch
	      currentpoint = rightpoint;
	    }
	}
	else if(goleft) {
	  // Current point cannot be in range, but points that
	  // are in range may lie along the left branch.
	  if(leftpoint>=0) {
	    currentpoint = leftpoint;
	  } else rightpoint=-1; // Dead end, make sure while-loop exits.
	} else if(goright) {
	  // Current point cannot be in range, but points that
	  // are in range may lie along the right branch.
	  if(rightpoint>=0) {
	    currentpoint = rightpoint;
	  } else leftpoint=-1;  // Dead end, make sure while-loop exits.
	} else {
	  // Program concluded it should go neither left nor right.
	  // The likely cause is that it encountered a NAN. Give up on this point.
	  leftpoint=rightpoint=-1;
	  cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
	  cerr << "Input point " << x << ", " << y <<", target point " << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ".\n";
	}
	// Close x-dimension case
      } else if(dim%2==0) {
	ydiff = kdvec[currentpoint].point.y - y;
	goright = (ydiff <= range); // possible hits lie to the left;
	goleft = (ydiff >= -range); // possible hits lie to the right;
	if(goleft && goright) {
	    // Current point might be within range.
	    xdiff = kdvec[currentpoint].point.x - x;
	    if(fabs(ydiff)<=range && (xdiff*xdiff + ydiff*ydiff)<=rng2) {
	      // Current point is within range. Add it to the output vector
	      indexvec.push_back(currentpoint);
	    }
	    if(leftpoint>=0) {
	       //Explore leftward first.
	      currentpoint = leftpoint;
	      if(rightpoint>=0) {
		// Rightward branch will also be explored later
		checknum++;
		if(checknum>checkit.size()) {
		  checkit.push_back(rightpoint);
		}
		else {
		  checkit[checknum-1] = rightpoint;
		}
	      }
	    } else if(rightpoint>=0) {
	      // Leftward branch is a dead end: explore rightward branch
	      currentpoint = rightpoint;
	    }
	}
	else if(goleft) {
	  // Current point cannot be in range, but points that
	  // are in range may lie along the left branch.
	  if(leftpoint>=0) {
	    currentpoint = leftpoint;
	  } else rightpoint = -1; // Dead end, make sure while loop exits.
	} else if(goright) {
	  // Current point cannot be in range, but points that
	  // are in range may lie along the right branch.
	  if(rightpoint>=0) {
	    currentpoint = rightpoint;
	  } else leftpoint=-1;  // Dead end, make sure while loop exits.
	} else {
	  // Program concluded it should go neither left nor right.
	  // The likely cause is that it encountered a NAN. Give up on this point.
	  leftpoint=rightpoint=-1;
	  cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
	  cerr << "Input point " << x << ", " << y <<", target point " << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ".\n";
	}
	// Note that we do not need to worry about the possiblity
	// that current point will get set to -1: i.e., we were
	// at a leaf or a one-sided branch. Such cases will
	// be caught by the while statement.
	// Close y-dimension case
      }
      // Close while-loop checking if we've hit a leaf.
    }
    // We have climbed up the tree to a leaf. Go backwards through
    // the checkit vector and see if there is anything to check.
    checknum=checkit.size();
    while(checknum>=1 && checkit[checknum-1]<0) checknum--;
    if(checknum<=0) {
      //There were no valid entries to check: we're done.
      notdone=0;
    } else {
      //Set currentpoint to the last valid entry in checkit
      currentpoint = checkit[checknum-1];
      //Mark this point as used.
      checkit[checknum-1]=-1;
      leftpoint=rightpoint=0;
    }
  }
  return(0);
}

long medind_6LDx2(const vector <point6LDx2> &pointvec, int dim)
{
  vector <point6LDx2> pvec = pointvec; //Mutable copy of immutable input vector
  for(long unsigned int i=0; i<pvec.size(); i++) pvec[i].i1=i; //Redefine indices
  long medpt = pvec.size()/2; // Central point of vector (it will be off by one half
                              // for a vector with even length, but we don't care).
  if(dim%6 == 1) sort(pvec.begin(), pvec.end(), lower_point6LDx2_x()); // Sort vector by x
  else if(dim%6 == 2) sort(pvec.begin(), pvec.end(), lower_point6LDx2_y()); // Sort vector by y
  else if(dim%6 == 3) sort(pvec.begin(), pvec.end(), lower_point6LDx2_z()); // Sort vector by z
  else if(dim%6 == 4) sort(pvec.begin(), pvec.end(), lower_point6LDx2_vx()); // Sort vector by vx
  else if(dim%6 == 5) sort(pvec.begin(), pvec.end(), lower_point6LDx2_vy()); // Sort vector by vy
  else if(dim%6 == 0) sort(pvec.begin(), pvec.end(), lower_point6LDx2_vz()); // Sort vector by vz
  else {
    cerr << "ERROR: medind_6LDx2 recieved invalid dimension " << dim << "\n";
    return(-1);
  }
  return(pvec[medpt].i1); // Output the index of the median point in
                             // the original, unsorted input vector.
}

// splitLDx2: January 05, 2022:
// Given a vector of type point6LDx2, split it into two halves,
// a left half with all the points lower than or equal to a specified
// split point along the chosen dimension (use dim = 1, 2, 3, 4, 5, or 6
// to split along x, y, z, vx, vy, or vz respectively).
int splitLDx2(const vector <point6LDx2> &pointvec, int dim, long unsigned int splitpoint, vector <point6LDx2> &left, vector <point6LDx2> &right)
{
  long unsigned int i=0;
  long double splitval = 0.0L;

  if(dim%6==1) {
    // split on x
    splitval = pointvec[splitpoint].x;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].x<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==2) {
    // split on y
    splitval = pointvec[splitpoint].y;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].y<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==3) {
    // split on z
    splitval = pointvec[splitpoint].z;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].z<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==4) {
    // split on vx
    splitval = pointvec[splitpoint].vx;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].vx<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==5) {
    // split on vy
    splitval = pointvec[splitpoint].vy;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].vy<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==0) {
    // split on vz
    splitval = pointvec[splitpoint].vz;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].vz<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else {
      cerr << "ERROR: splitLDx2 asked to split on undefined dimension " << dim << "\n";
      return(1);
  } 
  return(0);
}

// kdtree_6D01: January 05, 2022
// Given an input root point, presumed to have null
// right and left branches, load the branches and then
// call kdtree_6D01 on them recursively.
int kdtree_6D01(const vector <point6LDx2> &invec, int dim, long unsigned int splitpoint, long unsigned int kdroot, vector <KD_point6LDx2> &kdvec)
{
  int lmed=0;
  int rmed=0;
  int kdct = kdvec.size()-1;
  long leftrootkd=-1;
  long rightrootkd=-1;
  point6LDx2 point0 = point6LDx2(0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0,0);
  KD_point6LDx2 lp = KD_point6LDx2(point0,-1,-1,0,0);
  KD_point6LDx2 rp = KD_point6LDx2(point0,-1,-1,0,0);
  vector <point6LDx2> leftvec = {};
  vector <point6LDx2> rightvec = {};

  // Basic outline: split the input vector into a left and a right
  // half, where the left half is below (or level with) splitpoint
  // in the dimension specified by dim, and the right half is above
  // splitpoint. Find the median of the left and right vectors,
  // and make the left and right branches from kdroot point to
  // these medians. Then call kdtree_6D01 itself recursively on
  // each of these median points, to peform a new split along a
  // different dimension.
  splitLDx2(invec,dim,splitpoint,leftvec,rightvec);

  dim+=1;
  while(dim>6) dim-=6;

  if(leftvec.size()==1) {
    // Left branch is just a single leaf
    lp = KD_point6LDx2(leftvec[0],-1,-1,dim,0); // Define new point as a leaf: branches point nowhere
    kdvec.push_back(lp); // Add this new point to the KD tree.
    kdct++; // Keep track of how many point are in the tree
    kdvec[kdroot].left = kdct; // Stick the new point on the left branch of the input root.
  } else if(leftvec.size()<=0) {
    // There is no left branch
    kdvec[kdroot].left = -1;
  }
  if(rightvec.size()==1) {
    // Right branch is just a single leaf
    rp = KD_point6LDx2(rightvec[0],-1,-1,dim,0);
    kdvec.push_back(rp);
    kdct++;
    kdvec[kdroot].right = kdct;
  } else if(rightvec.size()<=0) {
    // There is no right branch
    kdvec[kdroot].right = -1;
  }
   
 if(leftvec.size()>1) {
    lmed = medind_6LDx2(leftvec,dim);
    lp = KD_point6LDx2(leftvec[lmed],-1,-1,dim,0);
    kdvec.push_back(lp);
    kdct++;
    kdvec[kdroot].left = kdct;
    leftrootkd = kdct;
 }
 
  if(rightvec.size()>1) {
    rmed = medind_6LDx2(rightvec,dim);
    rp = KD_point6LDx2(rightvec[rmed],-1,-1,dim,0);
    kdvec.push_back(rp);
    kdct++;
    kdvec[kdroot].right = kdct;
    rightrootkd = kdct;
  }
  // I moved these down out of the above loops, because I thought
  // that otherwise, a bunch of stuff might get pushed down by the
  // left loop that the right loop didn't know about.
  if(leftvec.size()>1 && leftrootkd>=0) kdtree_6D01(leftvec,dim,lmed,leftrootkd,kdvec);
  else if(leftvec.size()>1 && leftrootkd<0)
    {
      cerr << "Error, kdtree_6D01 finds leftroot less than zero with leftvec.size() = " << leftvec.size() << "\n";
    }
  if(rightvec.size()>1 && rightrootkd>=0) kdtree_6D01(rightvec,dim,rmed,rightrootkd,kdvec);
  else if(rightvec.size()>1 && rightrootkd<0)
    {
      cerr << "Error, kdtree_6D01 finds rightroot less than zero with rightvec.size() = " << rightvec.size() << "\n";
    }

  return(0);
}

// point6LDx2_dist: January 05, 2022:
// Calculate the distance in 6-dimensional parameter space bewteen
// two points of class point6LDx2.
long double point6LDx2_dist(const point6LDx2 &p1, const point6LDx2 &p2)
{
  return(sqrt(LDSQUARE(p1.x - p2.x) + LDSQUARE(p1.y - p2.y) + LDSQUARE(p1.z - p2.z) + LDSQUARE(p1.vx - p2.vx) + LDSQUARE(p1.vy - p2.vy) + LDSQUARE(p1.vz - p2.vz)));
}	 

// point6LDx2_dist2: January 05, 2022:
// Calculate the squared distance in 6-dimensional parameter space
// between two points of class point6LDx2.
long double point6LDx2_dist2(const point6LDx2 &p1, const point6LDx2 &p2)
{
  return(LDSQUARE(p1.x - p2.x) + LDSQUARE(p1.y - p2.y) + LDSQUARE(p1.z - p2.z) + LDSQUARE(p1.vx - p2.vx) + LDSQUARE(p1.vy - p2.vy) + LDSQUARE(p1.vz - p2.vz));
}	 

// kdrange_6D01: January 05, 2022:
// Given a k-d tree vector kdvec created by kdtree_6D01,
// perform a range-query about the specified point. Returns
// a vector indexing all of the points in the input k-d tree
// that lie within the specified range of the input coordinates.
// Assumes that kdvec[0] is the root of the k-d tree.
int kdrange_6D01(const vector <KD_point6LDx2> &kdvec, const point6LDx2 &querypoint, long double range, vector <long> &indexvec)
{
  long double rng2 = range*range;
  int notdone=1;
  int dim=1;
  int currentpoint=0;
  int leftpoint=0;
  int rightpoint=0;
  int goleft=0;
  int goright=0;
  long double pointdiff = 0.0L;
  long double pdist2 = 0.0L;
  vector <long> checkit={};
  long unsigned int checknum=0;

  while(notdone>0) {
    // Climb to the top of the k-d tree, keeping track
    // of potentially interesting unexplored branches
    // in the vector checkit.
    while(leftpoint>=0 || rightpoint>=0) {
      // Previous step did not end on a leaf.
      leftpoint = kdvec[currentpoint].left;
      rightpoint = kdvec[currentpoint].right;
      dim = kdvec[currentpoint].dim;
      if(dim%6==1) pointdiff = kdvec[currentpoint].point.x - querypoint.x;
      else if(dim%6==2) pointdiff = kdvec[currentpoint].point.y - querypoint.y;
      else if(dim%6==3) pointdiff = kdvec[currentpoint].point.z - querypoint.z;
      else if(dim%6==4) pointdiff = kdvec[currentpoint].point.vx - querypoint.vx;
      else if(dim%6==5) pointdiff = kdvec[currentpoint].point.vy - querypoint.vy;
      else if(dim%6==0) pointdiff = kdvec[currentpoint].point.vz - querypoint.vz;

      goright = (pointdiff <= range); // possible hits lie to the left;
      goleft = (pointdiff >= -range); // possible hits lie to the right;
      if(goleft && goright) {
	// Current point might be within range.
	pdist2 = point6LDx2_dist2(querypoint,kdvec[currentpoint].point);
	if(pdist2 <= rng2) {
	  // Current point is within range. Add it to the output vector
	  indexvec.push_back(currentpoint);
	}
	if(leftpoint>=0) {
	  //Explore leftward first.
	  currentpoint = leftpoint;
	  if(rightpoint>=0) {
	    // Rightward branch will also be explored later
	    checknum++;
	    if(checknum>checkit.size()) {
	      checkit.push_back(rightpoint);
	    }
	    else {
	      checkit[checknum-1] = rightpoint;
	    }
	  }
	}
	else if(rightpoint>=0) {
	  // Leftward branch is a dead end: explore rightward branch
	  currentpoint = rightpoint;
	}
      }
      else if(goleft) {
	// Current point cannot be in range, but points that
	// are in range may lie along the left branch.
	if(leftpoint>=0) {
	  currentpoint = leftpoint;
	} else rightpoint=-1; // Dead end, make sure while-loop exits.
      } else if(goright) {
	// Current point cannot be in range, but points that
	// are in range may lie along the right branch.
	if(rightpoint>=0) {
	  currentpoint = rightpoint;
	} else leftpoint=-1;  // Dead end, make sure while-loop exits.
      } else {
	// Program concluded it should go neither left nor right.
	// The likely cause is that it encountered a NAN. Give up on this point.
	leftpoint=rightpoint=-1;
	cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
	cerr << "Query point:\n";
	cerr << querypoint.x << ", " << querypoint.y << ", " << querypoint.z << ", " << querypoint.vx << ", " << querypoint.vy << ", " << querypoint.vz << "\n";
	cerr << "Target point:\n";
 	cerr << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ", " << kdvec[currentpoint].point.z << ", " << kdvec[currentpoint].point.vx << ", " << kdvec[currentpoint].point.vy << ", " << kdvec[currentpoint].point.vz << "\n";
     }
      // Close while-loop checking if we've hit a leaf.
    }
    // We have climbed up the tree to a leaf. Go backwards through
    // the checkit vector and see if there is anything to check.
    checknum=checkit.size();
    while(checknum>=1 && checkit[checknum-1]<0) checknum--;
    if(checknum<=0) {
      //There were no valid entries to check: we're done.
      notdone=0;
    } else {
      //Set currentpoint to the last valid entry in checkit
      currentpoint = checkit[checknum-1];
      //Mark this point as used.
      checkit[checknum-1]=-1;
      leftpoint=rightpoint=0;
    }
  }
  return(0);
}

long double cluster_stats6D01(const vector <KD_point6LDx2> &cluster, vector <long double> &meanvals, vector <long double> &rmsvals)
{
  if(cluster.size()<2) {
    cerr << "ERROR: cluster_stats6D01 called with only " << cluster.size() << " points.\n";
    return(-1.0L);
  }
  
  long double xmean, ymean, zmean, vxmean, vymean, vzmean;
  xmean = ymean = zmean = vxmean = vymean = vzmean = 0.0L;
  long double xrms, yrms, zrms, vxrms, vyrms, vzrms;
  xrms = yrms = zrms = vxrms = vyrms = vzrms = 0.0L;
  long double norm = cluster.size();
  long double posrms = 0.0L;
  long double velrms = 0.0L;
  long double totalrms = 0.0L;
  
  for(long i=0; i<long(cluster.size()); i++) {
    xmean += cluster[i].point.x;   
    ymean += cluster[i].point.y;   
    zmean += cluster[i].point.z;
    vxmean += cluster[i].point.vx;   
    vymean += cluster[i].point.vy;   
    vzmean += cluster[i].point.vz;
  }
  xmean /= norm;
  ymean /= norm;
  zmean /= norm;
  vxmean /= norm;
  vymean /= norm;
  vzmean /= norm;
  
  meanvals.push_back(xmean);
  meanvals.push_back(ymean);
  meanvals.push_back(zmean);
  meanvals.push_back(vxmean);
  meanvals.push_back(vymean);
  meanvals.push_back(vzmean);
  
  for(long i=0; i<long(cluster.size()); i++) {
    xrms += LDSQUARE(cluster[i].point.x - xmean);   
    yrms += LDSQUARE(cluster[i].point.y - ymean);   
    zrms += LDSQUARE(cluster[i].point.z - zmean);
    vxrms += LDSQUARE(cluster[i].point.vx - vxmean);   
    vyrms += LDSQUARE(cluster[i].point.vy - vymean);   
    vzrms += LDSQUARE(cluster[i].point.vz - vzmean);
  }

  xrms /= norm;
  yrms /= norm;
  zrms /= norm;
  vxrms /= norm;
  vyrms /= norm;
  vzrms /= norm;

  posrms = xrms + yrms + zrms;
  velrms = vxrms + vyrms + vzrms;
  totalrms = posrms + velrms;

  xrms = sqrt(xrms);
  yrms = sqrt(yrms);
  zrms = sqrt(zrms);
  vxrms = sqrt(vxrms);
  vyrms = sqrt(vyrms);
  vzrms = sqrt(vzrms);
  posrms = sqrt(posrms);
  velrms = sqrt(velrms);
  totalrms = sqrt(totalrms);

  rmsvals.push_back(xrms);
  rmsvals.push_back(yrms);
  rmsvals.push_back(zrms);
  rmsvals.push_back(vxrms);
  rmsvals.push_back(vyrms);
  rmsvals.push_back(vzrms);
  rmsvals.push_back(posrms);
  rmsvals.push_back(velrms);
  rmsvals.push_back(totalrms);

  return(totalrms);
}
  
  
// DBSCAN_6D01: January 06, 2022:
// Given an input 6-dimensional kdtree produced by kdtree_6D01,
// find clusters using the DBSCAN algorithm, with range querying
// enabled by kdrange_6D01.
#define MINSPAN 1.0 // Temporal span must be at least this large (in days) for a bona fide cluster
#define MINDAYSTEPS 2 // A bona fide cluster must have at least this many intra-point
                      // time intervals greater than INTRANIGHTSTEP days.
#define INTRANIGHTSTEP 0.3 // Minimum interval in days between successive points
                           // in a tracklet, to enable them to be counted as being
                           // on separate nights.
int DBSCAN_6D01(vector <KD_point6LDx2> &kdtree, long double clustrad, int npt, const vector <det_bsc> &detvec, const vector <string> &det_id_vec, vector <KD6_clust> &outclusters, string rmsfile)
{
  long kdnum = kdtree.size();
  long kdct=0;
  int clustptct=0;
  long clusternum=0;
  long fakeclusternum=0;
  vector <long> queryout;
  vector <long> subquery;
  vector <long> clusterind;
  point6LDx2 querypoint = point6LDx2(0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0, 0);
  vector <KD_point6LDx2> cluster;
  KD6_clust oneclust = KD6_clust(0,{},{},{});
  vector <long double> meanvec;
  vector <long double> rmsvec;
  long double trms = 0.0L;
  long  i=0;
  vector <long> pointind;
  vector <long> pointjunk;
  vector <long double> clustmjd;
  vector <long double> mjdstep;
  long double timespan = 0.0L;
  int numdaysteps=0;

  ofstream outstream1 {rmsfile};
  outstream1.precision(17);

  // Loop on points
  for(kdct=0; kdct<kdnum; kdct++) {
    if(kdtree[kdct].flag == 0) {
      // Current point has not yet been assigned.
      // Range-query current point.
      querypoint = kdtree[kdct].point;
      queryout = {};
      cluster = {};
      clusterind = {};
      kdrange_6D01(kdtree, querypoint, clustrad, queryout);
      // If it's alone, mark it as noise.
      if(queryout.size()<=1) {
	kdtree[kdct].flag = -1; // Noise point.
	// cout << "Point " << kdct << ": noise\n";
      }
      else if(long(queryout.size()) >= npt) {
	// This is a core point of a new cluster.
	cout << "Point " << kdct << ": cluster core with " << queryout.size() << " neighbors.\n";
	fakeclusternum++;
	kdtree[kdct].flag = fakeclusternum;
	// Begin loading cluster
	clusterind.push_back(kdct);
	cluster.push_back(kdtree[kdct]);
	// Loop on points in cluster.
	clustptct=0;
	while(clustptct<long(queryout.size())) {
	  if(kdtree[queryout[clustptct]].flag != 0) {
	    // Current point has already been considered: skip to the next.
	    clustptct++;
	  } else {
	    // Range-query current cluster point.
	    querypoint = kdtree[queryout[clustptct]].point;
	    subquery={};
	    kdrange_6D01(kdtree, querypoint, clustrad, subquery);
	    if(long(subquery.size())>=npt) {
	      // This point is a core point.
	      kdtree[queryout[clustptct]].flag = fakeclusternum;
	      clusterind.push_back(queryout[clustptct]);
	      cluster.push_back(kdtree[queryout[clustptct]]);
	      // Add additional points to queryout as appropriate
	      for(i=0;i<long(subquery.size());i++) {
		if(kdtree[subquery[i]].flag == 0) queryout.push_back(subquery[i]);
	      }
	    } else {
	      // This is a border point. Add it to the cluster, but
	      // do not add its neighbors to queryout.
	      clusterind.push_back(queryout[clustptct]);
	      kdtree[queryout[clustptct]].flag = fakeclusternum;
	      cluster.push_back(kdtree[queryout[clustptct]]);
	    }
	    // Move on to next point in queryout vector
	    clustptct++;
	    // Close statement testing points for core vs. border status
	  }
	  // Close loop over the whole cluster
	}
	// Just finished loading a cluster.
	// Calculate some cluster statistics.
	cout << "Found cluster with " << cluster.size() << " = " << clusterind.size() << "points.\n";

	// Map cluster to individual detections.
	// create vector of unique detection indices.
	pointind={};
	for(i=0;i<long(clusterind.size());i++) {
	  pointind.push_back(kdtree[clusterind[i]].point.i1);
	  pointind.push_back(kdtree[clusterind[i]].point.i2);
	}
	// Sort vector of detection indices
	sort(pointind.begin(), pointind.end());
	// Cull out duplicate entries
	pointjunk = pointind;
	pointind={};
	pointind.push_back(pointjunk[0]);
	for(i=1; i<long(pointjunk.size()); i++) {
	  if(pointjunk[i]!=pointjunk[i-1]) pointind.push_back(pointjunk[i]);
	}
	// Load vector of detection MJD's
	clustmjd = {};
	for(i=0; i<long(pointind.size()); i++) {
	  clustmjd.push_back(detvec[pointind[i]].MJD);
	}
	// Sort vector of MJD's
	sort(clustmjd.begin(), clustmjd.end());
	timespan = clustmjd[clustmjd.size()-1] - clustmjd[0];
	// Load vector of MJD steps
	mjdstep={};
	for(i=1; i<long(clustmjd.size()); i++) {
	  mjdstep.push_back(clustmjd[i]-clustmjd[i-1]);
	}
	// Count steps large enough to suggest a daytime period between nights.
	numdaysteps=0;	
	for(i=0; i<long(mjdstep.size()); i++) {
	  if(mjdstep[i]>INTRANIGHTSTEP) numdaysteps++;
	}
	cout << "Unique pts: " << pointind.size() << " span: " << timespan << " daysteps: " << numdaysteps << "\n";
	// Does cluster pass the criteria for a linked detection?
	if(timespan >= MINSPAN && numdaysteps >= MINDAYSTEPS) {
	  clusternum++;
	  cout << "Cluster passes discovery criteria: will be designated as cluster " << clusternum << "\n";
	  outstream1 << "Found cluster " << clusternum << " with " << cluster.size() << " = " << clusterind.size() << "points.\n";
	  outstream1 << "Unique pts: " << pointind.size() << " span: " << timespan << " daysteps: " << numdaysteps << "\n";
	  meanvec = rmsvec = {};
	  trms = cluster_stats6D01(cluster, meanvec, rmsvec);
	  cout << "Cluster pos RMS: " << rmsvec[0] << " " << rmsvec[1] << " " << rmsvec[2] <<  " total pos " << rmsvec[6] << "\n";
	  cout << "Cluster vel RMS: " << rmsvec[3] << " " << rmsvec[4] << " " << rmsvec[5] <<  " total vel " << rmsvec[7] << "\n";
	  cout << "Cluster total RMS: " << rmsvec[8] << " = " << trms << "\n";
	  outstream1 << "Cluster pos RMS: " << rmsvec[0] << " " << rmsvec[1] << " " << rmsvec[2] <<  " total pos " << rmsvec[6] << "\n";
	  outstream1 << "Cluster vel RMS: " << rmsvec[3] << " " << rmsvec[4] << " " << rmsvec[5] <<  " total vel " << rmsvec[7] << "\n";
	  outstream1 << "Cluster total RMS: " << rmsvec[8] << " = " << trms << "\n";
	  // Write individual detections to output file
	  for(i=0; i<long(pointind.size()); i++) {
	    outstream1 << i << " " << detvec[pointind[i]].MJD << " " << detvec[pointind[i]].RA << " " << detvec[pointind[i]].Dec << " " << det_id_vec[pointind[i]] << "\n";
	    cout << i << " " << detvec[pointind[i]].MJD << " " << detvec[pointind[i]].RA << " " << detvec[pointind[i]].Dec << " " << det_id_vec[pointind[i]] << "\n";
	  }
	  outstream1 << "\n";
	  cout << "\n";
	  // Load cluster into oneclust.
	  oneclust = KD6_clust(cluster.size(),clusterind,meanvec,rmsvec);
	  // Push oneclust onto output vector.
	  outclusters.push_back(oneclust);
	} else {
	  cout << "Cluster failed criteria for a bona fide discovery.\n\n";
	}
	// Close statement testing for cluster vs. noise points.
      }
      // Close statement finding the next un-tested point.
    }
    // Close loop over entire k-d tree.
  }
  return(clusternum);
}
#undef MINSPAN
#undef MINDAYSTEPS
#undef INTRANIGHTSTEP

// DBSCAN_6D02: January 06, 2022:
// Like DBSCAN_6D01, but without kludgy debugging stuff.
// Given an input 6-dimensional kdtree produced by kdtree_6D01,
// find clusters using the DBSCAN algorithm, with range querying
// enabled by kdrange_6D01.
int DBSCAN_6D02(vector <KD_point6LDx2> &kdtree, long double clustrad, int npt, vector <KD6_clust> &outclusters)
{
  long kdnum = kdtree.size();
  long kdct=0;
  int clustptct=0;
  long clusternum=0;
  vector <long> queryout;
  vector <long> subquery;
  vector <long> clusterind;
  point6LDx2 querypoint = point6LDx2(0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0, 0);
  vector <KD_point6LDx2> cluster;
  KD6_clust oneclust = KD6_clust(0,{},{},{});
  vector <long double> meanvec;
  vector <long double> rmsvec;
  long i=0;

  // Loop on points
  for(kdct=0; kdct<kdnum; kdct++) {
    if(kdtree[kdct].flag == 0) {
      // Current point has not yet been assigned.
      // Range-query current point.
      querypoint = kdtree[kdct].point;
      queryout = {};
      cluster = {};
      clusterind = {};
      kdrange_6D01(kdtree, querypoint, clustrad, queryout);
      // If it's alone, mark it as noise.
      if(queryout.size()<=1) {
	kdtree[kdct].flag = -1; // Noise point.
	// cout << "Point " << kdct << ": noise\n";
      }
      else if(long(queryout.size()) >= npt) {
	// This is a core point of a new cluster.
	cout << "Point " << kdct << ": cluster core with " << queryout.size() << " neighbors.\n";
	clusternum++;
	kdtree[kdct].flag = clusternum;
	// Begin loading cluster
	clusterind.push_back(kdct);
	cluster.push_back(kdtree[kdct]);
	// Loop on points in cluster.
	clustptct=0;
	while(clustptct<long(queryout.size())) {
	  if(kdtree[queryout[clustptct]].flag != 0) {
	    // Current point has already been considered: skip to the next.
	    clustptct++;
	  } else {
	    // Range-query current cluster point.
	    querypoint = kdtree[queryout[clustptct]].point;
	    subquery={};
	    kdrange_6D01(kdtree, querypoint, clustrad, subquery);
	    if(long(subquery.size())>=npt) {
	      // This point is a core point.
	      kdtree[queryout[clustptct]].flag = clusternum;
	      clusterind.push_back(queryout[clustptct]);
	      cluster.push_back(kdtree[queryout[clustptct]]);
	      // Add additional points to queryout as appropriate
	      for(i=0;i<long(subquery.size());i++) {
		if(kdtree[subquery[i]].flag == 0) queryout.push_back(subquery[i]);
	      }
	    } else {
	      // This is a border point. Add it to the cluster, but
	      // do not add its neighbors to queryout.
	      clusterind.push_back(queryout[clustptct]);
	      kdtree[queryout[clustptct]].flag = clusternum;
	      cluster.push_back(kdtree[queryout[clustptct]]);
	    }
	    // Move on to next point in queryout vector
	    clustptct++;
	    // Close statement testing points for core vs. border status
	  }
	  // Close loop over the whole cluster
	}
	// Just finished loading a cluster.
	// Calculate some cluster statistics.
	meanvec = rmsvec = {};
	cluster_stats6D01(cluster, meanvec, rmsvec);
	// Load cluster into oneclust.
	oneclust = KD6_clust(cluster.size(),clusterind,meanvec,rmsvec);
	// Push oneclust onto output vector.
	outclusters.push_back(oneclust);
	// Close statement testing for cluster vs. noise points.
      }
      // Close statement finding the next un-tested point.
    }
    // Close loop over entire k-d tree.
  }
  return(clusternum);
}

// There follows a suite of programs re-implementing the 6-dimensional
// k-d tree and DBSCAN algorithms with integerized 6-D points for speed.


// conv_6LD_to_6i: January 07, 2022:
// Integerize an object of class point6LDx2 (6-dimensional
// long double point with 2 long integer indices) into an
// object of class point6ix2, to enable clustering algorithms
// to run faster. NOTE WELL that since the same scale factor
// is used for position and velocity, it is implicitly
// assumed that the velocity has previously been converted
// to position units via multiplication by an appropriate
// characteristic timescale.
point6ix2 conv_6LD_to_6i(point6LDx2 p1, long double scale)
{
  point6ix2 p2 = point6ix2(0,0,0,0,0,0,0,0);
  p2.x = int(p1.x/scale + 0.5);
  p2.y = int(p1.y/scale + 0.5);
  p2.z = int(p1.z/scale + 0.5);
  p2.vx = int(p1.vx/scale + 0.5);
  p2.vy = int(p1.vy/scale + 0.5);
  p2.vz = int(p1.vz/scale + 0.5);
  p2.i1 = p1.i1;
  p2.i2 = p1.i2;
  return(p2);
}

// conv_6i_to_6LD: January 07, 2022:
// Reverse the effect of conv_6LD_to_6i: expand an integerized
// 6-dimensional state vector back out into long doubles.
// Warning: Massive loss of precision, use only for crude
// comparisons!
point6LDx2 conv_6i_to_6LD(point6ix2 p1, long double scale)
{
  point6LDx2 p2 = point6LDx2(0L,0L,0L,0L,0L,0L,0,0);
  p2.x = p1.x*scale;
  p2.y = p1.y*scale;
  p2.z = p1.z*scale;
  p2.vx = p1.vx*scale;
  p2.vy = p1.vy*scale;
  p2.vz = p1.vz*scale;
  p2.i1 = p1.i1;
  p2.i2 = p1.i2;
  return(p2);
}

// conv_6d_to_6i: March 28, 2023
// Integerize an object of class point6dx2 (6-dimensional
// double point with 2 long integer indices) into an
// object of class point6ix2, to enable clustering algorithms
// to run faster. NOTE WELL that since the same scale factor
// is used for position and velocity, it is implicitly
// assumed that the velocity has previously been converted
// to position units via multiplication by an appropriate
// characteristic timescale.
point6ix2 conv_6d_to_6i(point6dx2 p1, double scale)
{
  point6ix2 p2 = point6ix2(0,0,0,0,0,0,0,0);
  p2.x = int(p1.x/scale + 0.5);
  p2.y = int(p1.y/scale + 0.5);
  p2.z = int(p1.z/scale + 0.5);
  p2.vx = int(p1.vx/scale + 0.5);
  p2.vy = int(p1.vy/scale + 0.5);
  p2.vz = int(p1.vz/scale + 0.5);
  p2.i1 = p1.i1;
  p2.i2 = p1.i2;
  return(p2);
}

// conv_6i_to_6d: March 28, 2023:
// Reverse the effect of conv_6d_to_6i: expand an integerized
// 6-dimensional state vector back out into long doubles.
// Warning: Massive loss of precision, use only for crude
// comparisons!
point6dx2 conv_6i_to_6d(point6ix2 p1, double scale)
{
  point6dx2 p2 = point6dx2(0l,0l,0l,0l,0l,0l,0,0);
  p2.x = p1.x*scale;
  p2.y = p1.y*scale;
  p2.z = p1.z*scale;
  p2.vx = p1.vx*scale;
  p2.vy = p1.vy*scale;
  p2.vz = p1.vz*scale;
  p2.i1 = p1.i1;
  p2.i2 = p1.i2;
  return(p2);
}


long medind_6ix2(const vector <point6ix2> &pointvec, int dim)
{
  vector <point6ix2> pvec = pointvec; //Mutable copy of immutable input vector
  for(long unsigned int i=0; i<pvec.size(); i++) pvec[i].i1=i; //Redefine indices
  long medpt = pvec.size()/2; // Central point of vector (it will be off by one half
                              // for a vector with even length, but we don't care).
  if(dim%6 == 1) sort(pvec.begin(), pvec.end(), lower_point6ix2_x()); // Sort vector by x
  else if(dim%6 == 2) sort(pvec.begin(), pvec.end(), lower_point6ix2_y()); // Sort vector by y
  else if(dim%6 == 3) sort(pvec.begin(), pvec.end(), lower_point6ix2_z()); // Sort vector by z
  else if(dim%6 == 4) sort(pvec.begin(), pvec.end(), lower_point6ix2_vx()); // Sort vector by vx
  else if(dim%6 == 5) sort(pvec.begin(), pvec.end(), lower_point6ix2_vy()); // Sort vector by vy
  else if(dim%6 == 0) sort(pvec.begin(), pvec.end(), lower_point6ix2_vz()); // Sort vector by vz
  else {
    cerr << "ERROR: medind_6ix2 recieved invalid dimension " << dim << "\n";
    return(-1);
  }
  return(pvec[medpt].i1); // Output the index of the median point in
                             // the original, unsorted input vector.
}

// splitix2: January 07, 2022:
// Given a vector of type point6ix2, split it into two halves,
// a left half with all the points lower than or equal to a specified
// split point along the chosen dimension (use dim = 1, 2, 3, 4, 5, or 6
// to split along x, y, z, vx, vy, or vz respectively).
int splitix2(const vector <point6ix2> &pointvec, int dim, long unsigned int splitpoint, vector <point6ix2> &left, vector <point6ix2> &right)
{
  long unsigned int i=0;
  int splitval = 0;

  if(dim%6==1) {
    // split on x
    splitval = pointvec[splitpoint].x;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].x<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==2) {
    // split on y
    splitval = pointvec[splitpoint].y;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].y<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==3) {
    // split on z
    splitval = pointvec[splitpoint].z;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].z<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==4) {
    // split on vx
    splitval = pointvec[splitpoint].vx;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].vx<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==5) {
    // split on vy
    splitval = pointvec[splitpoint].vy;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].vy<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%6==0) {
    // split on vz
    splitval = pointvec[splitpoint].vz;
    for(i=0 ; i<pointvec.size(); i++) {
      if(i!=splitpoint && pointvec[i].vz<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else {
      cerr << "ERROR: splitix2 asked to split on undefined dimension " << dim << "\n";
      return(1);
  } 
  return(0);
}

// kdtree_6i01: January 05, 2022
// Given an input root point, presumed to have null
// right and left branches, load the branches and then
// call kdtree_6i01 on them recursively.
int kdtree_6i01(const vector <point6ix2> &invec, int dim, long unsigned int splitpoint, long unsigned int kdroot, vector <KD_point6ix2> &kdvec)
{
  int lmed=0;
  int rmed=0;
  int kdct = kdvec.size()-1;
  long leftrootkd=-1;
  long rightrootkd=-1;
  point6ix2 point0 = point6ix2(0,0,0,0,0,0,0,0);
  KD_point6ix2 lp = KD_point6ix2(point0,-1,-1,0,0);
  KD_point6ix2 rp = KD_point6ix2(point0,-1,-1,0,0);
  vector <point6ix2> leftvec = {};
  vector <point6ix2> rightvec = {};

  // Basic outline: split the input vector into a left and a right
  // half, where the left half is below (or level with) splitpoint
  // in the dimension specified by dim, and the right half is above
  // splitpoint. Find the median of the left and right vectors,
  // and make the left and right branches from kdroot point to
  // these medians. Then call kdtree_6i01 itself recursively on
  // each of these median points, to peform a new split along a
  // different dimension.
  splitix2(invec,dim,splitpoint,leftvec,rightvec);

  dim+=1;
  while(dim>6) dim-=6;

  if(leftvec.size()==1) {
    // Left branch is just a single leaf
    lp = KD_point6ix2(leftvec[0],-1,-1,dim,0); // Define new point as a leaf: branches point nowhere
    kdvec.push_back(lp); // Add this new point to the KD tree.
    kdct++; // Keep track of how many point are in the tree
    kdvec[kdroot].left = kdct; // Stick the new point on the left branch of the input root.
  } else if(leftvec.size()<=0) {
    // There is no left branch
    kdvec[kdroot].left = -1;
  }
  if(rightvec.size()==1) {
    // Right branch is just a single leaf
    rp = KD_point6ix2(rightvec[0],-1,-1,dim,0);
    kdvec.push_back(rp);
    kdct++;
    kdvec[kdroot].right = kdct;
  } else if(rightvec.size()<=0) {
    // There is no right branch
    kdvec[kdroot].right = -1;
  }
   
 if(leftvec.size()>1) {
    lmed = medind_6ix2(leftvec,dim);
    lp = KD_point6ix2(leftvec[lmed],-1,-1,dim,0);
    kdvec.push_back(lp);
    kdct++;
    kdvec[kdroot].left = kdct;
    leftrootkd = kdct;
 }
 
  if(rightvec.size()>1) {
    rmed = medind_6ix2(rightvec,dim);
    rp = KD_point6ix2(rightvec[rmed],-1,-1,dim,0);
    kdvec.push_back(rp);
    kdct++;
    kdvec[kdroot].right = kdct;
    rightrootkd = kdct;
  }
  // I moved these down out of the above loops, because I thought
  // that otherwise, a bunch of stuff might get pushed down by the
  // left loop that the right loop didn't know about.
  if(leftvec.size()>1 && leftrootkd>=0) kdtree_6i01(leftvec,dim,lmed,leftrootkd,kdvec);
  else if(leftvec.size()>1 && leftrootkd<0)
    {
      cerr << "Error, kdtree_6i01 finds leftroot less than zero with leftvec.size() = " << leftvec.size() << "\n";
    }
  if(rightvec.size()>1 && rightrootkd>=0) kdtree_6i01(rightvec,dim,rmed,rightrootkd,kdvec);
  else if(rightvec.size()>1 && rightrootkd<0)
    {
      cerr << "Error, kdtree_6i01 finds rightroot less than zero with rightvec.size() = " << rightvec.size() << "\n";
    }

  return(0);
}

// point6ix2_dist2: January 07, 2022:
// Calculate the squared distance in 6-dimensional parameter space
// between two points of class point6LDx2.
long point6ix2_dist2(const point6ix2 &p1, const point6ix2 &p2)
{
  return(LSQUARE(p1.x - p2.x) + LSQUARE(p1.y - p2.y) + LSQUARE(p1.z - p2.z) + LSQUARE(p1.vx - p2.vx) + LSQUARE(p1.vy - p2.vy) + LSQUARE(p1.vz - p2.vz));
}	 

// kdrange_6i01: January 07, 2022:
// Given a k-d tree vector kdvec created by kdtree_6i01,
// perform a range-query about the specified point. Returns
// a vector indexing all of the points in the input k-d tree
// that lie within the specified range of the input coordinates.
// Assumes that kdvec[0] is the root of the k-d tree.
int kdrange_6i01(const vector <KD_point6ix2> &kdvec, const point6ix2 &querypoint, long range, vector <long> &indexvec)
{
  long rng2 = range*range;
  int notdone=1;
  int dim=1;
  int currentpoint=0;
  int leftpoint=0;
  int rightpoint=0;
  int goleft=0;
  int goright=0;
  long pointdiff = 0;
  long pdist2 = 0;
  vector <long> checkit={};
  unsigned int checknum=0;

  while(notdone>0) {
    // Climb to the top of the k-d tree, keeping track
    // of potentially interesting unexplored branches
    // in the vector checkit.
    while(leftpoint>=0 || rightpoint>=0) {
      // Previous step did not end on a leaf.
      leftpoint = kdvec[currentpoint].left;
      rightpoint = kdvec[currentpoint].right;
      dim = kdvec[currentpoint].dim;
      if(dim%6==1) pointdiff = kdvec[currentpoint].point.x - querypoint.x;
      else if(dim%6==2) pointdiff = kdvec[currentpoint].point.y - querypoint.y;
      else if(dim%6==3) pointdiff = kdvec[currentpoint].point.z - querypoint.z;
      else if(dim%6==4) pointdiff = kdvec[currentpoint].point.vx - querypoint.vx;
      else if(dim%6==5) pointdiff = kdvec[currentpoint].point.vy - querypoint.vy;
      else if(dim%6==0) pointdiff = kdvec[currentpoint].point.vz - querypoint.vz;

      goright = (pointdiff <= range); // possible hits lie to the left;
      goleft = (pointdiff >= -range); // possible hits lie to the right;
      if(goleft && goright) {
	// Current point might be within range.
	pdist2 = point6ix2_dist2(querypoint,kdvec[currentpoint].point);
	if(pdist2 <= rng2) {
	  // Current point is within range. Add it to the output vector
	  indexvec.push_back(currentpoint);
	}
	if(leftpoint>=0) {
	  //Explore leftward first.
	  currentpoint = leftpoint;
	  if(rightpoint>=0) {
	    // Rightward branch will also be explored later
	    checknum++;
	    if(checknum>checkit.size()) {
	      checkit.push_back(rightpoint);
	    }
	    else {
	      checkit[checknum-1] = rightpoint;
	    }
	  }
	}
	else if(rightpoint>=0) {
	  // Leftward branch is a dead end: explore rightward branch
	  currentpoint = rightpoint;
	}
      }
      else if(goleft) {
	// Current point cannot be in range, but points that
	// are in range may lie along the left branch.
	if(leftpoint>=0) {
	  currentpoint = leftpoint;
	} else rightpoint=-1; // Dead end, make sure while-loop exits.
      } else if(goright) {
	// Current point cannot be in range, but points that
	// are in range may lie along the right branch.
	if(rightpoint>=0) {
	  currentpoint = rightpoint;
	} else leftpoint=-1;  // Dead end, make sure while-loop exits.
      } else {
	// Program concluded it should go neither left nor right.
	// The likely cause is that it encountered a NAN. Give up on this point.
	leftpoint=rightpoint=-1;
	cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
	cerr << "Query point:\n";
	cerr << querypoint.x << ", " << querypoint.y << ", " << querypoint.z << ", " << querypoint.vx << ", " << querypoint.vy << ", " << querypoint.vz << "\n";
	cerr << "Target point:\n";
 	cerr << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ", " << kdvec[currentpoint].point.z << ", " << kdvec[currentpoint].point.vx << ", " << kdvec[currentpoint].point.vy << ", " << kdvec[currentpoint].point.vz << "\n";
     }
      // Close while-loop checking if we've hit a leaf.
    }
    // We have climbed up the tree to a leaf. Go backwards through
    // the checkit vector and see if there is anything to check.
    checknum=checkit.size();
    while(checknum>=1 && checkit[checknum-1]<0) checknum--;
    if(checknum<=0) {
      //There were no valid entries to check: we're done.
      notdone=0;
    } else {
      //Set currentpoint to the last valid entry in checkit
      currentpoint = checkit[checknum-1];
      //Mark this point as used.
      checkit[checknum-1]=-1;
      leftpoint=rightpoint=0;
    }
  }
  return(0);
}

// cluster_stats6i01: January 07, 2022:
// Calculates some cluster statistics for integerized state
// vectors clusters. Reverses the integerization in an approximate
// sense by multiplying back through by intconvscale.
double cluster_stats6i01(const vector <KD_point6ix2> &cluster, double intconvscale, vector <double> &meanvals, vector <double> &rmsvals)
{
  if(cluster.size()<2) {
    cerr << "ERROR: cluster_stats6i01 called with only " << cluster.size() << " points.\n";
    return(-1.0L);
  }
  
  double xmean, ymean, zmean, vxmean, vymean, vzmean;
  xmean = ymean = zmean = vxmean = vymean = vzmean = 0.0l;
  double xrms, yrms, zrms, vxrms, vyrms, vzrms;
  xrms = yrms = zrms = vxrms = vyrms = vzrms = 0.0l;
  double norm = cluster.size();
  double posrms = 0.0l;
  double velrms = 0.0l;
  double totalrms = 0.0l;
  
  for(unsigned int i=0; i<cluster.size(); i++) {
    xmean += intconvscale * cluster[i].point.x;   
    ymean += intconvscale * cluster[i].point.y;   
    zmean += intconvscale * cluster[i].point.z;
    vxmean += intconvscale * cluster[i].point.vx;   
    vymean += intconvscale * cluster[i].point.vy;   
    vzmean += intconvscale * cluster[i].point.vz;
  }
  xmean /= norm;
  ymean /= norm;
  zmean /= norm;
  vxmean /= norm;
  vymean /= norm;
  vzmean /= norm;
  
  meanvals.push_back(xmean);
  meanvals.push_back(ymean);
  meanvals.push_back(zmean);
  meanvals.push_back(vxmean);
  meanvals.push_back(vymean);
  meanvals.push_back(vzmean);
  
  for(unsigned int i=0; i<cluster.size(); i++) {
    xrms += DSQUARE(intconvscale * cluster[i].point.x - xmean);   
    yrms += DSQUARE(intconvscale * cluster[i].point.y - ymean);   
    zrms += DSQUARE(intconvscale * cluster[i].point.z - zmean);
    vxrms += DSQUARE(intconvscale * cluster[i].point.vx - vxmean);   
    vyrms += DSQUARE(intconvscale * cluster[i].point.vy - vymean);   
    vzrms += DSQUARE(intconvscale * cluster[i].point.vz - vzmean);
  }

  xrms /= norm;
  yrms /= norm;
  zrms /= norm;
  vxrms /= norm;
  vyrms /= norm;
  vzrms /= norm;

  
  posrms = xrms + yrms + zrms;
  velrms = vxrms + vyrms + vzrms;
  totalrms = posrms + velrms;

  xrms = sqrt(xrms);
  yrms = sqrt(yrms);
  zrms = sqrt(zrms);
  vxrms = sqrt(vxrms);
  vyrms = sqrt(vyrms);
  vzrms = sqrt(vzrms);
  posrms = sqrt(posrms);
  velrms = sqrt(velrms);
  totalrms = sqrt(totalrms);

  rmsvals.push_back(xrms);
  rmsvals.push_back(yrms);
  rmsvals.push_back(zrms);
  rmsvals.push_back(vxrms);
  rmsvals.push_back(vyrms);
  rmsvals.push_back(vzrms);
  rmsvals.push_back(posrms);
  rmsvals.push_back(velrms);
  rmsvals.push_back(totalrms);

  return(totalrms);
}
  

// DBSCAN_6i01: January 07, 2022:
// Like DBSCAN_6D02 (NOT DBSCAN_6D01), but uses integerized state
// vectors for speed.
int DBSCAN_6i01(vector <KD_point6ix2> &kdtree, double clustrad, int npt, double intconvscale, vector <KD6i_clust> &outclusters, int verbose)
{
  long kdnum = kdtree.size();
  long kdct=0;
  int clustptct=0;
  long clusternum=0;
  vector <long> queryout;
  vector <long> subquery;
  vector <long> clusterind;
  point6ix2 querypoint = point6ix2(0, 0, 0, 0, 0, 0, 0, 0);
  vector <KD_point6ix2> cluster;
  KD6i_clust oneclust = KD6i_clust(0,{},{},{});
  vector <double> meanvec;
  vector <double> rmsvec;
  long i=0;

  // Loop on points
  for(kdct=0; kdct<kdnum; kdct++) {
    if(kdtree[kdct].flag == 0) {
      // Current point has not yet been assigned.
      // Range-query current point.
      querypoint = kdtree[kdct].point;
      queryout = {};
      cluster = {};
      clusterind = {};
      kdrange_6i01(kdtree, querypoint, clustrad, queryout);
      // If it's alone, mark it as noise.
      if(queryout.size()<=1) {
	kdtree[kdct].flag = -1; // Noise point.
      }
      else if(long(queryout.size()) >= npt) {
	// This is a core point of a new cluster.
	if(verbose>=1) cout << "Point " << kdct << ": cluster core with " << queryout.size() << " neighbors.\n";
	clusternum++;
	kdtree[kdct].flag = clusternum;
	// Begin loading cluster
	clusterind.push_back(kdct);
	cluster.push_back(kdtree[kdct]);
	// Loop on points in cluster.
	clustptct=0;
	while(clustptct<long(queryout.size())) {
	  if(kdtree[queryout[clustptct]].flag != 0) {
	    // Current point has already been considered: skip to the next.
	    clustptct++;
	  } else {
	    // Range-query current cluster point.
	    querypoint = kdtree[queryout[clustptct]].point;
	    subquery={};
	    kdrange_6i01(kdtree, querypoint, clustrad, subquery);
	    if(long(subquery.size())>=npt) {
	      // This point is a core point.
	      kdtree[queryout[clustptct]].flag = clusternum;
	      clusterind.push_back(queryout[clustptct]);
	      cluster.push_back(kdtree[queryout[clustptct]]);
	      // Add additional points to queryout as appropriate
	      for(i=0;i<long(subquery.size());i++) {
		if(kdtree[subquery[i]].flag == 0) queryout.push_back(subquery[i]);
	      }
	    } else {
	      // This is a border point. Add it to the cluster, but
	      // do not add its neighbors to queryout.
	      clusterind.push_back(queryout[clustptct]);
	      // Also, do not flag border points: leave them free to be
	      // claimed by multiple clusters.
	      // kdtree[queryout[clustptct]].flag = clusternum;
	      cluster.push_back(kdtree[queryout[clustptct]]);
	    }
	    // Move on to next point in queryout vector
	    clustptct++;
	    // Close statement testing points for core vs. border status
	  }
	  // Close loop over the whole cluster
	}
	// Just finished loading a cluster.
	if(long(cluster.size())>=npt) {	  
	  // This cluster has enough points to be considered.
	  // Calculate some cluster statistics.
	  meanvec = rmsvec = {};
	  cluster_stats6i01(cluster, intconvscale, meanvec, rmsvec);

	  // Load cluster into oneclust
	  oneclust = KD6i_clust(cluster.size(),clusterind,meanvec,rmsvec);
	  // Push oneclust onto output vector.
	  outclusters.push_back(oneclust);
	} else {
	  cerr << "WARNING: DBSCAN_6i01 internal cluster " << clusternum << " is a dud, with only " << cluster.size() << " points of " << npt << " required.\n";
	  clusternum --;
	}
	// Close statement testing for cluster vs. noise points.
      } else {
	// More than one point, but fewer than npt, lie within
	// the clustering radius. Hence, we cannot assign a definitive
	// status to any of the points yet. For now, do nothing.
	;
      }
      // Close statement finding the next un-tested point.
    }
    // Close loop over entire k-d tree.
  }
  return(clusternum);
}

		
int celestial_to_statevec(double RA, double Dec,double delta,point3d &baryvec)
{
  double x,y,z,theta,phi,thetapole,phipole;
  x = y = z = theta = phi = thetapole = phipole = 0.0;
  theta = Dec/DEGPRAD;
  phi = RA/DEGPRAD;
  thetapole = NEPDEC/DEGPRAD;
      
  x = -cos(theta)*sin(phi); // sin(phi) and cos(phi) are switched here
  y = cos(theta)*cos(phi);  // because we're rotating by 270 degrees: that's
  z = sin(theta);           // the RA of the Ecliptic Pole.
  baryvec.x = delta*y;
  baryvec.y = delta*(-x*sin(thetapole) + z*cos(thetapole));
  baryvec.z = delta*(x*cos(thetapole) + z*sin(thetapole));
  // -x and y are switched above becuase we are rotating by 90 degrees
  // after the pole-switch, to get the old North Celestial Pole
  // on the +y axis where it should be.
  return(0);
}

int celestial_to_statevecLD(long double RA, long double Dec,long double delta,point3LD &baryvec)
{
  long double x,y,z,theta,phi,thetapole,phipole;
  x = y = z = theta = phi = thetapole = phipole = 0.0;
  theta = Dec/DEGPRAD;
  phi = RA/DEGPRAD;
  thetapole = NEPDEC/DEGPRAD;
      
  x = -cos(theta)*sin(phi); // sin(phi) and cos(phi) are switched here
  y = cos(theta)*cos(phi);  // because we're rotating by 270 degrees: that's
  z = sin(theta);           // the RA of the Ecliptic Pole.
  baryvec.x = delta*y;
  baryvec.y = delta*(-x*sin(thetapole) + z*cos(thetapole));
  baryvec.z = delta*(x*cos(thetapole) + z*sin(thetapole));
  // -x and y are switched above becuase we are rotating by 90 degrees
  // after the pole-switch, to get the old North Celestial Pole
  // on the +y axis where it should be.
  return(0);
}

int celestial_to_stateunit(double RA, double Dec,point3d &baryvec)
{
  double x,y,z,theta,phi,thetapole,phipole;
  x = y = z = theta = phi = thetapole = phipole = 0.0;
  theta = Dec/DEGPRAD;
  phi = RA/DEGPRAD;
  thetapole = NEPDEC/DEGPRAD;
      
  x = -cos(theta)*sin(phi); // sin(phi) and cos(phi) are switched here
  y = cos(theta)*cos(phi);  // because we're rotating by 270 degrees: that's
  z = sin(theta);           // the RA of the Ecliptic Pole.
  baryvec.x = y;
  baryvec.y = -x*sin(thetapole) + z*cos(thetapole);
  baryvec.z = x*cos(thetapole) + z*sin(thetapole);
  // -x and y are switched above because we are rotating by 90 degrees
  // after the pole-switch, to get the old North Celestial Pole
  // on the +y axis where it should be.
  return(0);
}

int celestial_to_stateunitLD(long double RA, long double Dec, point3LD &baryvec)
{
  long double x,y,z,theta,phi,thetapole,phipole;
  x = y = z = theta = phi = thetapole = phipole = 0.0;
  theta = Dec/DEGPRAD;
  phi = RA/DEGPRAD;
  thetapole = NEPDEC/DEGPRAD;
      
  x = -cos(theta)*sin(phi); // sin(phi) and cos(phi) are switched here
  y = cos(theta)*cos(phi);  // because we're rotating by 270 degrees: that's
  z = sin(theta);           // the RA of the Ecliptic Pole.
  baryvec.x = y;
  baryvec.y = -x*sin(thetapole) + z*cos(thetapole);
  baryvec.z = x*cos(thetapole) + z*sin(thetapole);
  // -x and y are switched above because we are rotating by 90 degrees
  // after the pole-switch, to get the old North Celestial Pole
  // on the +y axis where it should be.
  return(0);
}

// get_csv_string01: Given a line read from a csv file, and an
// starting point along that line, read the next comma-separated value,
// and put it into the output string. If the read was successful, return
// the line index of the comma, newline, or EOF at the end of the value read.
// Otherwise, return -1 as an error code.
int get_csv_string01(const string &lnfromfile, string &outstring, int startpoint)
{
  unsigned int i=startpoint;
  char c='0';
  outstring="";
  while(i<lnfromfile.size() && c!=',' && c!='\n' && c!=EOF) {
    c=lnfromfile[i];
    if(c!=',' && c!='\n' && c!=EOF) outstring.push_back(c);
    i++;
  }
  if(outstring.size() > 0) return(i-1); // Worked fine.
  else return(-1); // Error code
}

// get_sv_string01: Given a line read from a file with values
// separated by one of the following: comma, space, tab, pipe, or ampersand;
// and a starting point along that line, read the next value,
// and put it into the output string. If the read was successful, return
// the line index of the comma, space, tab, pipe, ampersand, newline,
// or EOF at the end of the value read.
// Otherwise, return -1 as an error code.
int get_sv_string01(const string &lnfromfile, string &outstring, int startpoint)
{
  unsigned int i=startpoint;
  char c='0';
  outstring="";
  while(i<lnfromfile.size() && c!=',' && c!='&' && c!='|' && c!=' ' && c!='\t' && c!='\r' && c!='\n' && c!='\v' && c!='\f' && c!='\n' && c!=EOF) {
    c=lnfromfile[i];
    if(c!=',' && c!='\n' && c!=EOF) outstring.push_back(c);
    i++;
  }
  if(outstring.size() > 0) return(i-1); // Worked fine.
  else return(-1); // Error code
}

// read_horizons_file: November 2021:
// Given an input state-vector ephemeris file downloaded directly
// from JPL Horizons, read it into position and velocity vectors.
// Note that the default unit convention is km for positions and
// km/sec for velocities. Note also that JPL state-vector
// ephemerides use dynamical TT, which is ahead of UT1 by about
// 70 seconds in 2022. This program does NOT correct TT to UT1,
// but programs making use of the ouput mjd, position, and velocity
// vectors might need to.
int read_horizons_file(string infile, vector <double> &mjdvec, vector <point3d> &pos, vector <point3d> &vel)
{
  ifstream instream1 {infile};
  point3d pospoint = point3d(0.0,0.0,0.0);
  point3d velpoint = point3d(0.0,0.0,0.0);
  int reachedeof=0;
  int ondata=0;
  unsigned int i=0;
  char c = '0';
  int reachedend=0;
  string teststring, lnfromfile;
  double x,y,z,vx,vy,vz,MJD;
  x = y = z = vx = vy = vz = MJD = 0.0l;
  
  if(!instream1) {
    cerr << "ERROR: can't open input file " << infile << "\n";
    return(1);
  }
  while(reachedeof==0 && !reachedend) {
    while(!ondata && !reachedend) {
      // See if this line contains the code for start-of-data.
      lnfromfile = "";
      teststring = "";
      getline(instream1,lnfromfile);
      if(instream1.eof()) reachedeof=1; //End of file, fine.
      else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
      else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
 
      if(lnfromfile.size()>=5) {
	for(i=0;i<5;i++) {
	  teststring.push_back(lnfromfile[i]);
	}
	if(teststring == "$$SOE") ondata=1;
	else if(teststring == "$$EOE") reachedend=1;
      }
    }
    while(ondata && !reachedend && reachedeof==0) {
      lnfromfile = "";
      teststring = "";
      getline(instream1,lnfromfile);
      if(instream1.eof()) reachedeof=1; //End of file, fine.
      else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
      else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
      if(lnfromfile.size()>=5) {
	for(i=0;i<5;i++) {
	  teststring.push_back(lnfromfile[i]);
	}
	if(teststring == "$$EOE") reachedend=1;
      }
      if(!reachedend && reachedeof==0) {
	//Attempt to read entire four-line block.
	//First line has MJD
	teststring = "";
	c='0';
	i=0;
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!=' ' && c!='\n' && c!=EOF) { 
	  c=lnfromfile[i];
	  if(c!=' ' && c!='=' && c!='\n' && c!=EOF) teststring.push_back(c);
	  if(c==EOF) reachedeof=1;
	  i++;
	}
	MJD=stod(teststring);
	//Next line has x,y,z positions
	lnfromfile = "";
	teststring = "";
	getline(instream1,lnfromfile);
	if(instream1.eof()) reachedeof=1; //End of file, fine.
	else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
	else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
	//Read to first equals sign
	c='0';
	i=0;
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  i++;
	}
	//Read to next equals sign, loading into teststring to get X
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	x = stod(teststring);	    
	teststring = "";
	//Read to next equals sign, loading into teststring to get Y
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	y = stod(teststring);	    
	teststring = "";
	//Read to next equals sign, loading into teststring to get Z
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	z = stod(teststring);
	//Next line has x,y,z velocities
	lnfromfile = "";
	teststring = "";
	getline(instream1,lnfromfile);
	if(instream1.eof()) reachedeof=1; //End of file, fine.
	else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
	else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
	//Read to first equals sign
	c='0';
	i=0;
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  i++;
	}
	//Read to next equals sign, loading into teststring to get XV
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	vx = stod(teststring);
	teststring = "";
	//Read to next equals sign, loading into teststring to get VY
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	vy = stod(teststring);	    
	teststring = "";
	//Read to next equals sign, loading into teststring to get VZ
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	vz = stod(teststring);
	// Load output vectors
	pospoint = point3d(x,y,z);
	velpoint = point3d(vx,vy,vz);
	pos.push_back(pospoint);
	vel.push_back(velpoint);
	mjdvec.push_back(MJD-MJDOFF);
	// Next line is of no current interest: read and discard
	lnfromfile = "";
	teststring = "";
	getline(instream1,lnfromfile);
	if(instream1.eof()) reachedeof=1; //End of file, fine.
	else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
	else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
      }
    }
  }
  if(reachedeof==1 && ondata==1) {
    //Read file successfully to the end.
    return(0);
  }
  else if(reachedeof==1) {
    //Did not find any data
    return(1);
  }
  else return(reachedeof);
}

// read_horizons_fileLD: November 2021:
// Given an input state-vector ephemeris file downloaded directly
// from JPL Horizons, read it into position and velocity vectors.
// Note that the default unit convention is km for positions and
// km/sec for velocities. Note also that JPL state-vector
// ephemerides use dynamical TT, which is ahead of UT1 by about
// 70 seconds in 2022. This program does NOT correct TT to UT1,
// but programs making use of the ouput mjd, position, and velocity
// vectors might need to.
int read_horizons_fileLD(string infile, vector <long double> &mjdvec, vector <point3LD> &pos, vector <point3LD> &vel)
{
  ifstream instream1 {infile};
  point3LD pospoint = point3LD(0.0,0.0,0.0);
  point3LD velpoint = point3LD(0.0,0.0,0.0);
  int reachedeof=0;
  int ondata=0;
  unsigned int i=0;
  char c = '0';
  int reachedend=0;
  string teststring, lnfromfile;
  long double x,y,z,vx,vy,vz,MJD;
  x = y = z = vx = vy = vz = MJD = 0.0l;
  
  if(!instream1) {
    cerr << "ERROR: can't open input file " << infile << "\n";
    return(1);
  }
  while(reachedeof==0 && !reachedend) {
    while(!ondata && !reachedend) {
      // See if this line contains the code for start-of-data.
      lnfromfile = "";
      teststring = "";
      getline(instream1,lnfromfile);
      if(instream1.eof()) reachedeof=1; //End of file, fine.
      else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
      else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
 
      if(lnfromfile.size()>=5) {
	for(i=0;i<5;i++) {
	  teststring.push_back(lnfromfile[i]);
	}
	if(teststring == "$$SOE") ondata=1;
	else if(teststring == "$$EOE") reachedend=1;
      }
    }
    while(ondata && !reachedend && reachedeof==0) {
      lnfromfile = "";
      teststring = "";
      getline(instream1,lnfromfile);
      if(instream1.eof()) reachedeof=1; //End of file, fine.
      else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
      else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
      if(lnfromfile.size()>=5) {
	for(i=0;i<5;i++) {
	  teststring.push_back(lnfromfile[i]);
	}
	if(teststring == "$$EOE") reachedend=1;
      }
      if(!reachedend && reachedeof==0) {
	//Attempt to read entire four-line block.
	//First line has MJD
	teststring = "";
	c='0';
	i=0;
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!=' ' && c!='\n' && c!=EOF) { 
	  c=lnfromfile[i];
	  if(c!=' ' && c!='=' && c!='\n' && c!=EOF) teststring.push_back(c);
	  if(c==EOF) reachedeof=1;
	  i++;
	}
	MJD=stold(teststring);
	//Next line has x,y,z positions
	lnfromfile = "";
	teststring = "";
	getline(instream1,lnfromfile);
	if(instream1.eof()) reachedeof=1; //End of file, fine.
	else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
	else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
	//Read to first equals sign
	c='0';
	i=0;
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  i++;
	}
	//Read to next equals sign, loading into teststring to get X
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	x = stold(teststring);	    
	teststring = "";
	//Read to next equals sign, loading into teststring to get Y
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	y = stold(teststring);	    
	teststring = "";
	//Read to next equals sign, loading into teststring to get Z
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	z = stold(teststring);
	//Next line has x,y,z velocities
	lnfromfile = "";
	teststring = "";
	getline(instream1,lnfromfile);
	if(instream1.eof()) reachedeof=1; //End of file, fine.
	else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
	else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
	//Read to first equals sign
	c='0';
	i=0;
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  i++;
	}
	//Read to next equals sign, loading into teststring to get XV
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	vx = stold(teststring);
	teststring = "";
	//Read to next equals sign, loading into teststring to get VY
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	vy = stold(teststring);	    
	teststring = "";
	//Read to next equals sign, loading into teststring to get VZ
	if(i<lnfromfile.size()) c=lnfromfile[i];
	while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
	  c=lnfromfile[i];
	  if(c==EOF) reachedeof=1;
	  if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
	  i++;
	}
	vz = stold(teststring);
	// Load output vectors
	pospoint = point3LD(x,y,z);
	velpoint = point3LD(vx,vy,vz);
	pos.push_back(pospoint);
	vel.push_back(velpoint);
	mjdvec.push_back(MJD-MJDOFF);
	// Next line is of no current interest: read and discard
	lnfromfile = "";
	teststring = "";
	getline(instream1,lnfromfile);
	if(instream1.eof()) reachedeof=1; //End of file, fine.
	else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
	else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
      }
    }
  }
  if(reachedeof==1 && ondata==1) {
    //Read file successfully to the end.
    return(0);
  }
  else if(reachedeof==1) {
    //Did not find any data
    return(1);
  }
  else return(reachedeof);
}

// read_horizons_csv: April 19, 2023:
// Given an input state-vector ephemeris file downloaded directly
// from JPL Horizons, WITH THE OPTIONAL CSV FORMAT SELECTED,
// read it into position and velocity vectors.
// Note that the default unit convention is km for positions and
// km/sec for velocities. Note also that JPL state-vector
// ephemerides use dynamical TT, which is ahead of UT1 by about
// 70 seconds in 2022. This program does NOT correct TT to UT1,
// but programs making use of the ouput mjd, position, and velocity
// vectors might need to.
int read_horizons_csv(string infile, vector <double> &mjdvec, vector <point3d> &pos, vector <point3d> &vel)
{
  ifstream instream1 {infile};
  point3d pospoint = point3d(0.0,0.0,0.0);
  point3d velpoint = point3d(0.0,0.0,0.0);
  int reachedeof=0;
  int ondata=0;
  int badread=0;
  unsigned int i=0;
  int reachedend=0;
  string teststring, lnfromfile, stest;
  double x,y,z,vx,vy,vz,MJD;
  x = y = z = vx = vy = vz = MJD = 0.0l;
  long double JD = 0L;
  int startpoint=0;
  int endpoint=0;
  
  if(!instream1) {
    cerr << "ERROR: can't open input file " << infile << "\n";
    return(1);
  }
  while(reachedeof==0 && !reachedend) {
    while(!ondata && !reachedend) {
      // See if this line contains the code for start-of-data.
      lnfromfile = "";
      teststring = "";
      getline(instream1,lnfromfile);
      if(instream1.eof()) reachedeof=1; //End of file, fine.
      else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
      else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
 
      if(lnfromfile.size()>=5) {
	for(i=0;i<5;i++) {
	  teststring.push_back(lnfromfile[i]);
	}
	if(teststring == "$$SOE") ondata=1;
	else if(teststring == "$$EOE") reachedend=1;
      }
    }
    while(ondata && !reachedend && reachedeof==0) {
      lnfromfile = "";
      teststring = "";
      getline(instream1,lnfromfile);
      if(instream1.eof()) reachedeof=1; //End of file, fine.
      else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
      else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
      if(lnfromfile.size()>=5) {
	for(i=0;i<5;i++) {
	  teststring.push_back(lnfromfile[i]);
	}
	if(teststring == "$$EOE") reachedend=1;
      }
      if(!reachedend && reachedeof==0) {
	// Read JD, and subtract offset to obtain MJD
	startpoint=0;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { JD = stold(stest); }
	  catch(...) { cerr << "ERROR: cannot read MJD string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	MJD = JD-MJDOFF;
	// Read and discard the calendar date
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint<=0) badread=1;
	// Read the state-vector X position
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { x = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read x string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector Y position
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { y = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read y string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector Z position
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { z = stod(stest); } 
	  catch(...) { cerr << "ERROR: cannot read z string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector VX velocity
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { vx = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read vx string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector VY velocity
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { vy = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read vy string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector VZ velocity
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { vz = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read vz string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}	  
	else badread=1;
	// Load output vectors
	pospoint = point3d(x,y,z);
	velpoint = point3d(vx,vy,vz);
	pos.push_back(pospoint);
	vel.push_back(velpoint);
	mjdvec.push_back(MJD);
      }
    }
  }
  if(reachedeof==1 && ondata==1) {
    //Read file successfully to the end.
    return(0);
  }
  else if(reachedeof==1) {
    //Did not find any data
    return(1);
  }
  else return(reachedeof);
}

// read_horizons_csv: April 19, 2023:
// Given an input state-vector ephemeris file downloaded directly
// from JPL Horizons, WITH THE OPTIONAL CSV FORMAT SELECTED,
// read it into a vector of type earthstate.
// Note that the default unit convention is km for positions and
// km/sec for velocities. Note also that JPL state-vector
// ephemerides use dynamical TT, which is ahead of UT1 by about
// 70 seconds in 2022. This program does NOT correct TT to UT1,
// but programs making use of the ouput mjd, position, and velocity
// vectors might need to.
int read_horizons_csv(string infile, vector <EarthState> &earthpos)
{
  ifstream instream1 {infile};
  EarthState earthonce = EarthState(0l,0l,0l,0l,0l,0l,0l);
  int reachedeof=0;
  int ondata=0;
  int badread=0;
  unsigned int i=0;
  int reachedend=0;
  string teststring, lnfromfile, stest;
  double x,y,z,vx,vy,vz,MJD;
  x = y = z = vx = vy = vz = MJD = 0.0l;
  long double JD = 0L;
  int startpoint=0;
  int endpoint=0;
  
  if(!instream1) {
    cerr << "ERROR: can't open input file " << infile << "\n";
    return(1);
  }
  while(reachedeof==0 && !reachedend) {
    while(!ondata && !reachedend) {
      // See if this line contains the code for start-of-data.
      lnfromfile = "";
      teststring = "";
      getline(instream1,lnfromfile);
      if(instream1.eof()) reachedeof=1; //End of file, fine.
      else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
      else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
 
      if(lnfromfile.size()>=5) {
	for(i=0;i<5;i++) {
	  teststring.push_back(lnfromfile[i]);
	}
	if(teststring == "$$SOE") ondata=1;
	else if(teststring == "$$EOE") reachedend=1;
      }
    }
    while(ondata && !reachedend && reachedeof==0) {
      lnfromfile = "";
      teststring = "";
      getline(instream1,lnfromfile);
      if(instream1.eof()) reachedeof=1; //End of file, fine.
      else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
      else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
      if(lnfromfile.size()>=5) {
	for(i=0;i<5;i++) {
	  teststring.push_back(lnfromfile[i]);
	}
	if(teststring == "$$EOE") reachedend=1;
      }
      if(!reachedend && reachedeof==0) {
	// Read JD, and subtract offset to obtain MJD
	startpoint=0;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { JD = stold(stest); }
	  catch(...) { cerr << "ERROR: cannot read MJD string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	MJD = JD-MJDOFF;
	// Read and discard the calendar date
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint<=0) badread=1;
	// Read the state-vector X position
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { x = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read x string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector Y position
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { y = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read y string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector Z position
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { z = stod(stest); } 
	  catch(...) { cerr << "ERROR: cannot read z string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector VX velocity
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { vx = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read vx string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector VY velocity
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { vy = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read vy string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}
	else badread=1;
	// Read the state-vector VZ velocity
	startpoint = endpoint+1;
	if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	if(endpoint>0) {
	  try { vz = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read vz string " << stest << " from line " << lnfromfile << "\n";
	    badread = 1; }
	}	  
	else badread=1;
	// Load output vectors
	earthonce = EarthState(MJD,x,y,z,vx,vy,vz);
	earthpos.push_back(earthonce);
      }
    }
  }
  if(reachedeof==1 && ondata==1) {
    //Read file successfully to the end.
    return(0);
  }
  else if(reachedeof==1) {
    //Did not find any data
    return(1);
  } else if(reachedeof==0 && reachedend==0) {
    cerr << "ERROR: Stopped reading file " << infile << " before the end\n";
    cerr << "Last point " << earthpos.size() << ", last line " << lnfromfile << "\n";
    return(1);
  } else return(reachedeof);
}


// poleswitch01: December 9, 2021: given a double precision input
// celestial position IN RADIANS, and the celestial position of the pole
// of a new coordinate system, calculates the position
// of the point in the new coordinate system.  The desired
// RA for the old pole in the new coordinates is also required.
int poleswitch01(const double &inRA, const double &inDec, const double &poleRA, const double &poleDec, const double &oldpoleRA, double &newRA, double &newDec)
{
  double x,y,z,xp,yp,zp,thetap,phip;
  int badphip;
  x=y=z=xp=yp=zp=thetap=phip = 0l;
  
  z = sin(inDec);
  x = cos(inDec)*cos(inRA-poleRA);
  y = cos(inDec)*sin(inRA-poleRA);

  zp = z*sin(poleDec) + x*cos(poleDec);
  xp = x*sin(poleDec) - z*cos(poleDec);
  yp = y;

  if(zp>1.0)
    {
      printf("WEIRD ERROR: POLESWITCH HAS z prime > 1.0!\n");
      printf("THIS VIOLATES BASIC TRIGONOMETRY\n");
      printf("zp-1.0 = %le\n",zp-1.0);
      printf("inDec = %lf, inRA = %lf, poleDec=%lf, poleRA=%lf\n",inDec,inRA,poleDec,poleRA);
      printf("xyz, xp yp zp = %lf %lf %lf %lf %lf %lf\n",x,y,z,xp,yp,zp);
      printf("SETTING zp to exactly 1.0.\n");
      zp = 1.0l;
    }
  thetap = asin(zp);

  phip=0.0;
  if(y==0.0)
    {
      if(x>=0.0)
	{
	  phip = 0.0;
	}
      else if(x<0.0)
	{
	  phip = M_PI;
	}
    }
  else if(y>0.0)
    {
      phip = M_PI/2.0l - atan(xp/yp);
    }
  else if(y<0.0)
    {
      phip = 3.0l*M_PI/2.0l - atan(xp/yp);
    }

  fflush(stdout);

  phip+=(oldpoleRA-M_PI);

  badphip = 0;
  if(phip<0.0||phip>=2.0l*M_PI)
    {
      badphip = 1;
    }
  while(badphip==1)
    {
      if(phip<0.0) phip+=2.0l*M_PI;
      else if(phip>=2.0l*M_PI) phip-=2.0l*M_PI;
      badphip = 0;
      if(phip<0.0||phip>=2.0l*M_PI)
	{
	  badphip = 1;
	}
    }
  newRA = phip;
  newDec = thetap;
  return(0);
}


// poleswitch01LD: December 9, 2021: given a double precision input
// celestial position IN RADIANS, and the celestial position of the pole
// of a new coordinate system, calculates the position
// of the point in the new coordinate system.  The desired
// RA for the old pole in the new coordinates is also required.
int poleswitch01LD(const long double &inRA, const long double &inDec, const long double &poleRA, const long double &poleDec, const long double &oldpoleRA, long double &newRA, long double &newDec)
{
  long double x,y,z,xp,yp,zp,thetap,phip;
  int badphip;
  x=y=z=xp=yp=zp=thetap=phip = 0L;
  
  z = sin(inDec);
  x = cos(inDec)*cos(inRA-poleRA);
  y = cos(inDec)*sin(inRA-poleRA);

  zp = z*sin(poleDec) + x*cos(poleDec);
  xp = x*sin(poleDec) - z*cos(poleDec);
  yp = y;

  if(zp>1.0L)
    {
      printf("WEIRD ERROR: POLESWITCH HAS z prime > 1.0!\n");
      printf("THIS VIOLATES BASIC TRIGONOMETRY\n");
      printf("zp-1.0 = %Le\n",zp-1.0L);
      printf("inDec = %Lf, inRA = %Lf, poleDec=%Lf, poleRA=%Lf\n",inDec,inRA,poleDec,poleRA);
      printf("xyz, xp yp zp = %Lf %Lf %Lf %Lf %Lf %Lf\n",x,y,z,xp,yp,zp);
      printf("SETTING zp to exactly 1.0.\n");
      zp = 1.0L;
    }
  thetap = asin(zp);

  phip=0.0L;
  if(y==0.0L)
    {
      if(x>=0.0L)
	{
	  phip = 0.0L;
	}
      else if(x<0.0L)
	{
	  phip = M_PI;
	}
    }
  else if(y>0.0L)
    {
      phip = M_PI/2.0L - atan(xp/yp);
    }
  else if(y<0.0L)
    {
      phip = 3.0L*M_PI/2.0L - atan(xp/yp);
    }

  fflush(stdout);

  phip+=(oldpoleRA-M_PI);

  badphip = 0;
  if(phip<0.0L||phip>=2.0L*M_PI)
    {
      badphip = 1;
    }
  while(badphip==1)
    {
      if(phip<0.0L) phip+=2.0L*M_PI;
      else if(phip>=2.0L*M_PI) phip-=2.0L*M_PI;
      badphip = 0;
      if(phip<0.0L||phip>=2.0L*M_PI)
	{
	  badphip = 1;
	}
    }
  newRA = phip;
  newDec = thetap;
  return(0);
}

// poleswitch02: December 9, 2021: Exactly like poleswitch01, but all
// angles are in degrees.
int poleswitch02(const double &inRA, const double &inDec, const double &poleRA, const double &poleDec, const double &oldpoleRA, double &newRA, double &newDec)
{
  double x,y,z,xp,yp,zp,thetap,phip;
  int badphip;
  x=y=z=xp=yp=zp=thetap=phip = 0l;
  
  z = sin(inDec/DEGPRAD);
  x = cos(inDec/DEGPRAD)*cos(inRA/DEGPRAD-poleRA/DEGPRAD);
  y = cos(inDec/DEGPRAD)*sin(inRA/DEGPRAD-poleRA/DEGPRAD);

  zp = z*sin(poleDec/DEGPRAD) + x*cos(poleDec/DEGPRAD);
  xp = x*sin(poleDec/DEGPRAD) - z*cos(poleDec/DEGPRAD);
  yp = y;

  if(zp>1.0)
    {
      printf("WEIRD ERROR: POLESWITCH HAS z prime > 1.0!\n");
      printf("THIS VIOLATES BASIC TRIGONOMETRY\n");
      printf("zp-1.0 = %le\n",zp-1.0);
      printf("inDec = %lf, inRA = %lf, poleDec=%lf, poleRA=%lf\n",inDec,inRA,poleDec,poleRA);
      printf("xyz, xp yp zp = %lf %lf %lf %lf %lf %lf\n",x,y,z,xp,yp,zp);
      printf("SETTING zp to exactly 1.0.\n");
      zp = 1.0l;
    }
  thetap = asin(zp);

  phip=0.0;
  if(y==0.0)
    {
      if(x>=0.0)
	{
	  phip = 0.0;
	}
      else if(x<0.0)
	{
	  phip = M_PI;
	}
    }
  else if(y>0.0)
    {
      phip = M_PI/2.0l - atan(xp/yp);
    }
  else if(y<0.0)
    {
      phip = 3.0l*M_PI/2.0l - atan(xp/yp);
    }

  fflush(stdout);

  phip+=(oldpoleRA/DEGPRAD-M_PI);

  badphip = 0;
  if(phip<0.0||phip>=2.0l*M_PI)
    {
      badphip = 1;
    }
  while(badphip==1)
    {
      if(phip<0.0) phip+=2.0l*M_PI;
      else if(phip>=2.0l*M_PI) phip-=2.0l*M_PI;
      badphip = 0;
      if(phip<0.0||phip>=2.0l*M_PI)
	{
	  badphip = 1;
	}
    }
  newRA = phip*DEGPRAD;
  newDec = thetap*DEGPRAD;
  return(0);
}


// poleswitch02LD: May 03, 2022: Exactly like poleswitch01LD, but with
// all angles in degrees.
int poleswitch02LD(const long double &inRA, const long double &inDec, const long double &poleRA, const long double &poleDec, const long double &oldpoleRA, long double &newRA, long double &newDec)
{
  long double x,y,z,xp,yp,zp,thetap,phip;
  int badphip;
  x=y=z=xp=yp=zp=thetap=phip = 0L;
  
  z = sin(inDec/DEGPRAD);
  x = cos(inDec/DEGPRAD)*cos(inRA/DEGPRAD-poleRA/DEGPRAD);
  y = cos(inDec/DEGPRAD)*sin(inRA/DEGPRAD-poleRA/DEGPRAD);

  zp = z*sin(poleDec/DEGPRAD) + x*cos(poleDec/DEGPRAD);
  xp = x*sin(poleDec/DEGPRAD) - z*cos(poleDec/DEGPRAD);
  yp = y;

  if(zp>1.0L)
    {
      printf("WEIRD ERROR: POLESWITCH HAS z prime > 1.0!\n");
      printf("THIS VIOLATES BASIC TRIGONOMETRY\n");
      printf("zp-1.0 = %Le\n",zp-1.0L);
      printf("inDec = %Lf, inRA = %Lf, poleDec=%Lf, poleRA=%Lf\n",inDec,inRA,poleDec,poleRA);
      printf("xyz, xp yp zp = %Lf %Lf %Lf %Lf %Lf %Lf\n",x,y,z,xp,yp,zp);
      printf("SETTING zp to exactly 1.0.\n");
      zp = 1.0L;
    }
  thetap = asin(zp);

  phip=0.0L;
  if(y==0.0L)
    {
      if(x>=0.0L)
	{
	  phip = 0.0L;
	}
      else if(x<0.0L)
	{
	  phip = M_PI;
	}
    }
  else if(y>0.0L)
    {
      phip = M_PI/2.0L - atan(xp/yp);
    }
  else if(y<0.0L)
    {
      phip = 3.0L*M_PI/2.0L - atan(xp/yp);
    }

  fflush(stdout);

  phip+=(oldpoleRA/DEGPRAD-M_PI);

  badphip = 0;
  if(phip<0.0L||phip>=2.0L*M_PI)
    {
      badphip = 1;
    }
  while(badphip==1)
    {
      if(phip<0.0L) phip+=2.0L*M_PI;
      else if(phip>=2.0L*M_PI) phip-=2.0L*M_PI;
      badphip = 0;
      if(phip<0.0L||phip>=2.0L*M_PI)
	{
	  badphip = 1;
	}
    }
  newRA = phip*DEGPRAD;
  newDec = thetap*DEGPRAD;
  return(0);
}

/*November 24, 2021: precess01a: Given celestial coordinates ra1,dec1,
and Modified Julian Date mjd, if precesscon>=0, assume inputs are
J2000.0 and precess to epoch-of-date; otherwise assume inputs are
epoch-of-date and precess to J2000.0*/
int precess01a(double ra1,double dec1,double mjd,double *ra2,double *dec2,int precesscon)
{
  double ndays,tds,zetaa,thetaa,zaa,ra4,dec4,cosra,sinra;
  ndays = tds = zetaa = thetaa = zaa = ra4 = dec4 = cosra = sinra = 0.0l;
  
  /*time since standard epoch*/
  ndays = mjd-(double)51544; /*Number of days since Jan 1, 2000*/
  tds = ndays/(double)36525.0;

  /*cubic approximation to precession*/
  zetaa = ZET0 + ZET1*tds + ZET2*tds*tds + ZET3*tds*tds*tds + ZET4*tds*tds*tds*tds + ZET5*tds*tds*tds*tds*tds;
  zaa = Z0 + Z1*tds + Z2*tds*tds + Z3*tds*tds*tds + Z4*tds*tds*tds*tds + Z5*tds*tds*tds*tds*tds;
  thetaa = THET1*tds + THET2*tds*tds + THET3*tds*tds*tds + THET4*tds*tds*tds*tds + THET5*tds*tds*tds*tds*tds;

  /*transformation from arcseconds to radians*/
  zetaa*=(M_PI/648000.0);
  zaa*=(M_PI/648000.0);
  thetaa*=(M_PI/648000.0);

  if(precesscon>=0)
    {
      /*Precess given J2000.0 coords to epoch of date*/

      /*get new declination*/
      if(dec1!=M_PI/2.0)
	{
	  /*printf("precess01 has normal declination case\n");*/
	  dec4 = asin(cos(ra1+zetaa)*sin(thetaa)*cos(dec1) + cos(thetaa)*sin(dec1));
	}
      else
	{
	  /*printf("precess01 has polar declination case\n");*/
	  dec4 = asin(cos(thetaa));
	}
      /*if declination was obviously meant to be the pole, but
        it has gotten a little off by roundoff error, collapse
        it to the pole.*/
      if(fabs(dec4-M_PI/2.0)<SMALLANG)
	{
	  dec4 = M_PI/2.0;
	}
      /*get new right ascension*/
      if(dec1!=M_PI/2.0&&dec4!=M_PI/2.0)
	{
	  /*printf("precess01 has normal right ascension case\n");*/
	  cosra = (cos(ra1+zetaa)*cos(thetaa)*cos(dec1) - sin(thetaa)*sin(dec1))/cos(dec4);
	  sinra = (sin(ra1+zetaa)*cos(dec1))/cos(dec4);
	  if(sinra>=0.0)
	    {
	      ra4 = acos(cosra)+zaa;
	    }
	  else{
	    ra4 = 2.0*M_PI - acos(cosra)+zaa;
	  }
	}
      else if(dec1==M_PI/2.0&&dec4!=M_PI/2.0)
	{	  
	  /*printf("precess01 has polar input right ascension case\n");*/
	  ra4 = M_PI + zaa;
	}
      else if(dec4==M_PI/2.0)
	{
	  /*printf("precess01 has polar output right ascension case\n");*/
	  ra4 = 0.0;
	}
      else
	{
	  printf("IMPOSSIBLE CASE ERROR IN precess01\n");
	}
    }
  else
    {
      /*Deprecess given epoch of date coords to J2000.0*/

      /*get new declination*/
      if(dec1!=M_PI/2.0)
	{
	  /*printf("precess01 has normal declination case\n");*/
	  dec4 = asin(-cos(ra1-zaa)*sin(thetaa)*cos(dec1) + cos(thetaa)*sin(dec1));
	}
      else
	{
	  /*printf("precess01 has polar declination case\n");*/
	  dec4 = asin(cos(thetaa));
	}
      /*if declination was obviously meant to be the pole, but
        it has gotten a little off by roundoff error, collapse
        it to the pole.*/
      if(fabs(dec4-M_PI/2.0)<SMALLANG)
	{
	  dec4 = M_PI/2.0;
	}
      /*get new right ascension*/
      if(dec1!=M_PI/2.0&&dec4!=M_PI/2.0)
	{
	  /*printf("precess01 has normal right ascension case\n");*/
	  cosra = (cos(ra1-zaa)*cos(thetaa)*cos(dec1) + sin(thetaa)*sin(dec1))/cos(dec4);
	  sinra = (sin(ra1-zaa)*cos(dec1))/cos(dec4);
	  if(sinra>=0.0)
	    {
	      ra4 = acos(cosra)-zetaa;
	    }
	  else{
	    ra4 = 2.0*M_PI - acos(cosra)-zetaa;
	  }
	}
      else if(dec1==M_PI/2.0&&dec4!=M_PI/2.0)
	{
	  /*printf("precess01 has polar input right ascension case\n");*/
	  ra4 = 2.0*M_PI-zetaa; /*Note this could be wrong*/
                              /*There might be two solutions*/   
	}
      else if(dec4==M_PI/2.0)
	{
	  /*printf("precess01 has polar output right ascension case\n");*/
	  ra4 = 0.0;
	}
      else
	{
	  printf("IMPOSSIBLE CASE ERROR IN precess01\n");
	}
    }
  *ra2 = ra4;
  *dec2 = dec4;
  return(1);
}

int precess01aLD(long double ra1,long double dec1,long double mjd,long double *ra2,long double *dec2,int precesscon)
{
  long double ndays,tds,zetaa,thetaa,zaa,ra4,dec4,cosra,sinra;
  ndays = tds = zetaa = thetaa = zaa = ra4 = dec4 = cosra = sinra = 0.0L;
  
  /*time since standard epoch*/
  ndays = mjd-51544L; /*Number of days since Jan 1, 2000*/
  tds = ndays/36525.0L;

  /*cubic approximation to precession*/
  zetaa = ZET0 + ZET1*tds + ZET2*tds*tds + ZET3*tds*tds*tds + ZET4*tds*tds*tds*tds + ZET5*tds*tds*tds*tds*tds;
  zaa = Z0 + Z1*tds + Z2*tds*tds + Z3*tds*tds*tds + Z4*tds*tds*tds*tds + Z5*tds*tds*tds*tds*tds;
  thetaa = THET1*tds + THET2*tds*tds + THET3*tds*tds*tds + THET4*tds*tds*tds*tds + THET5*tds*tds*tds*tds*tds;

  /*transformation from arcseconds to radians*/
  zetaa*=(M_PI/648000.0L);
  zaa*=(M_PI/648000.0L);
  thetaa*=(M_PI/648000.0L);

  if(precesscon>=0)
    {
      /*Precess given J2000.0 coords to epoch of date*/

      /*get new declination*/
      if(dec1!=M_PI/2.0L)
	{
	  /*printf("precess01 has normal declination case\n");*/
	  dec4 = asin(cos(ra1+zetaa)*sin(thetaa)*cos(dec1) + cos(thetaa)*sin(dec1));
	}
      else
	{
	  /*printf("precess01 has polar declination case\n");*/
	  dec4 = asin(cos(thetaa));
	}
      /*if declination was obviously meant to be the pole, but
        it has gotten a little off by roundoff error, collapse
        it to the pole.*/
      if(fabs(dec4-M_PI/2.0L)<SMALLANG)
	{
	  dec4 = M_PI/2.0L;
	}
      /*get new right ascension*/
      if(dec1!=M_PI/2.0L && dec4!=M_PI/2.0L)
	{
	  /*printf("precess01 has normal right ascension case\n");*/
	  cosra = (cos(ra1+zetaa)*cos(thetaa)*cos(dec1) - sin(thetaa)*sin(dec1))/cos(dec4);
	  sinra = (sin(ra1+zetaa)*cos(dec1))/cos(dec4);
	  if(sinra>=0.0)
	    {
	      ra4 = acos(cosra)+zaa;
	    }
	  else{
	    ra4 = 2.0L*M_PI - acos(cosra)+zaa;
	  }
	}
      else if(dec1==M_PI/2.0L && dec4!=M_PI/2.0L)
	{	  
	  /*printf("precess01 has polar input right ascension case\n");*/
	  ra4 = M_PI + zaa;
	}
      else if(dec4==M_PI/2.0L)
	{
	  /*printf("precess01 has polar output right ascension case\n");*/
	  ra4 = 0.0;
	}
      else
	{
	  printf("IMPOSSIBLE CASE ERROR IN precess01\n");
	}
    }
  else
    {
      /*Deprecess given epoch of date coords to J2000.0*/

      /*get new declination*/
      if(dec1!=M_PI/2.0L)
	{
	  /*printf("precess01 has normal declination case\n");*/
	  dec4 = asin(-cos(ra1-zaa)*sin(thetaa)*cos(dec1) + cos(thetaa)*sin(dec1));
	}
      else
	{
	  /*printf("precess01 has polar declination case\n");*/
	  dec4 = asin(cos(thetaa));
	}
      /*if declination was obviously meant to be the pole, but
        it has gotten a little off by roundoff error, collapse
        it to the pole.*/
      if(fabs(dec4-M_PI/2.0L)<SMALLANG)
	{
	  dec4 = M_PI/2.0L;
	}
      /*get new right ascension*/
      if(dec1!=M_PI/2.0L && dec4!=M_PI/2.0L)
	{
	  /*printf("precess01 has normal right ascension case\n");*/
	  cosra = (cos(ra1-zaa)*cos(thetaa)*cos(dec1) + sin(thetaa)*sin(dec1))/cos(dec4);
	  sinra = (sin(ra1-zaa)*cos(dec1))/cos(dec4);
	  if(sinra>=0.0)
	    {
	      ra4 = acos(cosra)-zetaa;
	    }
	  else{
	    ra4 = 2.0L*M_PI - acos(cosra)-zetaa;
	  }
	}
      else if(dec1==M_PI/2.0L && dec4!=M_PI/2.0L)
	{
	  /*printf("precess01 has polar input right ascension case\n");*/
	  ra4 = 2.0L*M_PI-zetaa; /*Note this could be wrong*/
                              /*There might be two solutions*/   
	}
      else if(dec4==M_PI/2.0L)
	{
	  /*printf("precess01 has polar output right ascension case\n");*/
	  ra4 = 0.0;
	}
      else
	{
	  printf("IMPOSSIBLE CASE ERROR IN precess01\n");
	}
    }
  *ra2 = ra4;
  *dec2 = dec4;
  return(1);
}


// solvematrix01: November 23, 2021
// Given a matrix with dimensions eqnum,eqnum+1, interpret
// it as a system of eqnum linear equations in eqnum unknowns,
// with the first term in each equation being the constant term
// and the others being the coefficients of x1,x2,x3,etc;
// solve for the vector of x values or report the matrix to
// be singular.
int solvematrix01(const vector <vector <double>> &inmat, int eqnum, vector <double> &outvec, int verbose)
{
  int eqhi,termhi,eqct,termct,i,j;
  double max,pivot;
  vector <vector <double>> newmat;
  vector <double> coeffvec;
  vector <double> outvec2;

  eqhi=termhi=eqct=termct=i=j=0;
  max=pivot=0.0;
     
  if(eqnum==1)
    {
      if(inmat[0][1]!=0.0)
	{
	  outvec[0] = -inmat[0][0]/inmat[0][1];
	  return(0);
	}
      else
	{
	  /*The coefficient for x1 was zero, so it is
            impossible to solve*/
	  printf("ERROR: solvematrix01 fed a singular matrix!\n");
	  outvec[0]=0.0;
	  return(1);
	}
    }
  else
    {
      make_dmat(eqnum-1,eqnum,newmat);
      make_dvec(eqnum,coeffvec);
      make_dvec(eqnum-1,outvec2);
      /*REDUCE THE NUMBER OF EQUATIONS BY 1*/
      /*Find the coefficient with the largest absolute value*/
      eqhi=0;
      termhi=1;
      max = fabs(inmat[0][1]);
      for(eqct=0;eqct<eqnum;eqct++)
	{
	  for(termct=1;termct<eqnum+1;termct++)
	    {
	      if(max<=fabs(inmat[eqct][termct]))
		{
		  max = fabs(inmat[eqct][termct]);
		  eqhi=eqct;
		  termhi=termct;
		}
	    }
	}
      pivot=inmat[eqhi][termhi];
      if(verbose>=1) printf("At %lf, coefficent %d of equation %d was the largest\n",pivot,termhi-1,eqhi);
      if(max==0.0)
	{
	  printf("ERROR: solvematrix01 fed a singular matrix!\n");
	  for(eqct=0;eqct<eqnum;eqct++) outvec[eqct]=0.0;
	  return(1);
	}
      /*Solve equation eqhi for the x value corresponding to termhi*/
      j=0;
      coeffvec[0]=inmat[eqhi][0]/pivot;
      for(termct=1;termct<eqnum+1;termct++)
	{
	  if(termct!=termhi)
	    {
	      j+=1;
	      coeffvec[j]=inmat[eqhi][termct]/pivot;
	    }
	}
      if(verbose>=1) printf("Coefficient substitution vector:\n");
      if(verbose>=1) printf("%lf",coeffvec[0]);
      if(verbose>=1) for(j=1;j<eqnum;j++) printf(" %lf",coeffvec[j]);
      if(verbose>=1) printf("\n");
      /*Substitute this solution into the other equations,
        creating a new matrix with one fewer equations*/
      i=0;
      for(eqct=0;eqct<eqnum;eqct++)
	{
	  if(eqct!=eqhi)
	    {
	      j=0;
	      newmat[i][j]=inmat[eqct][0]-coeffvec[0]*inmat[eqct][termhi];
	      for(termct=1;termct<eqnum+1;termct++)
		{
		  if(termct!=termhi)
		    {
		      j+=1;
		      newmat[i][j]=inmat[eqct][termct]-coeffvec[j]*inmat[eqct][termhi];
		    }
		}
	      i+=1;
	    }
	}
      if(verbose>=1) printf("New reduced matrix:\n");
      for(i=0;i<eqnum-1;i++)
	{
	  if(verbose>=1) printf("%lf",newmat[i][0]);
	  if(verbose>=1) for(j=1;j<eqnum;j++) printf(" %lf",newmat[i][j]);
	  if(verbose>=1) printf("\n");
	}
      /*Call solvematrix01 recursively on this new matrix*/
      if(solvematrix01(newmat,eqnum-1,outvec2,verbose))
	{
	  printf("ERROR: recursive call of solvematrix01 failed\n");
	  for(eqct=0;eqct<eqnum;eqct++) outvec[eqct]=0.0;
	  return(1);
	}
      if(verbose>=1) printf("Recursive result\n");
      if(verbose>=1) printf("%lf",outvec2[0]);
      if(verbose>=1) for(i=1;i<eqnum-1;i++) printf(" %lf",outvec2[i]);
      if(verbose>=1) printf("\n");
      /*Load the solution for everything except the pivot*/
      i=0;
      for(eqct=0;eqct<eqnum;eqct++)
	{
	  if(eqct!=termhi-1)
	    {
	      outvec[eqct]=outvec2[i];
	      i+=1;
	    }
	}
      /*Load the solution for the pivot*/
      outvec[termhi-1]=-coeffvec[0];
      for(i=0;i<eqnum-1;i++) outvec[termhi-1]-=coeffvec[i+1]*outvec2[i];
    }
  return(0);
}

int solvematrix01LD(const vector <vector <long double>> &inmat, int eqnum, vector <long double> &outvec, int verbose)
{
  int eqhi,termhi,eqct,termct,i,j;
  long double max,pivot;
  vector <vector <long double>> newmat;
  vector <long double> coeffvec;
  vector <long double> outvec2;

  eqhi=termhi=eqct=termct=i=j=0;
  max=pivot=0.0;
     
  if(eqnum==1)
    {
      if(inmat[0][1]!=0.0)
	{
	  outvec[0] = -inmat[0][0]/inmat[0][1];
	  return(0);
	}
      else
	{
	  /*The coefficient for x1 was zero, so it is
            impossible to solve*/
	  printf("ERROR: solvematrix01 fed a singular matrix!\n");
	  outvec[0]=0.0;
	  return(1);
	}
    }
  else
    {
      make_LDmat(eqnum-1,eqnum,newmat);
      make_LDvec(eqnum,coeffvec);
      make_LDvec(eqnum-1,outvec2);
      /*REDUCE THE NUMBER OF EQUATIONS BY 1*/
      /*Find the coefficient with the largest absolute value*/
      eqhi=0;
      termhi=1;
      max = fabs(inmat[0][1]);
      for(eqct=0;eqct<eqnum;eqct++)
	{
	  for(termct=1;termct<eqnum+1;termct++)
	    {
	      if(max<=fabs(inmat[eqct][termct]))
		{
		  max = fabs(inmat[eqct][termct]);
		  eqhi=eqct;
		  termhi=termct;
		}
	    }
	}
      pivot=inmat[eqhi][termhi];
      if(verbose>=1) printf("At %Lf, coefficent %d of equation %d was the largest\n",pivot,termhi-1,eqhi);
      if(max==0.0)
	{
	  printf("ERROR: solvematrix01 fed a singular matrix!\n");
	  for(eqct=0;eqct<eqnum;eqct++) outvec[eqct]=0.0;
	  return(1);
	}
      /*Solve equation eqhi for the x value corresponding to termhi*/
      j=0;
      coeffvec[0]=inmat[eqhi][0]/pivot;
      for(termct=1;termct<eqnum+1;termct++)
	{
	  if(termct!=termhi)
	    {
	      j+=1;
	      coeffvec[j]=inmat[eqhi][termct]/pivot;
	    }
	}
      if(verbose>=1) printf("Coefficient substitution vector:\n");
      if(verbose>=1) printf("%Lf",coeffvec[0]);
      if(verbose>=1) for(j=1;j<eqnum;j++) printf(" %Lf",coeffvec[j]);
      if(verbose>=1) printf("\n");
      /*Substitute this solution into the other equations,
        creating a new matrix with one fewer equations*/
      i=0;
      for(eqct=0;eqct<eqnum;eqct++)
	{
	  if(eqct!=eqhi)
	    {
	      j=0;
	      newmat[i][j]=inmat[eqct][0]-coeffvec[0]*inmat[eqct][termhi];
	      for(termct=1;termct<eqnum+1;termct++)
		{
		  if(termct!=termhi)
		    {
		      j+=1;
		      newmat[i][j]=inmat[eqct][termct]-coeffvec[j]*inmat[eqct][termhi];
		    }
		}
	      i+=1;
	    }
	}
      if(verbose>=1) printf("New reduced matrix:\n");
      for(i=0;i<eqnum-1;i++)
	{
	  if(verbose>=1) printf("%Lf",newmat[i][0]);
	  if(verbose>=1) for(j=1;j<eqnum;j++) printf(" %Lf",newmat[i][j]);
	  if(verbose>=1) printf("\n");
	}
      /*Call solvematrix01 recursively on this new matrix*/
      if(solvematrix01LD(newmat,eqnum-1,outvec2,verbose))
	{
	  printf("ERROR: recursive call of solvematrix01 failed\n");
	  for(eqct=0;eqct<eqnum;eqct++) outvec[eqct]=0.0;
	  return(1);
	}
      if(verbose>=1) printf("Recursive result\n");
      if(verbose>=1) printf("%Lf",outvec2[0]);
      if(verbose>=1) for(i=1;i<eqnum-1;i++) printf(" %Lf",outvec2[i]);
      if(verbose>=1) printf("\n");
      /*Load the solution for everything except the pivot*/
      i=0;
      for(eqct=0;eqct<eqnum;eqct++)
	{
	  if(eqct!=termhi-1)
	    {
	      outvec[eqct]=outvec2[i];
	      i+=1;
	    }
	}
      /*Load the solution for the pivot*/
      outvec[termhi-1]=-coeffvec[0];
      for(i=0;i<eqnum-1;i++) outvec[termhi-1]-=coeffvec[i+1]*outvec2[i];
    }
  return(0);
}


// perfectpoly01: November 23, 2021:
// Given a y-vector and an x-vector with the same number N of points,
// calculate the coefficients of a polynomial function y=f(x), of order
// N-1, that precisely passes through all the points.
int perfectpoly01(const vector <double> &x, const vector <double> &y, vector <double> &fitvec)
{
  vector <vector <double>> dmatrix;
  vector <double> outvec;
  unsigned int i=0;
  unsigned int j=0;
  unsigned int k=0;
  int status=0;
  unsigned int npoints = x.size();
  if(y.size() != npoints) {
    cerr << "ERROR: x and y vectors in perfectpoly don't have the same number of points!\n";
    return(1);
  }
  if(npoints<=1) {
    cerr << "ERROR: perfectpoly cannot fit just a single point!\n";
    return(2);
  }
  make_dmat(npoints,npoints+1,dmatrix);
  make_dvec(npoints,outvec);
  //cout << "perfectpoly fitting matrix:\n";
  for(i=0;i<npoints;i++) {
    dmatrix[i][0] = -y[i];
    //cout << dmatrix[i][0] << " ";
    for(j=1;j<=npoints;j++) {
      dmatrix[i][j] = 1.0;
      for(k=2;k<=j;k++) dmatrix[i][j]*=x[i];
      //cout << dmatrix[i][j] << " ";
    }
    //cout << "\n";
  }
  status=solvematrix01(dmatrix,npoints,fitvec,0);
  return(status);
}

int perfectpoly01LD(const vector <long double> &x, const vector <long double> &y, vector <long double> &fitvec)
{
  vector <vector <long double>> dmatrix;
  vector <long double> outvec;
  unsigned int i=0;
  unsigned int j=0;
  unsigned int k=0;
  int status=0;
  unsigned int npoints = x.size();
  if(y.size() != npoints) {
    cerr << "ERROR: x and y vectors in perfectpoly don't have the same number of points!\n";
    return(1);
  }
  if(npoints<=1) {
    cerr << "ERROR: perfectpoly cannot fit just a single point!\n";
    return(2);
  }
  make_LDmat(npoints,npoints+1,dmatrix);
  make_LDvec(npoints,outvec);
  //cout << "perfectpoly fitting matrix:\n";
  for(i=0;i<npoints;i++) {
    dmatrix[i][0] = -y[i];
    //cout << dmatrix[i][0] << " ";
    for(j=1;j<=npoints;j++) {
      dmatrix[i][j] = 1.0;
      for(k=2;k<=j;k++) dmatrix[i][j]*=x[i];
      //cout << dmatrix[i][j] << " ";
    }
    //cout << "\n";
  }
  status=solvematrix01LD(dmatrix,npoints,fitvec,0);
  return(status);
}

// planetpos01: November 24, 2021:
// Given a vector of MJD values and a vector of 3-D planet positions,
// use polynomial interpolation to obtain a precise estimate of the
// 3-D planet position at the time detmjd. It is assumed that the
// input time detmjd is in UT1 or some reasonable approximation
// thereof, while the planet ephemeris vectors are in dynamical TT.
// Hence, a correction is applied to the input time before the
// interpolation. If the calling function actually has time in TT
// already, planetpos01 should be called with the correction
// pre-subtracted from detmjd, so it will cancel out internally.
int planetpos01(double detmjd, int polyorder, const vector <double> &posmjd, const vector <point3d> &planetpos, point3d &outpos)
{
  long fitnum = polyorder+1;
  long pointsbefore = fitnum - fitnum/2;
  long pbf=0;
  vector <double> xvec;
  vector <double> yvec;
  vector <double> fitvec;
  double tdelt=0;
  double sumvar=0;
  long i=0;
  long j=0;
  long k=0;
  make_dvec(fitnum,fitvec);

  // Convert input time from UT1 standard to the dynamical time (TT) used
  // by JPL Horizons for the state vectors.
  detmjd += TTDELTAT/SOLARDAY;
  // TT is ahead of UT1 because the Earth's rotation is slowing down.
  
  //Interpolate to find the planet's exact position at the time
  //of the detection.
  pbf=0;
  i=long(posmjd.size());
  if(long(planetpos.size())!=i) {
    cerr << "ERROR: planetpos01 finds time and position vectors\n";
    cerr << "to have different lengths\n";
    return(1);
  }
  while(i>0 && pbf<pointsbefore)
    {
      i--;
      if(posmjd[i]<detmjd) pbf++;
    }
  pbf=i;
  xvec={};
  yvec={};
  tdelt = detmjd-posmjd[pbf];
  // Load vectors to fit x-coordinate of Earth's position.
  for(i=pbf;i<pbf+fitnum;i++) {
    xvec.push_back(posmjd[i]-posmjd[pbf]);
    yvec.push_back(planetpos[i].x);
  }
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outpos.x = fitvec[0];
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outpos.x += sumvar;
  }
  // Load vector to fit y-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[i].y);
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outpos.y = fitvec[0];
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outpos.y += sumvar;
  }
  // Load vector to fit z-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[i].z);
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outpos.z = fitvec[0];
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outpos.z += sumvar;
  }
  return(0);
}

// planetpos01: April 19, 2023: Exactly like overloaded function
// above, but takes in an EarthState vector instead of the two
// separate vectors for MJD and position.
// Uses polynomial interpolation to obtain a precise estimate of the
// 3-D planet position at the time detmjd. It is assumed that the
// input time detmjd is in UT1 or some reasonable approximation
// thereof, while the planet ephemeris vectors are in dynamical TT.
// Hence, a correction is applied to the input time before the
// interpolation. If the calling function actually has time in TT
// already, planetpos01 should be called with the correction
// pre-subtracted from detmjd, so it will cancel out internally.
int planetpos01(double detmjd, int polyorder, const vector <EarthState> &planetpos, point3d &outpos)
{
  long fitnum = polyorder+1;
  long pointsbefore = fitnum - fitnum/2;
  long pbf=0;
  vector <double> xvec;
  vector <double> yvec;
  vector <double> fitvec;
  double tdelt=0;
  double sumvar=0;
  long i=0;
  long j=0;
  long k=0;
  make_dvec(fitnum,fitvec);

  // Convert input time from UT1 standard to the dynamical time (TT) used
  // by JPL Horizons for the state vectors.
  detmjd += TTDELTAT/SOLARDAY;
  // TT is ahead of UT1 because the Earth's rotation is slowing down.
  
  //Interpolate to find the planet's exact position at the time
  //of the detection.
  pbf=0;
  i=long(planetpos.size());
  
  while(i>0 && pbf<pointsbefore)
    {
      i--;
      if(planetpos[i].MJD < detmjd) pbf++;
    }
  pbf=i;
  xvec={};
  yvec={};
  tdelt = detmjd-planetpos[pbf].MJD;
  // Load vectors to fit x-coordinate of Earth's position.
  for(i=pbf;i<pbf+fitnum;i++) {
    xvec.push_back(planetpos[i].MJD-planetpos[pbf].MJD);
    yvec.push_back(planetpos[i].x);
  }
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outpos.x = fitvec[0];
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outpos.x += sumvar;
  }
  // Load vector to fit y-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[i].y);
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outpos.y = fitvec[0];
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outpos.y += sumvar;
  }
  // Load vector to fit z-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[i].z);
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outpos.z = fitvec[0];
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outpos.z += sumvar;
  }
  return(0);
}

// planetpos01LD: November 24, 2021:
// Given a vector of MJD values and a vector of 3-D planet positions,
// use polynomial interpolation to obtain a precise estimate of the
// 3-D planet position at the time detmjd. It is assumed that the
// input time detmjd is in UT1 or some reasonable approximation
// thereof, while the planet ephemeris vectors are in dynamical TT.
// Hence, a correction is applied to the input time before the
// interpolation. If the calling function actually has time in TT
// already, planetpos01 should be called with the correction
// pre-subtracted from detmjd, so it will cancel out internally.
int planetpos01LD(long double detmjd, int polyorder, const vector <long double> &posmjd, const vector <point3LD> &planetpos, point3LD &outpos)
{
  long fitnum = polyorder+1;
  long pointsbefore = fitnum - fitnum/2;
  long pbf=0;
  vector <long double> xvec;
  vector <long double> yvec;
  vector <long double> fitvec;
  long double tdelt=0;
  long double sumvar=0;
  long i=0;
  int j=0;
  int k=0;
  make_LDvec(fitnum,fitvec);

  // Convert input time from UT1 standard to the dynamical time (TT) used
  // by JPL Horizons for the state vectors.
  detmjd += TTDELTAT/SOLARDAY;
  // TT is ahead of UT1 because the Earth's rotation is slowing down.
  
  //Interpolate to find the planet's exact position at the time
  //of the detection.
  pbf=0;
  i=long(posmjd.size());
  if(long(planetpos.size())!=i) {
    cerr << "ERROR: planetpos01 finds time and position vectors\n";
    cerr << "to have different lengths\n";
    return(1);
  }
  while(i>0 && pbf<pointsbefore)
    {
      i--;
      if(posmjd[i]<detmjd) pbf++;
    }
  pbf=i;
  xvec={};
  yvec={};
  tdelt = detmjd-posmjd[pbf];
  // Load vectors to fit x-coordinate of the planet's position.
  for(i=pbf;i<pbf+fitnum;i++) {
    xvec.push_back(posmjd[i]-posmjd[pbf]);
    yvec.push_back(planetpos[i].x);
  }
  // Solve for polynomial interpolation
  perfectpoly01LD(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outpos.x = fitvec[0];
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outpos.x += sumvar;
  }
  // Load vector to fit y-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[i].y);
  // Solve for polynomial interpolation
  perfectpoly01LD(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outpos.y = fitvec[0];
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outpos.y += sumvar;
  }
  // Load vector to fit z-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[i].z);
  // Solve for polynomial interpolation
  perfectpoly01LD(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outpos.z = fitvec[0];
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outpos.z += sumvar;
  }
  return(0);
}

// planetposvel01: December 09, 2021:
// Given a vector of MJD values, a vector of 3-D planet positions,
// and a vector of 3-D planet velocities (expected usually to be for the Sun,
// but could be for another object), perform a polynomial fit to the velocity,
// and use it to calculate an interpolated velocity and an integrated position
// at the instant detmjd. It is assumed that the
// input time detmjd is in UT1 or some reasonable approximation
// thereof, while the planet ephemeris vectors are in dynamical TT.
// Hence, a correction is applied to the input time before the
// interpolation. If the calling function actually has time in TT
// already, planetpos01 should be called with the correction
// pre-subtracted from detmjd, so it will cancel out internally.
int planetposvel01(double detmjd, int polyorder, const vector <double> &posmjd, const vector <point3d> &planetpos, const vector <point3d> &planetvel, point3d &outpos, point3d &outvel)
{
  long fitnum = polyorder+1;
  long pointsbefore = fitnum - fitnum/2;
  long pbf=0;
  vector <double> xvec;
  vector <double> yvec;
  vector <double> fitvec;
  double tdelt=0;
  double sumvar=0;
  long i=0;
  int j=0;
  int k=0;
  make_dvec(fitnum,fitvec);
  
  // Convert input time from UT1 standard to the dynamical time (TT) used
  // by JPL Horizons for the state vectors.
  detmjd += TTDELTAT/SOLARDAY;
  // TT is ahead of UT1 because the Earth's rotation is slowing down.

  //Interpolate to find the planet's exact velocity at the time
  //of the detection.
  pbf=0;
  i=long(posmjd.size());
  if(long(planetpos.size())!=i) {
    cerr << "ERROR: planetpos01 finds time and position vectors\n";
    cerr << "to have different lengths\n";
    return(1);
  }
  while(i>0 && pbf<pointsbefore)
    {
      i--;
      if(posmjd[i]<detmjd) pbf++;
    }
  pbf=i;
  xvec={};
  yvec={};
  tdelt = detmjd-posmjd[pbf];
  // Load vectors to fit x-coordinate of the planet's velocity
  for(i=pbf;i<pbf+fitnum;i++) {
    xvec.push_back(posmjd[i]-posmjd[pbf]);
    yvec.push_back(planetvel[i].x);
  }
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated velocity and position
  outvel.x = fitvec[0];
  outpos.x = planetpos[pbf].x + fitvec[0]*tdelt*SOLARDAY;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.x += sumvar;
    sumvar *= tdelt*SOLARDAY/((long double)(j+1.0L)); // One more power of tdelt, for the position.
    outpos.x += sumvar;
  }
  // Load vector to fit y-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetvel[i].y);
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outvel.y = fitvec[0];
  outpos.y = planetpos[pbf].y + fitvec[0]*tdelt*SOLARDAY;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.y += sumvar;
    sumvar *= tdelt*SOLARDAY/((long double)(j+1.0L)); // One more power of tdelt, for the position.
    outpos.y += sumvar;
  }
  // Load vector to fit z-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetvel[i].z);
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outvel.z = fitvec[0];
  outpos.z = planetpos[pbf].z + fitvec[0]*tdelt*SOLARDAY;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.z += sumvar;
    sumvar *= tdelt*SOLARDAY/((long double)(j+1.0L)); // One more power of tdelt, for the position.
    outpos.z += sumvar;
  }
  return(0);
}

// planetposvel01: April 19, 2023:
// Like overloaded function immediately above, but takes a single
// input vector of type EarthState, in place of the three vectors
// posmjd, planetpos, and planetvel.
// Given a vector of MJD values, a vector of 3-D planet positions,
// and a vector of 3-D planet velocities (expected usually to be for the Sun,
// but could be for another object), perform a polynomial fit to the velocity,
// and use it to calculate an interpolated velocity and an integrated position
// at the instant detmjd. It is assumed that the
// input time detmjd is in UT1 or some reasonable approximation
// thereof, while the planet ephemeris vectors are in dynamical TT.
// Hence, a correction is applied to the input time before the
// interpolation. If the calling function actually has time in TT
// already, planetpos01 should be called with the correction
// pre-subtracted from detmjd, so it will cancel out internally.
int planetposvel01(double detmjd, int polyorder, const vector <EarthState> &planetpos, point3d &outpos, point3d &outvel)
{
  long fitnum = polyorder+1;
  long pointsbefore = fitnum - fitnum/2;
  long pbf=0;
  vector <double> xvec;
  vector <double> yvec;
  vector <double> fitvec;
  double tdelt=0;
  double sumvar=0;
  long i=0;
  int j=0;
  int k=0;
  make_dvec(fitnum,fitvec);
  
  // Convert input time from UT1 standard to the dynamical time (TT) used
  // by JPL Horizons for the state vectors.
  detmjd += TTDELTAT/SOLARDAY;
  // TT is ahead of UT1 because the Earth's rotation is slowing down.

  //Interpolate to find the planet's exact velocity at the time
  //of the detection.
  pbf=0;
  i=long(planetpos.size());
  while(i>0 && pbf<pointsbefore)
    {
      i--;
      if(planetpos[i].MJD < detmjd) pbf++;
    }
  pbf=i;
  xvec={};
  yvec={};
  tdelt = detmjd-planetpos[pbf].MJD;
  // Load vectors to fit x-coordinate of the planet's velocity
  for(i=pbf;i<pbf+fitnum;i++) {
    xvec.push_back(planetpos[i].MJD-planetpos[pbf].MJD);
    yvec.push_back(planetpos[i].vx);
  }
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated velocity and position
  outvel.x = fitvec[0];
  outpos.x = planetpos[pbf].x + fitvec[0]*tdelt*SOLARDAY;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.x += sumvar;
    sumvar *= tdelt*SOLARDAY/((long double)(j+1.0L)); // One more power of tdelt, for the position.
    outpos.x += sumvar;
  }
  // Load vector to fit y-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[i].vy);
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outvel.y = fitvec[0];
  outpos.y = planetpos[pbf].y + fitvec[0]*tdelt*SOLARDAY;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.y += sumvar;
    sumvar *= tdelt*SOLARDAY/((long double)(j+1.0L)); // One more power of tdelt, for the position.
    outpos.y += sumvar;
  }
  // Load vector to fit z-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[i].vz);
  // Solve for polynomial interpolation
  perfectpoly01(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outvel.z = fitvec[0];
  outpos.z = planetpos[pbf].z + fitvec[0]*tdelt*SOLARDAY;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.z += sumvar;
    sumvar *= tdelt*SOLARDAY/((long double)(j+1.0L)); // One more power of tdelt, for the position.
    outpos.z += sumvar;
  }
  return(0);
}

// planetposvel01LD: December 09, 2021:
// Given a vector of MJD values, a vector of 3-D planet positions,
// and a vector of 3-D planet velocities (expected usually to be for the Sun,
// but could be for another object), perform a polynomial fit to the velocity,
// and use it to calculate an interpolated velocity and an integrated position
// at the instant detmjd. It is assumed that the
// input time detmjd is in UT1 or some reasonable approximation
// thereof, while the planet ephemeris vectors are in dynamical TT.
// Hence, a correction is applied to the input time before the
// interpolation. If the calling function actually has time in TT
// already, planetpos01 should be called with the correction
// pre-subtracted from detmjd, so it will cancel out internally.
int planetposvel01LD(long double detmjd, int polyorder, const vector <long double> &posmjd, const vector <point3LD> &planetpos, const vector <point3LD> &planetvel, point3LD &outpos, point3LD &outvel)
{
  long fitnum = polyorder+1;
  long pointsbefore = fitnum - fitnum/2;
  long pbf=0;
  vector <long double> xvec;
  vector <long double> yvec;
  vector <long double> fitvec;
  long double tdelt=0;
  long double sumvar=0;
  long i=0;
  int j=0;
  int k=0;
  make_LDvec(fitnum,fitvec);
  
  // Convert input time from UT1 standard to the dynamical time (TT) used
  // by JPL Horizons for the state vectors.
  detmjd += TTDELTAT/SOLARDAY;
  // TT is ahead of UT1 because the Earth's rotation is slowing down.

  //Interpolate to find the planet's exact velocity at the time
  //of the detection.
  pbf=0;
  i=long(posmjd.size());
  if(long(planetpos.size())!=i) {
    cerr << "ERROR: planetpos01 finds time and position vectors\n";
    cerr << "to have different lengths\n";
    return(1);
  }
  while(i>0 && pbf<pointsbefore)
    {
      i--;
      if(posmjd[i]<detmjd) pbf++;
    }
  pbf=i;
  xvec={};
  yvec={};
  tdelt = detmjd-posmjd[pbf];
  // Load vectors to fit x-coordinate of the planet's velocity
  for(i=pbf;i<pbf+fitnum;i++) {
    xvec.push_back(posmjd[i]-posmjd[pbf]);
    yvec.push_back(planetvel[i].x);
  }
  // Solve for polynomial interpolation
  perfectpoly01LD(xvec,yvec,fitvec);
  // Calculate interpolated velocity and position
  outvel.x = fitvec[0];
  outpos.x = planetpos[pbf].x + fitvec[0]*tdelt*SOLARDAY;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.x += sumvar;
    sumvar *= tdelt*SOLARDAY/((long double)(j+1.0L)); // One more power of tdelt, for the position.
    outpos.x += sumvar;
  }
  // Load vector to fit y-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetvel[i].y);
  // Solve for polynomial interpolation
  perfectpoly01LD(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outvel.y = fitvec[0];
  outpos.y = planetpos[pbf].y + fitvec[0]*tdelt*SOLARDAY;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.y += sumvar;
    sumvar *= tdelt*SOLARDAY/((long double)(j+1.0L)); // One more power of tdelt, for the position.
    outpos.y += sumvar;
  }
  // Load vector to fit z-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetvel[i].z);
  // Solve for polynomial interpolation
  perfectpoly01LD(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outvel.z = fitvec[0];
  outpos.z = planetpos[pbf].z + fitvec[0]*tdelt*SOLARDAY;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.z += sumvar;
    sumvar *= tdelt*SOLARDAY/((long double)(j+1.0L)); // One more power of tdelt, for the position.
    outpos.z += sumvar;
  }
  return(0);
}


// nplanetpos01LD: November 30, 2021:
// like planetpos01LD, but calculates positions for
// not just one planet, but planetnum different planets.
// planetpos01LD assumes the position vectors of the
// respective planets have been pushed back, in order,
// into planetpos.  It is assumed that the
// input time detmjd is in UT1 or some reasonable approximation
// thereof, while the planet ephemeris vectors are in dynamical TT.
// Hence, a correction is applied to the input time before the
// interpolation. If the calling function actually has time in TT
// already, planetpos01 should be called with the correction
// pre-subtracted from detmjd, so it will cancel out internally.
int nplanetpos01LD(long double detmjd, int planetnum, int polyorder, const vector <long double> &posmjd, const vector <point3LD> &planetpos, vector <point3LD> &outpos)
{
  int fitnum = polyorder+1;
  int mjdsize = posmjd.size();
  int posvecsize = planetpos.size(); 
  int pointsbefore = fitnum - fitnum/2;
  int pbf=0;
  int planetct=0;
  vector <long double> xvec;
  vector <long double> yvec;
  vector <long double> fitvec;
  long double tdelt=0;
  long double sumvar=0;
  int i=0;
  int j=0;
  int k=0;
  point3LD postemp = point3LD(0,0,0);
  make_LDvec(fitnum,fitvec);

  // Convert input time from UT1 standard to the dynamical time (TT) used
  // by JPL Horizons for the state vectors.
  detmjd += TTDELTAT/SOLARDAY;
  // TT is ahead of UT1 because the Earth's rotation is slowing down.

  if(mjdsize*planetnum != posvecsize) {
    cout << "ERROR: vector sizes do not match in nplanetpos01LD!\n";
    cout << "Got a size of " << posvecsize << " for the position vector\n";
    cout << "with " << planetnum << " planets and an mjd vector size of " << mjdsize << "\n";
    cout << "Size of the position vector should be the product of these two:\n";
    cout << mjdsize << "x" << planetnum << " = " << mjdsize*planetnum << ", not " << posvecsize << "\n";
    return(-1);
  }
  
  // Identify the point pbf that is a bit before the specified
  // time detmjd, and is the appropriate point to start the interpolation.
  // The number of timesteps pbf should be before detmjd depends on
  // the order of the polynomial interpolation. For higher-order interpolations,
  // we need a larger number of total points, and detmjd should always be
  // near the center of the set of point being considered.
  pbf=0;
  i=mjdsize;
  while(i>0 && pbf<pointsbefore)
    {
      i--;
      if(posmjd[i]<detmjd) pbf++;
    }
  pbf=i;

  //Interpolate to find the exact position of each planet at detmjd.
  for(planetct=0; planetct<planetnum; planetct++) {
    xvec={};
    yvec={};
    tdelt = detmjd-posmjd[pbf];
    // Load vectors to fit x-coordinate of the planet's position.
    for(i=pbf;i<pbf+fitnum;i++) {
      xvec.push_back(posmjd[i]-posmjd[pbf]);
      yvec.push_back(planetpos[planetct*mjdsize+i].x);
    }
    // Solve for polynomial interpolation
    perfectpoly01LD(xvec,yvec,fitvec);
    // Calculate interpolated position.
    postemp.x = fitvec[0];
    for(j=1;j<fitnum;j++) {
      sumvar = fitvec[j]*tdelt;
      for(k=2;k<=j;k++) sumvar*=tdelt;
      postemp.x += sumvar;
    }
    // Load vector to fit y-coordinate
    yvec={};
    for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[planetct*mjdsize+i].y);
    // Solve for polynomial interpolation
    perfectpoly01LD(xvec,yvec,fitvec);
    // Calculate interpolated position.
    postemp.y = fitvec[0];
    for(j=1;j<fitnum;j++) {
      sumvar = fitvec[j]*tdelt;
      for(k=2;k<=j;k++) sumvar*=tdelt;
      postemp.y += sumvar;
    }
    // Load vector to fit z-coordinate
    yvec={};
    for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetpos[planetct*mjdsize+i].z);
    // Solve for polynomial interpolation
    perfectpoly01LD(xvec,yvec,fitvec);
    // Calculate interpolated position.
    postemp.z = fitvec[0];
    for(j=1;j<fitnum;j++) {
      sumvar = fitvec[j]*tdelt;
      for(k=2;k<=j;k++) sumvar*=tdelt;
      postemp.z += sumvar;
    }
    outpos.push_back(postemp);
  }
  return(0);
}

// nplanetgrab01LD: December01, 2021:
// Simply grab positions for every planet at a given, specified
// time step: that is, serve as nplanetpos01LD for the case where
// no interpolation is necessary because we want just a specific
// integer time step.
int nplanetgrab01LD(int pointrequest, int planetnum, const vector <long double> &posmjd, const vector <point3LD> &planetpos, vector <point3LD> &outpos)
{
  int mjdsize = posmjd.size();
  int posvecsize = planetpos.size(); 
  point3LD postemp = point3LD(0,0,0);
  int planetct=0;

  if(mjdsize*planetnum != posvecsize) {
    cout << "ERROR: vector sizes do not match in nplanetgrab01LD!\n";
    cout << "Got a size of " << posvecsize << " for the position vector\n";
    cout << "with " << planetnum << " planets and an mjd vector size of " << mjdsize << "\n";
    cout << "Size of the position vector should be the product of these two:\n";
    cout << mjdsize << "x" << planetnum << " = " << mjdsize*planetnum << ", not " << posvecsize << "\n";
    return(-1);
  }

  if(pointrequest<0 || pointrequest >= mjdsize)
    {
      cout << "ERROR: out-of-range point request in nplanetgrab01LD!\n";
      cout << "Requested point " << pointrequest << " when allowed range is 0-" << mjdsize-1 << "\n";
      return(2);
    }
   
  //Interpolate to find the exact position of each planet at detmjd.
  for(planetct=0; planetct<planetnum; planetct++) {
    postemp.x = planetpos[planetct*mjdsize+pointrequest].x;
    postemp.y = planetpos[planetct*mjdsize+pointrequest].y;
    postemp.z = planetpos[planetct*mjdsize+pointrequest].z;
    outpos.push_back(postemp);
  }
  return(0);
}

// observer_barycoords01: November 24, 2021:
// Given the MJD of an observation, and a file giving barycentric state-vector
// coordinates for the Earth, the longitude and MPC latitude sin and cos terms
// for an observatory, calculate the observer's topocentric position
// in barycentric state vector coordinates.
// Note that the handling of Earth's rotation assumes that the
// input MJD is UT1, while the ephemeris vectors posmjd
// and planetpos are in dynamical TT. Hence, after calculating
// aspects related to Earth's rotation with detmjd as input,
// planetpos01 is called which internally converts the input
// UT1 into TT.
int observer_barycoords01(double detmjd, int polyorder, double lon, double obscos, double obssine, const vector <double> &posmjd, const vector <point3d> &planetpos, point3d &outpos)
{
  double gmst=0;
  double djdoff = detmjd-double(51544.5);
  double zenithRA=0.0;
  double zenithDec=0.0;
  double junkRA=0.0;
  double junkDec=0.0;
  double crad = sqrt(obscos*obscos + obssine*obssine)*EARTHEQUATRAD;
  point3d obs_from_geocen = point3d(0,0,0);
  point3d geocen_from_barycen = point3d(0,0,0);
  
  gmst = double(18.697374558) + double(24.06570982441908)*djdoff;
  // Add the longitude, converted to hours.
  // Note: at this point it stops being gmst.
  gmst += lon/double(15.0);
  // Get a value between 0 and 24.0.
  while(gmst>=24.0) gmst-=double(24.0);
  while(gmst<0.0) gmst+=double(24.0);
  // Convert to degrees
  zenithRA = gmst * double(15.0);
  // Get zenithDec    
  if(obscos!=0.0) {
    zenithDec = atan(obssine/obscos)*DEGPRAD;
  } else if(obssine>=0.0) {
    zenithDec = 90.0;
  } else {
    zenithDec=-90.0;
  }
  // Now zenithRA and zenithDec are epoch-of-date coordinates.
  // If you want them in J2000.0, this is the place to convert them.
  int precesscon=-1; //Precess epoch-of-date to J2000.0
  junkRA = zenithRA/DEGPRAD;
  junkDec = zenithDec/DEGPRAD;
  precess01a(junkRA,junkDec,detmjd,&zenithRA,&zenithDec,precesscon);
  zenithRA*=DEGPRAD;
  zenithDec*=DEGPRAD;
  celestial_to_statevec(zenithRA,zenithDec,crad,obs_from_geocen);
  // crad is the distance from the geocenter to the observer, in AU.
  planetpos01(detmjd,polyorder,posmjd,planetpos,geocen_from_barycen);
  outpos.x = geocen_from_barycen.x + obs_from_geocen.x;
  outpos.y = geocen_from_barycen.y + obs_from_geocen.y;
  outpos.z = geocen_from_barycen.z + obs_from_geocen.z;
  return(0);
}

// observer_baryvel01: March 21, 2023:
// Exactly like observer_barycoords01, but also calculates the
// observer's barycentric velocity.
// Note that the handling of Earth's rotation assumes that the
// input MJD is UT1, while the ephemeris vectors posmjd
// and planetpos are in dynamical TT. Hence, after calculating
// aspects related to Earth's rotation with detmjd as input,
// planetpos01 is called which internally converts the input
// UT1 into TT.
int observer_baryvel01(double detmjd, int polyorder, double lon, double obscos, double obssine, const vector <double> &posmjd, const vector <point3d> &planetpos, const vector <point3d> &planetvel, point3d &outpos, point3d &outvel)
{
  double gmst=0;
  double djdoff = detmjd-51544.5l;
  double zenithRA=0.0l;
  double zenithDec=0.0l;
  double velRA=0.0l; // RA of the observer's geocentric velocity vector
  double velDec=0.0l; // Dec of the observer's geocentric velocity vector, always 0.
  double junkRA=0.0l;
  double junkDec=0.0l;
  double crad = sqrt(obscos*obscos + obssine*obssine)*EARTHEQUATRAD;
  double cvel = 2.0l*M_PI*obscos*EARTHEQUATRAD/SIDEREALDAY;
  point3d obs_from_geocen = point3d(0,0,0);
  point3d vel_from_geocen = point3d(0,0,0);
  point3d geocen_from_barycen = point3d(0,0,0);
  point3d vel_from_barycen = point3d(0,0,0);

  gmst = 18.697374558l + 24.06570982441908l*djdoff;
  // Add the longitude, converted to hours.
  // Note: at this point it stops being gmst.
  gmst += lon/15.0l;
  // Get a value between 0 and 24.0.
  while(gmst>=24.0l) gmst-=24.0l;
  while(gmst<0.0l) gmst+=24.0l;
  // Convert to degrees
  zenithRA = gmst * 15.0l;
  // Get zenithDec    
  if(obscos!=0.0l) {
    zenithDec = atan(obssine/obscos)*DEGPRAD;
  } else if(obssine>=0.0l) {
    zenithDec = 90.0l;
  } else {
    zenithDec=-90.0l;
  }
  // Calculate RA and Dec of the observer's geocentric velocity vector.
  velRA = zenithRA+90.0l;
  if(velRA >= 360.0l) velRA -= 360.0l;
  velDec = 0.0l;
  
  // Now zenithRA and zenithDec are epoch-of-date coordinates.
  // If you want them in J2000.0, this is the place to convert them.
  int precesscon=-1; //Precess epoch-of-date to J2000.0
  junkRA = zenithRA/DEGPRAD;
  junkDec = zenithDec/DEGPRAD;
  precess01a(junkRA,junkDec,detmjd,&zenithRA,&zenithDec,precesscon);
  zenithRA*=DEGPRAD;
  zenithDec*=DEGPRAD;
  celestial_to_statevec(zenithRA,zenithDec,crad,obs_from_geocen);
  // crad is the distance from the geocenter to the observer, in AU.
  // Now velRA and velDec are also epoch-of-date coordinates,
  // and hence should be converted to J2000.0.
  junkRA = velRA/DEGPRAD;
  junkDec = velDec/DEGPRAD;
  precess01a(junkRA,junkDec,detmjd,&velRA,&velDec,precesscon);
  velRA*=DEGPRAD;
  velDec*=DEGPRAD;
  celestial_to_statevec(velRA,velDec,cvel,vel_from_geocen);
  // cvel is the Earth's rotation velocity at the latitude of
  // the observer, in km/sec.

  planetposvel01(detmjd,polyorder,posmjd,planetpos,planetvel,geocen_from_barycen,vel_from_barycen);

  outpos.x = geocen_from_barycen.x + obs_from_geocen.x;
  outpos.y = geocen_from_barycen.y + obs_from_geocen.y;
  outpos.z = geocen_from_barycen.z + obs_from_geocen.z;
  outvel.x = vel_from_barycen.x + vel_from_geocen.x;
  outvel.y = vel_from_barycen.y + vel_from_geocen.y;
  outvel.z = vel_from_barycen.z + vel_from_geocen.z;
 
  cout << "spinvel: " << vel_from_geocen.x << " " << vel_from_geocen.y << " " << vel_from_geocen.z << "\n";
  cout << "orbvel: " << vel_from_barycen.x << " " << vel_from_barycen.y << " " << vel_from_barycen.z << "\n";
  cout << "total: " << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  return(0);
}

// observer_baryvel01: April 19, 2023:
// Like overloaded fuction just above, but takes as input a
// single vector of type EarthState, rather than three separate
// vectors posmjd, planetpos, and planetvel.
// Note that the handling of Earth's rotation assumes that the
// input MJD is UT1, while the ephemeris vectors posmjd
// and planetpos are in dynamical TT. Hence, after calculating
// aspects related to Earth's rotation with detmjd as input,
// planetpos01 is called which internally converts the input
// UT1 into TT.
int observer_baryvel01(double detmjd, int polyorder, double lon, double obscos, double obssine, const vector <EarthState> &earthpos, point3d &outpos, point3d &outvel)
{
  double gmst=0;
  double djdoff = detmjd-51544.5l;
  double zenithRA=0.0l;
  double zenithDec=0.0l;
  double velRA=0.0l; // RA of the observer's geocentric velocity vector
  double velDec=0.0l; // Dec of the observer's geocentric velocity vector, always 0.
  double junkRA=0.0l;
  double junkDec=0.0l;
  double crad = sqrt(obscos*obscos + obssine*obssine)*EARTHEQUATRAD;
  double cvel = 2.0l*M_PI*obscos*EARTHEQUATRAD/SIDEREALDAY;
  point3d obs_from_geocen = point3d(0,0,0);
  point3d vel_from_geocen = point3d(0,0,0);
  point3d geocen_from_barycen = point3d(0,0,0);
  point3d vel_from_barycen = point3d(0,0,0);

  gmst = 18.697374558l + 24.06570982441908l*djdoff;
  // Add the longitude, converted to hours.
  // Note: at this point it stops being gmst.
  gmst += lon/15.0l;
  // Get a value between 0 and 24.0.
  while(gmst>=24.0l) gmst-=24.0l;
  while(gmst<0.0l) gmst+=24.0l;
  // Convert to degrees
  zenithRA = gmst * 15.0l;
  // Get zenithDec    
  if(obscos!=0.0l) {
    zenithDec = atan(obssine/obscos)*DEGPRAD;
  } else if(obssine>=0.0l) {
    zenithDec = 90.0l;
  } else {
    zenithDec=-90.0l;
  }
  // Calculate RA and Dec of the observer's geocentric velocity vector.
  velRA = zenithRA+90.0l;
  if(velRA >= 360.0l) velRA -= 360.0l;
  velDec = 0.0l;
  
  // Now zenithRA and zenithDec are epoch-of-date coordinates.
  // If you want them in J2000.0, this is the place to convert them.
  int precesscon=-1; //Precess epoch-of-date to J2000.0
  junkRA = zenithRA/DEGPRAD;
  junkDec = zenithDec/DEGPRAD;
  precess01a(junkRA,junkDec,detmjd,&zenithRA,&zenithDec,precesscon);
  zenithRA*=DEGPRAD;
  zenithDec*=DEGPRAD;
  celestial_to_statevec(zenithRA,zenithDec,crad,obs_from_geocen);
  // crad is the distance from the geocenter to the observer, in AU.
  // Now velRA and velDec are also epoch-of-date coordinates,
  // and hence should be converted to J2000.0.
  junkRA = velRA/DEGPRAD;
  junkDec = velDec/DEGPRAD;
  precess01a(junkRA,junkDec,detmjd,&velRA,&velDec,precesscon);
  velRA*=DEGPRAD;
  velDec*=DEGPRAD;
  celestial_to_statevec(velRA,velDec,cvel,vel_from_geocen);
  // cvel is the Earth's rotation velocity at the latitude of
  // the observer, in km/sec.

  planetposvel01(detmjd,polyorder,earthpos,geocen_from_barycen,vel_from_barycen);

  outpos.x = geocen_from_barycen.x + obs_from_geocen.x;
  outpos.y = geocen_from_barycen.y + obs_from_geocen.y;
  outpos.z = geocen_from_barycen.z + obs_from_geocen.z;
  outvel.x = vel_from_barycen.x + vel_from_geocen.x;
  outvel.y = vel_from_barycen.y + vel_from_geocen.y;
  outvel.z = vel_from_barycen.z + vel_from_geocen.z;
 
  cout << "spinvel: " << vel_from_geocen.x << " " << vel_from_geocen.y << " " << vel_from_geocen.z << "\n";
  cout << "orbvel: " << vel_from_barycen.x << " " << vel_from_barycen.y << " " << vel_from_barycen.z << "\n";
  cout << "total: " << outvel.x << " " << outvel.y << " " << outvel.z << "\n";
  return(0);
}


// observer_barycoords01LD: November 24, 2021:
// Given the MJD of an observation, and a file giving barycentric state-vector
// coordinates for the Earth, the longitude and MPC latitude sin and cos terms
// for an observatory, calculate the observer's topocentric position
// in barycentric state vector coordinates.
// Note that the handling of Earth's rotation assumes that the
// input MJD is UT1, while the ephemeris vectors posmjd
// and planetpos are in dynamical TT. Hence, after calculating
// aspects related to Earth's rotation with detmjd as input,
// planetpos01 is called which internally converts the input
// UT1 into TT.
int observer_barycoords01LD(long double detmjd, int polyorder, long double lon, long double obscos, long double obssine, const vector <long double> &posmjd, const vector <point3LD> &planetpos, point3LD &outpos)
{
  long double gmst=0;
  long double djdoff = detmjd-51544.5L;
  long double zenithRA=0.0;
  long double zenithDec=0.0;
  long double junkRA=0.0;
  long double junkDec=0.0;
  long double crad = sqrt(obscos*obscos + obssine*obssine)*EARTHEQUATRAD;
  point3LD obs_from_geocen = point3LD(0,0,0);
  point3LD geocen_from_barycen = point3LD(0,0,0);
 
  gmst = 18.697374558L + 24.06570982441908L*djdoff;
  // Add the longitude, converted to hours.
  // Note: at this point it stops being gmst.
  gmst += lon/15.0L;
  // Get a value between 0 and 24.0.
  while(gmst>=24.0L) gmst -= 24.0L;
  while(gmst<0.0L) gmst += 24.0L;
  // Convert to degrees
  zenithRA = gmst * 15.0L;
  // Get zenithDec    
  if(obscos!=0.0L) {
    zenithDec = atan(obssine/obscos)*DEGPRAD;
  } else if(obssine>=0.0L) {
    zenithDec = 90.0L;
  } else {
    zenithDec=-90.0L;
  }
  // Now zenithRA and zenithDec are epoch-of-date coordinates.
  // If you want them in J2000.0, this is the place to convert them.
  int precesscon=-1; //Precess epoch-of-date to J2000.0
  junkRA = zenithRA/DEGPRAD;
  junkDec = zenithDec/DEGPRAD;
  precess01aLD(junkRA,junkDec,detmjd,&zenithRA,&zenithDec,precesscon);
  zenithRA*=DEGPRAD;
  zenithDec*=DEGPRAD;
  celestial_to_statevecLD(zenithRA,zenithDec,crad,obs_from_geocen);
  // crad is the distance from the geocenter to the observer, in AU.
  planetpos01LD(detmjd,polyorder,posmjd,planetpos,geocen_from_barycen);
  // cout << "obs_from_geocen: " << obs_from_geocen.x << " " << obs_from_geocen.y << " " << obs_from_geocen.z << " \n";
  // cout << "geocen_from_barycen: " << geocen_from_barycen.x << " " << geocen_from_barycen.y << " " << geocen_from_barycen.z << "\n";
  outpos.x = geocen_from_barycen.x + obs_from_geocen.x;
  outpos.y = geocen_from_barycen.y + obs_from_geocen.y;
  outpos.z = geocen_from_barycen.z + obs_from_geocen.z;
  return(0);
}


// helioproj01: November 26, 2021
int helioproj01(point3d unitbary, point3d obsbary,double heliodist,double &geodist, point3d &projbary)
{
  double a,b,c;
  double alphapos;
  // double opdistcos,sunelong,opelon;
  double barydist,obsdot;
  
  //cout << fixed << setprecision(9) << "helioproj01 input observer pos: " << obsbary.x << " " << obsbary.y << " " << obsbary.z << "\n";
  //cout << "barycentric unit vector: " << unitbary.x << " " << unitbary.y << " " << unitbary.z << "\n";

  barydist = sqrt(dotprod3d(obsbary,obsbary));
  obsdot = dotprod3d(unitbary,obsbary);
  // opdistcos = obsdot/barydist;
  // opelong = acos(opdistcos)*DEGPRAD;
  // if(opelong < 0.0) opelong = 90.0 - opelong;
  // sunelong = 180.0 - opelong;

  //cout << "barydist = " << barydist << ", obsdot = " << obsdot << ", opdistcos = " << opdistcos << ", opelong = " << opelong << ", sunelong = " << sunelong << "\n";
  
  // Default values, signify failure if returned.
  geodist = -1.0;
  projbary.x=projbary.y=projbary.z = 0.0;

  a = 1.0;
  b = 2.0*obsdot;
  c = barydist*barydist - heliodist*heliodist;
  
  if((b*b-4.0*a*c)<0.0) {
    // Quadratic equation for object's distance from the
    // observer has no solution: this unit vector extending
    // from Earth can never intersect this heliocentric
    // sphere. Most likely, the heliocentric sphere lies
    // entirely inside the vector from Earth.
    return(-1);
  } else
    {
      //cout << "a = " << a << ", b = " << b << ", c = " << c << "\n";
      alphapos = (-b + sqrt(b*b-4.0*a*c))/2.0/a;
      if(alphapos>0.0) {
	geodist=alphapos;
	projbary.x = obsbary.x + alphapos*unitbary.x;
	projbary.y = obsbary.y + alphapos*unitbary.y;
	projbary.z = obsbary.z + alphapos*unitbary.z;
	//cout << "Distance from Earth: " << alphapos << "\n";
	//cout << "Barycentric coordinates: " << projbary.x << " " << projbary.y << " " << projbary.z << "\n";
	return(0);
      }
      else return(-1);
    }
  return(-1);
}

// helioproj01LD: Given a unit vector unitbary giving the direction
// toward which an object was seen; a vector obsbary giving the full
// 3-D location of the observer; and a distance heliodist giving the
// object's distance from the origin (e.g., the sun or the barycenter),
// find out where the vector unitbary, extended indefinitely in the
// positive direction, intersects the sphere of radius heliodist.
// Provides the distance from the observer to the intersection point
// in geodist, and the full 3-D coordinates of the intersection point
// in projbary. Returns 0 if a valid solution was found, or -1 if
// there was no solution.
// NOTE: helioproj01LD considers only a solution where the vector is EXITING
// the sphere. Hence, for cases where the vector originates outside
// the sphere and intesects it twice, helioproj01LD will solve for
// and return only one of the two intersections.
int helioproj01LD(point3LD unitbary, point3LD obsbary, long double heliodist, long double &geodist, point3LD &projbary)
{
  long double a,b,c;
  long double alphapos,obsdot,barydist;
  // long double opdistcos,sunelong,opelong;
  
  //cout << fixed << setprecision(9) << "helioproj01 input observer pos: " << obsbary.x << " " << obsbary.y << " " << obsbary.z << "\n";
  //cout << "barycentric unit vector: " << unitbary.x << " " << unitbary.y << " " << unitbary.z << "\n";

  barydist = sqrt(dotprod3LD(obsbary,obsbary));
  obsdot = dotprod3LD(unitbary,obsbary);
  //opdistcos = obsdot/barydist;
  //opelong = acos(opdistcos)*DEGPRAD;
  //if(opelong < 0.0L) opelong = 90.0L - opelong;
  //sunelong = 180.0L - opelong;

  //cout << "barydist = " << barydist << ", obsdot = " << obsdot << ", opdistcos = " << opdistcos << ", opelong = " << opelong << ", sunelong = " << sunelong << "\n";
  
  // Default values, signify failure if returned.
  geodist = -1.0L;
  projbary.x=projbary.y=projbary.z = 0.0L;

  a = 1.0L;
  b = 2.0L*obsdot;
  c = barydist*barydist - heliodist*heliodist;
  
  if((b*b-4.0L*a*c)<0.0L) {
    // Quadratic equation for object's distance from the
    // observer has no solution: this unit vector extending
    // from Earth can never intersect this heliocentric
    // sphere. Most likely, the heliocentric sphere lies
    // entirely inside the vector from Earth.
    return(-1);
  } else
    {
      //cout << "a = " << a << ", b = " << b << ", c = " << c << "\n";
      alphapos = (-b + sqrt(b*b-4.0L*a*c))/2.0L/a;
      if(alphapos>0.0L) {
	geodist=alphapos;
	projbary.x = obsbary.x + alphapos*unitbary.x;
	projbary.y = obsbary.y + alphapos*unitbary.y;
	projbary.z = obsbary.z + alphapos*unitbary.z;
	//cout << "Distance from Earth: " << alphapos << "\n";
	//cout << "Barycentric coordinates: " << projbary.x << " " << projbary.y << " " << projbary.z << "\n";
	return(0);
      }
      else return(-1);
    }
  return(-1);
}

// helioproj02LD: June 22, 2022:
// Like helioproj01LD, but correctly solves cases where the vector
// originates outside the heliocentric sphere and interesects it twice.
// If valid solutions exist, returns the number of solutions found
// (either 1 or 2), with the solutions stored in the ouput vectors geodist
// and projbary. If no solution exists, returns -1 with output vectors
// empty.
int helioproj02LD(point3LD unitbary, point3LD obsbary, long double heliodist, vector <long double> &geodist, vector <point3LD> &projbary)
{
  long double a,b,c;
  long double alphapos,alphaneg,barydist,obsdot;
  //long double opdistcos,sunelong,opelong;
  point3LD barypos = point3LD(0.0L,0.0L,0.0L);
  
  barydist = sqrt(dotprod3LD(obsbary,obsbary));
  obsdot = dotprod3LD(unitbary,obsbary);
  //opdistcos = obsdot/barydist;
  //opelong = acos(opdistcos)*DEGPRAD;
  //if(opelong < 0.0L) opelong = 90.0L - opelong;
  //sunelong = 180.0L - opelong;

  // Make sure output vectors start out empty
  geodist = {};
  projbary = {};

  a = 1.0L;
  b = 2.0L*obsdot;
  c = barydist*barydist - heliodist*heliodist;
  
  if((b*b-4.0L*a*c)>=0.0L) {
    // Quadratic has two real solutions for geocentric distance.
    alphapos = (-b + sqrt(b*b-4.0L*a*c))/2.0L/a;
    if(alphapos>0.0L) {
      // The first solution is positive, and hence physically possible
      geodist.push_back(alphapos);
      barypos.x = obsbary.x + alphapos*unitbary.x;
      barypos.y = obsbary.y + alphapos*unitbary.y;
      barypos.z = obsbary.z + alphapos*unitbary.z;
      projbary.push_back(barypos);
      alphaneg = (-b - sqrt(b*b-4.0L*a*c))/2.0L/a;
      if(alphaneg>0.0L) {
	// The second solution is also positive, so they are both physically possible
	geodist.push_back(alphaneg);
	barypos.x = obsbary.x + alphaneg*unitbary.x;
	barypos.y = obsbary.y + alphaneg*unitbary.y;
	barypos.z = obsbary.z + alphaneg*unitbary.z;
	projbary.push_back(barypos);
	// Sanity check before we return
	if(geodist.size()==2 || projbary.size()==2) {
	  return(2);
	} else {
	  cerr << "ERROR: vector sizes don't match number of real solutions!\n";
	  cerr << "Sizes: " << geodist.size() << " " << projbary.size() << " should both be exactly 2\n";
	} 
      } else {
	// Only the first solution was positive.
	if(geodist.size()==1 || projbary.size()==1) {
	  return(1);
	} else {
	  cerr << "ERROR: vector sizes don't match number of real solutions!\n";
	  cerr << "Sizes: " << geodist.size() << " " << projbary.size() << " should both be exactly 1\n";
	} 
      }
    } else {
      // There were no positive real solutions for geocentric distance.
      // Since scalar distance cannot be negative, this means the quadratic
      // has no physical solution: the vector unitbary has no intersection
      // with the heliocentric sphere of radius heliodist in the positive direction.
      return(-1);
    }
  } else {
    // Quadratic equation for object's distance from the
    // observer has no real solution: this unit vector extending
    // from Earth can never intersect this heliocentric
    // sphere. Most likely, the heliocentric sphere lies
    // entirely inside the vector from Earth.
    return(-1);
  }
  // Should never reach this point, but return anyway just in case.
  return(-1);
}

// helioproj02: March 28, 2023
// Like helioproj02LD, but uses just double-precision, rather than
// long double. Correctly solves cases where the vector
// originates outside the heliocentric sphere and interesects it twice.
// If valid solutions exist, returns the number of solutions found
// (either 1 or 2), with the solutions stored in the ouput vectors geodist
// and projbary. If no solution exists, returns -1 with output vectors
// empty.
int helioproj02(point3d unitbary, point3d obsbary, double heliodist, vector <double> &geodist, vector <point3d> &projbary)
{
  double a,b,c;
  double alphapos,alphaneg,barydist,obsdot;
  point3d barypos = point3d(0.0l,0.0l,0.0l);
  
  barydist = sqrt(dotprod3d(obsbary,obsbary));
  obsdot = dotprod3d(unitbary,obsbary);

  // Make sure output vectors start out empty
  geodist = {};
  projbary = {};

  a = 1.0l;
  b = 2.0l*obsdot;
  c = barydist*barydist - heliodist*heliodist;
  
  if((b*b-4.0L*a*c)>=0.0l) {
    // Quadratic has two real solutions for geocentric distance.
    alphapos = (-b + sqrt(b*b-4.0L*a*c))/2.0L/a;
    if(alphapos>0.0L) {
      // The first solution is positive, and hence physically possible
      geodist.push_back(alphapos);
      barypos.x = obsbary.x + alphapos*unitbary.x;
      barypos.y = obsbary.y + alphapos*unitbary.y;
      barypos.z = obsbary.z + alphapos*unitbary.z;
      projbary.push_back(barypos);
      alphaneg = (-b - sqrt(b*b-4.0L*a*c))/2.0L/a;
      if(alphaneg>0.0L) {
	// The second solution is also positive, so they are both physically possible
	geodist.push_back(alphaneg);
	barypos.x = obsbary.x + alphaneg*unitbary.x;
	barypos.y = obsbary.y + alphaneg*unitbary.y;
	barypos.z = obsbary.z + alphaneg*unitbary.z;
	projbary.push_back(barypos);
	// Sanity check before we return
	if(geodist.size()==2 || projbary.size()==2) {
	  return(2);
	} else {
	  cerr << "ERROR: vector sizes don't match number of real solutions!\n";
	  cerr << "Sizes: " << geodist.size() << " " << projbary.size() << " should both be exactly 2\n";
	} 
      } else {
	// Only the first solution was positive.
	if(geodist.size()==1 || projbary.size()==1) {
	  return(1);
	} else {
	  cerr << "ERROR: vector sizes don't match number of real solutions!\n";
	  cerr << "Sizes: " << geodist.size() << " " << projbary.size() << " should both be exactly 1\n";
	} 
      }
    } else {
      // There were no positive real solutions for geocentric distance.
      // Since scalar distance cannot be negative, this means the quadratic
      // has no physical solution: the vector unitbary has no intersection
      // with the heliocentric sphere of radius heliodist in the positive direction.
      return(-1);
    }
  } else {
    // Quadratic equation for object's distance from the
    // observer has no real solution: this unit vector extending
    // from Earth can never intersect this heliocentric
    // sphere. Most likely, the heliocentric sphere lies
    // entirely inside the vector from Earth.
    return(-1);
  }
  // Should never reach this point, but return anyway just in case.
  return(-1);
}


// accelcalc01LD: December 01, 2021
// Given a vector of planet positions for a particular instant in time,
// calculate the resulting gravitational acceleration at the point targpos.
// The planet masses must be provided in the vector planetmasses, and must
// consist of GM products in km^3/sec^2, the default coordinate set used
// by JPL. The output acceleration will have units of km/sec^2.
int accelcalc01LD(int planetnum, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const point3LD &targpos, point3LD &accel)
{
  accel = point3LD(0L,0L,0L);
  point3LD relpos = point3LD(0L,0L,0L);
  long double distcubed=0L;
  int i=0;
  for(i=0;i<planetnum;i++)
    {
      // Calculate relative position vector pointing from target toward planet
      relpos.x = planetpos[i].x - targpos.x;
      relpos.y = planetpos[i].y - targpos.y;
      relpos.z = planetpos[i].z - targpos.z;
      distcubed = dotprod3LD(relpos,relpos);
      distcubed = distcubed*sqrt(distcubed); // planet-target distance raised to the third power.
      accel.x += planetmasses[i]*relpos.x/distcubed;
      accel.y += planetmasses[i]*relpos.y/distcubed;
      accel.z += planetmasses[i]*relpos.z/distcubed;
    }
  return(0);
}

// integrate_orbit_constac: December 07, 2021:
// Really crude and fast orbit integrator using constant
// acceleration between timesteps.
int integrate_orbit_constac(int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, long double mjdstart, point3LD startpos, point3LD startvel, long double mjdend, point3LD &endpos, point3LD &endvel)
{
  int startnear=0;
  int pointafter=0;
  int lastpoint=0;
  int i=0;
  int planetpointnum = planetmjd.size();
  point3LD singleaccel = point3LD(0L,0L,0L);
  point3LD singlevel = point3LD(0L,0L,0L);
  point3LD singlepos = point3LD(0L,0L,0L);
  vector <point3LD> planetnow;
  long double dt1=0L;
  int planetpointct=0;
  
  // Find the point just after mjdstart, and the point just before mjdend.
  for(i=0;i<planetpointnum-1;i++)
    {
      if(planetmjd[i]<=mjdstart && planetmjd[i+1]>mjdstart) pointafter=i+1;
      if(planetmjd[i]<mjdend && planetmjd[i+1]>=mjdend) lastpoint=i;
    }
  // Find the point nearest mjdstart, and the point nearest mjdend.
  if(fabs(planetmjd[pointafter]-mjdstart) < fabs(planetmjd[pointafter-1]-mjdstart)) startnear = pointafter;
  else startnear = pointafter-1;

  // Load starting position and velocity.
  singlepos = startpos;
  singlevel = startvel;
  
  // Calculate acceleration for the starting position and the nearest point.
  planetnow={};
  nplanetgrab01LD(startnear, planetnum, planetmjd, planetpos, planetnow);
  accelcalc01LD(planetnum, planetmasses, planetnow, singlepos, singleaccel);
  
  // Use that as the first constant acceleration, and integrate forward.
  if(mjdend <= planetmjd[pointafter]) {
    // This next planet file timestep is after mjdend:
    // Perform simple two-step integration.
    dt1 = (mjdend-mjdstart)*SOLARDAY;
    cout << "dt1 = " << dt1 << "\n";
    // Calculate position at mjdend.
    singlepos.x += singlevel.x*dt1 + singleaccel.x*0.5L*dt1*dt1;
    singlepos.y += singlevel.y*dt1 + singleaccel.y*0.5L*dt1*dt1;
    singlepos.z += singlevel.z*dt1 + singleaccel.z*0.5L*dt1*dt1;
    singlevel.x += singleaccel.x*dt1;
    singlevel.y += singleaccel.y*dt1;
    singlevel.z += singleaccel.z*dt1;
  } else {
    // Integration will span multiple timesteps.
    dt1 = (planetmjd[pointafter]-mjdstart)*SOLARDAY;
    for(planetpointct=pointafter; planetpointct<=lastpoint+1; planetpointct++) {
      // Calculate new position using current acceleration.
      singlepos.x += singlevel.x*dt1 + singleaccel.x*0.5L*dt1*dt1;
      singlepos.y += singlevel.y*dt1 + singleaccel.y*0.5L*dt1*dt1;
      singlepos.z += singlevel.z*dt1 + singleaccel.z*0.5L*dt1*dt1;
      singlevel.x += singleaccel.x*dt1;
      singlevel.y += singleaccel.y*dt1;
      singlevel.z += singleaccel.z*dt1;
      if(planetpointct<lastpoint) {
	// Calculate new acceleration.
	nplanetgrab01LD(planetpointct, planetnum, planetmjd, planetpos, planetnow);
	accelcalc01LD(planetnum, planetmasses, planetnow, singlepos, singleaccel);
	// Calculate new timestep normally
	dt1 = (planetmjd[planetpointct+1]-planetmjd[planetpointct])*SOLARDAY;
      } else if (planetpointct==lastpoint) {
	// Calculate new acceleration.
	nplanetgrab01LD(planetpointct, planetnum, planetmjd, planetpos, planetnow);
	accelcalc01LD(planetnum, planetmasses, planetnow, singlepos, singleaccel);
	// Calculate new timestep using the fact that we are on the last point.
	dt1 = (mjdend-planetmjd[planetpointct])*SOLARDAY;
      } else {
	// Do nothing: we are at the end of the loop.
	;
      }
    }
  }
  endpos = singlepos;
  endvel = singlevel;
  return(0);
}

// kep_transcendental: December 08 2021:
// Solve the trancendental Kepler Equation
// q = psi - e*sin(psi) for psi given q and e,
// returning a result guaranteed to be correct
// within tol, unless KEPTRANSITMAX iterations
// elapse without achieving this.
long double kep_transcendental(long double q, long double e, long double tol)
{
  int itct=0;
  long double psi_guess = M_PI; // Derivative is always positive: correct answer
                                // can always be reached starting from psi = pi.

  if(tol<=0L) {
    cerr << "ERROR: kep_trancendental called with non-positive tolerance " << tol << "\n";
    return(-99.9);
  }

  long double fpsi = psi_guess - e*sin(psi_guess) - q;
  long double fprime = 1.0L - e*cos(psi_guess);
  itct=0;
  cout.precision(17);
  while(itct<KEPTRANSITMAX && fabs(fpsi) > tol) {
    psi_guess += -fpsi/fprime;
    if(psi_guess >= 2.0L*M_PI) psi_guess = 2.0L*M_PI - tol;
    if(psi_guess < 0.0L) psi_guess = 0.0L + tol;
    fpsi = psi_guess - e*sin(psi_guess) - q;
    fprime = 1.0L - e*cos(psi_guess);
    // cout << "kep itct " << itct << "psi, fpsi, fprime: " << psi_guess << ", " << fpsi << ", " << fprime << "\n"; 
    itct++;
  }

  if(itct>=KEPTRANSITMAX) {
    cout.precision(21);
    cout << "Warning: kep_trancendental " << itct << " iters, still " << fpsi << " > tol = " << tol;
    cout << " Call was q = " << q << ", e = " << e << "\n";
  }
  // cout << "kep_transcendental obtained error of " << fpsi << " in only " << itct << " iterations\n";
  return(psi_guess);
}
    
double kep_transcendental(double q, double e, double tol)
{
  int itct=0;
  double psi_guess = M_PI; // Derivative is always positive: correct answer
                                // can always be reached starting from psi = pi.

  if(tol<=0L) {
    cerr << "ERROR: kep_trancendental called with non-positive tolerance " << tol << "\n";
    return(-99.9);
  }

  double fpsi = psi_guess - e*sin(psi_guess) - q;
  double fprime = 1.0L - e*cos(psi_guess);
  itct=0;
  cout.precision(17);
  while(itct<KEPTRANSITMAX && fabs(fpsi) > tol) {
    psi_guess += -fpsi/fprime;
    if(psi_guess >= 2.0L*M_PI) psi_guess = 2.0L*M_PI - tol;
    if(psi_guess < 0.0L) psi_guess = 0.0L + tol;
    fpsi = psi_guess - e*sin(psi_guess) - q;
    fprime = 1.0L - e*cos(psi_guess);
    // cout << "kep itct " << itct << "psi, fpsi, fprime: " << psi_guess << ", " << fpsi << ", " << fprime << "\n"; 
    itct++;
  }

  if(itct>=KEPTRANSITMAX) {
    cout.precision(21);
    cout << "Warning: kep_trancendental " << itct << " iters, still " << fpsi << " > tol = " << tol;
    cout << " Call was q = " << q << ", e = " << e << "\n";
  }
  // cout << "kep_transcendental obtained error of " << fpsi << " in only " << itct << " iterations\n";
  return(psi_guess);
}

// Keplerint: December 08, 2021:
// Integrate an orbit assuming we have a Keplerian 2-body problem
// with all the mass in the Sun, and the input position and velocity
// are relative to the Sun.
int Keplerint(const long double MGsun, const long double mjdstart, const point3LD &startpos, const point3LD &startvel, const long double mjdend, point3LD &endpos, point3LD &endvel)
{
  long double e,E,a,lscalar,r0,v0,r1,v1;
  point3LD lvec = point3LD(0L,0L,0L);
  point3LD r1unit = point3LD(0L,0L,0L);
  point3LD v1unit = point3LD(0L,0L,0L);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0L;
  long double costheta,theta0,theta1,radvel,cospsi,psi;
  costheta = theta0 = theta1 = radvel = cospsi = psi = 0L;
  long double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0L;
  long double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0L;
  long double sinev,thetav,v1ra,v1dec;
  sinev = thetav = v1ra = v1dec = 0L;
    
  // Calculate scalar input position
  r0 = sqrt(dotprod3LD(startpos,startpos));
  v0 = sqrt(dotprod3LD(startvel,startvel));
  
  // Calculate specific energy and angular momentum
  E = 0.5L*v0*v0 - MGsun/r0;
  lvec = crossprod3LD(startpos,startvel);
  lscalar = sqrt(dotprod3LD(lvec,lvec));
  if(E>=0L) {
    //cerr << "ERROR: Keplerint finds positive total energy: " << E << "\n";
    return(1);
  }
		 
  // Calculate a and e: orbital semimajor axis and eccentricity.
  cout.precision(17);
  a = -MGsun*0.5L/E;
  e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/MGsun/MGsun);
  if(e < 0L && e > -KEPTRANSTOL) {
    // eccentricity is formally negative, but very close to zero, so instead
    // of throwing an error, we simply set it to exactly zero.
    e = 0L;
  }
  if(e<0L || e>=1.0L) {
    cerr << "ERROR: Keplerint finds eccentricity out of range: " << e << "\n";
    return(1);
  }

  //cout << "semimajor axis = " << a/AU_KM << " and eccentricity = " << e << "\n";
	       
  // Calculate angle theta0 from perihelion using simple ellipse geometry.
  if(e>0L) costheta = ((a-a*e*e)/r0 - 1.0L)/e;
  else costheta = 1.0L;
  if(costheta>=-1.0L && costheta<=1.0L) theta0 = acos(costheta);
  else if(costheta<-1.0) {
    // cerr << "WARNING: Keplerint finds costheta+1.0 = " << costheta+1.0 << "\n";
    costheta = -1.0L;
    theta0 = M_PI;
  } else if(costheta>1.0) {
    // cerr << "WARNING: Keplerint finds costheta-1.0 = " << costheta-1.0 << "\n";
    costheta = 1.0L;
    theta0 = 0L;
  } else {
    cerr << "ERROR: Keplerint finds costheta = " << costheta << "\n";
    return(1);
  }
    
  radvel = dotprod3LD(startpos,startvel)/r0;
  //cout << "Radial velocity = " << radvel << " km/sec\n";
  
  if(radvel>=0) {
    // We are moving outward from perihelion: theta will be correct.
    //cout << "Moving outward from perihelion.\n";
    ;
  } else {
    // We are moving inward towards perihelion: theta needs adjustment.
    theta0 = 2.0L*M_PI - theta0;
    //cout << "Moving inward towards perihelion.\n";
  }
  // If e=0, we'll have costheta = 1.0 and theta0 = 0.0
  
  // Calculate Goldstein's psi variable from theta.
  cospsi = (costheta + e)/(1.0L + costheta*e);
  if(cospsi>=-1.0L && cospsi<=1.0L) psi = acos(cospsi);
  else if(cospsi<-1.0) {
    // cerr << "WARNING: Keplerint finds cospsi+1.0 = " << cospsi+1.0 << "\n";
    cospsi = -1.0L;
    psi = M_PI;
  } else if(cospsi>1.0) {
    // cerr << "WARNING: Keplerint finds cospsi-1.0 = " << cospsi-1.0 << "\n";
    cospsi = 1.0L;
    psi = 0L;  
  } else {
    cerr << "ERROR: Keplerint finds cospsi = " << cospsi << "\n";
    return(1);
  }
  if(radvel<0) {
    // We are moving inward towards perihelion: psi needs adjustment.
    psi = 2.0L*M_PI - psi;
  }
  // If e=0, we'll have cospsi = 1.0 and psi = 0.0

  // Calculate time since perihelion using psi.
  omega = sqrt(MGsun/(a*a*a));
  //cout << "Period = " << 2.0L*M_PI/omega/SOLARDAY/365.25 << " years\n";
  t0omega = psi - e*sin(psi);
  // If e=0, we'll have t0omega = 0.0
 
  // The new time t1 for which we want to re-evaluate psi is
  // given by t0 + mjdend-mjdstart.
  t1omega = t0omega + (mjdend-mjdstart)*SOLARDAY*omega;
  //cout << " t1omega = " << t1omega;
  while(t1omega >= 2.0L*M_PI) t1omega -= 2.0L*M_PI;
  while(t1omega < 0.0L) t1omega += 2.0L*M_PI;
  //cout << " t1omega = " << t1omega << "\n";
  // Solve Kepler's equation for psi(t1)
  psi = kep_transcendental(t1omega,e,KEPTRANSTOL);
  //cout << "New psi = " << psi*DEGPRAD;
  cospsi = cos(psi);
  //cout << " New cospsi = " << cospsi;
  // Calculate theta(t1) from psi(t1)
  if(1.0L - e*cospsi != 0.0L) {
    costheta = (cospsi - e)/(1.0L - e*cospsi);
    if(costheta >= -1.0L && costheta <= 1.0L) theta1 = acos(costheta);
    else if(costheta < -1.0L) {
      cout << "Warning: costheta = " << costheta << "\n";
      theta1 = M_PI;
    } else {
      cout << "Warning: costheta = " << costheta << "\n";
      theta1 = 0.0L;
    }
    if(psi>M_PI && theta1<=M_PI) theta1 = 2.0L*M_PI - theta1;
  } else {
    cerr << "Warning: e*cos(psi) = " << e*cospsi << " so 1 - e*cos(psi) = " << 1.0L - e*cospsi << "\n";
    theta1 = 0.0L;
  }
  while(theta1<0.0L) theta1 += 2.0L*M_PI;
  while(theta1>=2.0L*M_PI) theta1 -= 2.0L*M_PI;
  //cout << " New theta1 = " << theta1*DEGPRAD << " diff: " << theta1*DEGPRAD << " - " <<  theta0*DEGPRAD << " = " << (theta1-theta0)*DEGPRAD << "\n";

  // Calculate r(t1) from psi(t1)
  r1 = a*(1.0L - e*cospsi);
  // Calculate v1 from r1 and the known energy
  v1 = sqrt((E +  MGsun/r1)*2.0L);
  
  // Use vector algebra to find the full vector r(t1).
  // This vector is perpendicular to lvec, and is angled by theta1-theta0
  // relative to startpos.
  // Convert angular momentum vector to spherical coordinates
  celedeproj01LD(lvec,&lra,&ldec); // Note that output is in degrees.
  //cout << "Psuedo-celestial coords of angular momentum vector: " << lra << " " << ldec << "\n";
  celedeproj01LD(startpos,&r0ra,&r0dec); // Note that output is in degrees.
  //cout << "Psuedo-celestial coords of original position: " << r0ra << " " << r0dec << "\n";
  // Transform the starting unit vector into a coordinate system with
  // the angular momentum vector at the pole, and the old pole at RA=0
  poleswitch01LD(r0ra/DEGPRAD,r0dec/DEGPRAD,lra/DEGPRAD,ldec/DEGPRAD,0.0L,newra,newdec); // Output is radians
  //cout << "Orbital plane coords of original position: " << newra*DEGPRAD << " " << newdec*DEGPRAD << "\n";
  // Rotate starting unit vector around the angular momentum axis by
  // the calculated angle.
  newra += theta1-theta0;
  // cout << "Orbital plane coords of new position " << newra*DEGPRAD << " " << newdec*DEGPRAD << "\n";
  // The unit vector for the new position r1 is on the equator at this RA,
  // in the coordinate system that has the angular momentum vector at the pole.
  // Convert back to the original coordinate system.
  poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,r1ra,r1dec); // Output is radians
  // Now for the velocity. If the velocity is at right angle to the vector r1,
  // the product v1*r1 is the angular momentum. Otherwise, l/(v1*r1) is the sine
  // of the angle between v1 and r1.

  sinev = lscalar/v1/r1;
  if(sinev>=1.0L) thetav = 0.5L*M_PI;
  else if(sinev<0.0L) {
    cerr << "ERROR: negative angular momentum?\nv1,r1,v1*r1,lscalar,sinev = " << v1 << ", " << r1 << ", " << v1*r1 << ", " << lscalar << ", " << sinev << "\n";
    thetav = 0.0L;
  }
  else thetav = asin(sinev);
  if(theta1<=M_PI) {
    // Outward bound from perihelion.
    newra += thetav;
  } else {
    // Inward bound to perihelion
    newra += (M_PI - thetav);
  }
  poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,v1ra,v1dec); // Output is radians

  // cout << "Psuedo-celestial coords of new position: " << r1ra*DEGPRAD << " " << r1dec*DEGPRAD << "\n";
  r1unit = celeproj01LD(r1ra*DEGPRAD,r1dec*DEGPRAD);
  v1unit =celeproj01LD(v1ra*DEGPRAD,v1dec*DEGPRAD);
  
  endpos.x = r1unit.x*r1;
  endpos.y = r1unit.y*r1;
  endpos.z = r1unit.z*r1;
  endvel.x = v1unit.x*v1;
  endvel.y = v1unit.y*v1;
  endvel.z = v1unit.z*v1;
  
  return(0);
}

// Keplerint: March 28, 2023: Like Keplerint, but uses only
// double precision, not long double.
// Integrate an orbit assuming we have a Keplerian 2-body problem
// with all the mass in the Sun, and the input position and velocity
// are relative to the Sun.
int Keplerint(const double MGsun, const double mjdstart, const point3d &startpos, const point3d &startvel, const double mjdend, point3d &endpos, point3d &endvel)
{
  double e,E,a,lscalar,r0,v0,r1,v1;
  point3d lvec = point3d(0l,0l,0l);
  point3d r1unit = point3d(0l,0l,0l);
  point3d v1unit = point3d(0l,0l,0l);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0l;
  double costheta,theta0,theta1,radvel,cospsi,psi;
  costheta = theta0 = theta1 = radvel = cospsi = psi = 0l;
  double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0l;
  double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0l;
  double sinev,thetav,v1ra,v1dec;
  sinev = thetav = v1ra = v1dec = 0l;
    
  // Calculate scalar input position
  r0 = sqrt(dotprod3d(startpos,startpos));
  v0 = sqrt(dotprod3d(startvel,startvel));
  
  // Calculate specific energy and angular momentum
  E = 0.5l*v0*v0 - MGsun/r0;
  lvec = crossprod3d(startpos,startvel);
  lscalar = sqrt(dotprod3d(lvec,lvec));
  if(E>=0l) {
    //cerr << "ERROR: Keplerint finds positive total energy: " << E << "\n";
    return(1);
  }
		 
  // Calculate a and e: orbital semimajor axis and eccentricity.
  cout.precision(17);
  a = -MGsun*0.5l/E;
  e = sqrt(1.0l + 2.0l*E*lscalar*lscalar/MGsun/MGsun);
  if(e < 0l && e > -KEPTRANSTOL2) {
    // eccentricity is formally negative, but very close to zero, so instead
    // of throwing an error, we simply set it to exactly zero.
    e = 0l;
  }
  if(e<0l || e>=1.0l) {
    cerr << "ERROR: Keplerint finds eccentricity out of range: " << e << "\n";
    return(1);
  }

  // Calculate angle theta0 from perihelion using simple ellipse geometry.
  if(e>0l) costheta = ((a-a*e*e)/r0 - 1.0l)/e;
  else costheta = 1.0l;
  if(costheta>=-1.0l && costheta<=1.0l) theta0 = acos(costheta);
  else if(costheta<-1.0) {
    // cerr << "WARNING: Keplerint finds costheta+1.0 = " << costheta+1.0 << "\n";
    costheta = -1.0l;
    theta0 = M_PI;
  } else if(costheta>1.0) {
    // cerr << "WARNING: Keplerint finds costheta-1.0 = " << costheta-1.0 << "\n";
    costheta = 1.0l;
    theta0 = 0l;
  } else {
    cerr << "ERROR: Keplerint finds costheta = " << costheta << "\n";
    return(1);
  }
    
  radvel = dotprod3d(startpos,startvel)/r0;
  //cout << "Radial velocity = " << radvel << " km/sec\n";
  
  if(radvel>=0) {
    // We are moving outward from perihelion: theta will be correct.
    //cout << "Moving outward from perihelion.\n";
    ;
  } else {
    // We are moving inward towards perihelion: theta needs adjustment.
    theta0 = 2.0l*M_PI - theta0;
    //cout << "Moving inward towards perihelion.\n";
  }
  // If e=0, we'll have costheta = 1.0 and theta0 = 0.0
  
  // Calculate Goldstein's psi variable from theta.
  cospsi = (costheta + e)/(1.0l + costheta*e);
  if(cospsi>=-1.0l && cospsi<=1.0l) psi = acos(cospsi);
  else if(cospsi<-1.0) {
    // cerr << "WARNING: Keplerint finds cospsi+1.0 = " << cospsi+1.0 << "\n";
    cospsi = -1.0l;
    psi = M_PI;
  } else if(cospsi>1.0) {
    // cerr << "WARNING: Keplerint finds cospsi-1.0 = " << cospsi-1.0 << "\n";
    cospsi = 1.0l;
    psi = 0l;  
  } else {
    cerr << "ERROR: Keplerint finds cospsi = " << cospsi << "\n";
    return(1);
  }
  if(radvel<0) {
    // We are moving inward towards perihelion: psi needs adjustment.
    psi = 2.0l*M_PI - psi;
  }
  // If e=0, we'll have cospsi = 1.0 and psi = 0.0

  // Calculate time since perihelion using psi.
  omega = sqrt(MGsun/(a*a*a));
  //cout << "Period = " << 2.0L*M_PI/omega/SOLARDAY/365.25 << " years\n";
  t0omega = psi - e*sin(psi);
  // If e=0, we'll have t0omega = 0.0
 
  // The new time t1 for which we want to re-evaluate psi is
  // given by t0 + mjdend-mjdstart.
  t1omega = t0omega + (mjdend-mjdstart)*SOLARDAY*omega;
  //cout << " t1omega = " << t1omega;
  while(t1omega >= 2.0l*M_PI) t1omega -= 2.0l*M_PI;
  while(t1omega < 0.0l) t1omega += 2.0l*M_PI;
  //cout << " t1omega = " << t1omega << "\n";
  // Solve Kepler's equation for psi(t1)
  psi = kep_transcendental(t1omega,e,KEPTRANSTOL2);
  //cout << "New psi = " << psi*DEGPRAD;
  cospsi = cos(psi);
  //cout << " New cospsi = " << cospsi;
  // Calculate theta(t1) from psi(t1)
  if(1.0l - e*cospsi != 0.0l) {
    costheta = (cospsi - e)/(1.0l - e*cospsi);
    if(costheta >= -1.0l && costheta <= 1.0l) theta1 = acos(costheta);
    else if(costheta < -1.0l) {
      cout << "Warning: costheta = " << costheta << "\n";
      theta1 = M_PI;
    } else {
      cout << "Warning: costheta = " << costheta << "\n";
      theta1 = 0.0l;
    }
    if(psi>M_PI && theta1<=M_PI) theta1 = 2.0l*M_PI - theta1;
  } else {
    cerr << "Warning: e*cos(psi) = " << e*cospsi << " so 1 - e*cos(psi) = " << 1.0l - e*cospsi << "\n";
    theta1 = 0.0l;
  }
  while(theta1<0.0l) theta1 += 2.0l*M_PI;
  while(theta1>=2.0l*M_PI) theta1 -= 2.0l*M_PI;

  // Calculate r(t1) from psi(t1)
  r1 = a*(1.0L - e*cospsi);
  // Calculate v1 from r1 and the known energy
  v1 = sqrt((E +  MGsun/r1)*2.0L);
  
  // Use vector algebra to find the full vector r(t1).
  // This vector is perpendicular to lvec, and is angled by theta1-theta0
  // relative to startpos.
  // Convert angular momentum vector to spherical coordinates
  celedeproj01(lvec,&lra,&ldec); // Note that output is in degrees.
  //cout << "Psuedo-celestial coords of angular momentum vector: " << lra << " " << ldec << "\n";
  celedeproj01(startpos,&r0ra,&r0dec); // Note that output is in degrees.
  //cout << "Psuedo-celestial coords of original position: " << r0ra << " " << r0dec << "\n";
  // Transform the starting unit vector into a coordinate system with
  // the angular momentum vector at the pole, and the old pole at RA=0
  poleswitch01(r0ra/DEGPRAD,r0dec/DEGPRAD,lra/DEGPRAD,ldec/DEGPRAD,0.0L,newra,newdec); // Output is radians
  //cout << "Orbital plane coords of original position: " << newra*DEGPRAD << " " << newdec*DEGPRAD << "\n";
  // Rotate starting unit vector around the angular momentum axis by
  // the calculated angle.
  newra += theta1-theta0;
  // cout << "Orbital plane coords of new position " << newra*DEGPRAD << " " << newdec*DEGPRAD << "\n";
  // The unit vector for the new position r1 is on the equator at this RA,
  // in the coordinate system that has the angular momentum vector at the pole.
  // Convert back to the original coordinate system.
  poleswitch01(newra,0.0l,0.0l,ldec/DEGPRAD,lra/DEGPRAD,r1ra,r1dec); // Output is radians
  // Now for the velocity. If the velocity is at right angle to the vector r1,
  // the product v1*r1 is the angular momentum. Otherwise, l/(v1*r1) is the sine
  // of the angle between v1 and r1.

  sinev = lscalar/v1/r1;
  if(sinev>=1.0l) thetav = 0.5l*M_PI;
  else if(sinev<0.0l) {
    cerr << "ERROR: negative angular momentum?\nv1,r1,v1*r1,lscalar,sinev = " << v1 << ", " << r1 << ", " << v1*r1 << ", " << lscalar << ", " << sinev << "\n";
    thetav = 0.0l;
  }
  else thetav = asin(sinev);
  if(theta1<=M_PI) {
    // Outward bound from perihelion.
    newra += thetav;
  } else {
    // Inward bound to perihelion
    newra += (M_PI - thetav);
  }
  poleswitch01(newra,0.0l,0.0l,ldec/DEGPRAD,lra/DEGPRAD,v1ra,v1dec); // Output is radians

  // cout << "Psuedo-celestial coords of new position: " << r1ra*DEGPRAD << " " << r1dec*DEGPRAD << "\n";
  r1unit = celeproj01(r1ra*DEGPRAD,r1dec*DEGPRAD);
  v1unit =celeproj01(v1ra*DEGPRAD,v1dec*DEGPRAD);
  
  endpos.x = r1unit.x*r1;
  endpos.y = r1unit.y*r1;
  endpos.z = r1unit.z*r1;
  endvel.x = v1unit.x*v1;
  endvel.y = v1unit.y*v1;
  endvel.z = v1unit.z*v1;
  
  return(0);
}


// Kepler2dyn: May 31, 2022:
// Given Keplerian orbital parameters and a starting MJD,
// convert the Keplerian orbital parameters to barycentric
// Cartesian state vectors at the instant of the MJD.
int Kepler2dyn(const long double mjdnow, const keplerian_orbit &keporb, point3LD &outpos,  point3LD &outvel)
{
  long double meananom,theta,psi,cospsi,costheta;
  long double heliodist,ellipsearea,Period,sweeprate;
  long double radvel,angvel,tanvel,poleRA,poleDec,oldpoleRA;
  long double newRA,newDec,totalvel,thetaoc;
  long double vtheta1,vtheta2,posRA,posDec,velRA,velDec;
  vtheta1 = vtheta2 = posRA = posDec = velRA = velDec = 0.0L;
  // keporb.semimaj_axis        in AU
  // keporb.eccentricity        unitless
  // keporb.inclination         in degrees
  // keporb.long_ascend_node    Longitude of the ascending node, in degrees
  // keporb.arg_perihelion      Argument of perihelion, in degrees
  // keporb.mean_anom           Mean anomaly at the epoch, in degrees
  // keporb.mjd_epoch           Epoch for the orbit in MJD
  // keporb.mean_daily_motion   in degrees/day

  meananom = keporb.mean_anom + (mjdnow-keporb.mjd_epoch)*keporb.mean_daily_motion;
  //cout << "keporb.meananom = " << keporb.mean_anom << " meananom = " << meananom << "\n";
  
  // Solve Kepler's equation for psi (the true anomaly) given the mean anomaly.
  psi = kep_transcendental(meananom/DEGPRAD,keporb.eccentricity,KEPTRANSTOL);
  //cout << "New psi = " << psi*DEGPRAD;
  cospsi = cos(psi);
  //cout << " New cospsi = " << cospsi << "\n";
  // Calculate theta from psi
  if(1.0L - keporb.eccentricity*cospsi != 0.0L) {
    costheta = (cospsi - keporb.eccentricity)/(1.0L - keporb.eccentricity*cospsi);
    if(costheta >= -1.0L && costheta <= 1.0L) theta = acos(costheta);
    else if(costheta < -1.0L) {
      cout << "Warning: costheta = " << costheta << "\n";
      theta = M_PI;
    } else {
      cout << "Warning: costheta = " << costheta << "\n";
      theta = 0.0L;
    }
    if(psi>M_PI && theta<=M_PI) theta = 2.0L*M_PI - theta;
  } else {
    cerr << "Warning: e*cos(psi) = " << keporb.eccentricity*cospsi << " so 1 - e*cos(psi) = " << 1.0L - keporb.eccentricity*cospsi << "\n";
    theta = 0.0L;
  }
  while(theta<0.0L) theta += 2.0L*M_PI;
  while(theta>=2.0L*M_PI) theta -= 2.0L*M_PI;

  // Calculate heliodist from psi
  heliodist = keporb.semimaj_axis*(1.0L - keporb.eccentricity*cospsi);
  //cout << "keporb.semimaj_axis = " << keporb.semimaj_axis << " keporb.eccentricity = " << keporb.eccentricity << " heliodist = " << heliodist <<"\n";
  // Now effectively we have the asteroid's position fully
  //  specified in a coordinate system for which the asteroid's
  //  orbit defines the equatorial plane and its perihelion defines
  //  zero longitude. The current longitude is theta1, the 
  //  latitude is zero by definition, and heliodist
  //  is the radius. Two steps remain: (1) Calculate the velocity
  //  in this orbital coordinate system, and (2) tranform both
  //  position and velocity into heliocentric coordinates.

  // Calculate the Velocity
  ellipsearea = M_PI*keporb.semimaj_axis*keporb.semimaj_axis*sqrt(1.0L - LDSQUARE(keporb.eccentricity));
  // This area is in AU^2. The period of the orbit is:
  Period = 360.0L/keporb.mean_daily_motion;
  //cout << "Period = " << Period << " days\n";
  // This is the period in days, because keporb.mean_daily_motion is the
  //  mean daily motion in degrees/day
  sweeprate = ellipsearea/Period;
  // This is the rate at which area is swept out, in terms
  //  of AU^2/day. This will allow us to find the instantaneous
  //  angular velocity
  // Area of triangle swept out in time dt is r^2 * dt * angvel / 2.0
  //cout << "heliodist = " << heliodist << "\n";
  angvel = sweeprate*2.0L/LDSQUARE(heliodist);
  //cout << "angvel = " << angvel << " rad/day = " << angvel*DEGPRAD << " deg/day = " << angvel*DEGPRAD*150.0L << "arcsec/hr\n";
  // This is the angular velocity in radians/day 
  // All we need now is the radial component. 
  radvel = angvel*keporb.semimaj_axis*(1.0L - LDSQUARE(keporb.eccentricity))*(keporb.eccentricity*sin(theta)/LDSQUARE(1.0L + keporb.eccentricity*cos(theta)));
  // This is the radial velocity in AU/day. The formula is 
  //  derived in my October 05, 2016 notebook entry.
  tanvel = angvel*heliodist;
  //cout << "radvel, tanvel = " << radvel << " " << tanvel << "\n";
  // Calculate the angle of the velocity relative
  //  to the position vector
  if(tanvel>0.0) vtheta1 = M_PI/2.0L - atan(radvel/tanvel);
  else if(tanvel==0.0 && radvel>=0.0) vtheta1 = 0.0;
  else if(tanvel==0.0 && radvel<0.0) vtheta1 = M_PI;
  else if(tanvel<0.0) vtheta1 = 3.0L*M_PI/2.0L - atan(radvel/tanvel);

  // Convert to degrees
  vtheta1*=DEGPRAD;
  // Add in the angle of the position vector
  vtheta1 += theta*DEGPRAD;
  // Calculate the total velocity
  totalvel = sqrt(radvel*radvel + tanvel*tanvel);
  // Note that at this point, vtheta1 is in degrees
  //  and total vel is in AU/day.

  // Reckon from the line of nodes (intersection of the
  //  asteroidal orbit with the plane of the ecliptic)
  
  // Position vector:
  thetaoc = keporb.arg_perihelion+theta*DEGPRAD;
  while(thetaoc>=360.0L) thetaoc-=360.0L;
  // Velocity vector
  vtheta2 = keporb.arg_perihelion+vtheta1;
  while(vtheta2>=360.0L) vtheta2-=360.0L;

  //cout << "The angle from perihelion is " << theta*DEGPRAD << " degrees.\n";
  //cout << "The angle from the line of nodes is " << thetaoc << " degrees.\n";
  // Now in orbital coordinates, the asteroid's position
  //  has a 'declination' of zero (by definition) and a
  //  'right ascension' equal to thetaoc
  // RA in orbital coords
  posRA = thetaoc/DEGPRAD;
  velRA = vtheta2/DEGPRAD;
  // Dec in orbital coords
  posDec = 0.0L;
  velDec = 0.0L;
  // Now I need the orbital coordinates of the ecliptic pole.
  //  Since I have defined the line of nodes as the RA reference
  //  in orbital coordinates, and this line has to be perpendicular
  //  to the vector to the ecliptic pole, it follows that the
  //  orbital RA of the ecliptic pole is 90.0 degrees.
  //  The orbital declination of the ecliptic pole is ninety
  //  degrees minus the inclination.
  poleRA = M_PI/2.0L;
  poleDec = (90.0L - keporb.inclination)/DEGPRAD;
  // Now I need the right ascension of the old pole
  //  in the new coordinates.  This is equal to the
  //  longitude of the ascending node minus 90 degrees.
  oldpoleRA = (keporb.long_ascend_node - 90.0L)/DEGPRAD;
  poleswitch01LD(posRA,posDec,poleRA,poleDec,oldpoleRA,newRA,newDec); // Output is radians
  outpos.x = heliodist*cos(newDec)*cos(newRA);
  outpos.y = heliodist*cos(newDec)*sin(newRA);
  outpos.z = heliodist*sin(newDec);
  poleswitch01LD(velRA,velDec,poleRA,poleDec,oldpoleRA,newRA,newDec); // Output is radians  
  outvel.x = totalvel*cos(newDec)*cos(newRA);
  outvel.y = totalvel*cos(newDec)*sin(newRA);
  outvel.z = totalvel*sin(newDec);
  
  return(0);
}


// hyp_transcendental: April 25, 2022:
// Solve the hyperbolic form of the trancendental
// Kepler Equation q = e*sinh(psi) - psi for psi given q and e,
// returning a result guaranteed to be correct
// within tol, unless KEPTRANSITMAX iterations
// elapse without achieving this.
long double hyp_transcendental(long double q, long double e, long double tol)
{
  int itct=0;
  long double psi_guess = M_PI;
  
  if(tol<=0L) {
    cerr << "ERROR: hyp_trancendental called with non-positive tolerance " << tol << "\n";
    return(-99.9);
  }

  if(q>=0) psi_guess = 3.0;
  else psi_guess = -3.0;
  
  long double fpsi = e*sinh(psi_guess) - psi_guess - q;
  long double fprime = e*cosh(psi_guess) - 1.0L;
  itct=0;
  cout.precision(17);
  while(itct<KEPTRANSITMAX && fabs(fpsi) > tol) {
    psi_guess += -fpsi/fprime;
    fpsi = e*sinh(psi_guess) - psi_guess - q;
    fprime = e*cosh(psi_guess) - 1.0L;
    // cout << "kep itct " << itct << "psi, fpsi, fprime: " << psi_guess << ", " << fpsi << ", " << fprime << "\n"; 
    itct++;
  }

  if(itct>=KEPTRANSITMAX) {
    cout.precision(21);
    cout << "Warning: hyp_trancendental " << itct << " iters, still " << fpsi << " > tol = " << tol;
    cout << " Call was q = " << q << ", e = " << e << "\n";
  }
  // cout << "hyp_transcendental obtained error of " << fpsi << " in only " << itct << " iterations\n";
  return(psi_guess);
}
    
// Hyper_Kepint: April 25, 2022:
// Integrate a hyperbolic orbit assuming we have a Keplerian 2-body problem
// with all the mass in the Sun, and the input position and velocity
// are relative to the Sun.
int Hyper_Kepint(const long double MGsun, const long double mjdstart, const point3LD &startpos, const point3LD &startvel, const long double mjdend, point3LD &endpos, point3LD &endvel)
{
  long double e,E,a,lscalar,r0,v0,r1,v1;
  point3LD lvec = point3LD(0L,0L,0L);
  point3LD r1unit = point3LD(0L,0L,0L);
  point3LD v1unit = point3LD(0L,0L,0L);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0L;
  long double coshH,H0,theta0,theta1,radvel,H;
  coshH = H0 = theta0 = theta1 = radvel = H = 0L;
  long double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0L;
  long double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0L;
  long double x,y,junkra,junkdec,sinev,thetav,v1ra,v1dec;
  x = y = junkra = junkdec = sinev = thetav = v1ra = v1dec = 0L;
  int debug=0;
  
  // Calculate scalar input position
  r0 = sqrt(dotprod3LD(startpos,startpos));
  v0 = sqrt(dotprod3LD(startvel,startvel));
  
  // Calculate specific energy and angular momentum
  E = 0.5L*v0*v0 - MGsun/r0;
  lvec = crossprod3LD(startpos,startvel);
  lscalar = sqrt(dotprod3LD(lvec,lvec));
  if(debug>=2) cout << "Energy = " << E << ", angmom = " << lscalar << "\n";
  if(E<=0L) {
    //cerr << "ERROR: Hyper_Kepint finds negative total energy: " << E << "\n";
    return(1);
  }
		 
  // Calculate a and e: orbital semimajor axis and eccentricity.
  cout.precision(17);
  a = -MGsun*0.5L/E; // By convention, semimajor axis a is negative
                     // for hyperbolic orbits.
  e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/MGsun/MGsun);
  if(debug>=2) cout << "a = " << a/AU_KM << ", e = " << e << "\n";
  
  if(e<=1.0L) {
    cerr << "ERROR: Hyper_Kepint finds eccentricity out of range: " << e << "\n";
    return(1);
  }

  // Calculate the value of the hyperbolic anomaly H at the starting time
  if(e>1.0L) {
    coshH = (a-r0)/(a*e);
    if(coshH>=1.0L) H0 = acosh(coshH);
    else {
      cerr << "ERROR: Hyper_Kepint finds cosh(H) = " << coshH << "\n";
      return(1);
    }
  }
  radvel = dotprod3LD(startpos,startvel)/r0;
  if(debug>=2) cout << "H0 = " << H0 << ", radial velocity = " << radvel << " km/sec\n";
  
  if(radvel>=0) {
    // We are moving outward from perihelion: H0 will be correct
     ;
  } else {
    // We are moving inward towards perihelion: H0 needs adjustment.
    H0 = -H0;
   }
  
  if(debug>=2) cout << "H0 = " << H0 << ", radial velocity = " << radvel << " km/sec\n";

  // Calculate the angle theta0 from perihelion using H0.
  x = a*(cosh(H0)-e);
  y = -a*sqrt(e*e-1)*sinh(H0);
  if(y>0.0L) {
    theta0 = M_PI/2.0L - atan(x/y);
  } else if (y<0.0L) {
    theta0 = 3.0L*M_PI/2.0L - atan(x/y);
  } else {
    // Presumably y==0. There's also the possibility that
    // it could be a NAN, but we won't worry about that right here.
    if(x<0) theta0 = M_PI;
    else theta0 = 0.0L;
  }
  if(debug>=2) cout << "x = " << x/AU_KM << ", y = " << y/AU_KM << ", theta0 = " << theta0*M_PI/180.0 << "\n";

  // Calculate time since perihelion using H0.
  omega = sqrt(-MGsun/(a*a*a));
  //cout << "omega = " << omega << "\n";
  t0omega = e*sinh(H0) - H0;
  //cout << "t0omega = " << t0omega << "\n";
 
  // The new time t1 for which we want to re-evaluate psi is
  // given by t0 + mjdend-mjdstart.
  t1omega = t0omega + (mjdend-mjdstart)*SOLARDAY*omega;
  //cout << "t1omega = " << t1omega << "\n";
  // Solve the hyperbolic form of Kepler's equation for H(t1)
  H = hyp_transcendental(t1omega,e,KEPTRANSTOL);
  //cout << "H = " << H << "\n";

  // Calculate r(t1) from H(t1)
  r1 = a*(1.0L - e*cosh(H));
  // Calculate v1 from r1 and the known energy
  v1 = sqrt((E +  MGsun/r1)*2.0L);
  
  // Calculate the angle theta1 from perihelion from H(t1).
  x = a*(cosh(H)-e);
  y = -a*sqrt(e*e-1)*sinh(H);
  if(y>0.0L) {
    theta1 = M_PI/2.0L - atan(x/y);
  } else if (y<0.0L) {
    theta1 = 3.0L*M_PI/2.0L - atan(x/y);
  } else {
    // Presumably y==0. There's also the possibility that
    // it could be a NAN, but we won't worry about that right here.
    if(x<0) theta1 = M_PI;
    else theta1 = 0.0L;
  }

  // Use vector algebra to find the full vector r(t1).
  // This vector is perpendicular to lvec, and is angled by theta1-theta0
  // relative to startpos.
  // Convert angular momentum vector to spherical coordinates
  celedeproj01LD(lvec,&lra,&ldec); // Note that output is in degrees.
  //cout << "Psuedo-celestial coords of angular momentum vector: " << lra << " " << ldec << "\n";
  celedeproj01LD(startpos,&r0ra,&r0dec); // Note that output is in degrees.
  //cout << "Psuedo-celestial coords of original position: " << r0ra << " " << r0dec << "\n";
  // Transform the starting unit vector into a coordinate system with
  // the angular momentum vector at the pole, and the old pole at RA=0
  poleswitch01LD(r0ra/DEGPRAD,r0dec/DEGPRAD,lra/DEGPRAD,ldec/DEGPRAD,0.0L,newra,newdec); // Output is radians
  //cout << "Orbital plane coords of original position: " << newra*DEGPRAD << " " << newdec*DEGPRAD << "\n";
  // Rotate starting unit vector around the angular momentum axis by
  // the calculated angle.
  newra += theta1-theta0;
  // cout << "Orbital plane coords of new position " << newra*DEGPRAD << " " << newdec*DEGPRAD << "\n";
  // The unit vector for the new position r1 is on the equator at this RA,
  // in the coordinate system that has the angular momentum vector at the pole.
  // Convert back to the original coordinate system.
  poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,r1ra,r1dec); // Output is radians
  // Now for the velocity. If the velocity is at right angle to the vector r1,
  // the product v1*r1 is the angular momentum. Otherwise, l/(v1*r1) is the sine
  // of the angle between v1 and r1.

  sinev = lscalar/v1/r1;
  if(sinev>=1.0L) thetav = 0.5L*M_PI;
  else if(sinev<0.0L) {
    cerr << "ERROR: negative angular momentum?\nv1,r1,v1*r1,lscalar,sinev = " << v1 << ", " << r1 << ", " << v1*r1 << ", " << lscalar << ", " << sinev << "\n";
    thetav = 0.0L;
  }
  else thetav = asin(sinev);
  if(theta1<=M_PI) {
    // Outward bound from perihelion.
    newra += thetav;
  } else {
    // Inward bound to perihelion
    newra += (M_PI - thetav);
  }
  poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,v1ra,v1dec); // Output is radians

  // cout << "Psuedo-celestial coords of new position: " << r1ra*DEGPRAD << " " << r1dec*DEGPRAD << "\n";
  r1unit = celeproj01LD(r1ra*DEGPRAD,r1dec*DEGPRAD);
  v1unit =celeproj01LD(v1ra*DEGPRAD,v1dec*DEGPRAD);
  
  endpos.x = r1unit.x*r1;
  endpos.y = r1unit.y*r1;
  endpos.z = r1unit.z*r1;
  endvel.x = v1unit.x*v1;
  endvel.y = v1unit.y*v1;
  endvel.z = v1unit.z*v1;
  
  return(0);
}

// integrate_orbit01LD: December 01, 2021
// Uses the approximation of linearly varying acceleration
// to integrate the orbit of a massless test particle (e.g. asteroid)
// under the gravity of multiple 'planets'. It is assumed that
// in general these 'planets' will consist of the Sun, the
// eight major planets, and the Moon (possibly needed for
// cases of NEOs closely approaching the Earth). However,
// more or fewer planets may be used as desired.
int integrate_orbit01LD(int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, long double mjdstart, point3LD startpos, point3LD startvel, long double mjdend, point3LD &endpos, point3LD &endvel)
{
  vector <point3LD> planetstart;
  vector <point3LD> planetend;
  vector <point3LD> targaccel;
  vector <point3LD> targvel;
  vector <point3LD> targpos;
  point3LD singleaccel = point3LD(0L,0L,0L);
  point3LD singlevel = point3LD(0L,0L,0L);
  point3LD singlepos = point3LD(0L,0L,0L);
  int i=0;
  int planetpointnum = planetmjd.size();
  int planetpointct = 0;
  int pointafter=0;
  int lastpoint=0;
  long double dt1=0L;
  point3LD accelslope = point3LD(0L,0L,0L);

  if(mjdend<mjdstart) {
    cerr << "ERROR: integrate_orbit01LD called with end time (" << mjdend << ") before start time (" << mjdstart << ")\n";
    return(1);
  }
  else if(mjdstart<=planetmjd[1] || mjdend>=planetmjd[planetpointnum-1]) {
    cerr << "ERROR: integrate_orbit01LD called with time range (" << mjdstart << "-" << mjdend << ") outside range of planet vectors (" << planetmjd[1] << "-" << planetmjd[planetpointnum-1] << ")\n";
    return(1);
  }
  planetstart={};
  nplanetpos01LD(mjdstart,planetnum,5,planetmjd,planetpos,planetstart);
  
  // Calculate acceleration at starting point.
  accelcalc01LD(planetnum, planetmasses, planetstart, startpos, singleaccel);
  targaccel.push_back(singleaccel);
  targvel.push_back(startvel);
  targpos.push_back(startpos);

  if(targaccel.size()!=1) {
    cerr << "ERROR: nplanetpos01LD targaccel vector has " << targaccel.size() << " entries, should be 1\n";
    return(2);
  }
  if(targvel.size()!=1) {
    cerr << "ERROR: nplanetpos01LD targvel vector has " << targvel.size() << " entries, should be 1\n";
    return(2);
  }
  if(targpos.size()!=1) {
    cerr << "ERROR: nplanetpos01LD targpos vector has " << targpos.size() << " entries, should be 1\n";
    return(2);
  }

  // Find next entry in planetmjd vector.
  for(i=0;i<planetpointnum;i++)
    {
      if(planetmjd[i]>mjdstart) break;
    }
  pointafter=i; // This is the first planet file time step after mjdstart.
  if(mjdend <= planetmjd[pointafter]) {
    // This next planet file timestep is after mjdend:
    // Perform simple two-step integration.
    dt1 = (mjdend-mjdstart)*SOLARDAY;
    cout << "dt1 = " << dt1 << "\n";
    // First Approx: constant acceleration.
    singlepos.x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1;
    singlepos.y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1;
    singlepos.z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1;
    // Calculate accleration at this new position.
    planetend={};
    nplanetpos01LD(mjdend,planetnum,5,planetmjd,planetpos,planetend);
    accelcalc01LD(planetnum, planetmasses, planetend, singlepos, singleaccel);
    targaccel.push_back(singleaccel);
    if(targaccel.size()!=2) {
      cerr << "ERROR: nplanetpos01LD targaccel vector has " << targaccel.size() << " entries, should be 2\n";
      return(2);
    }
    // Better approx: linear variation of acceleration in time
    accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
    accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
    accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;
    cout << "Accel in m/sec^2: " << targaccel[0].x*1000.0 << " "<< targaccel[0].y*1000.0 << " " << targaccel[0].z*1000.0 << "\n";
    endpos.x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
    endpos.y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
    endpos.z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

    endvel.x = targvel[0].x + targaccel[0].x*dt1 + accelslope.x*0.5L*dt1*dt1;
    endvel.y = targvel[0].y + targaccel[0].y*dt1 + accelslope.y*0.5L*dt1*dt1;
    endvel.z = targvel[0].z + targaccel[0].z*dt1 + accelslope.z*0.5L*dt1*dt1;
  } else { 
    // The desired endpoint mjdend is more than a fractional timestep
    // away from the starting position.
    // Find the point just before the endpoint:
    // this is as far as our integration will go.
    for(i=0;i<planetpointnum-1;i++)
      {
      if(planetmjd[i+1]>=mjdend) break;
      }
    lastpoint = i;
    // cout << "Integrating from point " << pointafter << " to " << lastpoint << "\n";
    // cout << "MJD from " << planetmjd[pointafter] << " to " << planetmjd[lastpoint] << "\n"; 
    dt1 = (planetmjd[pointafter]-mjdstart)*SOLARDAY;
    for(planetpointct=pointafter;planetpointct<=lastpoint+1;planetpointct++)
      {
	// Regard input acceleration targaccel[0] as perfect and inviolable,
	// and produce the best possible value for the next time step.
	// First Approx: constant acceleration.
	singlepos.x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1;
	singlepos.y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1;
	singlepos.z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1;
	// Calculate acceleration at this new position.
	planetend={};
	if(planetpointct<=lastpoint) nplanetgrab01LD(planetpointct, planetnum, planetmjd, planetpos, planetend);
	else nplanetpos01LD(mjdend,planetnum,5,planetmjd,planetpos,planetend);
	accelcalc01LD(planetnum, planetmasses, planetend, singlepos, singleaccel);
	if(targaccel.size()<2) targaccel.push_back(singleaccel);
	else targaccel[1] = singleaccel;

	// Better approx: linear variation of acceleration in time
	accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
	accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
	accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;
    
	singlepos.x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
	singlepos.y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
	singlepos.z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

	// Re-calculate acceleration
	accelcalc01LD(planetnum, planetmasses, planetend, singlepos, singleaccel);
	targaccel[1] = singleaccel;

	// Re-calculate position and velocity.
	accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
	accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
	accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;

	singlepos.x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
	singlepos.y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
	singlepos.z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

	singlevel.x = targvel[0].x + targaccel[0].x*dt1 + accelslope.x*0.5L*dt1*dt1;
	singlevel.y = targvel[0].y + targaccel[0].y*dt1 + accelslope.y*0.5L*dt1*dt1;
	singlevel.z = targvel[0].z + targaccel[0].z*dt1 + accelslope.z*0.5L*dt1*dt1;

	if(targpos.size()<2) targpos.push_back(singlepos);
	else  targpos[1] = singlepos;
	if(targvel.size()<2) targvel.push_back(singlevel);
	else targvel[1] = singlevel;
	// Cycle target vectors for the next time step.
	targpos[0] = targpos[1];
	targvel[0] = targvel[1];
	targaccel[0] = targaccel[1];
	// Re-set dt1 for the next step.
	if(planetpointct<lastpoint) dt1 = (planetmjd[planetpointct+1] - planetmjd[planetpointct])*SOLARDAY;
	else dt1 = (mjdend - planetmjd[planetpointct])*SOLARDAY;
      }
    endpos = targpos[0];
    endvel = targvel[0];
  }

  return(0);
}

// integrate_orbit02LD: December 01, 2021
// Uses modeling of the acceleration as a polynomial of order n>1
// to integrate the orbit of a massless test particle (e.g. asteroid)
// under the gravity of multiple 'planets'. It is assumed that
// in general these 'planets' will consist of the Sun, the
// eight major planets, and the Moon (possibly needed for
// cases of NEOs closely approaching the Earth). However,
// more or fewer planets may be used as desired.
int integrate_orbit02LD(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, long double mjdstart, point3LD startpos, point3LD startvel, long double mjdend, point3LD &endpos, point3LD &endvel)
{
  vector <point3LD> planetsalltimes;
  vector <point3LD> planetsonce;
  vector <point3LD> targaccel;
  vector <point3LD> accelfit;
  vector <point3LD> targvel;
  vector <point3LD> targpos;
  vector <point3LD> accelmod;
  vector <long double> temptime;
  vector <long double> ppxvec;
  vector <long double> ppyvec;
  vector <long double> ppfitvec;
  point3LD singleaccel = point3LD(0L,0L,0L);
  point3LD singlevel = point3LD(0L,0L,0L);
  point3LD singlepos = point3LD(0L,0L,0L);
  int i=0;
  int j=0;
  int k=0;
  int endhere=-1;
  int planetpointnum = planetmjd.size();
  int pointafter=0;
  int latestpoint=0;
  int stepsin=0;
  long double dt0=0L;
  long double dt1=0L;
  long double dt2=0L;
  long double timemult=0L;
  point3LD accelslope = point3LD(0L,0L,0L);

  if(polyorder<2) {
    cerr << "ERROR: integrate_orbit02LD called with polyorder = " << polyorder << "\n";
    cerr << "polyorder must be at least 2!\n";
    return(1);
  }
  
  if(mjdend<mjdstart) {
    cerr << "ERROR: integrate_orbit02LD called with end time (" << mjdend << ") before start time (" << mjdstart << ")\n";
    return(1);
  }
  else if(mjdstart<=planetmjd[1] || mjdend>=planetmjd[planetpointnum-1]) {
    cerr << "ERROR: integrate_orbit02LD called with time range (" << mjdstart << "-" << mjdend << ") outside range of planet vectors (" << planetmjd[1] << "-" << planetmjd[planetpointnum-1] << ")\n";
    return(1);
  }

  // Make sure that relevant vectors for the polynomial fitting
  // are all large enough.
  for(i=0;i<=polyorder+1;i++) {
    targaccel.push_back(singleaccel);
    accelfit.push_back(singleaccel);
    targvel.push_back(singlevel);
    targpos.push_back(singlepos);
    accelmod.push_back(singleaccel);
    temptime.push_back(0L);
    ppfitvec.push_back(0L);
  }

  // Load the initial time and planet position vectors
  temptime[0] = mjdstart;
  planetsonce={};
  nplanetpos01LD(temptime[0],planetnum,5,planetmjd,planetpos,planetsonce);
  for(i=0;i<planetnum;i++) planetsalltimes.push_back(planetsonce[i]);
 
  for(i=0;i<planetpointnum;i++) {
    if(planetmjd[i]>mjdstart) break;
  }
  pointafter = i; // first point after mjdstart
  dt0 = (planetmjd[pointafter+1] - planetmjd[pointafter])*SOLARDAY/sqrt(M_PI);
  j=1;
  i=0;
  while(j<=polyorder+1)
    {
      if(planetmjd[pointafter+i]<mjdend || endhere>=0) {
	temptime[j] = planetmjd[pointafter+i];
	planetsonce={};	
	nplanetgrab01LD(pointafter+i, planetnum, planetmjd, planetpos, planetsonce);
	for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
	latestpoint = pointafter+i;
	j++;
	i++;
      } else if (planetmjd[pointafter+i] == mjdend) {
	// Weird case where the requested mjdend falls exactly on a timestep
	temptime[j] = mjdend;
	if(j<=polyorder) endhere=j;
	planetsonce={};
	nplanetpos01LD(temptime[j],planetnum,5,planetmjd,planetpos,planetsonce);
	for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
	j++;
	i++; // Must advance i to avoid two identical times in tempvec.
      } else {
	// More typical case where mjdend falls somewhere in betweeen two timesteps.
	temptime[j] = mjdend;
	if(j<=polyorder) endhere=j;
	planetsonce={};
	nplanetpos01LD(temptime[j],planetnum,5,planetmjd,planetpos,planetsonce);
	for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
	j++;
	// Don't advance i: we don't want to skip the next regular timestep.
      }
    }
  // Load starting position and velocity
  targvel[0] = startvel;
  targpos[0] = startpos;
  // Bootstrap up to a fit of order polyorder.
  // Calculate acceleration at starting point, loading planet positions from big vector.
  for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[i];
  accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[0], targaccel[0]);
  dt1 = (temptime[1]-temptime[0])*SOLARDAY;
  // First Approx: estimate next position, assuming constant acceleration.
  targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1;
  targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1;
  targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1;
  // Calculate acceleration at this new position.
  for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*1 + i];
  accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);
  
  // Second approx: linearly varying acceleration.
  accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
  accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
  accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;

  // Improved position for next time step.
  targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
  targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
  targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

  // Re-calculate acceleration at this improved position.
  accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);

  // Re-calculate improved acceleration slope.
  accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
  accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
  accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;
  
  // Improved velocity for next time step
  targvel[1].x = targvel[0].x + targaccel[0].x*dt1 + accelslope.x*0.5L*dt1*dt1;
  targvel[1].y = targvel[0].y + targaccel[0].y*dt1 + accelslope.y*0.5L*dt1*dt1;
  targvel[1].z = targvel[0].z + targaccel[0].z*dt1 + accelslope.z*0.5L*dt1*dt1;

  // Use linearly extrapolated acceleration to estimate position for
  // the next time step.
  dt1 = (temptime[2]-temptime[1])*SOLARDAY;
  targpos[2].x = targpos[1].x + targvel[1].x*dt1 + targaccel[1].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
  targpos[2].y = targpos[1].y + targvel[1].y*dt1 + targaccel[1].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
  targpos[2].z = targpos[1].z + targvel[1].z*dt1 + targaccel[1].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

  // Calculate acceleration for this extrapolated position.
  for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*2 + i];
  accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[2], targaccel[2]);

  // Now we have three acceleration points: can load for a full polynomial fit.
  for(stepsin=3;stepsin<=polyorder+1;stepsin++) {
    // Fit for x component of acceleration.
    ppxvec={};
    ppyvec={};
    for(i=0;i<stepsin;i++) {
      ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
      ppyvec.push_back(targaccel[i].x);
    }
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
    // Fit for y component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
    // Fit for z component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
    // Re-calculate all of the positions and velocities using this fit.
    for(j=1;j<=stepsin;j++) {
      dt2 = (temptime[j]-temptime[0])*SOLARDAY;
      // Positions
      targpos[j].x = targpos[0].x + targvel[0].x*dt2;
      targpos[j].y = targpos[0].y + targvel[0].y*dt2;
      targpos[j].z = targpos[0].z + targvel[0].z*dt2;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	targpos[j].x += accelfit[i].x*timemult;
	targpos[j].y += accelfit[i].y*timemult;
	targpos[j].z += accelfit[i].z*timemult;
      }
      // Velocities
      targvel[j].x = targvel[0].x;
      targvel[j].y = targvel[0].y;
      targvel[j].z = targvel[0].z;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	targvel[j].x += accelfit[i].x*timemult;
	targvel[j].y += accelfit[i].y*timemult;
	targvel[j].z += accelfit[i].z*timemult;
      }
      // Accelerations
      accelmod[j].x = 0L;
      accelmod[j].y = 0L;
      accelmod[j].z = 0L;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	accelmod[j].x += accelfit[i].x*timemult;
	accelmod[j].y += accelfit[i].y*timemult;
	accelmod[j].z += accelfit[i].z*timemult;
      }
    }
    // Re-calculate accelerations using these revised positions
    for(j=1;j<=stepsin;j++) {
      for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
      accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
    }
    cout.precision(17);
  
    // Perform new fits to revised accelerations
    // Fit for x component of acceleration.
    ppxvec={};
    ppyvec={};
    for(i=0;i<stepsin;i++) {
      ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
      ppyvec.push_back(targaccel[i].x);
    }
   // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
    // Fit for y component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
   // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
    // Fit for z component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
    for(i=0;i<stepsin;i++) {
    }
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
    // Re-calculate all of the positions and velocities using this fit.
    for(j=1;j<=stepsin;j++) {
      dt2 = (temptime[j]-temptime[0])*SOLARDAY;
      // Positions
      targpos[j].x = targpos[0].x + targvel[0].x*dt2;
      targpos[j].y = targpos[0].y + targvel[0].y*dt2;
      targpos[j].z = targpos[0].z + targvel[0].z*dt2;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	targpos[j].x += accelfit[i].x*timemult;
	targpos[j].y += accelfit[i].y*timemult;
	targpos[j].z += accelfit[i].z*timemult;
      }
      // Velocities
      targvel[j].x = targvel[0].x;
      targvel[j].y = targvel[0].y;
      targvel[j].z = targvel[0].z;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	targvel[j].x += accelfit[i].x*timemult;
	targvel[j].y += accelfit[i].y*timemult;
	targvel[j].z += accelfit[i].z*timemult;
      }
      // Accelerations
      accelmod[j].x = 0L;
      accelmod[j].y = 0L;
      accelmod[j].z = 0L;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	accelmod[j].x += accelfit[i].x*timemult;
	accelmod[j].y += accelfit[i].y*timemult;
	accelmod[j].z += accelfit[i].z*timemult;
      }
    }
    // Re-calculate accelerations using these revised positions
    for(j=1;j<=stepsin;j++) {
      for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
      accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
    }
  }
  // We are now set up for a full-order polynomial integration.
  // First, account for the possibility that the desired endpoint
  // has already been calculated.
  if(endhere>=0) {
      endpos = targpos[endhere];
      endvel = targvel[endhere];
      return(0);
  }
  // If we reach this point, proceed with the full polynomial integration.
  while(endhere<0) {
    // Cycle the dynamical vectors
    for(i=0;i<polyorder+1;i++) {
      temptime[i] = temptime[i+1];
      targaccel[i] = targaccel[i+1];
      targvel[i] = targvel[i+1];
      targpos[i] = targpos[i+1];
      for(j=0;j<planetnum;j++) planetsalltimes[planetnum*i + j] = planetsalltimes[planetnum*(i+1) + j];
    }
    // Load a new point into the planet and time vectors
    if(planetmjd[latestpoint+1]<mjdend) {
      // This is just a regular integration step.
      latestpoint+=1;
      temptime[polyorder+1] = planetmjd[latestpoint];
      planetsonce={};	
      nplanetgrab01LD(latestpoint, planetnum, planetmjd, planetpos, planetsonce);
      for(j=0;j<planetnum;j++) planetsalltimes[planetnum*(polyorder+1) + j] = planetsonce[j];
    }
    else {
      // We have arrived at the requested endpoint.
      temptime[polyorder+1] = mjdend;
      endhere = polyorder+1;
      planetsonce={};
      nplanetpos01LD(temptime[polyorder+1],planetnum,5,planetmjd,planetpos,planetsonce);
      for(j=0;j<planetnum;j++) planetsalltimes[planetnum*(polyorder+1) + j] = planetsonce[j];
    }
    // Fit for acceleration
    // x component of acceleration.
    ppxvec={};
    ppyvec={};
    for(i=0;i<polyorder+1;i++) {
      ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
      ppyvec.push_back(targaccel[i].x);
    }
   // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
    // Fit for y component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
    // Fit for z component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
    // Re-calculate all of the positions and velocities using this fit.
    for(j=1;j<=polyorder+1;j++) {
      dt2 = (temptime[j]-temptime[0])*SOLARDAY;
      // Positions
      targpos[j].x = targpos[0].x + targvel[0].x*dt2;
      targpos[j].y = targpos[0].y + targvel[0].y*dt2;
      targpos[j].z = targpos[0].z + targvel[0].z*dt2;
      for(i=0;i<polyorder+1;i++) {
	timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	targpos[j].x += accelfit[i].x*timemult;
	targpos[j].y += accelfit[i].y*timemult;
	targpos[j].z += accelfit[i].z*timemult;
      }
      // Velocities
      targvel[j].x = targvel[0].x;
      targvel[j].y = targvel[0].y;
      targvel[j].z = targvel[0].z;
      for(i=0;i<polyorder+1;i++) {
	timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	targvel[j].x += accelfit[i].x*timemult;
	targvel[j].y += accelfit[i].y*timemult;
	targvel[j].z += accelfit[i].z*timemult;
      }
      // Accelerations
      accelmod[j].x = 0L;
      accelmod[j].y = 0L;
      accelmod[j].z = 0L;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	accelmod[j].x += accelfit[i].x*timemult;
	accelmod[j].y += accelfit[i].y*timemult;
	accelmod[j].z += accelfit[i].z*timemult;
      }
    }
    // Re-calculate accelerations using these revised positions
    for(j=1;j<=polyorder+1;j++) {
      cout << "RECALCULATING ACCELERATIONS!\n";
      for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
      accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
    }

    // Use these revised accelerations to re-do the fits
    // x component of acceleration.
    ppxvec={};
    ppyvec={};
    for(i=0;i<polyorder+1;i++) {
      ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
      ppyvec.push_back(targaccel[i].x);
    }
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
    // Fit for y component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
   // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
    // Fit for z component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
    // Re-calculate all of the positions and velocities using this fit.
    for(j=1;j<=polyorder+1;j++) {
      dt2 = (temptime[j]-temptime[0])*SOLARDAY;
      // Positions
      targpos[j].x = targpos[0].x + targvel[0].x*dt2;
      targpos[j].y = targpos[0].y + targvel[0].y*dt2;
      targpos[j].z = targpos[0].z + targvel[0].z*dt2;
      for(i=0;i<polyorder+1;i++) {
	timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	targpos[j].x += accelfit[i].x*timemult;
	targpos[j].y += accelfit[i].y*timemult;
	targpos[j].z += accelfit[i].z*timemult;
      }
      // Velocities
      targvel[j].x = targvel[0].x;
      targvel[j].y = targvel[0].y;
      targvel[j].z = targvel[0].z;
      for(i=0;i<polyorder+1;i++) {
	timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	targvel[j].x += accelfit[i].x*timemult;
	targvel[j].y += accelfit[i].y*timemult;
	targvel[j].z += accelfit[i].z*timemult;
      }
      // Accelerations
      accelmod[j].x = 0L;
      accelmod[j].y = 0L;
      accelmod[j].z = 0L;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	accelmod[j].x += accelfit[i].x*timemult;
	accelmod[j].y += accelfit[i].y*timemult;
	accelmod[j].z += accelfit[i].z*timemult;
      }
    }
    // Re-calculate accelerations using these revised positions
    for(j=1;j<=polyorder+1;j++) {
      for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
      accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
    }
    // We have now gone through two iterations of extrapolation
    // to predict the next acceleration point as accurately as possible.
    // The next step of the loop will move the extrapolated point back
    // by one step, and use it to start extrapolating a new point,
    // at the same time refining the former extrapolated points.
    
    if(endhere>=0) {
      endpos = targpos[endhere];
      endvel = targvel[endhere];
    }
  }
  return(0);
}

// iswhitespace: December 2021
// Does the input integer c correspond to a whitespace character?
int iswhitespace(int c)
{
  if(c==' ' || c=='\t' || c=='\r' || c=='\n' || c=='\v' || c=='\f') return(1);
  else return(0);
}

// readconfigLD: December 2021
// Read a single long double parameter from a file stream, where
// the calling function guarantees that every line in the input file
// stream is either a comment line with # as the first character,
// or else the line we mean to read, which begins with the single
// long double parameter to be read. The point is to enable reading
// a configuration file where any number of explanatory comment lines
// may precede the desired parameter.
int readconfigLD(ifstream &instream1, long double *ldval)
{
  string lnfromfile;
  string stest;
  unsigned int i=0;
  int c = '0';
  *ldval = 0L;
  getline(instream1,lnfromfile);
  if(instream1.eof()) {
    // EOF: should not happen because the file should adhere to the
    // specifications given in the calling function, and hence the
    // calling function should stop before it hits the end of the file.
    return(-1);
  } else if(instream1.fail()) {
    return(-2); // Worse problem.
  } else if(instream1.bad()) {
    return(-3); // Even worse.
  } else if(lnfromfile[0]=='#') {
    // It was a comment line. Not an error, but signal
    // the detection of a comment in the return.
    return(1);
  } else {
    // Apparently a valid line.
    c='0';
    i=0;
    while(i<lnfromfile.size() && !iswhitespace(c) && c!=EOF) {
      c=lnfromfile[i];
      if(!iswhitespace(c) && c!=EOF) stest.push_back(c);
      i++;
    }
    *ldval = stold(stest);
    return(0);
  }
}

// readconfigd: December 2021
// Read a single double-precision parameter from a file stream, where
// the calling function guarantees that every line in the input file
// stream is either a comment line with # as the first character,
// or else the line we mean to read, which begins with the single
// double-precision parameter to be read. The point is to enable reading
// a configuration file where any number of explanatory comment lines
// may precede the desired parameter.
int readconfigd(ifstream &instream1, double *dval)
{
  string lnfromfile;
  string stest;
  unsigned int i=0;
  int c = '0';
  *dval = 0L;
  getline(instream1,lnfromfile);
  if(instream1.eof()) {
    // EOF: should not happen because the file should adhere to the
    // specifications given in the calling function, and hence the
    // calling function should stop before it hits the end of the file.
    return(-1);
  } else if(instream1.fail()) {
    return(-2); // Worse problem.
  } else if(instream1.bad()) {
    return(-3); // Even worse.
  } else if(lnfromfile[0]=='#') {
    // It was a comment line. Not an error, but signal
    // the detection of a comment in the return.
    return(1);
  } else {
    // Apparently a valid line.
    c='0';
    i=0;
    while(i<lnfromfile.size() && !iswhitespace(c) && c!=EOF) {
      c=lnfromfile[i];
      if(!iswhitespace(c) && c!=EOF) stest.push_back(c);
      i++;
    }
    *dval = stod(stest);
    return(0);
  }
}

// readconfigd: December 2021
// Read a single integer parameter from a file stream, where
// the calling function guarantees that every line in the input file
// stream is either a comment line with # as the first character,
// or else the line we mean to read, which begins with the single
// integer parameter to be read. The point is to enable reading
// a configuration file where any number of explanatory comment lines
// may precede the desired parameter.
int readconfigint(ifstream &instream1, int *ival)
{
  string lnfromfile;
  string stest;
  unsigned int i=0;
  int c = '0';
  *ival = 0;
  getline(instream1,lnfromfile);
  if(instream1.eof()) {
    // EOF: should not happen because the file should adhere to the
    // specifications given in the calling function, and hence the
    // calling function should stop before it hits the end of the file.
    return(-1);
  } else if(instream1.fail()) {
    return(-2); // Worse problem.
  } else if(instream1.bad()) {
    return(-3); // Even worse.
  } else if(lnfromfile[0]=='#') {
    // It was a comment line. Not an error, but signal
    // the detection of a comment in the return.
    return(1);
  } else {
    // Apparently a valid line.
    c='0';
    i=0;
    while(i<lnfromfile.size() && !iswhitespace(c) && c!=EOF) {
      c=lnfromfile[i];
      if(!iswhitespace(c) && c!=EOF) stest.push_back(c);
      i++;
    }
    *ival = stoi(stest);
    return(0);
  }
}

// Read a single string parameter from a file stream, where
// the calling function guarantees that every line in the input file
// stream is either a comment line with # as the first character,
// or else the line we mean to read, which begins with the single
// string parameter to be read. The point is to enable reading
// a configuration file where any number of explanatory comment lines
// may precede the desired parameter.
int readconfigstring(ifstream &instream1, string &sval)
{
  string lnfromfile;
  string stest;
  unsigned int i=0;
  int c = '0';
  sval={};
  getline(instream1,lnfromfile);
  if(instream1.eof()) {
    // EOF: should not happen because the file should adhere to the
    // specifications given in the calling function, and hence the
    // calling function should stop before it hits the end of the file.
    return(-1);
  } else if(instream1.fail()) {
    return(-2); // Worse problem.
  } else if(instream1.bad()) {
    return(-3); // Even worse.
  } else if(lnfromfile[0]=='#') {
    // It was a comment line. Not an error, but signal
    // the detection of a comment in the return.
    return(1);
  } else {
    // Apparently a valid line.
    c='0';
    i=0;
    while(i<lnfromfile.size() && !iswhitespace(c) && c!=EOF) {
      c=lnfromfile[i];
      if(!iswhitespace(c) && c!=EOF) stest.push_back(c);
      i++;
    }
    sval = stest;
    return(0);
  }
}

// read_accel_fileLD: December 16, 2021:
// Read a file giving heliocentric distance, radial velocity,
// normalization, and acceleration. This file is expected to
// have a one-line header, marked as a non-data line by the fact
// that it begins with #. There follow any number of lines
// whose first for columns are (1) heliocentric distance,
// (2) heliocentric radial velocity, (3) normalization, and
// (4) heliocentric acceleration. Additional columns are expected
// but are not read or used. Lines with distance or
// normalization equal to (or less than) zero are ignored
// as invalid but do not produce errors. The distances,
// velocities, and accelerations for all valid lines are
// output in the vectors heliodist, heliovel, and helioacc,
// which are expected to be empty when the function is called.
//
// Note that the first input files for this function were produced
// by the CCode program Kepler_dyn11.c. The acceleration values
// were obtained by averaging the actual acceleration for
// every instance where any known asteroid 
// was found in the specified bin of distance and radial velocity
// over a several-year integration with 1-day sampling.
int read_accel_fileLD(string accelfile, vector <long double> &heliodist, vector <long double> &heliovel, vector <long double> &helioacc)
{
  string lnfromfile;
  string stest;
  char c = '0';
  long double dist,vel,norm,acc;
  dist = vel = norm = acc = 0L;

  ifstream instream1 {accelfile};
  if(!instream1) {
    cerr << "ERROR: can't open input acceleration file " << accelfile << "\n";
    return(1);
  }
 
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    // Read first character in the current line.
    instream1 >> c;
    if(c=='#') {
      // Comment line in file: skip to the end of the line.
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) getline(instream1,lnfromfile);
    } else if (!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      // Put the character back
      instream1.unget();
      // Read distance, velocity, normalization, and acceleration.
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) instream1 >> dist;
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) instream1 >> vel;
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) instream1 >> norm;
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) instream1 >> acc;
      // Skip the rest of the line.
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) getline(instream1,lnfromfile);
       // Finished reading. velocity and acceleration
      // are allowed to be zero or negative, but distance
      // and normalization must be strictly positive.
      if(dist>0.0L && norm >0.0L)
	{
	  heliodist.push_back(dist);
	  heliovel.push_back(vel);
	  helioacc.push_back(acc);
	}
    }
  }
  if(instream1.eof()) {
    return(0); // Reached end of file, fine.
  } else if(instream1.fail()) {
    return(2); // Some problem.
  } else if(instream1.bad()) {
    return(3); // Worse problem.
  }
  // if we get here, we probably finished reading OK.
  return(0);
}

// read_longitude_fileLD: August 25, 2022:
// Read a file giving the angular velocity and acceleration in
// terms of the heliocentric ecliptic longitude.
// This file is expected to have a one-line header, marked as
// a non-data line by the fact that it begins with #. There
// follow any number of lines whose first two columns are
// (1) velocity in ecliptic longitude (deg/day), and (2)
// accleration in ecliptic longitude (deg/day^2).
int read_longitude_fileLD(string accelfile, vector <long double> &longitude_vel, vector <long double> &longitude_acc)
{
  string lnfromfile;
  string stest;
  char c = '0';
  long double vel,acc;
  vel = acc = 0L;

  ifstream instream1 {accelfile};
  if(!instream1) {
    cerr << "ERROR: can't open input acceleration file " << accelfile << "\n";
    return(1);
  }
 
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    // Read first character in the current line.
    instream1 >> c;
    if(c=='#') {
      // Comment line in file: skip to the end of the line.
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) getline(instream1,lnfromfile);
    } else if (!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      // Put the character back
      instream1.unget();
      // Read velocity and acceleration.
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) instream1 >> vel;
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) instream1 >> acc;
      // Skip the rest of the line.
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) getline(instream1,lnfromfile);
       // Finished reading.
      longitude_vel.push_back(vel);
      longitude_acc.push_back(acc);
    }
  }
  if(instream1.eof()) {
    return(0); // Reached end of file, fine.
  } else if(instream1.fail()) {
    return(2); // Some problem.
  } else if(instream1.bad()) {
    return(3); // Worse problem.
  }
  // if we get here, we probably finished reading OK.
  return(0);
}

// weight_posvel_rms: December 16, 2021:
// Given vectors containing the positions and velocities
// for a set of 6-D state vectors (e.g. a cluster of heliocentric
// 'arrows' produced by heliolinc), calculate the RMS relative
// to the cluster mean (i.e., the STDEV) of x, y, z, vz, vy, and vz.
// Convert the velocity STDEV values to match the positional ones
// in terms of units and scale by multiplying them by a characteristic
// time dtime. Store the 3 positional and 3 scaled velocity
// STDEV values in the vector rmsvec. Return the overall STDEV:
// that is, the quadrature sum of the six individual values
// with the velocity terms pre-weighted via multiplication by dtime.
long double weight_posvel_rms(const vector <point3LD> &poscluster,const vector <point3LD> &velcluster,const long double dtime, vector <long double> &rmsvec)
{
  unsigned int pnum = poscluster.size();
  unsigned int i=0;
  long double norm,x,y,z,vx,vy,vz;
  x = y = z = vx = vy = vz = 0L;
  norm = pnum;
  long double xrms, yrms, zrms, vxrms, vyrms, vzrms;
  xrms = yrms = zrms = vxrms = vyrms = vzrms = 0L;

  cout << "weight_posvel_rms pnum = " << pnum << " " << norm << "\n"; 
  
  if(pnum <= 0) {
    cerr << "ERROR: weight_posvel_rms called with no input points!\n";
    return(-1.0L);
  } else if(velcluster.size() != pnum) {
    cerr << "ERROR: weight_posvel_rms finds input position and velocity vectors\n";
    cerr << "do not agree in size!\n";
    return(-1.0L);
  }
  for(i=0;i<pnum;i++) {
    x += poscluster[i].x;
    y += poscluster[i].y;
    z += poscluster[i].z;
    vx += velcluster[i].x;
    vy += velcluster[i].y;
    vz += velcluster[i].z;
  }

  x/=norm;
  y/=norm;
  z/=norm;
  vx/=norm;
  vy/=norm;
  vz/=norm;
  
  for(i=0;i<pnum;i++) {
    xrms += LDSQUARE(poscluster[i].x - x);
    yrms += LDSQUARE(poscluster[i].y - y);
    zrms += LDSQUARE(poscluster[i].z - z);
    vxrms += LDSQUARE(velcluster[i].x - vx);
    vyrms += LDSQUARE(velcluster[i].y - vy);
    vzrms += LDSQUARE(velcluster[i].z - vz);
  }
  
  cout << "weight_posvel_rms pnum = " << pnum << " " << norm <<  " " << xrms <<  " " << yrms <<  " " << zrms <<  " " << vxrms <<  " " << vyrms <<  " " << vzrms << "\n"; 
  
  xrms = sqrt(xrms/norm);
  yrms = sqrt(yrms/norm);
  zrms = sqrt(zrms/norm);
  vxrms = sqrt(vxrms/norm)*dtime; // Convert velocity to position
  vyrms = sqrt(vyrms/norm)*dtime; // through multiplication by a
  vzrms = sqrt(vzrms/norm)*dtime; // characteristic timescale.

  
  rmsvec.push_back(xrms);
  rmsvec.push_back(yrms);
  rmsvec.push_back(zrms);
  rmsvec.push_back(vxrms);
  rmsvec.push_back(vyrms);
  rmsvec.push_back(vzrms);

  return(sqrt(xrms*xrms + yrms*yrms + zrms*zrms + vxrms*vxrms + vyrms*vyrms + vzrms*vzrms));
}


// linfituw01: January 20, 2022
// Simple and crude utility program, does an unweighted
// linear fit of the form y = mx * b, for m = slope, b = intercept
int linfituw01(const vector <double> &x, const vector <double> &y, double &slope, double &intercept)
{
  int i;
  int pointnum = x.size();
  double delta,xal,yal,xty,xsq,nsum;

  if(pointnum<=1) {
    cerr << "ERROR: linfituw01 CALLED WITH ONLY ONE POINT\n";
    return(1);
  }

  xal = yal = xty = xsq = nsum = 0.0;
  for(i=0;i<pointnum;i++) {
    xal += x[i];
    yal += y[i];
    xsq += x[i]*x[i];
    xty += x[i]*y[i];
    nsum += 1.0l;
  }
  delta = nsum*xsq - xal*xal;
  if(delta==0.0) {
    cerr << "ERROR: linfituw01 has non-finite slope\n";
    return(1);
  }
  intercept = (xsq*yal - xal*xty)/delta;
  slope = (nsum*xty - xal*yal)/delta;

  return(0);
}

// linfit01: January 03, 2023
// Simple utility program, does an weighted
// linear fit of the form y = mx * b, for m = slope, b = intercept
int linfit01(const vector <double> &x, const vector <double> &y, const vector <double> &yerr, double &slope, double &intercept)
{
  int i;
  int pointnum = x.size();
  double delta,xal,yal,xty,xsq,nsum;

  if(pointnum<=1) {
    cerr << "ERROR: linfit01 CALLED WITH ONLY ONE POINT\n";
    return(1);
  }

  xal = yal = xty = xsq = nsum = 0.0;
  for(i=0;i<pointnum;i++) {
    xal += x[i]/DSQUARE(yerr[i]);
    yal += y[i]/DSQUARE(yerr[i]);
    xsq += x[i]*x[i]/DSQUARE(yerr[i]);
    xty += x[i]*y[i]/DSQUARE(yerr[i]);
    nsum += 1.0l/DSQUARE(yerr[i]);
  }
  delta = nsum*xsq - xal*xal;
  if(delta==0.0) {
    cerr << "ERROR: linfituw01 has non-finite slope\n";
    return(1);
  }
  intercept = (xsq*yal - xal*xty)/delta;
  slope = (nsum*xty - xal*yal)/delta;

  return(0);
}

// multilinfit01: October 15, 2019, translated from C on Jan 03, 2023
// Finds the unweighted least-squares fit modeling the input vector
// yvec (length pnum) as a linear combination of fitnum other
// vectors supplied in the matrix xmat (size fitnum x pnum). The
// vector of best-fit coefficients for xmat is given in avec.*/
int multilinfit01(const vector <double> &yvec, const vector <vector <double>> &xmat, int pnum, int fitnum, vector <double> &avec, int verbose)
{
  double fitpar=0l;
  vector <double> fitvec;
  vector <vector <double>> fitmat;
  int pct,fitct,k;
  
  // Load fitmat for input into solvematrix01
  for(fitct=0;fitct<fitnum;fitct++)
    {
      // First the constant term -- that is, the term that does not
      //  multiply any of the fitting coefficients -- which is also
      // the only term that involves yvec
      fitpar=0l;
      fitvec={};
      for(pct=0;pct<pnum;pct++) fitpar -= yvec[pct]*xmat[fitct][pct];
      fitvec.push_back(fitpar);
      /*Now the actual coefficients*/
      for(k=0;k<fitnum;k++) { 
	fitpar = 0.0;
	for(pct=0;pct<pnum;pct++) fitpar += xmat[fitct][pct]*xmat[k][pct];
	fitvec.push_back(fitpar);
      }
      fitmat.push_back(fitvec);
    }
  solvematrix01(fitmat,fitnum,avec,verbose);
  return(0);
}

// multilinfit02: June 10, 2022, translated from C on Jan 03, 2023
// Finds the WEIGHTED least-squares fit modeling the input vector
// yvec (length pnum) as a linear combination of fitnum other
// vectors supplied in the matrix xmat (size fitnum x pnum). The
// vector of best-fit coefficients for xmat is given in avec.*/
int multilinfit02(const vector <double> &yvec, const vector <double> &sigvec, const vector <vector <double>> &xmat, int pnum, int fitnum, vector <double> &avec, int verbose)
{
  double fitpar;
  vector <double> fitvec;
  vector <vector <double>> fitmat;
  int pct,fitct,k;
  
  /*Load fitmat for input into solvematrix01*/
  for(fitct=0;fitct<fitnum;fitct++)
    {
      /*First the constant term -- that is, the term that does not
        multiply any of the fitting coefficients -- which is also
        the only term that involves yvec*/
      fitpar = 0l;
      fitvec={};
      for(pct=0;pct<pnum;pct++) {
	if(isnormal(sigvec[pct])) fitpar -= yvec[pct]*xmat[fitct][pct]/DSQUARE(sigvec[pct]);
      }
      fitvec.push_back(fitpar);
      /*Now the actual coefficients*/
      for(k=0;k<fitnum;k++) {
	fitpar = 0l;
	for(pct=0;pct<pnum;pct++) {
	  if(isnormal(sigvec[pct])) fitpar += xmat[fitct][pct]*xmat[k][pct]/DSQUARE(sigvec[pct]);
	}
	fitvec.push_back(fitpar);
      }
      fitmat.push_back(fitvec);
    }
  solvematrix01(fitmat,fitnum,avec,verbose);

  return(0);
}

// multilinfit02b: January 20, 2023, supposed to be a faster version
// of multilinfit02, because it does not re-calculate redundant
// terms in the matrix, and because it takes an input variance
// vector so as not to have to do any internal squaring.
// Also for speed, does not call isnormal(). Hence, the calling function
// is responsible for not inputting NANs or points with zero error.
// Finds the WEIGHTED least-squares fit modeling the input vector
// yvec (length pnum) as a linear combination of fitnum other
// vectors supplied in the matrix xmat (size fitnum x pnum). The
// vector of best-fit coefficients for xmat is given in avec.*/
int multilinfit02b(const vector <double> &yvec, const vector <double> &varvec, const vector <vector <double>> &xmat, int pnum, int fitnum, vector <double> &avec, int verbose)
{
  double fitpar;
  vector <vector <double>> fitmat;
  int pct,fitct,k;
  
  make_dmat(fitnum,fitnum+1,fitmat);
  
  /*Load fitmat for input into solvematrix01*/
  for(fitct=0;fitct<fitnum;fitct++)
    {
      /*First the constant term -- that is, the term that does not
        multiply any of the fitting coefficients -- which is also
        the only term that involves yvec*/
      fitpar = 0l;
       for(pct=0;pct<pnum;pct++) {
	 fitpar -= yvec[pct]*xmat[fitct][pct]/varvec[pct];
      }
       fitmat[fitct][0] = fitpar;
      /*Now the actual coefficients*/
      for(k=fitct;k<fitnum;k++) {
	fitpar = 0l;
	for(pct=0;pct<pnum;pct++) {
	  fitpar += xmat[fitct][pct]*xmat[k][pct]/varvec[pct];
	}
	fitmat[fitct][k+1] = fitpar;
      }
    }
  // Fill in the lower diagonal
  for(fitct=0;fitct<fitnum;fitct++) {
    for(k=0;k<fitct;k++) {
      fitmat[fitct][k+1] = fitmat[k][fitct+1];
    }
  }
   
  solvematrix01(fitmat,fitnum,avec,verbose);
  return(0);
}

  
// arc2cel01: September 09, 2020
// Given a central point and the arc distance and celestial
// position angle to second point, calculate the celestial 
// coordinates of this second point. All input and output
// quantities are in degrees. Note that this is the reverse
// process of, e.g. distradec02, which finds
// the position angle and arc distance between two points
// on the celestial sphere. The current program finds the
// celestial coordinates of the second point, given the
// first point, and the arc distance and position angle.
// NOTE WELL: here the arc dist is in degrees.
int arc2cel01(double racenter,double deccenter,double dist,double pa,double &outra,double &outdec)
{
  double colat1,tpa,rpa,arc,coscolat,colat2;
  double cosdra,deltaRA;

  tpa=pa;
  while(tpa>=360.0l) tpa-=360.0l;
  while(tpa<0.0l) tpa+=360.0l;

  // Handle trivial cases
  if(dist==0.0l)
    {
      outra = racenter;
      outdec = deccenter;
      return(0);
    }
  else if(dist==180.0l)
    {
      outra = racenter + 180.0l;
      if(outra >= 360.0l) outra -= 360.0l;
      outdec = -deccenter;
      return(0);
    }
  else if(deccenter==90.0l)
    {
      cerr << "WARNING: arc2cel01 called with starting point at north pole!\n";
      outra = tpa;
      outdec = 90.0l - dist;
      return(0);
    }
  else if(deccenter==-90.0l)
    {
      cerr << "WARNING: arc2cel01 called with starting point at south pole!\n";
      outra = tpa;
      outdec = -90.0l + dist;
      return(0);
    }

  colat1 = M_PI/(double)2.0l - deccenter/DEGPRAD;
  rpa = tpa/DEGPRAD;
  arc = dist/DEGPRAD;

  coscolat = cos(arc)*cos(colat1) + sin(arc)*sin(colat1)*cos(rpa);
  if(coscolat>1.0l) {
    if(WARN_INVERSE_TRIG>0) cerr << "WARNING: arc2cel01 attempting to take arccos of 1 + " << coscolat-1.0l << "\n";
    colat2 = 0.0l;
  } else colat2 = acos(coscolat);
  outdec = 90.0l - colat2*DEGPRAD;
  if(sin(colat2)<=0.0l)
    {
      outra = 0.0l;
      return(0);
    }
  cosdra = (cos(arc) - cos(colat1)*cos(colat2)) / (sin(colat1)*sin(colat2));
  if(cosdra>1.0l) {
    if(WARN_INVERSE_TRIG>0) cerr  << "WARNING: arc2cel01 attempting to take arccos of 1 + " << cosdra-1.0l << "\n";
    deltaRA = 0.0l;
  } else deltaRA = acos(cosdra)*DEGPRAD;

  // Direction of RA change
  if(tpa<=180.0l)
    {
      // Change is to the east
      outra = racenter + deltaRA;
      while(outra>=360.0l) outra-=360.0l;
      return(0);
    }
  else
    {
      // Change is to the west
      outra = racenter - deltaRA;
      while(outra<0.0l) outra+=360.0l;
      return(0);
    }
  return(0);
}
    

// obscode_lookup: March 09, 2022:
// Look up an observatory code from a list, and copy the
// coordinates to obslon, plxcos, and plxsin.
int obscode_lookup(const vector <observatory> &observatory_list, const char* obscode, double &obslon, double &plxcos,double &plxsin)
{
  int i=0;
  int nobs = observatory_list.size();
  for(i=0; i<nobs; i++) {
    // cout << "obscode_lookup comparing " << obscode << " with " << observatory_list[i].obscode << ":";
    if(stringnmatch01(observatory_list[i].obscode,obscode,3)==0) {
      // cout << " it\'s a match!\n";
      obslon = observatory_list[i].obslon;
      plxcos = observatory_list[i].plxcos;
      plxsin = observatory_list[i].plxsin;
      return(0);
    } // else cout << " not a match\n";
  }
  cerr << "ERROR: observatory " << obscode << " not found in list\n";
  return(1);
}

// intzero01i: March 11, 2022:
// Given an input integer i, add leading zeros as
// needed to fill out a string of length n. For example,
// intzero01i(9,4) will produce the string "0009".
// Works on negative integers, producing one fewer leading
// zeros than with an otherwise-identical positive number,
// with a negative sign in place of the first zero.
string intzero01i(const int i, const int n)
{
  int itemp=i;
  string outstring;
  int leadzero=n-1;
  int isneg=0;
  int j=1;
  if(i<0) {
    isneg=1;
    itemp=-i;
    leadzero=n-2;
  }

  // Load a string with the positive expression of i,
  // with no leading zeros.
  stringstream ss;
  ss << itemp;
  string str = ss.str();
  j=1;
  while(j<itemp && leadzero>0) {
    j*=10;
    if(j<=itemp) leadzero--;
  }

  if(isneg==1) outstring.push_back('-');
  for(j=0;j<leadzero;j++) outstring.push_back('0');
  for(j=0;j<int(str.size());j++) outstring.push_back(str[j]);
  return(outstring);
}

// get_col_vector01: Given a line from a file, containing fields
// delimited by spaces, tabs, or commas, load a string vector with
// the strings corresponding to each individual column, and return
// the total number of columns that were found. Note that only one
// space, comma, or tab should be used between successive fields.
// If, e.g., two spaces are used, get_col_vector01 will interpret
// them as delimiting a field of zero length. The reason it is
// designed to work this way is that in csv files, consecutive commas
// are in fact frequently designed to indicate null fields.
int get_col_vector01(const string &lnfromfile, vector <string> &outvec)
{
  unsigned int i=0;
  int stringct=0;
  char c='0';
  string outstring;

  if(lnfromfile.size()<=0) return(-1);

  outvec={};
  c = lnfromfile[i];
  while(i<lnfromfile.size()) {
    outstring=""; //Setup to read the string from the next column
    while(i<lnfromfile.size() && c!=',' && c!=' ' && c!='\t' && c!='\n' && c!=EOF) {
      if(c!=',' && c!=' ' && c!='\t' && c!='\n' && c!=EOF) outstring.push_back(c);
      i++;
      c=lnfromfile[i];
    }
    stringct++;
    outvec.push_back(outstring);
    i++;
    c=lnfromfile[i];
  }
  return(stringct);
}

// MJD2mpcdate: March 16, 2022: C++ version of old C code for
// converting an MJD to MPC date format.
//
// MJDtoMPCdate01: August 10, 2015: convert MJD to MPC-formatted
// date: that is, year, month, and then decimal day.
// Date must be after the year 1900. January 1, 1900 appears to
// be MJD 15020, so December 31, 1899 must have been MJD 15019.
// This date will be our reference.*/
int mjd2mpcdate(double MJD,int &year,int &month,double &day)
{
  double daynum;
  int dayct,yearct,leapcheck,cencheck,fourcheck;
  int isleap,daymonth[13],daymonthleap[13],i;

  daynum = MJD - (double)15019; /*Days since Dec 31, 1899*/
  if(daynum<=0.0)
    {
      cerr << "ERROR: date before 1900! ABORTING!\n";
      return(0);
    }

  daymonth[1] = 31;
  daymonth[2] = 28;
  daymonth[3] = 31;
  daymonth[4] = 30;
  daymonth[5] = 31;
  daymonth[6] = 30;
  daymonth[7] = 31;
  daymonth[8] = 31;
  daymonth[9] = 30;
  daymonth[10] = 31;
  daymonth[11] = 30;
  daymonth[12] = 31;

  for(i=1;i<=12;i++) daymonthleap[i] = daymonth[i];
  daymonthleap[2] = 29;

  dayct=daynum;
  day = daynum-(double)dayct;

  /*Now count up from 1900 to see what year it is*/
  yearct=1900;
  isleap=0;
  while((isleap==0&&dayct>365)||(isleap==1&&dayct>366))
    {
      /*printf("year = %d, dayct = %d, isleap = %d\n",yearct,dayct,isleap);*/
      /*Subtract the days for year number yearct*/
      if(isleap==0) dayct-=365;
      else if(isleap==1) dayct-=366;
      /*Go on to the next year*/
      yearct+=1;
      /*Find out if it is a leap year*/
      isleap = 0;
      leapcheck = yearct/4;
      if(yearct==(leapcheck*4))
	{
	  /*The year is evenly divisible by four: it's a leap year.*/
	  isleap = 1;
	  cencheck = yearct/100;
	  if(yearct==(cencheck*100))
	    {
	      /*No, wait: it's also the turn of a century: it's
                not a leap year.*/
	      isleap = 0;
	      fourcheck = yearct/400;
	      /*Except if it's divisible by 400, it is a leap year
                after all.*/
	      if(yearct==(fourcheck*400)) isleap = 1;
	    }
	}
    }
  /*OK, yearct is now equal to the correct year, and
    isleap correctly indicates whether or not this is
    a leap year.*/
  year = yearct;

  /*Find out what month it is*/
  i=1;
  if(isleap==0)
    {
      while(dayct>daymonth[i])
	{
	  dayct-=daymonth[i];
	  i+=1;
	}
    }
  else if(isleap==1)
    {
      while(dayct>daymonthleap[i])
	{
	  dayct-=daymonthleap[i];
	  i+=1;
	}
    }
  month = i;
  day+=(double)dayct; 
 return(1);
}

// stringline01: March 22, 2022: Given an input string read
// from a file, pull out individual pieces delimited by spaces,
// commas, or tabs, and load them in an input string vector.
// Return the number of distinct strings that were read.
int stringline01(const string &lnfromfile, vector <string> &outstrings) {
  char c='\0';
  string onestring;
  int i=0;
  int linelen = lnfromfile.size();
  int stringnum=0;
  
  outstrings={};
  if(linelen<=0) return(0);

  while(i<linelen) {
    onestring="";
    c=lnfromfile[i];
    i++;
    while(i<=linelen && c!=' ' && c!=',' && c!='\t' && c!='\n' && c!=EOF) {
      onestring.push_back(c);
      if(i<linelen) c=lnfromfile[i];
      i++;
    }
    if(onestring.size()>0) {
      stringnum++;
      outstrings.push_back(onestring);
    }
  }
  return(stringnum);
}

long double ldmedian(const vector <long double> &invec) {
  vector <long double> invec2 = invec; // Mutable copy of immutable input.
  int npts = invec2.size();
  int pt;
  pt=0;
  if(npts<=0) {
    cerr << "ERROR: ldmedian called on an input vector of length zero\n";
    return(0.0L);
  } else if(npts==1) {
    return(invec[0]);
  } else if(npts%2 == 0) {
    // We have an even number of points
    // Sort invec.
    sort(invec2.begin(),invec2.end());
    // Median is average of two middle points
    pt = npts/2;
    return((invec2[pt-1] + invec2[pt])/2.0L);
  } else {
    // We have an odd number of points.
    // Sort invec.
    sort(invec2.begin(),invec2.end());
    // Median is value of the middle point
    pt = (npts-1)/2;
    return(invec2[pt]);
  }
}
    
int ldmedian_minmax(const vector <long double> &invec, long double &median, long double &min, long double &max) {
  vector <long double> invec2 = invec; // Mutable copy of immutable input.
  int npts = invec2.size();
  int pt;
  pt=0;
  if(npts<=0) {
    cerr << "ERROR: ldmedian_minmax called on an input vector of length zero\n";
    median = min = max = 0.0L;
    return(1);
  } else if(npts==1) {
    median = min = max = invec[0];
    return(0);
  } else {
    // Sort invec.
    sort(invec2.begin(),invec2.end());
    min = invec2[0];
    max = invec2[npts-1];
    if(npts%2 == 0) {
      // We have an even number of points
      // Median is average of two middle points
      pt = npts/2;
      median = (invec2[pt-1] + invec2[pt])/2.0L;
    } else {
      // We have an odd number of points.
      // Median is value of the middle point
      pt = (npts-1)/2;
      median = invec2[pt];
    }
    return(0);
  }
}
    
double dmedian(const vector <double> &invec) {
  vector <double> invec2 = invec; // Mutable copy of immutable input.
  int npts = invec2.size();
  int pt;
  pt=0;
  if(npts<=0) {
    cerr << "ERROR: dmedian called on an input vector of length zero\n";
    return(0.0l);
  } else if(npts==1) {
    return(invec[0]);
  } else if(npts%2 == 0) {
    // We have an even number of points
    // Sort invec.
    sort(invec2.begin(),invec2.end());
    // Median is average of two middle points
    pt = npts/2;
    return((invec2[pt-1] + invec2[pt])/2.0l);
  } else {
    // We have an odd number of points.
    // Sort invec.
    sort(invec2.begin(),invec2.end());
    // Median is value of the middle point
    pt = (npts-1)/2;
    return(invec2[pt]);
  }
}
    
int dmedian_minmax(const vector <double> &invec, double &median, double &min, long double &max) {
  vector <double> invec2 = invec; // Mutable copy of immutable input.
  int npts = invec2.size();
  int pt;
  pt=0;
  if(npts<=0) {
    cerr << "ERROR: dmedian_minmax called on an input vector of length zero\n";
    median = min = max = 0.0l;
    return(1);
  } else if(npts==1) {
    median = min = max = invec[0];
    return(0);
  } else {
    // Sort invec.
    sort(invec2.begin(),invec2.end());
    min = invec2[0];
    max = invec2[npts-1];
    if(npts%2 == 0) {
      // We have an even number of points
      // Median is average of two middle points
      pt = npts/2;
      median = (invec2[pt-1] + invec2[pt])/2.0l;
    } else {
      // We have an odd number of points.
      // Median is value of the middle point
      pt = (npts-1)/2;
      median = invec2[pt];
    }
    return(0);
  }
}
    
float fmedian(const vector <float> &invec) {
  vector <float> invec2 = invec; // Mutable copy of immutable input.
  int npts = invec2.size();
  int pt;
  pt=0;
  if(npts<=0) {
    cerr << "ERROR: fmedian called on an input vector of length zero\n";
    return(0.0);
  } else if(npts==1) {
    return(invec[0]);
  } else if(npts%2 == 0) {
    // We have an even number of points
    // Sort invec.
    sort(invec2.begin(),invec2.end());
    // Median is average of two middle points
    pt = npts/2;
    return((invec2[pt-1] + invec2[pt])/2.0l);
  } else {
    // We have an odd number of points.
    // Sort invec.
    sort(invec2.begin(),invec2.end());
    // Median is value of the middle point
    pt = (npts-1)/2;
    return(invec2[pt]);
  }
}
    
int fmedian_minmax(const vector <float> &invec, float &median, float &min, float &max) {
  vector <float> invec2 = invec; // Mutable copy of immutable input.
  int npts = invec2.size();
  int pt;
  pt=0;
  if(npts<=0) {
    cerr << "ERROR: fmedian_minmax called on an input vector of length zero\n";
    median = min = max = 0.0;
    return(1);
  } else if(npts==1) {
    median = min = max = invec[0];
    return(0);
  } else {
    // Sort invec.
    sort(invec2.begin(),invec2.end());
    min = invec2[0];
    max = invec2[npts-1];
    if(npts%2 == 0) {
      // We have an even number of points
      // Median is average of two middle points
      pt = npts/2;
      median = (invec2[pt-1] + invec2[pt])/2.0l;
    } else {
      // We have an odd number of points.
      // Median is value of the middle point
      pt = (npts-1)/2;
      median = invec2[pt];
    }
    return(0);
  }
}
    
// stateunit_to_celestial: April 08, 2022:
// Given an input unit vector in the solar system barycentric coordinate
// system, calculate the RA, Dec coordinates to which it points.
// Example: if the input unit vector is the relative position of
// an asteroid relative to an observer in the barycentric system,
// the output will be the asteroid's RA and Dec in the sky as seen
// by that observer.
// Note: this is an improved version that replaces an older (April 8, 2022)
// routine of the same name that was slower because it used an unnecessarily
// clumsy sequence of mathematical operations, including a call to poleswitch.
int stateunit_to_celestial(point3d &baryvec, double &RA, double &Dec)
{
  double poleDec,yszc,yczs,xe;

  poleDec = NEPDEC/DEGPRAD;
  yszc = baryvec.y*sin(poleDec) - baryvec.z*cos(poleDec);
  yczs = baryvec.y*cos(poleDec) + baryvec.z*sin(poleDec);
  xe = baryvec.x;

  if(yszc==0 && xe>=0) RA = 0.0;
  else if(yszc==0 && xe<0) RA = M_PI;
  else if(yszc>0) RA = M_PI/2.0L - atan(xe/yszc);
  else if(yszc<0) RA = 3.0L*M_PI/2.0L - atan(xe/yszc);
  else {
    //Logically excluded case
    cerr << "Logically excluded case in stateunit_to_celestial\n";
    cerr << "xe = " << xe << " yszc = " << yszc << "\n";
    return(1);
  }

  if(yczs>1.0) {
    if(WARN_INVERSE_TRIG>0) cerr << "Warning: stateunit_to_celestial attempting to take arcsin of 1 + " << yczs-1.0L << "\n";
    Dec = M_PI/2.0L;
  } else if(yczs<-1.0) {
    if(WARN_INVERSE_TRIG>0) cerr << "Warning: stateunit_to_celestial attempting to take arcsin of -1 - " << yczs+1.0L << "\n";
    Dec = -M_PI/2.0L;
  } else {
    Dec = asin(yczs);
  }

  RA *= DEGPRAD;
  Dec *= DEGPRAD;
  return(0);
}

// stateunitLD_to_celestial: September 01, 2022:
// Given an input unit vector in the solar system barycentric coordinate
// system, calculate the RA, Dec coordinates to which it points.
// Example: if the input unit vector is the relative position of
// an asteroid relative to an observer in the barycentric system,
// the output will be the asteroid's RA and Dec in the sky as seen
// by that observer.
// Note: this is an improved version that replaces an older routine of
// the same name that was slower because it used an unnecessarily clumsy
// sequence of mathematical operations, including a call to poleswitch.
int stateunitLD_to_celestial(point3LD &baryvec, long double &RA, long double &Dec)
{
  long double poleDec,yszc,yczs,xe;

  poleDec = NEPDEC/DEGPRAD;
  yszc = baryvec.y*sin(poleDec) - baryvec.z*cos(poleDec);
  yczs = baryvec.y*cos(poleDec) + baryvec.z*sin(poleDec);
  xe = baryvec.x;

  if(yszc==0 && xe>=0) RA = 0.0;
  else if(yszc==0 && xe<0) RA = M_PI;
  else if(yszc>0) RA = M_PI/2.0L - atan(xe/yszc);
  else if(yszc<0) RA = 3.0L*M_PI/2.0L - atan(xe/yszc);
  else {
    //Logically excluded case
    cerr << "Logically excluded case in stateunitLD_to_celestial\n";
    cerr << "xe = " << xe << " yszc = " << yszc << "\n";
    return(1);
  }

  if(yczs>1.0) {
    if(WARN_INVERSE_TRIG>0) cerr << "Warning: stateunitLD_to_celestial attempting to take arcsin of 1 + " << yczs-1.0L << "\n";
    Dec = M_PI/2.0L;
  } else if(yczs<-1.0) {
    if(WARN_INVERSE_TRIG>0) cerr << "Warning: stateunitLD_to_celestial attempting to take arcsin of -1 - " << yczs+1.0L << "\n";
    Dec = -M_PI/2.0L;
  } else {
    Dec = asin(yczs);
  }

  RA *= DEGPRAD;
  Dec *= DEGPRAD;
  return(0);
}


#define DEBUG 0

// integrate_orbit03LD: April 08, 2022
// Uses modeling of the acceleration as a polynomial of order n>1
// to integrate the orbit of a massless test particle (e.g. asteroid)
// under the gravity of multiple 'planets'. It is assumed that
// in general these 'planets' will consist of the Sun, the
// eight major planets, and the Moon (possibly needed for
// cases of NEOs closely approaching the Earth). However,
// more or fewer planets may be used as desired.
// Note that the vector obsMJD is assumed to be time-sorted, and
// serious failures will result if it is not.
int integrate_orbit03LD(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const vector <long double> &obsMJD, point3LD startpos, point3LD startvel, long double mjdstart, vector <point3LD> &obspos, vector <point3LD> &obsvel)
{
  vector <point3LD> planetsalltimes;
  vector <point3LD> planetsonce;
  vector <point3LD> targaccel;
  vector <point3LD> accelfit;
  vector <point3LD> targvel;
  vector <point3LD> targpos;
  vector <point3LD> accelmod;
  vector <long double> temptime;
  vector <long double> ppxvec;
  vector <long double> ppyvec;
  vector <long double> ppfitvec;
  vector <long double> obsMJD2 = obsMJD; // Mutable copy of immutable input vector.
  point3LD singleaccel = point3LD(0L,0L,0L);
  point3LD singlevel = point3LD(0L,0L,0L);
  point3LD singlepos = point3LD(0L,0L,0L);
  int i=0;
  int j=0;
  int k=0;
  int obsnum=obsMJD.size();
  int obsct=0;
  int planetpointnum = planetmjd.size();
  vector <long double> forwardmjd;
  vector <long double> backwardmjd;
  int pointafter=0;
  int pointbefore=0;
  int latestpoint=0;
  int stepsin=0;
  long double dt0=0L;
  long double dt1=0L;
  long double dt2=0L;
  long double timemult=0L;
  point3LD accelslope = point3LD(0L,0L,0L);
  int obsbefore=0;
  int obsafter=0;
  long double lastobstime=0L;
  long double firstobstime=0L;
  int obspoint=0;

  // Correct the input observed times, assumed to be UT1, by the
  // corretion to TT (terrestrial time) used, e.g. by JPL Horizons.
  for(obsct=0;obsct<obsnum;obsct++) {
    obsMJD2[obsct] += TTDELTAT/SOLARDAY; // UT is always behind TT, hence UT must be adjusted forward.
  }
    
  firstobstime = obsMJD2[0];
  lastobstime = obsMJD2[obsnum-1];
  if(DEBUG>1) cout << "Integration will span times from MJD " << firstobstime << " to " << lastobstime << "\n";

  if(polyorder<2) {
    cerr << "ERROR: integrate_orbit03LD called with polyorder = " << polyorder << "\n";
    cerr << "polyorder must be at least 2!\n";
    return(1);
  }
  
  if(lastobstime<firstobstime) {
    cerr << "ERROR: integrate_orbit03LD called with end time (" << lastobstime << ") before start time (" << firstobstime << ")\n";
    return(1);
  } else if(mjdstart<=planetmjd[1] || firstobstime<=planetmjd[1] || lastobstime>=planetmjd[planetpointnum-1]) {
    cerr << "ERROR: integrate_orbit03LD called with start time " << mjdstart << " or time range )" << firstobstime << "-" << lastobstime << ") outside range of planet vectors (" << planetmjd[1] << "-" << planetmjd[planetpointnum-1] << ")\n";
    return(1);
  }
  // Make sure that relevant vectors for the polynomial fitting
  // are all large enough.
  for(i=0;i<=polyorder+1;i++) {
    targaccel.push_back(singleaccel);
    accelfit.push_back(singleaccel);
    targvel.push_back(singlevel);
    targpos.push_back(singlepos);
    accelmod.push_back(singleaccel);
    temptime.push_back(0L);
    ppfitvec.push_back(0L);
  }

  // Make sure the output vectors are large enough
  obspos={};
  obsvel={};
  for(obsct=0;obsct<=obsnum;obsct++) {
    obsvel.push_back(singlevel);
    obspos.push_back(singlepos);
  }
  // Are the observations before or after the starting time, or both?
  obsbefore=obsafter=0;
  for(obsct=0;obsct<=obsnum;obsct++) {
    if(obsMJD2[obsct]<mjdstart) obsbefore=1;
    if(obsMJD2[obsct]>=mjdstart) obsafter=1;
  }

  if(DEBUG>1) cout << "Checking for forward integration\n";
  if(obsafter==1) {
  if(DEBUG>1) cout << "Forward integration will be performed\n";
    // Integrate forward in time to find the position of
    // the target asteroid at all positions after starttime.

    // Load the initial time and planet position vectors
    planetsonce={};
    forwardmjd={};
    planetsalltimes={};
    temptime[0] = mjdstart;
    forwardmjd.push_back(mjdstart);
    nplanetpos01LD(temptime[0]-TTDELTAT/SOLARDAY,planetnum,5,planetmjd,planetpos,planetsonce);
    for(i=0;i<planetnum;i++) planetsalltimes.push_back(planetsonce[i]);
    if(DEBUG>1) cout  << fixed << setprecision(6) << "mjdstart = " << mjdstart << " loaded: " << forwardmjd[0] << "\n";
 
    // Find the first observation time simultaneous with or after starttime
    obsct=0;
    while(obsMJD2[obsct]<mjdstart) obsct++;
    obspoint=obsct;
    
    // Find the first point in the planet files that is after mjdstart
    for(i=0;i<planetpointnum;i++) {
      if(planetmjd[i]>mjdstart) break;
    }
    pointafter = i; // first point after mjdstart

    if(DEBUG>1) cout  << fixed << setprecision(6) << "Starting points for forward integration:\n";
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Observations, point " << obspoint << ", time " << obsMJD2[obspoint] << "\n";
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Planet file, point " << pointafter << ", time " << planetmjd[pointafter] << "\n";
    
    // Load a vector with times and planet positions for all the planet file
    // points spanning the times of the observations.
    
    dt0 = (planetmjd[pointafter+1] - planetmjd[pointafter])*SOLARDAY/sqrt(M_PI);
    
    j=1;
    i=0;
    obsct=obspoint;
    // j counts elements actually loaded into forwardmjd and planetsalltimes
    // i counts steps in the planetfile past pointafter
    // obsct indexes the observation vector.
    for(obsct=obspoint;obsct<obsnum;obsct++) {
      while(planetmjd[pointafter+i]<obsMJD2[obsct]) {
	// Load any planet file points until we get past obsMJD2[obsct]
	forwardmjd.push_back(planetmjd[pointafter+i]);
	if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " > planetmjd[" << pointafter+i << "] = " << planetmjd[pointafter+i] << ": loaded forwardmjd[" << j << "] = " << forwardmjd[j] << "\n";
	planetsonce={};
	nplanetgrab01LD(pointafter+i, planetnum, planetmjd, planetpos, planetsonce);
	for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
	i++;
	j++;
      }
      // Load the next observation point.
      forwardmjd.push_back(obsMJD2[obsct]);
      planetsonce={};
      nplanetpos01LD(obsMJD2[obsct]-TTDELTAT/SOLARDAY,planetnum,5,planetmjd,planetpos,planetsonce);
      for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
      j++;
      if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " < planetmjd[" << pointafter+i << "] = " << planetmjd[pointafter+i] << ": loaded forwardmjd[" << j-1 << "] = " << forwardmjd[j-1] << "\n";
    }
    // Add extra planet points until we have polyorder+2
    while(j<polyorder+2) {
      forwardmjd.push_back(planetmjd[pointafter+i]);
      planetsonce={};
      nplanetgrab01LD(pointafter+i, planetnum, planetmjd, planetpos, planetsonce);
      for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
      if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " < planetmjd[" << pointafter+i << "] = " << planetmjd[pointafter+i] << ": loaded forwardmjd[" << j-1 << "] = " << forwardmjd[j-1] << "\n";
      i++;
      j++;
    }
    if(DEBUG>1) {
      cout << fixed << setprecision(6)  << "Loaded " << forwardmjd.size() << " points in forwardmjd\n";
      for(i=0;i<long(forwardmjd.size());i++) {
	cout  << fixed << setprecision(6) << "Forward MJD = " << forwardmjd[i] << "\n";
      }
    }
    // Now we've loaded all the times we want for the forward integration
    // Load temporary time vector
    for(j=0;j<=polyorder+1;j++) temptime[j] = forwardmjd[j];

    //cout << fixed << setprecision(6)  << "Loaded " << forwardmjd.size() << " points in forwardmjd\n";
    //cout << fixed << setprecision(6)  << "Loaded " << planetsalltimes.size() << " points in planetsalltimes\n";
    
    // Load starting position and velocity
    targvel[0] = startvel;
    targpos[0] = startpos;
    j=0;
    //cout << temptime[j] << " " << targpos[j].x << " " << targpos[j].y << " " << targpos[j].z << " " << targvel[j].x << " " << targvel[j].y << " " << targvel[j].z << "\n";
    // Bootstrap up to a fit of order polyorder.
    // Calculate acceleration at starting point, loading planet positions from big vector.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[0], targaccel[0]);
    dt1 = (temptime[1]-temptime[0])*SOLARDAY;
    // First Approx: estimate next position, assuming constant acceleration.
    targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1;
    targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1;
    targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1;
    // Calculate acceleration at this new position.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*1 + i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);
  
    // Second approx: linearly varying acceleration.
    accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
    accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
    accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;

    // Improved position for next time step.
    targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
    targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
    targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

    // Re-calculate acceleration at this improved position.
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);

    // Re-calculate improved acceleration slope.
    accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
    accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
    accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;
  
    // Improved velocity for next time step
    targvel[1].x = targvel[0].x + targaccel[0].x*dt1 + accelslope.x*0.5L*dt1*dt1;
    targvel[1].y = targvel[0].y + targaccel[0].y*dt1 + accelslope.y*0.5L*dt1*dt1;
    targvel[1].z = targvel[0].z + targaccel[0].z*dt1 + accelslope.z*0.5L*dt1*dt1;

    // Use linearly extrapolated acceleration to estimate position for
    // the next time step.
    dt1 = (temptime[2]-temptime[1])*SOLARDAY;
    targpos[2].x = targpos[1].x + targvel[1].x*dt1 + targaccel[1].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
    targpos[2].y = targpos[1].y + targvel[1].y*dt1 + targaccel[1].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
    targpos[2].z = targpos[1].z + targvel[1].z*dt1 + targaccel[1].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

    // Calculate acceleration for this extrapolated position.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*2 + i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[2], targaccel[2]);

    // Now we have three acceleration points: can load for a full polynomial fit.
    for(stepsin=3;stepsin<=polyorder+1;stepsin++) {
      // Fit for x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<stepsin;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=stepsin;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=stepsin;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
      cout.precision(17);
  
      // Perform new fits to revised accelerations
      // Fit for x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<stepsin;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
      for(i=0;i<stepsin;i++) {
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=stepsin;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=stepsin;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
    }

    // We are now set up for a full-order polynomial integration.
    // We have valid positions in targpos, targvel, and temptime
    // for indices from 0 to polyorder+1.
    // See if any of the target points have already been calculated
    j=0;
    obsct=obspoint;
    while(j<=polyorder+1 && obsct<obsnum) {
      if(obsMJD2[obsct]==temptime[j]) {
	obspos[obsct]=targpos[j];
	obsvel[obsct]=targvel[j];
	j++;
	obsct++;
      } else if(obsMJD2[obsct]<temptime[j]) obsct++;
      else if(temptime[j]<obsMJD2[obsct]) j++;
      else {
	cerr << "Impossible time comparison case: " << temptime[j] << " " << obsMJD2[obsct] << "\n";
	return(4);
      }
    }
    latestpoint=polyorder+1;
    // Proceed with the full polynomial integration.
    while(latestpoint<long(forwardmjd.size())-1) {
      latestpoint++;
      // Cycle the dynamical vectors
      for(i=0;i<polyorder+1;i++) {
	temptime[i] = temptime[i+1];
	targaccel[i] = targaccel[i+1];
	targvel[i] = targvel[i+1];
	targpos[i] = targpos[i+1];
      }
      // Load a new point into temptime
      temptime[polyorder+1] = forwardmjd[latestpoint];
      // Fit for acceleration
      // x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<polyorder+1;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=polyorder+1;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=polyorder+1;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + i];
	//cout << "Using planet point " << latestpoint+j-polyorder-1 << " Earth at x = " << planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + 3].x << "\n";
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }

      // Use these revised accelerations to re-do the fits
      // x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<polyorder+1;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=polyorder+1;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=polyorder+1;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
	//cout << "Using planet point " << latestpoint+j-polyorder-1 << " Earth at x = " << planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + 3].x << "\n";
      }
      j=polyorder+1;
      //cout << temptime[j] << " " << targpos[j].x << " " << targpos[j].y << " " << targpos[j].z << " " << targvel[j].x << " " << targvel[j].y << " " << targvel[j].z << "\n";
      // Load any target points that have newly been calculated,
      // or re-calculated
     j=0;
      obsct=obspoint;
       while(j<=polyorder+1 && obsct<obsnum) {
	if(obsMJD2[obsct]==temptime[j]) {
	  obspos[obsct]=targpos[j];
	  obsvel[obsct]=targvel[j];
	  j++;
	  obsct++;
	} else if(obsMJD2[obsct]<temptime[j]) obsct++;
	else if(temptime[j]<obsMJD2[obsct]) j++;
	else {
	  cerr << "Impossible time comparison case: " << temptime[j] << " " << obsMJD2[obsct] << "\n";
	  return(4);
	}
      }
      // We have now gone through two iterations of extrapolation
      // to predict the next acceleration point as accurately as possible.
      // The next step of the loop will move the extrapolated point back
      // by one step, and use it to start extrapolating a new point,
      // at the same time refining the former extrapolated points.
    }
  } // END OF FORWARD INTEGRATION
  if(DEBUG>1) {
    cout << "Results of forward integration:\n";
    for(j=0;j<long(temptime.size());j++) {
      cout  << fixed << setprecision(6) << temptime[j] << " " << targpos[j].x/AU_KM << " " << targpos[j].y/AU_KM << " " << targpos[j].z/AU_KM << " " << targvel[j].x << " " << targvel[j].y << " " << targvel[j].z << "\n";
    }
  }

  // NOW PEFORM BACKWARD INTEGRATION, IF NECESSARY
  if(DEBUG>1) cout << "Checking for backward integration\n";
  if(obsbefore==1) {
  if(DEBUG>1) cout << "Backward integration will be performed\n";
    // Integrate backward in time to find the position of
    // the target asteroid at all positions before starttime.

    // Load the initial time and planet position vectors
    planetsonce={};
    backwardmjd={};
    planetsalltimes={};
    temptime[0] = -mjdstart;
    backwardmjd.push_back(-mjdstart);
    nplanetpos01LD(mjdstart-TTDELTAT/SOLARDAY,planetnum,5,planetmjd,planetpos,planetsonce);
    for(i=0;i<planetnum;i++) planetsalltimes.push_back(planetsonce[i]);
    if(DEBUG>1) cout  << fixed << setprecision(6) << "mjdstart = " << mjdstart << " loaded: " << backwardmjd[0] << "\n";

    // Find the first observation time before starttime
    obsct=obsnum-1;
    while(obsMJD2[obsct]>=mjdstart) obsct--;
    obspoint=obsct;
    
    // Find the first point in the planet files that is before mjdstart
    for(i=planetpointnum-1;i>0;i--) {
      if(planetmjd[i]<mjdstart) break;
    }
    pointbefore = i; // first point before mjdstart

    if(DEBUG>1) cout  << fixed << setprecision(6) << "Starting points for backward integration:\n";
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Observations, point " << obspoint << ", time " << obsMJD2[obspoint] << "\n";
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Planet file, point " << pointbefore << ", time " << planetmjd[pointbefore] << "\n";
    
    // Load a vector with times and planet positions for all the planet file
    // points spanning the times of the observations.
    
    // Time scaling factor, designed to avoid integers or near-zero values
    dt0 = (planetmjd[pointbefore+1] - planetmjd[pointbefore])*SOLARDAY/sqrt(M_PI);
    
    j=1;
    i=0;
    obsct=obspoint;
    // j counts elements actually loaded into forwardmjd and planetsalltimes
    // i counts steps in the planetfile before pointbefore
    // obsct indexes the observation vector.
    for(obsct=obspoint; obsct>=0; obsct--) {
      while(planetmjd[pointbefore-i]>obsMJD2[obsct]) {
	// Load any planet file points until we get to one before obsMJD2[obsct]
      	backwardmjd.push_back(-planetmjd[pointbefore-i]);
	if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " < planetmjd[" << pointbefore-i << "] = " << planetmjd[pointbefore-i] << ": loaded backwardmjd[" << j << "] = " << backwardmjd[j] << "\n";
	planetsonce={};
	nplanetgrab01LD(pointbefore-i, planetnum, planetmjd, planetpos, planetsonce);
	for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
	i++;
	j++;
      }
      // Load the next observation point.
      backwardmjd.push_back(-obsMJD2[obsct]);
      planetsonce={};
      nplanetpos01LD(obsMJD2[obsct]-TTDELTAT/SOLARDAY,planetnum,5,planetmjd,planetpos,planetsonce);
      for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
      j++;
      if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " > planetmjd[" << pointbefore-i << "] = " << planetmjd[pointbefore-i] << ": loaded backwardmjd[" << j-1 << "] = " << backwardmjd[j-1] << "\n";
    }
    // Add extra planet points until we have polyorder+2
    while(j<polyorder+2) {
      backwardmjd.push_back(-planetmjd[pointbefore-i]);
      planetsonce={};
      nplanetgrab01LD(pointbefore-i, planetnum, planetmjd, planetpos, planetsonce);
      for(k=0;k<planetnum;k++) planetsalltimes.push_back(planetsonce[k]);
      if(DEBUG>1) cout << "obsMJD2[" << obsct << "] = " << obsMJD2[obsct] << " > planetmjd[" << pointbefore-i << "] = " << planetmjd[pointbefore-i] << ": loaded backwardmjd[" << j-1 << "] = " << backwardmjd[j-1] << "\n";
      i++;
      j++;
    }
    if(DEBUG>1) {
      cout << fixed << setprecision(6)  << "Loaded " << backwardmjd.size() << " points in backwardmjd\n";
      for(i=0;i<long(backwardmjd.size());i++) {
	cout  << fixed << setprecision(6) << "Backward MJD = " << backwardmjd[i] << "\n";
      }
    }

    // Now we've loaded all the times we want for the backward integration
    // Load temporary time vector
    for(j=0;j<=polyorder+1;j++) temptime[j] = backwardmjd[j];

    // Load starting position and velocity, negative velocity for backward integration
    targpos[0] = startpos;
    targvel[0].x = -startvel.x;
    targvel[0].y = -startvel.y;
    targvel[0].z = -startvel.z;
    // Bootstrap up to a fit of order polyorder.
    // Calculate acceleration at starting point, loading planet positions from big vector.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[0], targaccel[0]);
    dt1 = (temptime[1]-temptime[0])*SOLARDAY;
    // First Approx: estimate next position, assuming constant acceleration.
    targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1;
    targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1;
    targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1;
    // Calculate acceleration at this new position.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*1 + i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);
  
    // Second approx: linearly varying acceleration.
    accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
    accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
    accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;

    // Improved position for next time step.
    targpos[1].x = targpos[0].x + targvel[0].x*dt1 + targaccel[0].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
    targpos[1].y = targpos[0].y + targvel[0].y*dt1 + targaccel[0].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
    targpos[1].z = targpos[0].z + targvel[0].z*dt1 + targaccel[0].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

    // Re-calculate acceleration at this improved position.
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);

    // Re-calculate improved acceleration slope.
    accelslope.x = (targaccel[1].x-targaccel[0].x)/dt1;
    accelslope.y = (targaccel[1].y-targaccel[0].y)/dt1;
    accelslope.z = (targaccel[1].z-targaccel[0].z)/dt1;
  
    // Improved velocity for next time step
    targvel[1].x = targvel[0].x + targaccel[0].x*dt1 + accelslope.x*0.5L*dt1*dt1;
    targvel[1].y = targvel[0].y + targaccel[0].y*dt1 + accelslope.y*0.5L*dt1*dt1;
    targvel[1].z = targvel[0].z + targaccel[0].z*dt1 + accelslope.z*0.5L*dt1*dt1;

    // Use linearly extrapolated acceleration to estimate position for
    // the next time step.
    dt1 = (temptime[2]-temptime[1])*SOLARDAY;
    targpos[2].x = targpos[1].x + targvel[1].x*dt1 + targaccel[1].x*0.5L*dt1*dt1 + accelslope.x*dt1*dt1*dt1/6.0L;
    targpos[2].y = targpos[1].y + targvel[1].y*dt1 + targaccel[1].y*0.5L*dt1*dt1 + accelslope.y*dt1*dt1*dt1/6.0L;
    targpos[2].z = targpos[1].z + targvel[1].z*dt1 + targaccel[1].z*0.5L*dt1*dt1 + accelslope.z*dt1*dt1*dt1/6.0L;

    // Calculate acceleration for this extrapolated position.
    for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*2 + i];
    accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[2], targaccel[2]);

    // Now we have three acceleration points: can load for a full polynomial fit.
    for(stepsin=3;stepsin<=polyorder+1;stepsin++) {
      // Fit for x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<stepsin;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=stepsin;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=stepsin;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
      cout.precision(17);
  
      // Perform new fits to revised accelerations
      // Fit for x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<stepsin;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
      for(i=0;i<stepsin;i++) {
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=stepsin;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=stepsin;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*j + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
      
    }

    // We are now set up for a full-order polynomial integration.
    // We have valid positions in targpos, targvel, and temptime
    // for indices from 0 to polyorder+1.
    // See if any of the target points have already been calculated
    j=0;
    obsct=obspoint;
    while(j<=polyorder+1 && obsct>=0) {
      if(obsMJD2[obsct]==-temptime[j]) {
	obspos[obsct]=targpos[j];
	obsvel[obsct].x=-targvel[j].x;
	obsvel[obsct].y=-targvel[j].y;
	obsvel[obsct].z=-targvel[j].z;
	j++;
	obsct--;
      } else if(obsMJD2[obsct]>-temptime[j]) obsct--;
      else if(-temptime[j]>obsMJD2[obsct]) j++;
      else {
	cerr << "Impossible time comparison case: " << -temptime[j] << " " << obsMJD2[obsct] << "\n";
	return(4);
      }
    }
    latestpoint=polyorder+1;
    // Proceed with the full polynomial integration.
    while(latestpoint<long(backwardmjd.size())-1) {
      latestpoint++;
      // Cycle the dynamical vectors
      for(i=0;i<polyorder+1;i++) {
	temptime[i] = temptime[i+1];
	targaccel[i] = targaccel[i+1];
	targvel[i] = targvel[i+1];
	targpos[i] = targpos[i+1];
      }
      // Load a new point into temptime
      temptime[polyorder+1] = backwardmjd[latestpoint];
      // Fit for acceleration
      // x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<polyorder+1;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=polyorder+1;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=polyorder+1;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }

      // Use these revised accelerations to re-do the fits
      // x component of acceleration.
      ppxvec={};
      ppyvec={};
      for(i=0;i<polyorder+1;i++) {
	ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
	ppyvec.push_back(targaccel[i].x);
      }
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
      // Fit for y component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
      // Fit for z component of acceleration. Note that we have
      // already loaded the time vector ppxvec.
      ppyvec={};
      for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
      // Perform fit, and store in accelfit.
      perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
      for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
      // Re-calculate all of the positions and velocities using this fit.
      for(j=1;j<=polyorder+1;j++) {
	dt2 = (temptime[j]-temptime[0])*SOLARDAY;
	// Positions
	targpos[j].x = targpos[0].x + targvel[0].x*dt2;
	targpos[j].y = targpos[0].y + targvel[0].y*dt2;
	targpos[j].z = targpos[0].z + targvel[0].z*dt2;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	  targpos[j].x += accelfit[i].x*timemult;
	  targpos[j].y += accelfit[i].y*timemult;
	  targpos[j].z += accelfit[i].z*timemult;
	}
	// Velocities
	targvel[j].x = targvel[0].x;
	targvel[j].y = targvel[0].y;
	targvel[j].z = targvel[0].z;
	for(i=0;i<polyorder+1;i++) {
	  timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	  targvel[j].x += accelfit[i].x*timemult;
	  targvel[j].y += accelfit[i].y*timemult;
	  targvel[j].z += accelfit[i].z*timemult;
	}
	// Accelerations
	accelmod[j].x = 0L;
	accelmod[j].y = 0L;
	accelmod[j].z = 0L;
	for(i=0;i<stepsin;i++) {
	  timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	  accelmod[j].x += accelfit[i].x*timemult;
	  accelmod[j].y += accelfit[i].y*timemult;
	  accelmod[j].z += accelfit[i].z*timemult;
	}
      }
      // Re-calculate accelerations using these revised positions
      for(j=1;j<=polyorder+1;j++) {
	for(i=0;i<planetnum;i++) planetsonce[i] = planetsalltimes[planetnum*(latestpoint+j-polyorder-1) + i];
	accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      }
      j=polyorder+1;
      //cout << temptime[j] << " " << targpos[j].x << " " << targpos[j].y << " " << targpos[j].z << " " << targvel[j].x << " " << targvel[j].y << " " << targvel[j].z << "\n";
      // Load any target points that have newly been calculated,
      // or re-calculated
      j=0;
      obsct=obspoint;
      while(j<=polyorder+1 && obsct>=0) {
	if(obsMJD2[obsct]==-temptime[j]) {
	  obspos[obsct] = targpos[j];
	  obsvel[obsct].x = -targvel[j].x;
	  obsvel[obsct].y = -targvel[j].y;
	  obsvel[obsct].z = -targvel[j].z;
	  j++;
	  obsct--;
	} else if(obsMJD2[obsct]>-temptime[j]) obsct--;
	else if(-temptime[j]>obsMJD2[obsct]) j++;
	else {
	  cerr << "Impossible time comparison case: " << -temptime[j] << " " << obsMJD2[obsct] << "\n";
	  return(4);
	}
      }
      // We have now gone through two iterations of extrapolation
      // to predict the next acceleration point as accurately as possible.
      // The next step of the loop will move the extrapolated point back
      // by one step, and use it to start extrapolating a new point,
      // at the same time refining the former extrapolated points.
    }
  } // Closes the statement doing backward integration.
 
  return(0);
}

// tortoisechi01: April 11, 2022:
// Get chi-square value based on an input simplex point.
long double tortoisechi01(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, const vector <long double> &scalestate, long double timescale, long double mjdstart, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid)
{
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  int obsct;
  int obsnum = obsMJD.size();
  long double light_travel_time;
  point3LD outpos = point3LD(0,0,0);
  long double outRA=0L;
  long double outDec=0L;
  long double ldval=0L;
  long double chisq=0L;
  double dval;
  resid = fitRA = fitDec = {};
  
  // Input point scalestate is supposed to hold a 6-D state vector
  // in units of AU and AU/timescale, where timescale is a value in
  // days used to convert velocity to distance units

  // Convert to km and km/sec
  point3LD startpos = point3LD(scalestate[0]*AU_KM,scalestate[1]*AU_KM,scalestate[2]*AU_KM);
  point3LD startvel = point3LD(scalestate[3]*AU_KM/SOLARDAY/timescale,scalestate[4]*AU_KM/SOLARDAY/timescale,scalestate[5]*AU_KM/SOLARDAY/timescale);

  // Integrate orbit.
  integrate_orbit03LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, obsMJD, startpos, startvel, mjdstart, obspos, obsvel);
		  
  for(obsct=0;obsct<obsnum;obsct++) {
    // Initial approximation of the coordinates relative to the observer
    outpos.x = obspos[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - observerpos[obsct].z;
    // Initial approximation of the observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Convert to meters and divide by the speed of light to get the light travel time.
    light_travel_time = ldval*1000.0/CLIGHT;
    // Light-travel-time corrected version of coordinates relative to the observer
    outpos.x = obspos[obsct].x - light_travel_time*obsvel[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - light_travel_time*obsvel[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - light_travel_time*obsvel[obsct].z - observerpos[obsct].z;
    // Light-travel-time corrected observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    // Project onto the celestial sphere.
    stateunitLD_to_celestial(outpos, outRA, outDec);
    if(DEBUG>1) cout  << fixed << setprecision(6) << "Input MJD " << obsMJD[obsct] << ": " << obsRA[obsct] << " "  << obsDec[obsct] << " "  << " Output: " << outRA << ": " << outDec <<  "\n";
    dval = distradec01(obsRA[obsct],obsDec[obsct],outRA,outDec);
    dval *= 3600.0L; // Convert to arcsec
    fitRA.push_back(outRA);
    fitDec.push_back(outDec);
    resid.push_back(dval);
  }
  chisq=0.0L;
  for(obsct=0;obsct<obsnum;obsct++) {
    chisq += LDSQUARE(resid[obsct]/sigastrom[obsct]);
    if(DEBUG>0) cout << "Residual for point " << obsct << " is " << resid[obsct] << "\n";
  }
  return(chisq);
}

// integrate_orbit04LD: May 06, 2022: Like integrate_orbit03LD,
// but less complex in that it outputs positions only at times
// corresponding to the timesteps of the input planet files,
// not at a customized input list of 'observations times',
// and the start and end times MUST fall exactly on timesteps
// of the planet files (indeed, they are specified not by floating-point
// MJD values, but by integer indexes to the planet files, which must
// be determined by the calling function.
// The current program is simpler and faster
// than integrate_orbit03LD, as long as customized times are not needed.
//
// Description of ancestor program integrate_orbit03LD:
// Uses modeling of the acceleration as a polynomial of order n>1
// to integrate the orbit of a massless test particle (e.g. asteroid)
// under the gravity of multiple 'planets'. It is assumed that
// in general these 'planets' will consist of the Sun, the
// eight major planets, and the Moon (possibly needed for
// cases of NEOs closely approaching the Earth). However,
// more or fewer planets may be used as desired.
// Note that the vector obsMJD is assumed to be time-sorted, and
// serious failures will result if it is not.
int integrate_orbit04LD(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, point3LD startpos, point3LD startvel, int startpoint, int endpoint, vector <long double> &outMJD, vector <point3LD> &outpos, vector <point3LD> &outvel)
{
  vector <point3LD> planetsalltimes;
  vector <point3LD> planetsonce;
  vector <point3LD> targaccel;
  vector <point3LD> accelfit;
  vector <point3LD> targvel;
  vector <point3LD> targpos;
  vector <point3LD> accelmod;
  vector <long double> temptime;
  vector <long double> ppxvec;
  vector <long double> ppyvec;
  vector <long double> ppfitvec;
  point3LD singleaccel = point3LD(0L,0L,0L);
  point3LD singlevel = point3LD(0L,0L,0L);
  point3LD singlepos = point3LD(0L,0L,0L);
  int i=0;
  int j=0;
  int outnum = endpoint-startpoint+1;
  int outct=0;
  int latestpoint=0;
  int stepsin=0;
  long double dt0=0L;
  long double dt2=0L;
  long double timemult=0L;
  point3LD accelslope = point3LD(0L,0L,0L);

  if(polyorder<2) {
    cerr << "ERROR: integrate_orbit04LD called with polyorder = " << polyorder << "\n";
    cerr << "polyorder must be at least 4!\n";
    return(1);
  }
  
  if(endpoint<startpoint) {
    cerr << "ERROR: integrate_orbit04LD called with end point (" << endpoint << ") before starting point (" << startpoint << ")\n";
    return(1);
  } else if(startpoint<0 || endpoint>=long(planetmjd.size())) {
    cerr << "ERROR: integrate_orbit04LD called with starting point " << startpoint << " or endpoint" << endpoint << " outside range of planet vectors (0 - " << planetmjd.size() << ")\n";
    return(1);
  }
  // Make sure that relevant vectors for the polynomial fitting
  // are all large enough.
  for(i=0;i<=polyorder+1;i++) {
    targaccel.push_back(singleaccel);
    accelfit.push_back(singleaccel);
    targvel.push_back(singlevel);
    targpos.push_back(singlepos);
    accelmod.push_back(singleaccel);
    temptime.push_back(0L);
    ppfitvec.push_back(0L);
  }

  // Make sure the output vectors are large enough, and load
  // outMJD with the actual output times.
  outpos={};
  outvel={};
  outMJD={};
  for(outct=0;outct<=outnum;outct++) {
    outvel.push_back(singlevel);
    outpos.push_back(singlepos);
    outMJD.push_back(planetmjd[startpoint+outct]);
  }

  // Load the initial time vector
  for(j=0;j<=polyorder+1;j++) temptime[j] = planetmjd[startpoint+j];
    
  // Load starting position and velocity
  targvel[0] = startvel;
  targpos[0] = startpos;
  j=0;

  // Bootstrap up to a fit of order polyorder.
  // Calculate acceleration at starting point, loading planet positions from big vector.
  planetsonce={};
  nplanetgrab01LD(startpoint, planetnum, planetmjd, planetpos, planetsonce);
  accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[0], targaccel[0]);
  dt0 = (temptime[1]-temptime[0])*SOLARDAY;

  // First Approx: estimate next position, assuming constant acceleration.
  targpos[1].x = targpos[0].x + targvel[0].x*dt0 + targaccel[0].x*0.5L*dt0*dt0;
  targpos[1].y = targpos[0].y + targvel[0].y*dt0 + targaccel[0].y*0.5L*dt0*dt0;
  targpos[1].z = targpos[0].z + targvel[0].z*dt0 + targaccel[0].z*0.5L*dt0*dt0;

  // Calculate acceleration at this new position.
  planetsonce={};
  nplanetgrab01LD(startpoint+1, planetnum, planetmjd, planetpos, planetsonce);
  accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);
  
  // Second approx: linearly varying acceleration.
  accelslope.x = (targaccel[1].x-targaccel[0].x)/dt0;
  accelslope.y = (targaccel[1].y-targaccel[0].y)/dt0;
  accelslope.z = (targaccel[1].z-targaccel[0].z)/dt0;

  // Improved position for next time step.
  targpos[1].x = targpos[0].x + targvel[0].x*dt0 + targaccel[0].x*0.5L*dt0*dt0 + accelslope.x*dt0*dt0*dt0/6.0L;
  targpos[1].y = targpos[0].y + targvel[0].y*dt0 + targaccel[0].y*0.5L*dt0*dt0 + accelslope.y*dt0*dt0*dt0/6.0L;
  targpos[1].z = targpos[0].z + targvel[0].z*dt0 + targaccel[0].z*0.5L*dt0*dt0 + accelslope.z*dt0*dt0*dt0/6.0L;

  // Re-calculate acceleration at this improved position.
  accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[1], targaccel[1]);

  // Re-calculate improved acceleration slope.
  accelslope.x = (targaccel[1].x-targaccel[0].x)/dt0;
  accelslope.y = (targaccel[1].y-targaccel[0].y)/dt0;
  accelslope.z = (targaccel[1].z-targaccel[0].z)/dt0;
  
  // Improved velocity for next time step
  targvel[1].x = targvel[0].x + targaccel[0].x*dt0 + accelslope.x*0.5L*dt0*dt0;
  targvel[1].y = targvel[0].y + targaccel[0].y*dt0 + accelslope.y*0.5L*dt0*dt0;
  targvel[1].z = targvel[0].z + targaccel[0].z*dt0 + accelslope.z*0.5L*dt0*dt0;

  // Use linearly extrapolated acceleration to estimate position for
  // the next time step.
  dt0 = (temptime[2]-temptime[1])*SOLARDAY;
  targpos[2].x = targpos[1].x + targvel[1].x*dt0 + targaccel[1].x*0.5L*dt0*dt0 + accelslope.x*dt0*dt0*dt0/6.0L;
  targpos[2].y = targpos[1].y + targvel[1].y*dt0 + targaccel[1].y*0.5L*dt0*dt0 + accelslope.y*dt0*dt0*dt0/6.0L;
  targpos[2].z = targpos[1].z + targvel[1].z*dt0 + targaccel[1].z*0.5L*dt0*dt0 + accelslope.z*dt0*dt0*dt0/6.0L;
  // Calculate acceleration for this extrapolated position.
  planetsonce={};
  nplanetgrab01LD(startpoint+2, planetnum, planetmjd, planetpos, planetsonce);
  accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[2], targaccel[2]);

  // Now we have three acceleration points: can load for a full polynomial fit.
  for(stepsin=3;stepsin<=polyorder+1;stepsin++) {
    // Fit for x component of acceleration.
    ppxvec={};
    ppyvec={};
    for(i=0;i<stepsin;i++) {
      ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
      ppyvec.push_back(targaccel[i].x);
    }
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
    // Fit for y component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
    // Fit for z component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
    // Re-calculate all of the positions and velocities using this fit.
    for(j=1;j<=stepsin;j++) {
      dt2 = (temptime[j]-temptime[0])*SOLARDAY;
      // Positions
      targpos[j].x = targpos[0].x + targvel[0].x*dt2;
      targpos[j].y = targpos[0].y + targvel[0].y*dt2;
      targpos[j].z = targpos[0].z + targvel[0].z*dt2;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	targpos[j].x += accelfit[i].x*timemult;
	targpos[j].y += accelfit[i].y*timemult;
	targpos[j].z += accelfit[i].z*timemult;
      }
      // Velocities
      targvel[j].x = targvel[0].x;
      targvel[j].y = targvel[0].y;
      targvel[j].z = targvel[0].z;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	targvel[j].x += accelfit[i].x*timemult;
	targvel[j].y += accelfit[i].y*timemult;
	targvel[j].z += accelfit[i].z*timemult;
      }
      // Accelerations
      accelmod[j].x = 0L;
      accelmod[j].y = 0L;
      accelmod[j].z = 0L;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	accelmod[j].x += accelfit[i].x*timemult;
	accelmod[j].y += accelfit[i].y*timemult;
	accelmod[j].z += accelfit[i].z*timemult;
      }
    }
    // Re-calculate accelerations using these revised positions
    for(j=1;j<=stepsin;j++) {
      planetsonce={};
      nplanetgrab01LD(startpoint+j, planetnum, planetmjd, planetpos, planetsonce);
      accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
    }
    cout.precision(17);
  
    // Perform new fits to revised accelerations
    // Fit for x component of acceleration.
    ppxvec={};
    ppyvec={};
    for(i=0;i<stepsin;i++) {
      ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
      ppyvec.push_back(targaccel[i].x);
    }
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].x = ppfitvec[i];
    // Fit for y component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].y);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].y = ppfitvec[i];
    // Fit for z component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<stepsin;i++) ppyvec.push_back(targaccel[i].z);
    for(i=0;i<stepsin;i++) {
    }
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<stepsin;i++) accelfit[i].z = ppfitvec[i];
    // Re-calculate all of the positions and velocities using this fit.
    for(j=1;j<=stepsin;j++) {
      dt2 = (temptime[j]-temptime[0])*SOLARDAY;
      // Positions
      targpos[j].x = targpos[0].x + targvel[0].x*dt2;
      targpos[j].y = targpos[0].y + targvel[0].y*dt2;
      targpos[j].z = targpos[0].z + targvel[0].z*dt2;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	targpos[j].x += accelfit[i].x*timemult;
	targpos[j].y += accelfit[i].y*timemult;
	targpos[j].z += accelfit[i].z*timemult;
      }
      // Velocities
      targvel[j].x = targvel[0].x;
      targvel[j].y = targvel[0].y;
      targvel[j].z = targvel[0].z;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	targvel[j].x += accelfit[i].x*timemult;
	targvel[j].y += accelfit[i].y*timemult;
	targvel[j].z += accelfit[i].z*timemult;
      }
      // Accelerations
      accelmod[j].x = 0L;
      accelmod[j].y = 0L;
      accelmod[j].z = 0L;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	accelmod[j].x += accelfit[i].x*timemult;
	accelmod[j].y += accelfit[i].y*timemult;
	accelmod[j].z += accelfit[i].z*timemult;
      }
    }
    // Re-calculate accelerations using these revised positions
    for(j=1;j<=stepsin;j++) {
      planetsonce={};
      nplanetgrab01LD(startpoint+j, planetnum, planetmjd, planetpos, planetsonce);
      accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
    }
  }

  // We are now set up for a full-order polynomial integration.
  // We have valid positions in targpos, targvel, and temptime
  // for indices from 0 to polyorder+1.
  // Load already calculated points into the output vectors
  for(i=0;i<=polyorder+1;i++) {
    outvel[i] = targvel[i];
    outpos[i] = targpos[i];
  }
  // Define the current reference point.
  latestpoint=polyorder+1;
  // Proceed with the full polynomial integration.
  while(latestpoint<outnum) {
    latestpoint++;
    // Cycle the dynamical vectors
    for(i=0;i<polyorder+1;i++) {
      temptime[i] = temptime[i+1];
      targaccel[i] = targaccel[i+1];
      targvel[i] = targvel[i+1];
      targpos[i] = targpos[i+1];
    }
    // Load a new point into temptime
    temptime[polyorder+1] = outMJD[latestpoint];
    // Fit for acceleration
    // x component of acceleration.
    ppxvec={};
    ppyvec={};
    for(i=0;i<polyorder+1;i++) {
      ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
      ppyvec.push_back(targaccel[i].x);
    }
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
    // Fit for y component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
    // Fit for z component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
    // Re-calculate all of the positions and velocities using this fit.
    for(j=1;j<=polyorder+1;j++) {
      dt2 = (temptime[j]-temptime[0])*SOLARDAY;
      // Positions
      targpos[j].x = targpos[0].x + targvel[0].x*dt2;
      targpos[j].y = targpos[0].y + targvel[0].y*dt2;
      targpos[j].z = targpos[0].z + targvel[0].z*dt2;
      for(i=0;i<polyorder+1;i++) {
	timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	targpos[j].x += accelfit[i].x*timemult;
	targpos[j].y += accelfit[i].y*timemult;
	targpos[j].z += accelfit[i].z*timemult;
      }
      // Velocities
      targvel[j].x = targvel[0].x;
      targvel[j].y = targvel[0].y;
      targvel[j].z = targvel[0].z;
      for(i=0;i<polyorder+1;i++) {
	timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	targvel[j].x += accelfit[i].x*timemult;
	targvel[j].y += accelfit[i].y*timemult;
	targvel[j].z += accelfit[i].z*timemult;
      }
      // Accelerations
      accelmod[j].x = 0L;
      accelmod[j].y = 0L;
      accelmod[j].z = 0L;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	accelmod[j].x += accelfit[i].x*timemult;
	accelmod[j].y += accelfit[i].y*timemult;
	accelmod[j].z += accelfit[i].z*timemult;
      }
    }
    // Re-calculate accelerations using these revised positions
    for(j=1;j<=polyorder+1;j++) {
      planetsonce={};
      nplanetgrab01LD(latestpoint+j-polyorder-1, planetnum, planetmjd, planetpos, planetsonce);
      accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
    }

    // Use these revised accelerations to re-do the fits
    // x component of acceleration.
    ppxvec={};
    ppyvec={};
    for(i=0;i<polyorder+1;i++) {
      ppxvec.push_back((temptime[i] - temptime[0])*SOLARDAY/dt0);
      ppyvec.push_back(targaccel[i].x);
    }
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].x = ppfitvec[i];
    // Fit for y component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].y);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].y = ppfitvec[i];
    // Fit for z component of acceleration. Note that we have
    // already loaded the time vector ppxvec.
    ppyvec={};
    for(i=0;i<polyorder+1;i++) ppyvec.push_back(targaccel[i].z);
    // Perform fit, and store in accelfit.
    perfectpoly01LD(ppxvec,ppyvec,ppfitvec);
    for(i=0;i<polyorder+1;i++) accelfit[i].z = ppfitvec[i];
    // Re-calculate all of the positions and velocities using this fit.
    for(j=1;j<=polyorder+1;j++) {
      dt2 = (temptime[j]-temptime[0])*SOLARDAY;
      // Positions
      targpos[j].x = targpos[0].x + targvel[0].x*dt2;
      targpos[j].y = targpos[0].y + targvel[0].y*dt2;
      targpos[j].z = targpos[0].z + targvel[0].z*dt2;
      for(i=0;i<polyorder+1;i++) {
	timemult = intpowLD(dt2,2+i)*factorialLD(i)/factorialLD(2+i)/intpowLD(dt0,i);
	targpos[j].x += accelfit[i].x*timemult;
	targpos[j].y += accelfit[i].y*timemult;
	targpos[j].z += accelfit[i].z*timemult;
      }
      // Velocities
      targvel[j].x = targvel[0].x;
      targvel[j].y = targvel[0].y;
      targvel[j].z = targvel[0].z;
      for(i=0;i<polyorder+1;i++) {
	timemult = intpowLD(dt2,1+i)/intpowLD(dt0,i)/((long double)(1+i));
	targvel[j].x += accelfit[i].x*timemult;
	targvel[j].y += accelfit[i].y*timemult;
	targvel[j].z += accelfit[i].z*timemult;
      }
      // Accelerations
      accelmod[j].x = 0L;
      accelmod[j].y = 0L;
      accelmod[j].z = 0L;
      for(i=0;i<stepsin;i++) {
	timemult = intpowLD(dt2,i)/intpowLD(dt0,i);
	accelmod[j].x += accelfit[i].x*timemult;
	accelmod[j].y += accelfit[i].y*timemult;
	accelmod[j].z += accelfit[i].z*timemult;
      }
    }
    // Re-calculate accelerations using these revised positions
    for(j=1;j<=polyorder+1;j++) {
      planetsonce={};
      nplanetgrab01LD(latestpoint+j-polyorder-1, planetnum, planetmjd, planetpos, planetsonce);
      accelcalc01LD(planetnum, planetmasses, planetsonce, targpos[j], targaccel[j]);
      outvel[latestpoint+j-polyorder-1] = targvel[j];
      outpos[latestpoint+j-polyorder-1] = targpos[j];
    }

    // We have now gone through two iterations of extrapolation
    // to predict the next acceleration point as accurately as possible.
    // The next step of the loop will move the extrapolated point back
    // by one step, and use it to start extrapolating a new point,
    // at the same time refining the former extrapolated points.
  }
 
  return(0);
}


#undef DEBUG

double gaussian_deviate()
{
  long double x,y,rsq;
  double g1;
  rsq = 0.0L;
  do {
    x = 2.0L*rand()/RAND_MAX - 1.0L;
    y = 2.0L*rand()/RAND_MAX - 1.0L;
    rsq = x*x + y*y;
  } while(rsq>=1.0L || rsq<=0.0L);
  g1 = sqrt(-2.0l*log(rsq)/rsq);
  return(x*g1);
}

int uvw_to_galcoord(const double &u, const double &v, const double &w, double &RA, double &Dec)
{
  double vtot = sqrt(u*u + v*v + w*w);
  double sinedec = 0l;

  if(vtot==0) {
    RA=Dec=0l;
    return(1);
  }
  sinedec = w/vtot;
  if(sinedec>1.0l) {
    if(WARN_INVERSE_TRIG>0) cout << "Warning: attempting to take arcsine of w/vtot = 1.0 + " << sinedec-1.0l << "\n";
    Dec = 90.0l;
  } else if(sinedec<-1.0l) {
    if(WARN_INVERSE_TRIG>0) cout << "Warning: attempting to take arcsine of w/vtot = -1.0 - " << -sinedec-1.0l << "\n";
    Dec = -90.0l;
  } else if(isnormal(sinedec) || sinedec==0.0l) {
    Dec = asin(sinedec)*180.0l/M_PI;
  } else {
    cerr << "ERROR: bad value for sinedec: w/vtot = " << w << "/" << vtot << " = " << sinedec << "\n";
    RA=Dec=0l;
    return(-1);
  }
  
  if(v==0.0l && u>=0.0l) {
    RA = 0.0l;
  } else if(v==0.0l && u<0.0l) {
    RA = 180.0l;
  } else if(v>0.0l) {
    RA = 90.0l - atan(u/v)*180.0l/M_PI;
  } else if(v<0.0l) {
    RA = 270.0l - atan(u/v)*180.0l/M_PI;
  } else {
    cerr << "ERROR: bad value for u or v, u = " << u << ", v = " << v << "\n";
    RA=0l;
    return(-1);
  }
  return(0);
}

// unitvar: May 05, 2022, generator of random variable
// uniformly distributed from 0 to 1. Must be initialized
// by calling function from a string seed as follows:
//   seed_seq seed (stringseed.begin(),stringseed.end());
//   mt19937_64 generator (seed);
// Must be initialized only once: then any number of calls
// may be made, not only to this function but also to any
// others that use the mt19937_64 generator (e.g.,
// gaussian_deviate_mt. It is an error to initialize the
// generator more than once (e.g., in the entirely of main()).
long double unitvar(mt19937_64 &generator)
{
  long double uv = (long double)(generator())/RAND_MAX_64;
  return(uv);
}

double gaussian_deviate_mt(mt19937_64 &generator)
{
  long double x,y,rsq;
  double g1;
  rsq = 0.0L;
  do {
    x = 2.0L*(long double)(generator())/RAND_MAX_64 - 1.0L;
    y = 2.0L*(long double)(generator())/RAND_MAX_64 - 1.0L;
    rsq = x*x + y*y;
  } while(rsq>=1.0L || rsq<=0.0L);
  g1 = sqrt(-2.0l*log(rsq)/rsq);
  return(x*g1);
}

// multilinfit01: June 24, 2022
// Finds the WEIGHTED least-squares fit modeling the input vector
// yvec (length pnum) as a linear combination of fitnum other
// vectors supplied in the matrix xmat (size fitnum x pnum). The
// vector of best-fit coefficients for xmat is given in avec.*/
int multilinfit01(const vector <double> &yvec, const vector <double> &sigvec, const vector <vector <double>> &xmat, int pnum, int fitnum, vector <double> &avec)
{
  vector <vector <double>> fitmat;
  int pct,fitct,k;
  int verbose=0;
  fflush(stdout);

  make_dmat(fitnum,fitnum+1,fitmat);
  avec={};
  make_dvec(fitnum,avec);
  
  // Load fitmat for input into solvematrix01
  for(fitct=0;fitct<fitnum;fitct++) {
    // First the constant term -- that is, the term that does not
    // multiply any of the fitting coefficients -- which is also
    // the only term that involves yvec
    fitmat[fitct][0]=0.0;
    for(pct=0;pct<pnum;pct++) {
      if(isnormal(sigvec[pct])) fitmat[fitct][0] -= yvec[pct]*xmat[fitct][pct]/DSQUARE(sigvec[pct]);
    }
    /*Now the actual coefficients*/
    for(k=0;k<fitnum;k++) {
      fitmat[fitct][k+1] = 0.0;
      for(pct=0;pct<pnum;pct++) {
	if(isnormal(sigvec[pct])) fitmat[fitct][k+1] += xmat[fitct][pct]*xmat[k][pct]/DSQUARE(sigvec[pct]);
      }
    }
  }
  solvematrix01(fitmat,fitnum,avec,verbose);
  return(0);
}

// polyfit01: June 24, 2022:
// Finds the WEIGHTED least-squares fit modeling the input vector
// yvec (length pnum) as a polynomial of order polyorder in the
// input vector xvec (length pnum), and outputs the coefficients
// in order from constant to highest-order in the vector avec.
// vectors supplied in the matrix xmat (size fitnum x pnum). The
// vector of best-fit coefficients for xmat is given in avec.*/
int polyfit01(const vector <double> &yvec, const vector <double> &sigvec, const vector <double> &xvec, int pnum, int polyorder, vector <double> &avec)
{
  vector <vector <double>> xmat;
  int pct,fitct;
  avec = {};
  make_dmat(polyorder+1,pnum,xmat);
  
  for(pct=0; pct<pnum; pct++) {
    xmat[0][pct] = 1.0;
    for(fitct=1; fitct<=polyorder; fitct++) {
      xmat[fitct][pct] = intpowD(xvec[pct],fitct);
    }
  }
  multilinfit01(yvec, sigvec, xmat, pnum, polyorder+1, avec);
  return(0);
}

// vaneproj01LD: Given a unit vector unitbary giving the direction
// toward which an object was seen, find its intersection with a
// heliocentric vane of constant ecliptic longitude.
int vaneproj01LD(point3LD unitbary, point3LD obsbary, long double ecliplon, long double &geodist, point3LD &projbary)
{
  long double normdot1,normaldist;
  point3LD plane_normvec = point3LD(0L,0L,0L);
  point3LD plane_to_obs = point3LD(0L,0L,0L);
  
  // 1. The input heliocentric ecliptic longitude defines a plane.
  //    Calculate the unit vector normal from the sun.
  plane_normvec = point3LD(-sin(ecliplon/DEGPRAD),cos(ecliplon/DEGPRAD),0L);

  // 2. We already have the instantaneous position of the observer in obsbary
  // 3. Find the point on the plane closest to the observer, as follows:
  // 3a. Calculate dot-product of the sun-observer vector and the plane normal.
  //     This is the distance from the observer to the nearest point on the plane
  normaldist = dotprod3LD(obsbary,plane_normvec);
  if(!isnormal(normaldist)) return(-1); // Mainly this is to catch the case that the observer
                                        // is already in the plane, in which case the
                                        // dot product is exactly zero and fails the
                                        // isnormal test.
    
  // 3c. Multiply the plane-normal unit vector by the resulting physical length
  plane_to_obs.x = normaldist*plane_normvec.x;
  plane_to_obs.y = normaldist*plane_normvec.y;
  plane_to_obs.z = normaldist*plane_normvec.z;
    
  // 4. Calculate the normalized dot product of the vector from the observer to periobs
  //     and the observation unit vector. Reject the point if the normalized dot product
  //     is too small.
  normdot1 = -dotprod3LD(plane_to_obs,unitbary)/fabs(normaldist);
  if(normdot1<=0L) return(-1); // Observer was looking away from the plane
  
  // 5. Divide the length of the vector from the observer to periobs by the normalized
  //    dot product.
  geodist = fabs(normaldist)/normdot1;
  
  // 6. Mutiply the observation unit vector by the resulting physical length.
  // 7. Add the resulting physical vector to the instantaneous position of the observer.
  projbary.x = obsbary.x + unitbary.x*geodist;
  projbary.y = obsbary.y + unitbary.y*geodist;
  projbary.z = obsbary.z + unitbary.z*geodist;
  
  return(0);
}

// 2pt_KepQ: October 21, 2022:
// Calculate and return the value of the function Q used in Section 6.11
// of J. M. A. Danby's Foundations of Celestial Mechanics, in the context
// of solving the two-point boundary value problem for a Kepler orbit.
// This function is not particularly profound or magical, as can be
// seen in the source code below.
long double Twopoint_KepQ(long double x)
{
  if(x>0L && x<=1.0L) return(1.5L/x/sqrt(x)*(0.5L*asin(2.0L*x-1.0L)*0.5L - sqrt(x-x*x) + M_PI/4.0L));
  else if(x<=0L) return(0.75L/(-x*sqrt(-x))*(2.0L*sqrt(x*x-x) - log(1.0L - 2.0L*x + 2.0L*sqrt(x*x-x))));
  else {
    cerr << "ERROR:  Twopoint_KepQ called with out-of-range argument " << x << "\n";
    return(-99);
  }
  //return(0.75L*(theta - sin(theta))/intpowLD(sin(0.5L*theta),3));
}

#define DEBUG_2PTBVP 0

// 2pt_Kepler_vel: October 21, 2022:
// Given two points in an object's orbit (as 3-D Cartesian
// vectors relative to the sun), and the time it takes to
// move from the first point to the second, solve for the object's
// Keplerian orbit, and in particular find its vector velocity
// when at the first point. Input positions are in units of km,
// timediff is in units of days, and the output velocity will
// be in km/sec.
// This code closely follows the derivation in Section 6.11 of
// J. M. A. Danby's Foundations of Celestial Mechanics.
int Twopoint_Kepler_vel(const long double MGsun, const point3LD startpoint, const point3LD endpoint, const long double timediff, point3LD &startvel, int itmax)
{
  long double r1 = vecabs3LD(startpoint);
  long double r2 = vecabs3LD(endpoint);
  point3LD pdiff = point3LD(endpoint.x - startpoint.x, endpoint.y - startpoint.y, endpoint.z - startpoint.z);
  long double c = vecabs3LD(pdiff);
  long double lambda1 = sqrt(r1+r2+c);
  long double lambda2 = sqrt(r1+r2-c);
  long double k = sqrt(MGsun);

  // Determine the sign-specifier X
  long double ac = (r1+r2+c)/4.0L;
  long double nc = k/ac/sqrt(ac); // This is in radians per second.
  long double dc = 2.0L*asin(sqrt(lambda2/lambda1));
  long double dtc = (M_PI - dc + sin(dc))/nc;
  long double X=1.0L;
  long double Y=1.0L;
  
  if(timediff*SOLARDAY>dtc) X=-1.0L;

  cout << "Q(0.5) = " << Twopoint_KepQ(0.5) << ", Q(-0.5) = " << Twopoint_KepQ(-0.5) << "\n";
  
  cout << "Initial setup stuff:\n";
  cout << "r1 = " << r1/AU_KM << ", r2 = " << r2/AU_KM << ", c = " << c/AU_KM << ", lambdas = " << lambda1 << ", " << lambda2 << ", k = " << k << ", X = " << X << "\n";
  
  long double z = 2.0L/lambda1/lambda1; // Initial guess for z = 1/a
  // Calculate original value of function to be minimized.
  long double Q1 = Twopoint_KepQ(z*lambda1*lambda1/4.0L);
  long double Q2 = Twopoint_KepQ(z*lambda2*lambda2/4.0L);
  long double f = (1.0L/k)*((1.0L/6.0L)*(X*Q1*lambda1*lambda1*lambda1 - Y*Q2*lambda2*lambda2*lambda2) + (1.0L-X)*M_PI/(z*sqrt(z))) - timediff*SOLARDAY;
  long double fprime = (1.0L/(4.0L*k*z))*(X*lambda1*lambda1*lambda1*(1.0L/sqrt(1.0L - lambda1*lambda1*z/4.0L) - Q1) - Y*lambda2*lambda2*lambda2*(1.0L/sqrt(1.0L - lambda2*lambda2*z/4.0L) - Q2)) - (1.0L-X)*3.0L*M_PI/(2.0L*k)/(z*z*sqrt(z));

  int itnum=0;
  long double deltaz = -f/fprime;

  cout << "0th iteration:\n";
  cout << "z = " << z << " = 1/" << 1.0L/z/AU_KM << " AU, Q1 = " << Q1 << ", Q2 = " << Q2 << ", f = " << f << ", fprime = " << fprime << ", deltaz = " << deltaz << "\n";
  
  while(fabs(f) > KEPTRANSTOL && itnum<itmax) {
    cout << itnum << " f = " << f << ": z = " << z << ", a = " << 1.0L/z << "\n";
    z += deltaz;
    Q1 = Twopoint_KepQ(z*lambda1*lambda1/4.0L);
    Q2 = Twopoint_KepQ(z*lambda2*lambda2/4.0L);
    f = (1.0L/k)*((1.0L/6.0L)*(X*Q1*lambda1*lambda1*lambda1 - Y*Q2*lambda2*lambda2*lambda2) + (1.0L-X)*M_PI/z/sqrt(z)) - timediff*86400.0L;
    fprime = (1.0L/(4.0L*k*z))*(X*lambda1*lambda1*lambda1*(1.0L/sqrt(1.0L - lambda1*lambda1*z/4.0L) - Q1) - Y*lambda2*lambda2*lambda2*(1.0L/sqrt(1.0L - lambda2*lambda2*z/4.0L) - Q2)) - (1.0L-X)*3.0L*M_PI/(2.0L*k)/(z*z*sqrt(z));

    deltaz = -f/fprime;
    cout << "l12z = " << lambda1*lambda1*z << ", l22z = " << lambda2*lambda2*z << "\n";
    cout << "z = " << z << " = 1/" << 1.0L/z/AU_KM << " AU, Q1 = " << Q1 << ", Q2 = " << Q2 << ", f = " << f << ", fprime = " << fprime << ", deltaz = " << deltaz << "\n";
    itnum++;
  }
  cout << itnum << " f = " << f << ": z = " << z << ", a = " << 1.0L/z << "\n";
  if(itnum==itmax) cerr << "WARNING: Twopoint_Kepler_vel reached iteration limit " << itmax << " with deltaz = " << deltaz << "\n";
  cout << "z = " << z << ", a = " << 1.0L/z << "\n";
  return(1);
}

// Keplerint_multipoint01: November 02, 2022: Like Keplerint, but does the
// calculation for a bunch of points simultaneously. Note that
// we assume the observation times and mjdstart are in UT1, which
// means that JPL Horizons state vectors cannot be used directly
// for mjdstart: one would have to correct the nominal value of
// mjdstart corresponding to the JPL Horizons state vectors.
// Description of ancestor program Keplerint:
// Integrate an orbit assuming we have a Keplerian 2-body problem
// with all the mass in the Sun, and the input position and velocity
// are relative to the Sun.
int Keplerint_multipoint01(const long double MGsun, const long double mjdstart, const vector <long double> &obsMJD, const point3LD &startpos, const point3LD &startvel, vector <point3LD> &obspos, vector <point3LD> &obsvel)
{
  long double e,E,a,lscalar,r0,v0,r1,v1;
  point3LD lvec = point3LD(0L,0L,0L);
  point3LD r1unit = point3LD(0L,0L,0L);
  point3LD v1unit = point3LD(0L,0L,0L);
  point3LD targpos = point3LD(0L,0L,0L);
  point3LD targvel = point3LD(0L,0L,0L);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0L;
  long double costheta,theta0,theta1,radvel,cospsi,psi;
  costheta = theta0 = theta1 = radvel = cospsi = psi = 0L;
  long double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0L;
  long double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0L;
  long double sinev,thetav,v1ra,v1dec;
  sinev = thetav = v1ra = v1dec = 0L;
  int obsct=0;
  int obsnum = obsMJD.size();
 
  // Calculate scalar input position
  r0 = vecabs3LD(startpos);
  v0 = vecabs3LD(startvel);
  
  // Calculate specific energy and angular momentum
  E = 0.5L*v0*v0 - MGsun/r0;
  lvec = crossprod3LD(startpos,startvel);
  lscalar = sqrt(dotprod3LD(lvec,lvec));
  if(E>=0L) {
    //cerr << "ERROR: Keplerint finds positive total energy: " << E << "\n";
    return(1);
  }
		 
  // Calculate a and e: orbital semimajor axis and eccentricity.
  cout.precision(17);
  a = -MGsun*0.5L/E;
  e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/MGsun/MGsun);
  if(e<0L || e>=1.0L) {
    cerr << "ERROR: Keplerint finds eccentricity out of range: " << e << "\n";
    return(1);
  }

  //cout << "semimajor axis = " << a/AU_KM << " and eccentricity = " << e << "\n";
	       
  // Calculate angle theta0 from perihelion using simple ellipse geometry.
  if(e>0L) {
    costheta = ((a-a*e*e)/r0 - 1.0L)/e;
    if(costheta>=-1.0L && costheta<=1.0L) theta0 = acos(costheta);
    else {
      cerr << "ERROR: Keplerint finds costheta = " << costheta << "\n";
      return(1);
    }
  }
  radvel = dotprod3LD(startpos,startvel)/r0;
  //cout << "Radial velocity = " << radvel << " km/sec\n";
  
  if(radvel>=0) {
    // We are moving outward from perihelion: theta will be correct.
    //cout << "Moving outward from perihelion.\n";
    ;
  } else {
    // We are moving inward towards perihelion: theta needs adjustment.
    theta0 = 2.0L*M_PI - theta0;
    //cout << "Moving inward towards perihelion.\n";
  }
  //cout << "theta0 = " << theta0*DEGPRAD << "\n";
  
  // Calculate Goldstein's psi variable from theta.
  cospsi = (costheta + e)/(1.0L + costheta*e);
  if(cospsi>=-1.0L && cospsi<=1.0L) psi = acos(cospsi);
  else {
    cerr << "ERROR: Keplerint finds cospsi = " << cospsi << "\n";
    return(1);
  }
  if(radvel<0) {
    // We are moving inward towards perihelion: psi needs adjustment.
    psi = 2.0L*M_PI - psi;
  }
  //cout << "psi = " << psi*DEGPRAD << "\n";
 
  // Calculate time since perihelion using psi.
  omega = sqrt(MGsun/(a*a*a));
  //cout << "Period = " << 2.0L*M_PI/omega/SOLARDAY/365.25 << " years\n";
  t0omega = psi - e*sin(psi);
  //cout << "t0omega = " << t0omega;

  // Loop on all times-of-observation, and calculate the target position at those times
  obspos = obsvel = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    // The new time t1 for which we want to re-evaluate psi is
    // given by t0 + obsMJD[obsct]-mjdstart.
    t1omega = t0omega + (obsMJD[obsct]-mjdstart)*SOLARDAY*omega;
    while(t1omega > 2.0L*M_PI) t1omega -= 2.0L*M_PI;
    while(t1omega < 0.0L) t1omega += 2.0L*M_PI;
    // Solve Kepler's equation for psi(t1)
    psi = kep_transcendental(t1omega,e,KEPTRANSTOL);
    cospsi = cos(psi);
    // Calculate theta(t1) from psi(t1)
    if(1.0L - e*cospsi != 0.0L) {
      costheta = (cospsi - e)/(1.0L - e*cospsi);
      if(costheta >= -1.0L && costheta <= 1.0L) theta1 = acos(costheta);
      else if(costheta < -1.0L) {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = M_PI;
      } else {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = 0.0L;
      }
      if(psi>M_PI && theta1<=M_PI) theta1 = 2.0L*M_PI - theta1;
    } else {
      cerr << "Warning: e*cos(psi) = " << e*cospsi << " so 1 - e*cos(psi) = " << 1.0L - e*cospsi << "\n";
      theta1 = 0.0L;
    }
    while(theta1<0.0L) theta1 += 2.0L*M_PI;
    while(theta1>=2.0L*M_PI) theta1 -= 2.0L*M_PI;

    // Calculate r(t1) from psi(t1)
    r1 = a*(1.0L - e*cospsi);
    // Calculate v1 from r1 and the known energy
    v1 = sqrt((E +  MGsun/r1)*2.0L);
  
    // Use vector algebra to find the full vector r(t1).
    // This vector is perpendicular to lvec, and is angled by theta1-theta0
    // relative to startpos.
    // Convert angular momentum vector to spherical coordinates
    celedeproj01LD(lvec,&lra,&ldec); // Note that output is in degrees.
    celedeproj01LD(startpos,&r0ra,&r0dec); // Note that output is in degrees.
    // Transform the starting unit vector into a coordinate system with
    // the angular momentum vector at the pole, and the old pole at RA=0
    poleswitch01LD(r0ra/DEGPRAD,r0dec/DEGPRAD,lra/DEGPRAD,ldec/DEGPRAD,0.0L,newra,newdec); // Output is radians
    // Rotate starting unit vector around the angular momentum axis by
    // the calculated angle.
    newra += theta1-theta0;
    // The unit vector for the new position r1 is on the equator at this RA,
    // in the coordinate system that has the angular momentum vector at the pole.
    // Convert back to the original coordinate system.
    poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,r1ra,r1dec); // Output is radians
    // Now for the velocity. If the velocity is at right angle to the vector r1,
    // the product v1*r1 is the angular momentum. Otherwise, l/(v1*r1) is the sine
    // of the angle between v1 and r1.

    sinev = lscalar/v1/r1;
    if(sinev>=1.0L) thetav = 0.5L*M_PI;
    else if(sinev<0.0L) {
      cerr << "ERROR: negative angular momentum?\nv1,r1,v1*r1,lscalar,sinev = " << v1 << ", " << r1 << ", " << v1*r1 << ", " << lscalar << ", " << sinev << "\n";
      thetav = 0.0L;
    }
    else thetav = asin(sinev);
    if(theta1<=M_PI) {
      // Outward bound from perihelion.
      newra += thetav;
    } else {
      // Inward bound to perihelion
      newra += (M_PI - thetav);
    }
    poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,v1ra,v1dec); // Output is radians

    r1unit = celeproj01LD(r1ra*DEGPRAD,r1dec*DEGPRAD);
    v1unit =celeproj01LD(v1ra*DEGPRAD,v1dec*DEGPRAD);
  
    targpos.x = r1unit.x*r1;
    targpos.y = r1unit.y*r1;
    targpos.z = r1unit.z*r1;
    targvel.x = v1unit.x*v1;
    targvel.y = v1unit.y*v1;
    targvel.z = v1unit.z*v1;
    obspos.push_back(targpos);
    obsvel.push_back(targvel);
  }
    
  return(0);
}

// Keplerint_multipoint02: November 08, 2022:
// Like Keplerint_multipoint01, but relays several of its
// internally-calculated orbital parameters to the calling function. 
// These include the semimajor axis (in km), the eccentricity, 
// and the initial angle-from-perihelion (in radians).
//
// Description of related program Keplerint_multipoint01:
// Like Keplerint, but does the
// calculation for a bunch of points simultaneously. Note that
// we assume the observation times and mjdstart are in UT1, which
// means that JPL Horizons state vectors cannot be used directly
// for mjdstart: one would have to correct the nominal value of
// mjdstart corresponding to the JPL Horizons state vectors.
// Description of ancestor program Keplerint:
// Integrate an orbit assuming we have a Keplerian 2-body problem
// with all the mass in the Sun, and the input position and velocity
// are relative to the Sun.
int Keplerint_multipoint02(const long double MGsun, const long double mjdstart, const vector <long double> &obsMJD, const point3LD &startpos, const point3LD &startvel, vector <point3LD> &obspos, vector <point3LD> &obsvel, long double *semimajor_axis, long double *eccen, long double *angperi)
{
  long double e,E,a,lscalar,r0,v0,r1,v1;
  point3LD lvec = point3LD(0L,0L,0L);
  point3LD r1unit = point3LD(0L,0L,0L);
  point3LD v1unit = point3LD(0L,0L,0L);
  point3LD targpos = point3LD(0L,0L,0L);
  point3LD targvel = point3LD(0L,0L,0L);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0L;
  long double costheta,theta0,theta1,radvel,cospsi,psi;
  costheta = theta0 = theta1 = radvel = cospsi = psi = 0L;
  long double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0L;
  long double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0L;
  long double sinev,thetav,v1ra,v1dec;
  sinev = thetav = v1ra = v1dec = 0L;
  int obsct=0;
  int obsnum = obsMJD.size();
 
  // Calculate scalar input position
  r0 = vecabs3LD(startpos);
  v0 = vecabs3LD(startvel);
  
  // Calculate specific energy and angular momentum
  E = 0.5L*v0*v0 - MGsun/r0;
  lvec = crossprod3LD(startpos,startvel);
  lscalar = sqrt(dotprod3LD(lvec,lvec));
  if(E>=0L) {
    //cerr << "ERROR: Keplerint finds positive total energy: " << E << "\n";
    return(1);
  }
		 
  // Calculate a and e: orbital semimajor axis and eccentricity.
  cout.precision(17);
  a = -MGsun*0.5L/E;
  e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/MGsun/MGsun);
  if(e<0L || e>=1.0L) {
    cerr << "ERROR: Keplerint finds eccentricity out of range: " << e << "\n";
    return(1);
  }

  //cout << "semimajor axis = " << a/AU_KM << " and eccentricity = " << e << "\n";
	       
  // Calculate angle theta0 from perihelion using simple ellipse geometry.
  if(e>0L) {
    costheta = ((a-a*e*e)/r0 - 1.0L)/e;
    if(costheta>=-1.0L && costheta<=1.0L) theta0 = acos(costheta);
    else {
      cerr << "ERROR: Keplerint finds costheta = " << costheta << "\n";
      return(1);
    }
  }
  radvel = dotprod3LD(startpos,startvel)/r0;
  //cout << "Radial velocity = " << radvel << " km/sec\n";
  
  if(radvel>=0) {
    // We are moving outward from perihelion: theta will be correct.
    //cout << "Moving outward from perihelion.\n";
    ;
  } else {
    // We are moving inward towards perihelion: theta needs adjustment.
    theta0 = 2.0L*M_PI - theta0;
    //cout << "Moving inward towards perihelion.\n";
  }
  //cout << "theta0 = " << theta0*DEGPRAD << "\n";
  *semimajor_axis = a;
  *eccen = e;
  *angperi = theta0;
  
  // Calculate Goldstein's psi variable from theta.
  cospsi = (costheta + e)/(1.0L + costheta*e);
  if(cospsi>=-1.0L && cospsi<=1.0L) psi = acos(cospsi);
  else {
    cerr << "ERROR: Keplerint finds cospsi = " << cospsi << "\n";
    return(1);
  }
  if(radvel<0) {
    // We are moving inward towards perihelion: psi needs adjustment.
    psi = 2.0L*M_PI - psi;
  }
  //cout << "psi = " << psi*DEGPRAD << "\n";
 
  // Calculate time since perihelion using psi.
  omega = sqrt(MGsun/(a*a*a));
  //cout << "Period = " << 2.0L*M_PI/omega/SOLARDAY/365.25 << " years\n";
  t0omega = psi - e*sin(psi);
  //cout << "t0omega = " << t0omega;

  // Loop on all times-of-observation, and calculate the target position at those times
  obspos = obsvel = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    // The new time t1 for which we want to re-evaluate psi is
    // given by t0 + obsMJD[obsct]-mjdstart.
    t1omega = t0omega + (obsMJD[obsct]-mjdstart)*SOLARDAY*omega;
    while(t1omega > 2.0L*M_PI) t1omega -= 2.0L*M_PI;
    while(t1omega < 0.0L) t1omega += 2.0L*M_PI;
    // Solve Kepler's equation for psi(t1)
    psi = kep_transcendental(t1omega,e,KEPTRANSTOL);
    cospsi = cos(psi);
    // Calculate theta(t1) from psi(t1)
    if(1.0L - e*cospsi != 0.0L) {
      costheta = (cospsi - e)/(1.0L - e*cospsi);
      if(costheta >= -1.0L && costheta <= 1.0L) theta1 = acos(costheta);
      else if(costheta < -1.0L) {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = M_PI;
      } else {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = 0.0L;
      }
      if(psi>M_PI && theta1<=M_PI) theta1 = 2.0L*M_PI - theta1;
    } else {
      cerr << "Warning: e*cos(psi) = " << e*cospsi << " so 1 - e*cos(psi) = " << 1.0L - e*cospsi << "\n";
      theta1 = 0.0L;
    }
    while(theta1<0.0L) theta1 += 2.0L*M_PI;
    while(theta1>=2.0L*M_PI) theta1 -= 2.0L*M_PI;

    // Calculate r(t1) from psi(t1)
    r1 = a*(1.0L - e*cospsi);
    // Calculate v1 from r1 and the known energy
    v1 = sqrt((E +  MGsun/r1)*2.0L);
  
    // Use vector algebra to find the full vector r(t1).
    // This vector is perpendicular to lvec, and is angled by theta1-theta0
    // relative to startpos.
    // Convert angular momentum vector to spherical coordinates
    celedeproj01LD(lvec,&lra,&ldec); // Note that output is in degrees.
    celedeproj01LD(startpos,&r0ra,&r0dec); // Note that output is in degrees.
    // Transform the starting unit vector into a coordinate system with
    // the angular momentum vector at the pole, and the old pole at RA=0
    poleswitch01LD(r0ra/DEGPRAD,r0dec/DEGPRAD,lra/DEGPRAD,ldec/DEGPRAD,0.0L,newra,newdec); // Output is radians
    // Rotate starting unit vector around the angular momentum axis by
    // the calculated angle.
    newra += theta1-theta0;
    // The unit vector for the new position r1 is on the equator at this RA,
    // in the coordinate system that has the angular momentum vector at the pole.
    // Convert back to the original coordinate system.
    poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,r1ra,r1dec); // Output is radians
    // Now for the velocity. If the velocity is at right angle to the vector r1,
    // the product v1*r1 is the angular momentum. Otherwise, l/(v1*r1) is the sine
    // of the angle between v1 and r1.

    sinev = lscalar/v1/r1;
    if(sinev>=1.0L) thetav = 0.5L*M_PI;
    else if(sinev<0.0L) {
      cerr << "ERROR: negative angular momentum?\nv1,r1,v1*r1,lscalar,sinev = " << v1 << ", " << r1 << ", " << v1*r1 << ", " << lscalar << ", " << sinev << "\n";
      thetav = 0.0L;
    }
    else thetav = asin(sinev);
    if(theta1<=M_PI) {
      // Outward bound from perihelion.
      newra += thetav;
    } else {
      // Inward bound to perihelion
      newra += (M_PI - thetav);
    }
    poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,v1ra,v1dec); // Output is radians

    r1unit = celeproj01LD(r1ra*DEGPRAD,r1dec*DEGPRAD);
    v1unit =celeproj01LD(v1ra*DEGPRAD,v1dec*DEGPRAD);
  
    targpos.x = r1unit.x*r1;
    targpos.y = r1unit.y*r1;
    targpos.z = r1unit.z*r1;
    targvel.x = v1unit.x*v1;
    targvel.y = v1unit.y*v1;
    targvel.z = v1unit.z*v1;
    obspos.push_back(targpos);
    obsvel.push_back(targvel);
  }
    
  return(0);
}

int Keplerint_multipoint02(const double MGsun, const double mjdstart, const vector <double> &obsMJD, const point3d &startpos, const point3d &startvel, vector <point3d> &obspos, vector <point3d> &obsvel, double *semimajor_axis, double *eccen, double *angperi)
{
  double e,E,a,lscalar,r0,v0,r1,v1;
  point3d lvec = point3d(0L,0L,0L);
  point3d r1unit = point3d(0L,0L,0L);
  point3d v1unit = point3d(0L,0L,0L);
  point3d targpos = point3d(0L,0L,0L);
  point3d targvel = point3d(0L,0L,0L);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0L;
  double costheta,theta0,theta1,radvel,cospsi,psi;
  costheta = theta0 = theta1 = radvel = cospsi = psi = 0L;
  double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0L;
  double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0L;
  double sinev,thetav,v1ra,v1dec;
  sinev = thetav = v1ra = v1dec = 0L;
  int obsct=0;
  int obsnum = obsMJD.size();
 
  // Calculate scalar input position
  r0 = vecabs3d(startpos);
  v0 = vecabs3d(startvel);
  
  // Calculate specific energy and angular momentum
  E = 0.5L*v0*v0 - MGsun/r0;
  lvec = crossprod3d(startpos,startvel);
  lscalar = sqrt(dotprod3d(lvec,lvec));
  if(E>=0L) {
    //cerr << "ERROR: Keplerint finds positive total energy: " << E << "\n";
    return(1);
  }
		 
  // Calculate a and e: orbital semimajor axis and eccentricity.
  cout.precision(17);
  a = -MGsun*0.5L/E;
  e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/MGsun/MGsun);
  if(e<0L || e>=1.0L) {
    cerr << "ERROR: Keplerint finds eccentricity out of range: " << e << "\n";
    return(1);
  }

  //cout << "semimajor axis = " << a/AU_KM << " and eccentricity = " << e << "\n";
	       
  // Calculate angle theta0 from perihelion using simple ellipse geometry.
  if(e>0L) {
    costheta = ((a-a*e*e)/r0 - 1.0L)/e;
    if(costheta>=-1.0L && costheta<=1.0L) theta0 = acos(costheta);
    else {
      cerr << "ERROR: Keplerint finds costheta = " << costheta << "\n";
      return(1);
    }
  }
  radvel = dotprod3d(startpos,startvel)/r0;
  //cout << "Radial velocity = " << radvel << " km/sec\n";
  
  if(radvel>=0) {
    // We are moving outward from perihelion: theta will be correct.
    //cout << "Moving outward from perihelion.\n";
    ;
  } else {
    // We are moving inward towards perihelion: theta needs adjustment.
    theta0 = 2.0L*M_PI - theta0;
    //cout << "Moving inward towards perihelion.\n";
  }
  //cout << "theta0 = " << theta0*DEGPRAD << "\n";
  *semimajor_axis = a;
  *eccen = e;
  *angperi = theta0;
  
  // Calculate Goldstein's psi variable from theta.
  cospsi = (costheta + e)/(1.0L + costheta*e);
  if(cospsi>=-1.0L && cospsi<=1.0L) psi = acos(cospsi);
  else {
    cerr << "ERROR: Keplerint finds cospsi = " << cospsi << "\n";
    return(1);
  }
  if(radvel<0) {
    // We are moving inward towards perihelion: psi needs adjustment.
    psi = 2.0L*M_PI - psi;
  }
  //cout << "psi = " << psi*DEGPRAD << "\n";
 
  // Calculate time since perihelion using psi.
  omega = sqrt(MGsun/(a*a*a));
  //cout << "Period = " << 2.0L*M_PI/omega/SOLARDAY/365.25 << " years\n";
  t0omega = psi - e*sin(psi);
  //cout << "t0omega = " << t0omega;

  // Loop on all times-of-observation, and calculate the target position at those times
  obspos = obsvel = {};
  for(obsct=0;obsct<obsnum;obsct++) {
    // The new time t1 for which we want to re-evaluate psi is
    // given by t0 + obsMJD[obsct]-mjdstart.
    t1omega = t0omega + (obsMJD[obsct]-mjdstart)*SOLARDAY*omega;
    while(t1omega > 2.0L*M_PI) t1omega -= 2.0L*M_PI;
    while(t1omega < 0.0L) t1omega += 2.0L*M_PI;
    // Solve Kepler's equation for psi(t1)
    psi = kep_transcendental(t1omega,e,KEPTRANSTOL2);
    cospsi = cos(psi);
    // Calculate theta(t1) from psi(t1)
    if(1.0L - e*cospsi != 0.0L) {
      costheta = (cospsi - e)/(1.0L - e*cospsi);
      if(costheta >= -1.0L && costheta <= 1.0L) theta1 = acos(costheta);
      else if(costheta < -1.0L) {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = M_PI;
      } else {
	cout << "Warning: costheta = " << costheta << "\n";
	theta1 = 0.0L;
      }
      if(psi>M_PI && theta1<=M_PI) theta1 = 2.0L*M_PI - theta1;
    } else {
      cerr << "Warning: e*cos(psi) = " << e*cospsi << " so 1 - e*cos(psi) = " << 1.0L - e*cospsi << "\n";
      theta1 = 0.0L;
    }
    while(theta1<0.0L) theta1 += 2.0L*M_PI;
    while(theta1>=2.0L*M_PI) theta1 -= 2.0L*M_PI;

    // Calculate r(t1) from psi(t1)
    r1 = a*(1.0L - e*cospsi);
    // Calculate v1 from r1 and the known energy
    v1 = sqrt((E +  MGsun/r1)*2.0L);
  
    // Use vector algebra to find the full vector r(t1).
    // This vector is perpendicular to lvec, and is angled by theta1-theta0
    // relative to startpos.
    // Convert angular momentum vector to spherical coordinates
    celedeproj01(lvec,&lra,&ldec); // Note that output is in degrees.
    celedeproj01(startpos,&r0ra,&r0dec); // Note that output is in degrees.
    // Transform the starting unit vector into a coordinate system with
    // the angular momentum vector at the pole, and the old pole at RA=0
    poleswitch01(r0ra/DEGPRAD,r0dec/DEGPRAD,lra/DEGPRAD,ldec/DEGPRAD,0.0L,newra,newdec); // Output is radians
    // Rotate starting unit vector around the angular momentum axis by
    // the calculated angle.
    newra += theta1-theta0;
    // The unit vector for the new position r1 is on the equator at this RA,
    // in the coordinate system that has the angular momentum vector at the pole.
    // Convert back to the original coordinate system.
    poleswitch01(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,r1ra,r1dec); // Output is radians
    // Now for the velocity. If the velocity is at right angle to the vector r1,
    // the product v1*r1 is the angular momentum. Otherwise, l/(v1*r1) is the sine
    // of the angle between v1 and r1.

    sinev = lscalar/v1/r1;
    if(sinev>=1.0L) thetav = 0.5L*M_PI;
    else if(sinev<0.0L) {
      cerr << "ERROR: negative angular momentum?\nv1,r1,v1*r1,lscalar,sinev = " << v1 << ", " << r1 << ", " << v1*r1 << ", " << lscalar << ", " << sinev << "\n";
      thetav = 0.0L;
    }
    else thetav = asin(sinev);
    if(theta1<=M_PI) {
      // Outward bound from perihelion.
      newra += thetav;
    } else {
      // Inward bound to perihelion
      newra += (M_PI - thetav);
    }
    poleswitch01(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,v1ra,v1dec); // Output is radians

    r1unit = celeproj01(r1ra*DEGPRAD,r1dec*DEGPRAD);
    v1unit =celeproj01(v1ra*DEGPRAD,v1dec*DEGPRAD);
  
    targpos.x = r1unit.x*r1;
    targpos.y = r1unit.y*r1;
    targpos.z = r1unit.z*r1;
    targvel.x = v1unit.x*v1;
    targvel.y = v1unit.y*v1;
    targvel.z = v1unit.z*v1;
    obspos.push_back(targpos);
    obsvel.push_back(targvel);
  }
    
  return(0);
}


// orbitchi01: November 02, 2022:
// Get chi-square value based on input state vectors,
// using 2-body Keplerian integration, rather than n-body.
// Input state vectors are expected to be in km and km/sec.
// Note that the output vectors fitRA, fitDec, and resid are
// null-wiped inside orbitchi01, so it isn't necessary for the
// calling function to wipe them.
long double orbitchi01(const point3LD &objectpos, const point3LD &objectvel, const long double mjdstart, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid)
{
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  int obsct;
  int obsnum = obsMJD.size();
  long double light_travel_time;
  point3LD outpos = point3LD(0,0,0);
  long double outRA=0L;
  long double outDec=0L;
  long double ldval=0L;
  long double chisq=0L;
  double dval;
  resid = fitRA = fitDec = {};
  int status=0;

  if(DEBUG_2PTBVP>1) cout << "Input start pos: " << objectpos.x << " "  << objectpos.y << " "  << objectpos.z << "\n";
  
  // Integrate orbit.
  status=0;
  status = Keplerint_multipoint01(GMSUN_KM3_SEC2,mjdstart,obsMJD,objectpos,objectvel,obspos,obsvel);
  if(status!=0) {
    // Keplerint_multipoint01 failed, likely because input state vectors lead
    // to an unbound orbit.
    return(LARGERR);
  }
  if(DEBUG_2PTBVP>1) cout << "Recovered start pos: " << obspos[0].x << " "  << obspos[0].y << " "  << obspos[0].z << "\n";


  for(obsct=0;obsct<obsnum;obsct++) {
    // Initial approximation of the coordinates relative to the observer
    outpos.x = obspos[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - observerpos[obsct].z;
    // Initial approximation of the observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Convert to meters and divide by the speed of light to get the light travel time.
    light_travel_time = ldval*1000.0/CLIGHT;
    // Light-travel-time corrected version of coordinates relative to the observer
    outpos.x = obspos[obsct].x - light_travel_time*obsvel[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - light_travel_time*obsvel[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - light_travel_time*obsvel[obsct].z - observerpos[obsct].z;
    // Output diagnostic stuff
    if(DEBUG_2PTBVP>1) {
      if(obsct == 0) cout << "Calculated start pos:\n" << obspos[obsct].x  - light_travel_time*obsvel[obsct].x << " " << obspos[obsct].y  - light_travel_time*obsvel[obsct].y << " " << obspos[obsct].z  - light_travel_time*obsvel[obsct].z << "\n";
      if(obsct == obsnum-1)  cout << "Calculated end pos:\n" << obspos[obsct].x  - light_travel_time*obsvel[obsct].x << " " << obspos[obsct].y  - light_travel_time*obsvel[obsct].y << " " << obspos[obsct].z  - light_travel_time*obsvel[obsct].z << "\n";
    }
    // Light-travel-time corrected observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    // Project onto the celestial sphere.
    stateunitLD_to_celestial(outpos, outRA, outDec);
    dval = distradec01(obsRA[obsct],obsDec[obsct],outRA,outDec);
    dval *= 3600.0L; // Convert to arcsec
    fitRA.push_back(outRA);
    fitDec.push_back(outDec);
    resid.push_back(dval);
  }
  chisq=0.0L;
  for(obsct=0;obsct<obsnum;obsct++) {
    chisq += LDSQUARE(resid[obsct]/sigastrom[obsct]);
  }
  return(chisq);
}

// orbitchi02: November 08, 2022:
// Like orbitchi01, but uses Keplerint_multipoint02
// rather than  Keplerint_multipoint01, in order to relay
// more of the internally calculated orbital parameters to
// the calling function. These include the orbital semimajor
// axis a (km), the eccentricity e, and the initial angle from
// perihelion (radians).
//
// Description of ancestor program orbitchi01:
// Get chi-square value based on input state vectors,
// using 2-body Keplerian integration, rather than n-body.
// Input state vectors are expected to be in km and km/sec.
// Note that the output vectors fitRA, fitDec, and resid are
// null-wiped inside orbitchi01, so it isn't necessary for the
// calling function to wipe them.
long double orbitchi02(const point3LD &objectpos, const point3LD &objectvel, const long double mjdstart, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid, long double *semimajor_axis, long double *eccen, long double *angperi)
{
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  int obsct;
  int obsnum = obsMJD.size();
  long double light_travel_time;
  point3LD outpos = point3LD(0,0,0);
  long double outRA=0L;
  long double outDec=0L;
  long double ldval=0L;
  long double chisq=0L;
  double dval;
  resid = fitRA = fitDec = {};
  int status=0;

  if(DEBUG_2PTBVP>1) cout << "Input start pos: " << objectpos.x << " "  << objectpos.y << " "  << objectpos.z << "\n";
  
  // Integrate orbit.
  status=0;
  status = Keplerint_multipoint02(GMSUN_KM3_SEC2,mjdstart,obsMJD,objectpos,objectvel,obspos,obsvel,semimajor_axis,eccen,angperi);
  if(status!=0) {
    // Keplerint_multipoint02 failed, likely because input state vectors lead
    // to an unbound orbit.
    return(LARGERR);
  }
  if(DEBUG_2PTBVP>1) cout << "Recovered start pos: " << obspos[0].x << " "  << obspos[0].y << " "  << obspos[0].z << "\n";


  for(obsct=0;obsct<obsnum;obsct++) {
    // Initial approximation of the coordinates relative to the observer
    outpos.x = obspos[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - observerpos[obsct].z;
    // Initial approximation of the observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Convert to meters and divide by the speed of light to get the light travel time.
    light_travel_time = ldval*1000.0/CLIGHT;
    // Light-travel-time corrected version of coordinates relative to the observer
    outpos.x = obspos[obsct].x - light_travel_time*obsvel[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - light_travel_time*obsvel[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - light_travel_time*obsvel[obsct].z - observerpos[obsct].z;
    // Output diagnostic stuff
    if(DEBUG_2PTBVP>1) {
      if(obsct == 0) cout << "Calculated start pos:\n" << obspos[obsct].x  - light_travel_time*obsvel[obsct].x << " " << obspos[obsct].y  - light_travel_time*obsvel[obsct].y << " " << obspos[obsct].z  - light_travel_time*obsvel[obsct].z << "\n";
      if(obsct == obsnum-1)  cout << "Calculated end pos:\n" << obspos[obsct].x  - light_travel_time*obsvel[obsct].x << " " << obspos[obsct].y  - light_travel_time*obsvel[obsct].y << " " << obspos[obsct].z  - light_travel_time*obsvel[obsct].z << "\n";
    }
    // Light-travel-time corrected observer-target distance
    ldval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= ldval;
    outpos.y /= ldval;
    outpos.z /= ldval;
    // Project onto the celestial sphere.
    stateunitLD_to_celestial(outpos, outRA, outDec);
    dval = distradec01(obsRA[obsct],obsDec[obsct],outRA,outDec);
    dval *= 3600.0L; // Convert to arcsec
    fitRA.push_back(outRA);
    fitDec.push_back(outDec);
    resid.push_back(dval);
  }
  chisq=0.0L;
  for(obsct=0;obsct<obsnum;obsct++) {
    chisq += LDSQUARE(resid[obsct]/sigastrom[obsct]);
  }
  return(chisq);
}

double orbitchi02(const point3d &objectpos, const point3d &objectvel, const double mjdstart, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid, double *semimajor_axis, double *eccen, double *angperi)
{
  vector <point3d> obspos;
  vector <point3d> obsvel;
  int obsct;
  int obsnum = obsMJD.size();
  double light_travel_time;
  point3d outpos = point3d(0,0,0);
  double outRA=0l;
  double outDec=0l;
  double dval=0l;
  double chisq=0l;
  resid = fitRA = fitDec = {};
  int status=0;

  if(DEBUG_2PTBVP>1) cout << "Input start pos: " << objectpos.x << " "  << objectpos.y << " "  << objectpos.z << "\n";
  
  // Integrate orbit.
  status=0;
  status = Keplerint_multipoint02(GMSUN_KM3_SEC2,mjdstart,obsMJD,objectpos,objectvel,obspos,obsvel,semimajor_axis,eccen,angperi);
  if(status!=0) {
    // Keplerint_multipoint02 failed, likely because input state vectors lead
    // to an unbound orbit.
    return(LARGERR);
  }
  if(DEBUG_2PTBVP>1) cout << "Recovered start pos: " << obspos[0].x << " "  << obspos[0].y << " "  << obspos[0].z << "\n";


  for(obsct=0;obsct<obsnum;obsct++) {
    // Initial approximation of the coordinates relative to the observer
    outpos.x = obspos[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - observerpos[obsct].z;
    // Initial approximation of the observer-target distance
    dval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Convert to meters and divide by the speed of light to get the light travel time.
    light_travel_time = dval*1000.0/CLIGHT;
    // Light-travel-time corrected version of coordinates relative to the observer
    outpos.x = obspos[obsct].x - light_travel_time*obsvel[obsct].x - observerpos[obsct].x;
    outpos.y = obspos[obsct].y - light_travel_time*obsvel[obsct].y - observerpos[obsct].y;
    outpos.z = obspos[obsct].z - light_travel_time*obsvel[obsct].z - observerpos[obsct].z;
    // Output diagnostic stuff
    if(DEBUG_2PTBVP>1) {
      if(obsct == 0) cout << "Calculated start pos:\n" << obspos[obsct].x  - light_travel_time*obsvel[obsct].x << " " << obspos[obsct].y  - light_travel_time*obsvel[obsct].y << " " << obspos[obsct].z  - light_travel_time*obsvel[obsct].z << "\n";
      if(obsct == obsnum-1)  cout << "Calculated end pos:\n" << obspos[obsct].x  - light_travel_time*obsvel[obsct].x << " " << obspos[obsct].y  - light_travel_time*obsvel[obsct].y << " " << obspos[obsct].z  - light_travel_time*obsvel[obsct].z << "\n";
    }
    // Light-travel-time corrected observer-target distance
    dval = sqrt(outpos.x*outpos.x + outpos.y*outpos.y + outpos.z*outpos.z);
    // Calculate unit vector
    outpos.x /= dval;
    outpos.y /= dval;
    outpos.z /= dval;
    // Project onto the celestial sphere.
    stateunit_to_celestial(outpos, outRA, outDec);
    dval = distradec01(obsRA[obsct],obsDec[obsct],outRA,outDec);
    dval *= 3600.0L; // Convert to arcsec
    fitRA.push_back(outRA);
    fitDec.push_back(outDec);
    resid.push_back(dval);
  }
  chisq=0.0L;
  for(obsct=0;obsct<obsnum;obsct++) {
    chisq += DSQUARE(resid[obsct]/sigastrom[obsct]);
  }
  return(chisq);
}


// TwopointF: October 26, 2022:
// Given input values for k = sqrt(GMsun), delta-t, lambda1 = sqrt(r1+r2+c),
// lambda2 = sqrt(r1+r2-c), and the semimajor axis a in km, evaluate a function
// of the semimajor axis a that, if minimized, solves the Kepler two-point
// boundary value problem. Distances are in units of km, but delta-t is in units
// of days.
long double TwopointF(long double a, long double k, long double lambda1, long double lambda2, long double deltat, long double Xsign, long double Ysign)
{
  long double kterm, asin1, sq1, asin2, sq2, Xterm;
  //  long double dkterm, dasin1, dsq1, dasin2, dsq2,a2;
  // Evaluate term-by-term, for clarity
  kterm = k*deltat*SOLARDAY/a/sqrt(a);
  asin1 = Xsign*2.0L*asin(lambda1/2.0L/sqrt(a));
  sq1 = Xsign*lambda1*sqrt(4.0L*a-lambda1*lambda1)/2.0L/a;
  asin2 = Ysign*2.0L*asin(lambda2/2.0L/sqrt(a));
  sq2 = Ysign*lambda2*sqrt(4.0L*a-lambda2*lambda2)/2.0L/a;
  if(Xsign==-1.0L) Xterm = -2.0L*M_PI;
  else Xterm = 0.0L;
  
  //a2 = a+1000.0L;
  //dkterm = k*deltat*SOLARDAY/a2/sqrt(a2) - k*deltat*SOLARDAY/a/sqrt(a);
  //dasin1 = 2.0L*asin(lambda1/2.0L/sqrt(a2)) - 2.0L*asin(lambda1/2.0L/sqrt(a)); 
  //dsq1 = lambda1*sqrt(4.0L*a2-lambda1*lambda1)/2.0L/a2 - lambda1*sqrt(4.0L*a-lambda1*lambda1)/2.0L/a;
  //dasin2 = 2.0L*asin(lambda2/2.0L/sqrt(a2)) - 2.0L*asin(lambda2/2.0L/sqrt(a));
  //dsq2 = lambda2*sqrt(4.0L*a2-lambda2*lambda2)/2.0L/a2 - lambda2*sqrt(4.0L*a-lambda2*lambda2)/2.0L/a;

  //dkterm /= 1000.0L;
  //dasin1 /= 1000.0L;
  //dsq1 /= 1000.0L;
  //dasin2 /= 1000.0L;
  //dsq2 /= 1000.0L;

  //cout << "Numder by term: " << dkterm << " " << dasin1 << " " << dsq1 << " " << dasin2 << " " << dsq2 << "\n";
  
  return(kterm - asin1 + sq1 + asin2 - sq2 +  Xterm);
}

double TwopointF(double a, double k, double lambda1, double lambda2, double deltat, double Xsign, double Ysign)
{
  double kterm, asin1, sq1, asin2, sq2, Xterm;
  //  long double dkterm, dasin1, dsq1, dasin2, dsq2,a2;
  // Evaluate term-by-term, for clarity
  kterm = k*deltat*SOLARDAY/a/sqrt(a);
  asin1 = Xsign*2.0L*asin(lambda1/2.0L/sqrt(a));
  sq1 = Xsign*lambda1*sqrt(4.0L*a-lambda1*lambda1)/2.0L/a;
  asin2 = Ysign*2.0L*asin(lambda2/2.0L/sqrt(a));
  sq2 = Ysign*lambda2*sqrt(4.0L*a-lambda2*lambda2)/2.0L/a;
  if(Xsign==-1.0L) Xterm = -2.0L*M_PI;
  else Xterm = 0.0L;
  
  //a2 = a+1000.0L;
  //dkterm = k*deltat*SOLARDAY/a2/sqrt(a2) - k*deltat*SOLARDAY/a/sqrt(a);
  //dasin1 = 2.0L*asin(lambda1/2.0L/sqrt(a2)) - 2.0L*asin(lambda1/2.0L/sqrt(a)); 
  //dsq1 = lambda1*sqrt(4.0L*a2-lambda1*lambda1)/2.0L/a2 - lambda1*sqrt(4.0L*a-lambda1*lambda1)/2.0L/a;
  //dasin2 = 2.0L*asin(lambda2/2.0L/sqrt(a2)) - 2.0L*asin(lambda2/2.0L/sqrt(a));
  //dsq2 = lambda2*sqrt(4.0L*a2-lambda2*lambda2)/2.0L/a2 - lambda2*sqrt(4.0L*a-lambda2*lambda2)/2.0L/a;

  //dkterm /= 1000.0L;
  //dasin1 /= 1000.0L;
  //dsq1 /= 1000.0L;
  //dasin2 /= 1000.0L;
  //dsq2 /= 1000.0L;

  //cout << "Numder by term: " << dkterm << " " << dasin1 << " " << dsq1 << " " << dasin2 << " " << dsq2 << "\n";
  
  return(kterm - asin1 + sq1 + asin2 - sq2 +  Xterm);
}


 // TwopointFprime: October 26, 2022:
// Given input values for k = sqrt(GMsun), delta-t, lambda1 = sqrt(r1+r2+c),
// lambda2 = sqrt(r1+r2-c), and the semimajor axis a in km, evaluate the
// derivative with respect to the semimajor axis 'a' a function that,
// if minimized, solves the Kepler two-point
// boundary value problem. Distances are in units of km, but delta-t is in units
// of days.
long double TwopointFprime(long double a, long double k, long double lambda1, long double lambda2, long double deltat, long double Xsign, long double Ysign)
{
  long double kterm, asin1, sq1, asin2, sq2, foural1, foural2;
  // Terms that appear a lot, pre-calculated for simplicity
  foural1 = 4.0L*a - lambda1*lambda1;
  foural2 = 4.0L*a - lambda2*lambda2;
  
  // Evaluate term-by-term, for clarity
  kterm = -3.0L/2.0L*k*deltat*SOLARDAY/a/a/sqrt(a);
  asin1 = Xsign*lambda1/a/sqrt(foural1);
  sq1 = Xsign*(lambda1*lambda1*lambda1 - 2.0L*a*lambda1)/(2.0L*a*a*sqrt(foural1));
  asin2 = -Ysign*lambda2/a/sqrt(foural2);
  sq2 = Ysign*(2.0L*a*lambda2 - lambda2*lambda2*lambda2)/(2.0L*a*a*sqrt(foural2));
  
  //cout << "Fprime by terms: " << kterm << " " << asin1 << " " << sq1 << " " << asin2 << " " << sq2 << " " << kterm + asin1 + sq1 + asin2 + sq2 << "\n";
  
  return(kterm + asin1 + sq1 + asin2 + sq2);
}

double TwopointFprime(double a, double k, double lambda1, double lambda2, double deltat, double Xsign, double Ysign)
{
  double kterm, asin1, sq1, asin2, sq2, foural1, foural2;
  // Terms that appear a lot, pre-calculated for simplicity
  foural1 = 4.0L*a - lambda1*lambda1;
  foural2 = 4.0L*a - lambda2*lambda2;
  
  // Evaluate term-by-term, for clarity
  kterm = -3.0L/2.0L*k*deltat*SOLARDAY/a/a/sqrt(a);
  asin1 = Xsign*lambda1/a/sqrt(foural1);
  sq1 = Xsign*(lambda1*lambda1*lambda1 - 2.0L*a*lambda1)/(2.0L*a*a*sqrt(foural1));
  asin2 = -Ysign*lambda2/a/sqrt(foural2);
  sq2 = Ysign*(2.0L*a*lambda2 - lambda2*lambda2*lambda2)/(2.0L*a*a*sqrt(foural2));
  
  //cout << "Fprime by terms: " << kterm << " " << asin1 << " " << sq1 << " " << asin2 << " " << sq2 << " " << kterm + asin1 + sq1 + asin2 + sq2 << "\n";
  
  return(kterm + asin1 + sq1 + asin2 + sq2);
}

  
// eccen_calc_fast: November 01, 2022:
// Given the semimajor axis of an ellipse, and two 3-D points
// on the ellipse, calculate the eccentricity of
// the ellipse, or -1.0 if there is no valid elliptical solution
// This version uses only the eccentric anomaly (psi), and not
// the true anomaly (theta), resulting in a single solution).
// The eccentric anomaly solution seemed on development tests
// to be slightly less accurate, but the differences appeared
// far too small to matter, and not worth the reduction in
// speed.
int eccen_calc_fast(long double a, point3LD rvec1, point3LD rvec2, long double *e, long double *theta, long double Xsign, long double Ysign)
{
  long double r1 = vecabs3LD(rvec1);
  long double r2 = vecabs3LD(rvec2);
  long double costheta,thetatest,r2test;
  long double cos_delta_theta = dotprod3LD(rvec1,rvec2)/r1/r2;
  long double delta_theta = acos(cos_delta_theta);
 
  // Parameters for solution using the eccentric anomaly
  long double c = sqrt(LDSQUARE(rvec2.x-rvec1.x) + LDSQUARE(rvec2.y-rvec1.y) + LDSQUARE(rvec2.z-rvec1.z));
  long double lambda1 = r1+r2+c;
  long double lambda2 = r1+r2-c;
  long double eps_star = 2.0L*asin(sqrt(lambda1/4.0L/a));
  long double delt_star = 2.0L*asin(sqrt(lambda2/4.0L/a));
  long double eps = Ysign*eps_star;
  long double delt = M_PI*(1.0L - Xsign) + Xsign*delt_star;
  long double delta_psi = eps - delt;
  long double amr1 = a-r1;
  long double amr2 = a-r2;
  long double cosdpsi = cos(delta_psi);
  long double cos2dpsi = LDSQUARE(cosdpsi);
  long double sindpsi = sin(delta_psi);
  long double sin2dpsi = LDSQUARE(sindpsi);
  
  long double eea = sqrt(cos2dpsi*amr1*amr1 - 2.0L*cosdpsi*amr1*amr2 + amr2*amr2 + sin2dpsi*amr1*amr1)/sindpsi/a;
  costheta = (a*(1.0L-eea*eea) - r1)/eea/r1;
  thetatest = acos(costheta);
  r2test = a*(1-eea*eea)/(1.0L +  eea*cos(thetatest + delta_theta));
  long double delta1 = fabs(r2test-r2);
  r2test = a*(1-eea*eea)/(1.0L +  eea*cos(2.0L*M_PI - thetatest + delta_theta));
  long double delta2 = fabs(r2test-r2);
  //cout << "eea deltas: " << delta1 << " " << delta2 << "\n";
  if(delta2<delta1) thetatest = 2.0L*M_PI - thetatest;
  //cout << "eea theta = " << thetatest*DEGPRAD << "\n";
  *e = eea;
  *theta = thetatest;

  return(0);
}

// Twopoint_Kepler_v1: November 01, 2022:
// Given two points in an object's orbit (as 3-D Cartesian
// vectors relative to the sun), and the time it takes to
// move from the first point to the second, solve for the 
// semimajor axis a of the object's Keplerian orbit, and
// from these derive the initial velocity in units of km/sec.
// Input positions are in units of km,
// timediff is in units of days.
// This code closely follows the derivation in Section 6.11 of
// J. M. A. Danby's Foundations of Celestial Mechanics.
//
// November 8: commented out calculation of eccentricity, which
// previously had been included in the argument list as
// long double *e, to increase speed.
point3LD Twopoint_Kepler_v1(const long double GMsun, const point3LD startpos, const point3LD endpos, const long double timediff, const long double Ysign, long double *a, int itmax, int verbose)
{
  // General quantities, and those having to do with solving for the semimajor axis.
  long double r1 = vecabs3LD(startpos);
  long double r2 = vecabs3LD(endpos);
  point3LD pdiff = point3LD(endpos.x - startpos.x, endpos.y - startpos.y, endpos.z - startpos.z);
  point3LD v1 = point3LD(0.0L,0.0L,0.0L);
  long double c = vecabs3LD(pdiff);
  long double lambda1 = sqrt(r1+r2+c);
  long double lambda2 = sqrt(r1+r2-c);
  long double k = sqrt(GMsun);

  // Catch the hyperbolic case
  long double dtp = (intpowLD(lambda1,3) - Ysign*intpowLD(lambda2,3))/6.0L/k;
  if(DEBUG_2PTBVP > 1) cout << "Checking for hyperbolic case: " << timediff << " <= " << dtp/SOLARDAY << "\n";
  if(timediff*SOLARDAY <= dtp) {
    if(verbose>=2) cerr << "ERROR: Twopoint_Kepler_v1 has hyperbolic case\n";
    if(verbose>=2) cerr << "Returning with velocity set to zero\n";
    return(v1);
  }
  
  // Determine the sign-specifier X
  long double ac = (r1+r2+c)/4.0L;
  long double nc = k/ac/sqrt(ac); // This is in radians per second.
  long double dc = 2.0L*asin(sqrt(lambda2/lambda1));
  long double dtc = (M_PI - dc + sin(dc))/nc;
  long double X=1.0L;
  long double Y=Ysign;
  long double aorb;
  long double delta_aorb;
  long double f;
  long double fprime;
  //long double eccen,thetaperi;
  
  if(timediff*SOLARDAY>dtc) X=-1.0L;

  if(DEBUG_2PTBVP>1) cout << "Initial setup stuff:\n";
  if(DEBUG_2PTBVP>1) cout << "r1 = " << r1/AU_KM << ", r2 = " << r2/AU_KM << ", c = " << c/AU_KM << ", lambdas = " << lambda1 << ", " << lambda2 << ", k = " << k << ", X = " << X << ", Y = " << Y << "\n";

  //for(ai=10;ai<=200;ai++) {
  //  aorb = (long double)ai*0.01L*AU_KM;
  //  if(aorb > lambda1*lambda1/4.0L) {
  //    cout << ai << " " << aorb/AU_KM << " " << TwopointF(aorb, k, lambda1, lambda2, timediff) << "\n";
  //  }
  //}

  aorb = 0.5L*r1 + 0.5L*r2;
  f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
  fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);

  int itnum=0;
  delta_aorb = -f/fprime;
  while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
  while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
    // Don't let it drop into the NAN range
    delta_aorb /= 2.0L;
  }
 
  if(DEBUG_2PTBVP>1) cout << "0th iteration:\n";
  if(DEBUG_2PTBVP>1) cout << "a = " << aorb << " = " << aorb/AU_KM << " AU, f = " << f << ", fprime = " << fprime << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
  
  while((fabs(f) > KEP2PBVPTOL || fabs(delta_aorb) > KEP2PBVPTOL*aorb) && aorb<1000000.0L*AU_KM && itnum<itmax) {
    aorb += delta_aorb;
    f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    f2 = TwopointF(aorb+1000.0L, k, lambda1, lambda2, timediff);
    fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    fprime2 = (f2-f)/1000.0L;
    delta_aorb = -f/fprime;
    while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
    while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
      // Don't let it drop into the NAN range
      delta_aorb /= 2.0L;
    }
    // Calculate the eccentricity
    // eccen_calc_precise(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    // cout << "iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = " << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    //eccen_calc_fast(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    //if(DEBUG_2PTBVP==1) {
    //  cout << "iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = ";
    //  cout << scientific << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    // }
    itnum++;
  }
  if(fabs(f) > KEP2PBVPTOL || fabs(delta_aorb) > KEP2PBVPTOL*aorb) {
    // We never found a solution. Re-start near the minimum allowable:
    aorb = lambda1*lambda1/3.99L;
    f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
    fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);
    itnum=0;
    delta_aorb = -f/fprime;
    while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
    while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
      // Don't let it drop into the NAN range
      delta_aorb /= 2.0L;
    }
 
    if(DEBUG_2PTBVP>1) cout << "2nd try, 0th iteration, :\n";
    if(DEBUG_2PTBVP>1) cout << "a = " << aorb << " = " << aorb/AU_KM << " AU, f = " << f << ", fprime = " << fprime << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    while((fabs(f) > KEP2PBVPTOL || fabs(delta_aorb) > KEP2PBVPTOL*aorb)  && aorb<1000000.0L*AU_KM && itnum<itmax) {
    aorb += delta_aorb;
    f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    f2 = TwopointF(aorb+1000.0L, k, lambda1, lambda2, timediff);
    fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    fprime2 = (f2-f)/1000.0L;
    delta_aorb = -f/fprime;
    while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
    while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
      // Don't let it drop into the NAN range
      delta_aorb /= 2.0L;
    }
    // Calculate the eccentricity
    // eccen_calc_precise(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    // cout << "2nd try, iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = " << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    //eccen_calc_fast(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    //if(DEBUG_2PTBVP==1) {
    //  cout << "2nd try, iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = ";
    //  cout << scientific << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    //}
    itnum++;
    }
    if(fabs(f) > KEP2PBVPTOL || fabs(delta_aorb) > KEP2PBVPTOL*aorb) {
      if(verbose>=2) cerr << "ERROR: Twopoint_Kepler_v1 failed to converge\n";
      if(verbose>=2) cerr << "Returning with velocity set to zero\n";
      *a=-1.0L;
      return(v1);
    }
  }
  *a = aorb;
  //*e = eccen;

  // Quantities having to do with solving for the velocity
  long double alpha,beta,gamma,ffunc,gfunc;
  
  // Calculate the velocity
  alpha = 0.5L*lambda1;
  beta = 0.5L*lambda2;
  gamma = alpha*sqrt(1.0L - lambda2*lambda2/4.0L/aorb) -  X*Y*beta*sqrt(1.0L - lambda1*lambda1/4.0L/aorb);
  ffunc = 1.0L - 2*gamma*gamma/r1;
  gfunc = (4.0L/k)*Y*alpha*beta*gamma;

  v1.x = (1.0L/gfunc)*(endpos.x - ffunc*startpos.x);
  v1.y = (1.0L/gfunc)*(endpos.y - ffunc*startpos.y);
  v1.z = (1.0L/gfunc)*(endpos.z - ffunc*startpos.z);
  return(v1);
}

point3d Twopoint_Kepler_v1(const double GMsun, const point3d startpos, const point3d endpos, const double timediff, const double Ysign, double *a, int itmax, int verbose)
{
  // General quantities, and those having to do with solving for the semimajor axis.
  double r1 = vecabs3d(startpos);
  double r2 = vecabs3d(endpos);
  point3d pdiff = point3d(endpos.x - startpos.x, endpos.y - startpos.y, endpos.z - startpos.z);
  point3d v1 = point3d(0.0L,0.0L,0.0L);
  double c = vecabs3d(pdiff);
  double lambda1 = sqrt(r1+r2+c);
  double lambda2 = sqrt(r1+r2-c);
  double k = sqrt(GMsun);

  // Catch the hyperbolic case
  double dtp = (intpowD(lambda1,3) - Ysign*intpowD(lambda2,3))/6.0L/k;
  if(DEBUG_2PTBVP > 1) cout << "Checking for hyperbolic case: " << timediff << " <= " << dtp/SOLARDAY << "\n";
  if(timediff*SOLARDAY <= dtp) {
    if(verbose>=2) cerr << "ERROR: Twopoint_Kepler_v1 has hyperbolic case\n";
    if(verbose>=2) cerr << "Returning with velocity set to zero\n";
    return(v1);
  }
  
  // Determine the sign-specifier X
  double ac = (r1+r2+c)/4.0L;
  double nc = k/ac/sqrt(ac); // This is in radians per second.
  double dc = 2.0L*asin(sqrt(lambda2/lambda1));
  double dtc = (M_PI - dc + sin(dc))/nc;
  double X=1.0L;
  double Y=Ysign;
  double aorb;
  double delta_aorb;
  double f;
  double fprime;
  //double eccen,thetaperi;
  
  if(timediff*SOLARDAY>dtc) X=-1.0L;

  if(DEBUG_2PTBVP>1) cout << "Initial setup stuff:\n";
  if(DEBUG_2PTBVP>1) cout << "r1 = " << r1/AU_KM << ", r2 = " << r2/AU_KM << ", c = " << c/AU_KM << ", lambdas = " << lambda1 << ", " << lambda2 << ", k = " << k << ", X = " << X << ", Y = " << Y << "\n";

  //for(ai=10;ai<=200;ai++) {
  //  aorb = (double)ai*0.01L*AU_KM;
  //  if(aorb > lambda1*lambda1/4.0L) {
  //    cout << ai << " " << aorb/AU_KM << " " << TwopointF(aorb, k, lambda1, lambda2, timediff) << "\n";
  //  }
  //}

  aorb = 0.5L*r1 + 0.5L*r2;
  f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
  fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);

  int itnum=0;
  delta_aorb = -f/fprime;
  while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
  while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
    // Don't let it drop into the NAN range
    delta_aorb /= 2.0L;
  }
 
  if(DEBUG_2PTBVP>1) cout << "0th iteration:\n";
  if(DEBUG_2PTBVP>1) cout << "a = " << aorb << " = " << aorb/AU_KM << " AU, f = " << f << ", fprime = " << fprime << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
  
  while((fabs(f) > KEP2PBVPTOL2 || fabs(delta_aorb) > KEP2PBVPTOL2*aorb) && aorb<1000000.0L*AU_KM && itnum<itmax) {
    aorb += delta_aorb;
    f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    f2 = TwopointF(aorb+1000.0L, k, lambda1, lambda2, timediff);
    fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    fprime2 = (f2-f)/1000.0L;
    delta_aorb = -f/fprime;
    while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
    while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
      // Don't let it drop into the NAN range
      delta_aorb /= 2.0L;
    }
    // Calculate the eccentricity
    // eccen_calc_precise(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    // cout << "iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = " << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    //eccen_calc_fast(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    //if(DEBUG_2PTBVP==1) {
    //  cout << "iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = ";
    //  cout << scientific << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    // }
    itnum++;
  }
  if(fabs(f) > KEP2PBVPTOL2 || fabs(delta_aorb) > KEP2PBVPTOL2*aorb) {
    // We never found a solution. Re-start near the minimum allowable:
    aorb = lambda1*lambda1/3.99L;
    f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
    fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);
    itnum=0;
    delta_aorb = -f/fprime;
    while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
    while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
      // Don't let it drop into the NAN range
      delta_aorb /= 2.0L;
    }
 
    if(DEBUG_2PTBVP>1) cout << "2nd try, 0th iteration, :\n";
    if(DEBUG_2PTBVP>1) cout << "a = " << aorb << " = " << aorb/AU_KM << " AU, f = " << f << ", fprime = " << fprime << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    while((fabs(f) > KEP2PBVPTOL2 || fabs(delta_aorb) > KEP2PBVPTOL2*aorb)  && aorb<1000000.0L*AU_KM && itnum<itmax) {
    aorb += delta_aorb;
    f = TwopointF(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    f2 = TwopointF(aorb+1000.0L, k, lambda1, lambda2, timediff);
    fprime = TwopointFprime(aorb, k, lambda1, lambda2, timediff, X, Y);
    //    fprime2 = (f2-f)/1000.0L;
    delta_aorb = -f/fprime;
    while(fabs(delta_aorb/aorb)>0.1) delta_aorb/=10.0L; // Don't change too much, might go negative.
    while(delta_aorb < 0.0L && aorb+delta_aorb < lambda1*lambda1/4.0L) {
      // Don't let it drop into the NAN range
      delta_aorb /= 2.0L;
    }
    // Calculate the eccentricity
    // eccen_calc_precise(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    // cout << "2nd try, iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = " << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    //eccen_calc_fast(aorb, startpos, endpos, &eccen, &thetaperi, X, Y);
    //if(DEBUG_2PTBVP==1) {
    //  cout << "2nd try, iteration " << itnum << ": a = " << aorb/AU_KM << " AU, e,theta = " << eccen << " " << thetaperi*DEGPRAD << ", f = ";
    //  cout << scientific << f << ", fprime = " << fprime << " = " << fprime2 << ", delta_aorb = " << delta_aorb/AU_KM << "\n";
    //}
    itnum++;
    }
    if(fabs(f) > KEP2PBVPTOL2 || fabs(delta_aorb) > KEP2PBVPTOL2*aorb) {
      if(verbose>=2) cerr << "ERROR: Twopoint_Kepler_v1 failed to converge\n";
      if(verbose>=2) cerr << "Returning with velocity set to zero\n";
      *a=-1.0L;
      return(v1);
    }
  }
  *a = aorb;
  //*e = eccen;

  // Quantities having to do with solving for the velocity
  double alpha,beta,gamma,ffunc,gfunc;
  
  // Calculate the velocity
  alpha = 0.5L*lambda1;
  beta = 0.5L*lambda2;
  gamma = alpha*sqrt(1.0L - lambda2*lambda2/4.0L/aorb) -  X*Y*beta*sqrt(1.0L - lambda1*lambda1/4.0L/aorb);
  ffunc = 1.0L - 2*gamma*gamma/r1;
  gfunc = (4.0L/k)*Y*alpha*beta*gamma;

  v1.x = (1.0L/gfunc)*(endpos.x - ffunc*startpos.x);
  v1.y = (1.0L/gfunc)*(endpos.y - ffunc*startpos.y);
  v1.z = (1.0L/gfunc)*(endpos.z - ffunc*startpos.z);
  return(v1);
}

// geodist_to_3Dpos01: November 02, 2022:
// Given and input RA, Dec position, observer's position, and
// geocentric distance, output a 3D position for the object.
// The geocentric distance is expected in AU, but the output
// position vector has units of km.
point3LD geodist_to_3Dpos01(long double RA, long double Dec, point3LD observerpos, long double geodist)
{
  point3LD baryvec = point3LD(0.0L,0.0L,0.0L);
  
  celestial_to_stateunitLD(RA, Dec, baryvec);
  baryvec.x = observerpos.x + baryvec.x*geodist*AU_KM;
  baryvec.y = observerpos.y + baryvec.y*geodist*AU_KM;
  baryvec.z = observerpos.z + baryvec.z*geodist*AU_KM;
  return(baryvec);
}

point3d geodist_to_3dpos01(double RA, double Dec, point3d observerpos, double geodist)
{
  point3d baryvec = point3d(0.0l,0.0l,0.0l);
  
  celestial_to_stateunit(RA, Dec, baryvec);
  baryvec.x = observerpos.x + baryvec.x*geodist*AU_KM;
  baryvec.y = observerpos.y + baryvec.y*geodist*AU_KM;
  baryvec.z = observerpos.z + baryvec.z*geodist*AU_KM;
  return(baryvec);
}

// Herget_unboundcheck01: November 03, 2022:
// Quickly check if input parameters for a Method of Herget orbit fit
// imply an unbound (hyperbolic) orbit.
int Herget_unboundcheck01(long double geodist1, long double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec)
{
  long numobs = long(obsMJD.size());
  if(long(obsRA.size()) != numobs || long(obsDec.size()) != numobs || long(observerpos.size()) != numobs) {
    cerr << "ERROR: Hergetchi01 finds unequal lenths among input vectors:\n";
    cerr << "observed MJD, RA, Dec, sigastrom, and observerpos have lengths " << numobs << " " << obsRA.size() << " " << obsDec.size() << " " << observerpos.size() << "\n";
    return(-1);
  }
  if(Hergetpoint2<=Hergetpoint1 || Hergetpoint1<0 || Hergetpoint2>=numobs) {
    cerr << "ERROR: Hergetchi01 has invalid input reference points:\n";
    cerr << "Starting point " << Hergetpoint1 << " and ending point " << Hergetpoint2 << ", where allowed range is 0 to " << numobs-1 << "\n";
    return(-1);
  }
  
  point3LD startpos = geodist_to_3Dpos01(obsRA[Hergetpoint1], obsDec[Hergetpoint1], observerpos[Hergetpoint1], geodist1);
  point3LD endpos = geodist_to_3Dpos01(obsRA[Hergetpoint2], obsDec[Hergetpoint2], observerpos[Hergetpoint2], geodist2);
  // Time difference should include a light-travel-time correction. The sign is determined
  // by the fact that if the object gets further away, the object time moves backward
  // relative to the observer time. Hence, if the object gets further away (i.e., geodist2>geodist1),
  // the object experiences less time than the observer beween the two observations, because
  // the observer is looking further back in time at the second observation.
  long double deltat = obsMJD[Hergetpoint2] - obsMJD[Hergetpoint1] - (geodist2-geodist1)/CLIGHT_AUDAY;

  long double r1 = vecabs3LD(startpos);
  long double r2 = vecabs3LD(endpos);
  point3LD pdiff = point3LD(endpos.x - startpos.x, endpos.y - startpos.y, endpos.z - startpos.z);
  long double c = vecabs3LD(pdiff);
  long double lambda1 = sqrt(r1+r2+c);
  long double lambda2 = sqrt(r1+r2-c);
  long double k = sqrt(GMSUN_KM3_SEC2);
  long double Ysign = 1.0L;
  
  // Catch the hyperbolic case
  long double dtp = (intpowLD(lambda1,3) - Ysign*intpowLD(lambda2,3))/6.0L/k;
  if(deltat*SOLARDAY <= dtp) {
    return(1);
  } else return(0);
}

int Herget_unboundcheck01(double geodist1, double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec)
{
  long numobs = long(obsMJD.size());
  if(long(obsRA.size()) != numobs || long(obsDec.size()) != numobs || long(observerpos.size()) != numobs) {
    cerr << "ERROR: Hergetchi01 finds unequal lenths among input vectors:\n";
    cerr << "observed MJD, RA, Dec, sigastrom, and observerpos have lengths " << numobs << " " << obsRA.size() << " " << obsDec.size() << " " << observerpos.size() << "\n";
    return(-1);
  }
  if(Hergetpoint2<=Hergetpoint1 || Hergetpoint1<0 || Hergetpoint2>=numobs) {
    cerr << "ERROR: Hergetchi01 has invalid input reference points:\n";
    cerr << "Starting point " << Hergetpoint1 << " and ending point " << Hergetpoint2 << ", where allowed range is 0 to " << numobs-1 << "\n";
    return(-1);
  }
  
  point3d startpos = geodist_to_3dpos01(obsRA[Hergetpoint1], obsDec[Hergetpoint1], observerpos[Hergetpoint1], geodist1);
  point3d endpos = geodist_to_3dpos01(obsRA[Hergetpoint2], obsDec[Hergetpoint2], observerpos[Hergetpoint2], geodist2);
  // Time difference should include a light-travel-time correction. The sign is determined
  // by the fact that if the object gets further away, the object time moves backward
  // relative to the observer time. Hence, if the object gets further away (i.e., geodist2>geodist1),
  // the object experiences less time than the observer beween the two observations, because
  // the observer is looking further back in time at the second observation.
  double deltat = obsMJD[Hergetpoint2] - obsMJD[Hergetpoint1] - (geodist2-geodist1)/CLIGHT_AUDAY;

  double r1 = vecabs3d(startpos);
  double r2 = vecabs3d(endpos);
  point3d pdiff = point3d(endpos.x - startpos.x, endpos.y - startpos.y, endpos.z - startpos.z);
  double c = vecabs3d(pdiff);
  double lambda1 = sqrt(r1+r2+c);
  double lambda2 = sqrt(r1+r2-c);
  double k = sqrt(GMSUN_KM3_SEC2);
  double Ysign = 1.0L;
  
  // Catch the hyperbolic case
  long double dtp = (intpowD(lambda1,3) - Ysign*intpowD(lambda2,3))/6.0L/k;
  if(deltat*SOLARDAY <= dtp) {
    return(1);
  } else return(0);
}

// Hergetchi01: November 02, 2022:
// Calculate the chi-square value between input observations and
// the result of orbit determination using the Method of Herget,
// where the Method of Herget is performed between Hergetpoint1 and
// Hergetpoint2 of the input observations, under the assumption that
// the object was a distance geodist1 from the observer when
// the observation corresponding to Hergetpoint1 was made, and
// a distance geodist2 from the observer when the observation
// corresponding to Hergetpoint2 was made. The vectors are all
// in units of km and km/sec, but geodist1 and geodist2 are
// in units of AU.
// Note that the output vectors fitRA, fitDec, and resid are
// null-wiped inside Hergetchi01, so it isn't necessary for the
// calling function to wipe them.
// The vector orbit holds a, e, mjd, and the state vectors, for now.
long double Hergetchi01(long double geodist1, long double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid, vector <long double> &orbit, int verbose)
{
  int i=0;
  long numobs = long(obsMJD.size());
  if(long(obsRA.size()) != numobs || long(obsDec.size()) != numobs || long(sigastrom.size()) != numobs || long(observerpos.size()) != numobs) {
    cerr << "ERROR: Hergetchi01 finds unequal lenths among input vectors:\n";
    cerr << "observed MJD, RA, Dec, sigastrom, and observerpos have lengths " << numobs << " " << obsRA.size() << " " << obsDec.size() << " " <<  sigastrom.size() << " " << observerpos.size() << "\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0L);
    return(LARGERR);
  }
  if(Hergetpoint2<=Hergetpoint1 || Hergetpoint1<0 || Hergetpoint2>=numobs) {
    cerr << "ERROR: Hergetchi01 has invalid input reference points:\n";
    cerr << "Starting point " << Hergetpoint1 << " and ending point " << Hergetpoint2 << ", where allowed range is 0 to " << numobs-1 << "\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0L);
    return(LARGERR);
  }
  
  point3LD startpos = geodist_to_3Dpos01(obsRA[Hergetpoint1], obsDec[Hergetpoint1], observerpos[Hergetpoint1], geodist1);
  point3LD endpos = geodist_to_3Dpos01(obsRA[Hergetpoint2], obsDec[Hergetpoint2], observerpos[Hergetpoint2], geodist2);
  // Time difference should include a light-travel-time correction. The sign is determined
  // by the fact that if the object gets further away, the object time moves backward
  // relative to the observer time. Hence, if the object gets further away (i.e., geodist2>geodist1),
  // the object experiences less time than the observer beween the two observations, because
  // the observer is looking further back in time at the second observation.
  long double deltat = obsMJD[Hergetpoint2] - obsMJD[Hergetpoint1] - (geodist2-geodist1)/CLIGHT_AUDAY;
  long double a,e,angperi;
  point3LD startvel = Twopoint_Kepler_v1(GMSUN_KM3_SEC2, startpos, endpos, deltat, 1.0L, &a, 100, verbose);
  if(startvel.x == 0.0L && startvel.y == 0.0L && startvel.z == 0.0L) {
    // This is a failure code for Twopoint_Kepler_v1, which is returned if, e.g.,
    // the implied orbit is hyperbolic. Nothing to be done here, but pass
    // the failure code up the call chain.
    if(verbose>=2) cerr << "ERROR: Hergetchi01 received failure code from Twopoint_Kepler_v1\n";
    if(verbose>=2) cerr << "On input distances " << geodist1 << " and " << geodist2 << "\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0L);
    return(LARGERR);
  }
    
  long double orbchi;
  // Note that the output vectors fitRA, fitDec, and resid are null-wiped
  // internally in orbitchi01, so it isn't necessary to wipe them here.
  orbchi = orbitchi02(startpos, startvel, obsMJD[Hergetpoint1]-geodist1/CLIGHT_AUDAY, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, &a, &e, &angperi);

  orbit={};
  orbit.push_back(a);
  orbit.push_back(e);
  orbit.push_back(obsMJD[Hergetpoint1]-geodist1/CLIGHT_AUDAY);
  orbit.push_back(startpos.x);
  orbit.push_back(startpos.y);
  orbit.push_back(startpos.z);
  orbit.push_back(startvel.x);
  orbit.push_back(startvel.y);
  orbit.push_back(startvel.z);
  
  if(DEBUG_2PTBVP>1) cout << "Target startpos:\n" << startpos.x << " " << startpos.y << " " << startpos.z << "\n";
  if(DEBUG_2PTBVP>1) cout << "Target endpos:\n" << endpos.x << " " << endpos.y << " " << endpos.z << "\n";
  
  return(orbchi);
}

double Hergetchi01(double geodist1, double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid, vector <double> &orbit, int verbose)
{
  int i=0;
  long numobs = long(obsMJD.size());
  if(long(obsRA.size()) != numobs || long(obsDec.size()) != numobs || long(sigastrom.size()) != numobs || long(observerpos.size()) != numobs) {
    cerr << "ERROR: Hergetchi01 finds unequal lenths among input vectors:\n";
    cerr << "observed MJD, RA, Dec, sigastrom, and observerpos have lengths " << numobs << " " << obsRA.size() << " " << obsDec.size() << " " <<  sigastrom.size() << " " << observerpos.size() << "\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0l);
    return(LARGERR2);
  }
  if(Hergetpoint2<=Hergetpoint1 || Hergetpoint1<0 || Hergetpoint2>=numobs) {
    cerr << "ERROR: Hergetchi01 has invalid input reference points:\n";
    cerr << "Starting point " << Hergetpoint1 << " and ending point " << Hergetpoint2 << ", where allowed range is 0 to " << numobs-1 << "\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0l);
    return(LARGERR2);
  }
  
  point3d startpos = geodist_to_3dpos01(obsRA[Hergetpoint1], obsDec[Hergetpoint1], observerpos[Hergetpoint1], geodist1);
  point3d endpos = geodist_to_3dpos01(obsRA[Hergetpoint2], obsDec[Hergetpoint2], observerpos[Hergetpoint2], geodist2);
  // Time difference should include a light-travel-time correction. The sign is determined
  // by the fact that if the object gets further away, the object time moves backward
  // relative to the observer time. Hence, if the object gets further away (i.e., geodist2>geodist1),
  // the object experiences less time than the observer beween the two observations, because
  // the observer is looking further back in time at the second observation.
  double deltat = obsMJD[Hergetpoint2] - obsMJD[Hergetpoint1] - (geodist2-geodist1)/CLIGHT_AUDAY;
  double a,e,angperi;
  point3d startvel = Twopoint_Kepler_v1(GMSUN_KM3_SEC2, startpos, endpos, deltat, 1.0l, &a, 100, verbose);
  if(startvel.x == 0.0L && startvel.y == 0.0L && startvel.z == 0.0L) {
    // This is a failure code for Twopoint_Kepler_v1, which is returned if, e.g.,
    // the implied orbit is hyperbolic. Nothing to be done here, but pass
    // the failure code up the call chain.
    if(verbose>=2) cerr << "ERROR: Hergetchi01 received failure code from Twopoint_Kepler_v1\n";
    if(verbose>=2) cerr << "On input distances " << geodist1 << " and " << geodist2 << "\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0l);
    return(LARGERR2);
  }
    
  long double orbchi;
  // Note that the output vectors fitRA, fitDec, and resid are null-wiped
  // internally in orbitchi01, so it isn't necessary to wipe them here.
  orbchi = orbitchi02(startpos, startvel, obsMJD[Hergetpoint1]-geodist1/CLIGHT_AUDAY, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, &a, &e, &angperi);

  orbit={};
  orbit.push_back(a);
  orbit.push_back(e);
  orbit.push_back(obsMJD[Hergetpoint1]-geodist1/CLIGHT_AUDAY);
  orbit.push_back(startpos.x);
  orbit.push_back(startpos.y);
  orbit.push_back(startpos.z);
  orbit.push_back(startvel.x);
  orbit.push_back(startvel.y);
  orbit.push_back(startvel.z);
  
  if(DEBUG_2PTBVP>1) cout << "Target startpos:\n" << startpos.x << " " << startpos.y << " " << startpos.z << "\n";
  if(DEBUG_2PTBVP>1) cout << "Target endpos:\n" << endpos.x << " " << endpos.y << " " << endpos.z << "\n";
  
  return(orbchi);
}


// Herget_simplex_init: November 09, 2022:
// Initialize a 2-D simplex suitable for downhill simplex orbit
// fitting using the Method of Herget. The parameter space has
// two dimensions: geodist1 and geodist2 -- that is, the distance
// from the Earth at the instant of the first and last observations
// of the object. The initialization uses three types, each based
// on a parameter simpscale which is expected to be strictly
// positive but less than (sqrt(5)-1)/2. The calling function is
// responsible to ensure an acceptable value for the inputs, though
// Herget_simplex_int will fix some types of bad input. simptype=0 uses multiplicative
// scaling to create an approximately equilateral triangle:
// (geodist1,geodist2), (geodist1*(1-simpscale),geodist2), (geodist1,geodist2*(1-simpscale)).
// simptype=1 creates a simplex elongated along the direction defined
// by geodist1=geodist2:
// (geodist1,geodist2),
// (geodist1*(1 + simpscale - simpscale^2),geodist2*(1 + simpscale + simpscale^2)),
// (geodist1*(1 - simpscale - simpscale^2),geodist2*(1 - simpscale + simpscale^2)).
// simptype=2 uses subtraction to create a precisely equilateral triangle:
// (geodist1,geodist2), (geodist1-simpscale,geodist2), (geodist2-simpscale,geodist1).
// simptype=3 creates an elongated simplex like simptype=1, but with
// sign-switching that apparently makes it work badly: this option exists
// only to probe details of the resulting bad behavior.
int Herget_simplex_int(long double geodist1, long double geodist2, long double simpscale, long double simplex[3][2], int simptype)
{
  if(simptype==0) {
    // Define initial simplex
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1*(1.0L - simpscale);
    // Make sure we don't have a negative distance.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5L*(simplex[1][0] + geodist1);
    simplex[1][1] = geodist2;
    simplex[2][0] = geodist1;
    simplex[2][1] = geodist2*(1.0L - simpscale);
    // Make sure we don't have a negative distance.
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5L*(simplex[2][1] + geodist2);
  } else if(simptype==1) {
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1*(1.0L + simpscale - LDSQUARE(simpscale));
    simplex[1][1] = geodist2*(1.0L + simpscale + LDSQUARE(simpscale));
    // Make sure we don't have negative distances.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5L*(simplex[1][0] + geodist1);
    while(simplex[1][1]<=0.0) simplex[1][1] = 0.5L*(simplex[1][1] + geodist2);
    simplex[2][0] = geodist1*(1.0L - simpscale - LDSQUARE(simpscale));
    simplex[2][1] = geodist2*(1.0L - simpscale + LDSQUARE(simpscale));
    // Make sure we don't have negative distances.
    while(simplex[2][0]<=0.0) simplex[2][0] = 0.5L*(simplex[2][0] + geodist1);
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5L*(simplex[2][1] + geodist2);
  } else if(simptype==2) {
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1 - simpscale;
    // Make sure we don't have a negative distance.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5L*(simplex[1][0] + geodist1);
    simplex[1][1] = geodist2;
    simplex[2][0] = geodist1;
    simplex[2][1] = geodist2- simpscale;
    // Make sure we don't have a negative distance.
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5L*(simplex[2][1] + geodist2);
  } else if(simptype==3) {
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1*(1.0L + simpscale + LDSQUARE(simpscale));
    simplex[1][1] = geodist2*(1.0L + simpscale - LDSQUARE(simpscale));
    // Make sure we don't have negative distances.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5L*(simplex[1][0] + geodist1);
    while(simplex[1][1]<=0.0) simplex[1][1] = 0.5L*(simplex[1][1] + geodist2);
    simplex[2][0] = geodist1*(1.0L - simpscale - LDSQUARE(simpscale));
    simplex[2][1] = geodist2*(1.0L - simpscale + LDSQUARE(simpscale));
    // Make sure we don't have negative distances.
    while(simplex[2][0]<=0.0) simplex[2][0] = 0.5L*(simplex[2][0] + geodist1);
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5L*(simplex[2][1] + geodist2);
  } else {
    // default to case 0, multiplicative approximately equilateral:
    // Define initial simplex
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1*(1.0L - simpscale);
    // Make sure we don't have a negative distance.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5L*(simplex[1][0] + geodist1);
    simplex[1][1] = geodist2;
    simplex[2][0] = geodist1;
    simplex[2][1] = geodist2*(1.0L - simpscale);
    // Make sure we don't have a negative distance.
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5L*(simplex[2][1] + geodist2);
  }
  return(0);
}

int Herget_simplex_int(double geodist1, double geodist2, double simpscale, double simplex[3][2], int simptype)
{
  if(simptype==0) {
    // Define initial simplex
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1*(1.0l - simpscale);
    // Make sure we don't have a negative distance.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5l*(simplex[1][0] + geodist1);
    simplex[1][1] = geodist2;
    simplex[2][0] = geodist1;
    simplex[2][1] = geodist2*(1.0l - simpscale);
    // Make sure we don't have a negative distance.
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5l*(simplex[2][1] + geodist2);
  } else if(simptype==1) {
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1*(1.0l + simpscale - DSQUARE(simpscale));
    simplex[1][1] = geodist2*(1.0l + simpscale + DSQUARE(simpscale));
    // Make sure we don't have negative distances.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5l*(simplex[1][0] + geodist1);
    while(simplex[1][1]<=0.0) simplex[1][1] = 0.5l*(simplex[1][1] + geodist2);
    simplex[2][0] = geodist1*(1.0L - simpscale - DSQUARE(simpscale));
    simplex[2][1] = geodist2*(1.0L - simpscale + DSQUARE(simpscale));
    // Make sure we don't have negative distances.
    while(simplex[2][0]<=0.0) simplex[2][0] = 0.5l*(simplex[2][0] + geodist1);
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5l*(simplex[2][1] + geodist2);
  } else if(simptype==2) {
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1 - simpscale;
    // Make sure we don't have a negative distance.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5l*(simplex[1][0] + geodist1);
    simplex[1][1] = geodist2;
    simplex[2][0] = geodist1;
    simplex[2][1] = geodist2- simpscale;
    // Make sure we don't have a negative distance.
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5l*(simplex[2][1] + geodist2);
  } else if(simptype==3) {
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1*(1.0l + simpscale + DSQUARE(simpscale));
    simplex[1][1] = geodist2*(1.0l + simpscale - DSQUARE(simpscale));
    // Make sure we don't have negative distances.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5l*(simplex[1][0] + geodist1);
    while(simplex[1][1]<=0.0) simplex[1][1] = 0.5l*(simplex[1][1] + geodist2);
    simplex[2][0] = geodist1*(1.0L - simpscale - DSQUARE(simpscale));
    simplex[2][1] = geodist2*(1.0L - simpscale + DSQUARE(simpscale));
    // Make sure we don't have negative distances.
    while(simplex[2][0]<=0.0) simplex[2][0] = 0.5l*(simplex[2][0] + geodist1);
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5l*(simplex[2][1] + geodist2);
  } else {
    // default to case 0, multiplicative approximately equilateral:
    // Define initial simplex
    simplex[0][0] = geodist1;
    simplex[0][1] = geodist2;
    simplex[1][0] = geodist1*(1.0L - simpscale);
    // Make sure we don't have a negative distance.
    while(simplex[1][0]<=0.0) simplex[1][0] = 0.5l*(simplex[1][0] + geodist1);
    simplex[1][1] = geodist2;
    simplex[2][0] = geodist1;
    simplex[2][1] = geodist2*(1.0l - simpscale);
    // Make sure we don't have a negative distance.
    while(simplex[2][1]<=0.0) simplex[2][1] = 0.5l*(simplex[2][1] + geodist2);
  }
  return(0);
}
  

#define SIMP_EXPAND_NUM 200
#define SIMP_EXPAND_FAC 20.0L
#define SIMP_MAXCT_EXPAND 2000
#define SIMP_MAXCT_TOTAL 10000

// Hergetfit01: November 03, 2022:
// Use Hergetchi01 to perform orbit fitting using the Method of Herget,
// and a downhill simplex method applied to the 2-dimensional space of
// geodist1 and geodist2.
// The vector orbit holds a [0], e [1], mjd [2], and the state vectors [3-8] on
// return of Hergerchi01(). Hergetfit01 pushes back one additional
// datum: the number of orbit evaluations (~iterations) required
// to reach convergence [9].
long double Hergetfit01(long double geodist1, long double geodist2, long double simplex_scale, int simptype, long double ftol, int point1, int point2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid, vector <long double> &orbit, int verbose)
{
  int Hergetpoint1, Hergetpoint2;
  long double simprange;
  long double simplex[3][2];
  long double simpchi[3];
  long double refdist[2],trialdist[2];
  long double chisq, bestchi, worstchi, newchi;
  long double global_bestchi = LARGERR;
  long double global_bestd1 = geodist1;
  long double global_bestd2 = geodist2;
  int i,j,worstpoint, bestpoint;
  int unboundsimplex[3];
  int simp_eval_ct=0;
  int simp_total_ct=0;
  if(simplex_scale<=0.0L || simplex_scale>=SIMPLEX_SCALE_LIMIT) {
    cerr << "WARNING: simplex scale must be between 0 and " << SIMPLEX_SCALE_LIMIT << "\n";
    cerr << "Input out-of-range value " << simplex_scale << " will be reseset to ";
    simplex_scale = SIMPLEX_SCALEFAC;
    cerr << simplex_scale << "\n";
  }
  
  // Input points are indexed from 1; apply offset
  Hergetpoint1 = point1-1;
  Hergetpoint2 = point2-1;

  if(DEBUG_2PTBVP>0) {
    cout << "Herget points: " << Hergetpoint1 << " " << Hergetpoint2 << "\n";
    for(i=0;i<long(obsMJD.size());i++) {
      cout << "Input observerpos " << i << ": " << obsMJD[i] << " " << observerpos[i].x << " " << observerpos[i].y << " " << observerpos[i].z << "\n";
    }
  }

  // SETUP FOR DOWNHILL SIMPLEX SEARCH
  Herget_simplex_int(geodist1, geodist2, simplex_scale, simplex, simptype);
  
  // See if the simplex leads to hyperbolic orbits
  for(i=0;i<3;i++) {
    unboundsimplex[i] = Herget_unboundcheck01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec);
    if(unboundsimplex[i]!=0 && unboundsimplex[i]!=1) {
      // Input is fatally flawed, exit.
      cerr << "ERROR: fatally flawed input to downhill simplex, dists " << simplex[i][0] << " " << simplex[i][1] << "\n";
      cerr << "points " << Hergetpoint1 << " and " << Hergetpoint2 << " out of allowed range 0 to " << obsMJD.size() << "\n";
      for(i=0;i<10;i++) orbit.push_back(-1.0L);
      return(LARGERR);
    }
  }
  while(unboundsimplex[0]==1 && unboundsimplex[1]==1 && unboundsimplex[2]==1) {
    // All the points are bad, shrink all the distances.
    if(verbose>=2) cout << "All points are hyperbolic with simplex:\n";
    for(i=0;i<3;i++) {
      if(verbose>=2) cout << simplex[i][0] << " " << simplex[i][1] << "\n";
      simplex[i][0]*=HERGET_DOWNSCALE;
      simplex[i][1]*=HERGET_DOWNSCALE;
      unboundsimplex[i] = Herget_unboundcheck01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec);
    }
  }
  // If we get here, there must be at least one good point.
  bestpoint = worstpoint = -1;
  for(i=2;i>=0;i--) {
    if(unboundsimplex[i]==0) bestpoint=i;
    if(unboundsimplex[i]==1) worstpoint=i;
  }
  if(bestpoint<0) {
    cerr << "Logically impossible case involving hyperbolic simplex points\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0L);
    return(LARGERR);
  }
  if(worstpoint>=0) {
    if(verbose>=2) cout << "Good simplex point " << bestpoint << ": " << simplex[bestpoint][0] << " " << simplex[bestpoint][1] << "\n";
    // There is at least one bad point.
    for(i=0;i<3;i++) {
      while(unboundsimplex[i]==1 && sqrt(LDSQUARE(simplex[i][0]-simplex[bestpoint][0]) + LDSQUARE(simplex[i][1]-simplex[bestpoint][1])) > MINHERGETDIST) {
	// Bring the bad point closer to a good point until it stops being bad.
	if(verbose>=2) cout << "Modifying bad simplex point " << i << ": " << simplex[i][0] << " " << simplex[i][1] << "\n";
	simplex[i][0] = 0.5L*(simplex[i][0]+simplex[bestpoint][0]);
	simplex[i][1] = 0.5L*(simplex[i][1]+simplex[bestpoint][1]);
	unboundsimplex[i] = Herget_unboundcheck01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec);
      }
    }
  }
  if(unboundsimplex[0]!=0 || unboundsimplex[1]!=0 || unboundsimplex[2]!=0) {
    // We tried everything and still couldn't come up with a good simplex.
    cerr << "Attempt to get viable simplex failed:\n";
    for(i=0;i<3;i++) {
      cerr << simplex[i][0] << " " << simplex[i][1] << " unbound = " << unboundsimplex[i] << "\n";
    }
    cerr << "Aborting\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0L);
    return(LARGERR);
  } else {
    if(verbose>=2) cout << "Good input simplex:\n";
    for(i=0;i<3;i++) {
      if(verbose>=2) cout << simplex[i][0] << " " << simplex[i][1] << " unbound = " << unboundsimplex[i] << "\n";
    }
  }
  
  for(i=0;i<3;i++) simpchi[i]=LARGERR;
  // Normally, the while loop immediately below should execute only once,
  // but if the input distances don't lead to a bound orbit, it will
  // reduce them and try again.
  while((simpchi[0] == LARGERR || simpchi[1] == LARGERR || simpchi[2] == LARGERR) && simplex[0][0]>MINHERGETDIST) {
    // Calculate chi-square values for each point in the initial simplex
    // Note that the output vectors fitRA, fitDec, and resid are null-wiped
    // internally, so it isn't necessary to wipe them here.
    for(i=0;i<3;i++) {
      if(verbose>=2) cout << "Calling Hergetchi01 with distances " << simplex[i][0] << " " << simplex[i][1] << " : ";
      simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
      if(verbose>=2) cout << "reduced chi-square value is " << simpchi[i]/obsMJD.size() << "\n";
      simp_eval_ct++;
      simp_total_ct++;
    }
    if(simpchi[0] == LARGERR || simpchi[1] == LARGERR || simpchi[2] == LARGERR) {
      if(verbose>=2) cout << "Hergetchi01 returned failure code with simplex:\n";
      for(i=0;i<3;i++) {
	if(verbose>=2) cout << simplex[i][0] << " " << simplex[i][1] << "\n";
	simplex[i][0]*=HERGET_DOWNSCALE;
	simplex[i][1]*=HERGET_DOWNSCALE;
      }
    }
  }
  if(simplex[0][0]<=MINHERGETDIST) {
    if(verbose>=2) cerr << "ERROR: no acceptable solutions found for the Kepler two-point boundary value problem:\n";
    if(verbose>=2) cerr << "Method of Herget cannot proceed with these data\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0L);
    return(LARGERR);
  }
  if(verbose>=2) cout << "Reduced chi-square value for input distances is " << simpchi[0]/obsMJD.size() << "\n";
  
  // Find best and worst points
  worstpoint=bestpoint=0;
  bestchi = worstchi = simpchi[0];
  for(i=1;i<3;i++) {
    if(simpchi[i]<bestchi) {
      bestchi = simpchi[i];
      bestpoint=i;
    }
    if(simpchi[i]>worstchi) {
      worstchi = simpchi[i];
      worstpoint=i;
    }
  }
  simprange = (worstchi-bestchi)/bestchi;

  // LAUNCH DOWNHILL SIMPLEX SEARCH
  while(simprange>ftol && simp_total_ct <= SIMP_MAXCT_TOTAL) {
    if(verbose>=2) cout << fixed << setprecision(6) << "Eval " << simp_total_ct << ": Best reduced chi-square value is " << bestchi/obsMJD.size() << ", range is " << simprange << ", vector is " << simplex[bestpoint][0] << " "  << simplex[bestpoint][1] << "\n";
    
    // Try to reflect away from worst point
    // Find mean over all the points except the worst one
    refdist[0] = refdist[1] = 0.0L;
    for(i=0;i<3;i++) {
      if(i!=worstpoint) {
	refdist[0] += simplex[i][0]/2.0L;
	refdist[1] += simplex[i][1]/2.0L;
      }
    }
    // Calculate new trial point
    trialdist[0] = refdist[0] - (simplex[worstpoint][0] - refdist[0]);
    trialdist[1] = refdist[1] - (simplex[worstpoint][1] - refdist[1]);
    // Calculate chi-square value at this new point
    chisq = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
    simp_eval_ct++;
    simp_total_ct++;
    if(chisq<bestchi) {
      // Very good result. Let this point replace worstpoint in the simplex
      for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
      simpchi[worstpoint]=chisq;
     // Extrapolate further in this direction: maybe we can do even better
      trialdist[0] = refdist[0] - 2.0L*(simplex[worstpoint][0] - refdist[0]);
      trialdist[1] = refdist[1] - 2.0L*(simplex[worstpoint][1] - refdist[1]);
      newchi = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
      simp_eval_ct++;
      simp_total_ct++;
      if(newchi<chisq) {
	// Let this even better point replace worstpoint in the simplex
	for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
	simpchi[worstpoint]=newchi;
      }
      // This closes the case where reflecting away from the
      // worst point was a big success.
    } else {
      // Reflecting away from the worst point wasn't great, but
      // we'll see what we can manage.
      if(chisq<worstchi) { 
	// The new point was at least better than the previous worst.
	// Add it to the simplex in place of the worst point
	for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
	simpchi[worstpoint]=chisq;
      } else {
	// The new point was really no good.
	// This is the part of the story where we give up on
	// reflecting away from the worst point, and we try
	// something else.
	// First, try contracting away from the bad point,
	// instead of reflecting away from it.
	trialdist[0] = 0.5L*(simplex[worstpoint][0] + refdist[0]);
	trialdist[1] = 0.5L*(simplex[worstpoint][1] + refdist[1]);
	// Calculate chi-square value at this new point
	chisq = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
	simp_eval_ct++;
	simp_total_ct++;
	if(chisq<worstchi) {
	  // The new point is better than the previous worst point
	  // Add it to the simplex in place of the worst point
	  for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
	  simpchi[worstpoint]=chisq;
	} else {
	  // Even contracting away from the bad point didn't help.
	  // Only one thing left to try: contract toward the best point.
	  // This means each point will become an average of the best
	  // point and its former self.
	  for(i=0;i<3;i++) {
	    if(i!=bestpoint) {
	      simplex[i][0] = 0.5L*(simplex[i][0] + simplex[bestpoint][0]);
	      simplex[i][1] = 0.5L*(simplex[i][1] + simplex[bestpoint][1]);
	      simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
	      simp_eval_ct++;
	      simp_total_ct++;
	    }
	  }
	  // Close case where nothing worked but contracting around the best point.
	}
	// Close case where reflecting away from the worst point did not work. 
      }
      // Close case where reflecting away from the worst point was not a big success.
    }
    // Expand the simplex if we've been running for a long time
    if(simp_eval_ct>SIMP_EXPAND_NUM && simp_total_ct <= SIMP_MAXCT_EXPAND) {
      // Zero the counter
      simp_eval_ct=0;
      // Find center of the simplex
      refdist[0] = (simplex[0][0] + simplex[1][0] + simplex[2][0])/3.0L;
      refdist[1] = (simplex[0][1] + simplex[1][1] + simplex[2][1])/3.0L;
      // Expand the simplex
      for(i=0;i<3;i++) {
	simplex[i][0] = refdist[0] + (simplex[i][0]-refdist[0])*SIMP_EXPAND_FAC;
	simplex[i][1] = refdist[1] + (simplex[i][1]-refdist[1])*SIMP_EXPAND_FAC;
      }
      // Re-evaluate the chi-square values
      for(i=0;i<3;i++) {
	simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
      }
    }

    // Identify best and worst points for next iteration.
    worstpoint=bestpoint=0;
    bestchi = worstchi = simpchi[0];
    for(i=1;i<3;i++) {
      if(simpchi[i]<bestchi) {
	bestchi = simpchi[i];
	bestpoint=i;
	} 
      if(simpchi[i]>worstchi) {
	worstchi = simpchi[i];
	worstpoint=i;
      }
    }
    if(bestchi<global_bestchi) {
      global_bestchi = bestchi;
      global_bestd1 = simplex[bestpoint][0];
      global_bestd2 = simplex[bestpoint][1];
    }

    if(bestchi<LARGERR) simprange = (worstchi-bestchi)/bestchi;
    else {
      if(verbose>=2) cout << "WARNING: probing a simplex with no valid points!\n";
      // We have problems: expanding the simplex resulted in no
      // acceptable points at all.
      simprange = LARGERR;
      // Try to find something reasonable.
      while((simpchi[0] == LARGERR || simpchi[1] == LARGERR || simpchi[2] == LARGERR) && simplex[0][0]>MINHERGETDIST) {
	// Calculate chi-square values for each point in the initial simplex
	// Note that the output vectors fitRA, fitDec, and resid are null-wiped
	// internally, so it isn't necessary to wipe them here.
	for(i=0;i<3;i++) {
	  if(verbose>=2) cout << "Calling Hergetchi01 with distances " << simplex[i][0] << " " << simplex[i][1] << "\n";
	  simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
	  simp_eval_ct++;
	  simp_total_ct++;
	}
	if(simpchi[0] == LARGERR || simpchi[1] == LARGERR || simpchi[2] == LARGERR) {
	  if(verbose>=2) cout << "Hergetchi01 returned failure code with simplex:\n";
	  for(i=0;i<3;i++) {
	    if(verbose>=2) cout << simplex[i][0] << " " << simplex[i][1] << "\n";
	    simplex[i][0]*=HERGET_DOWNSCALE;
	    simplex[i][1]*=HERGET_DOWNSCALE;
	  }
	}
      }
      if(simplex[0][0]<=MINHERGETDIST) {
	if(verbose>=2) cerr << "ERROR: no acceptable solutions found for the Kepler two-point boundary value problem:\n";
	if(verbose>=2) cerr << "Method of Herget cannot proceed with these data\n";
	for(i=0;i<10;i++) orbit.push_back(-1.0L);
	return(LARGERR);
      } else {
	// We did eventually find an acceptable simplex
	// Find best and worst points
	worstpoint=bestpoint=0;
	bestchi = worstchi = simpchi[0];
	for(i=1;i<3;i++) {
	  if(simpchi[i]<bestchi) {
	    bestchi = simpchi[i];
	    bestpoint=i;
	  }
	  if(simpchi[i]>worstchi) {
	    worstchi = simpchi[i];
	    worstpoint=i;
	  }
	}
	simprange = (worstchi-bestchi)/bestchi;
	if(bestchi<global_bestchi) {
	  global_bestchi = bestchi;
	  global_bestd1 = simplex[bestpoint][0];
	  global_bestd2 = simplex[bestpoint][1];
	}
	// Close case where we eventually found a viable simplex
      }
      // Close case where we had an unviable simplex and had to try to fix it.
    }
    // Close main optimization loop.
  }
  
  cout << fixed << setprecision(6) << "Best reduced chi-square value was " << global_bestchi/obsMJD.size() << ", with geocentric distances " << global_bestd1 << " and " << global_bestd2 << "\n";
  
  // Perform fit with final best parameters
  chisq = Hergetchi01(global_bestd1, global_bestd2, Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
  orbit.push_back((long double)simp_total_ct);
  return(chisq);
}

// Hergetfit01: April 11, 2023:
// Like Hergetfit01, but avoids long doubles.
// Use Hergetchi01 to perform orbit fitting using the Method of Herget,
// and a downhill simplex method applied to the 2-dimensional space of
// geodist1 and geodist2.
// The vector orbit holds a [0], e [1], mjd [2], and the state vectors [3-8] on
// return of Hergerchi01(). Hergetfit01 pushes back one additional
// datum: the number of orbit evaluations (~iterations) required
// to reach convergence [9].
double Hergetfit01(double geodist1, double geodist2, double simplex_scale, int simptype, double ftol, int point1, int point2, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid, vector <double> &orbit, int verbose)
{
  int Hergetpoint1, Hergetpoint2;
  double simprange;
  double simplex[3][2];
  double simpchi[3];
  double refdist[2],trialdist[2];
  double chisq, bestchi, worstchi, newchi;
  double global_bestchi = LARGERR2;
  double global_bestd1 = geodist1;
  double global_bestd2 = geodist2;
  int i,j,worstpoint, bestpoint;
  int unboundsimplex[3];
  int simp_eval_ct=0;
  int simp_total_ct=0;
  if(simplex_scale<=0.0L || simplex_scale>=SIMPLEX_SCALE_LIMIT) {
    cerr << "WARNING: simplex scale must be between 0 and " << SIMPLEX_SCALE_LIMIT << "\n";
    cerr << "Input out-of-range value " << simplex_scale << " will be reseset to ";
    simplex_scale = SIMPLEX_SCALEFAC;
    cerr << simplex_scale << "\n";
  }
  
  // Input points are indexed from 1; apply offset
  Hergetpoint1 = point1-1;
  Hergetpoint2 = point2-1;

  if(verbose>=2) {
    cout << "Herget points: " << Hergetpoint1 << " " << Hergetpoint2 << "\n";
    for(i=0;i<long(obsMJD.size());i++) {
      cout << "Input observerpos " << i << ": " << obsMJD[i] << " " << observerpos[i].x << " " << observerpos[i].y << " " << observerpos[i].z << "\n";
    }
  }

  // SETUP FOR DOWNHILL SIMPLEX SEARCH
  Herget_simplex_int(geodist1, geodist2, simplex_scale, simplex, simptype);
  
  // See if the simplex leads to hyperbolic orbits
  for(i=0;i<3;i++) {
    unboundsimplex[i] = Herget_unboundcheck01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec);
    if(unboundsimplex[i]!=0 && unboundsimplex[i]!=1) {
      // Input is fatally flawed, exit.
      cerr << "ERROR: fatally flawed input to downhill simplex, dists " << simplex[i][0] << " " << simplex[i][1] << "\n";
      cerr << "points " << Hergetpoint1 << " and " << Hergetpoint2 << " out of allowed range 0 to " << obsMJD.size() << "\n";
      for(i=0;i<10;i++) orbit.push_back(-1.0l);
      return(LARGERR2);
    }
  }
  while(unboundsimplex[0]==1 && unboundsimplex[1]==1 && unboundsimplex[2]==1) {
    // All the points are bad, shrink all the distances.
    if(verbose>=2) cout << "All points are hyperbolic with simplex:\n";
    for(i=0;i<3;i++) {
      if(verbose>=2) cout << simplex[i][0] << " " << simplex[i][1] << "\n";
      simplex[i][0]*=HERGET_DOWNSCALE;
      simplex[i][1]*=HERGET_DOWNSCALE;
      unboundsimplex[i] = Herget_unboundcheck01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec);
    }
  }
  // If we get here, there must be at least one good point.
  bestpoint = worstpoint = -1;
  for(i=2;i>=0;i--) {
    if(unboundsimplex[i]==0) bestpoint=i;
    if(unboundsimplex[i]==1) worstpoint=i;
  }
  if(bestpoint<0) {
    cerr << "Logically impossible case involving hyperbolic simplex points\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0l);
    return(LARGERR2);
  }
  if(worstpoint>=0) {
    if(verbose>=2) cout << "Good simplex point " << bestpoint << ": " << simplex[bestpoint][0] << " " << simplex[bestpoint][1] << "\n";
    // There is at least one bad point.
    for(i=0;i<3;i++) {
      while(unboundsimplex[i]==1 && sqrt(LDSQUARE(simplex[i][0]-simplex[bestpoint][0]) + LDSQUARE(simplex[i][1]-simplex[bestpoint][1])) > MINHERGETDIST) {
	// Bring the bad point closer to a good point until it stops being bad.
	if(verbose>=2) cout << "Modifying bad simplex point " << i << ": " << simplex[i][0] << " " << simplex[i][1] << "\n";
	simplex[i][0] = 0.5L*(simplex[i][0]+simplex[bestpoint][0]);
	simplex[i][1] = 0.5L*(simplex[i][1]+simplex[bestpoint][1]);
	unboundsimplex[i] = Herget_unboundcheck01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec);
      }
    }
  }
  if(unboundsimplex[0]!=0 || unboundsimplex[1]!=0 || unboundsimplex[2]!=0) {
    // We tried everything and still couldn't come up with a good simplex.
    cerr << "Attempt to get viable simplex failed:\n";
    for(i=0;i<3;i++) {
      cerr << simplex[i][0] << " " << simplex[i][1] << " unbound = " << unboundsimplex[i] << "\n";
    }
    cerr << "Aborting\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0l);
    return(LARGERR2);
  } else {
    if(verbose>=2) cout << "Good input simplex:\n";
    for(i=0;i<3;i++) {
      if(verbose>=2) cout << simplex[i][0] << " " << simplex[i][1] << " unbound = " << unboundsimplex[i] << "\n";
    }
  }
  
  for(i=0;i<3;i++) simpchi[i]=LARGERR2;
  // Normally, the while loop immediately below should execute only once,
  // but if the input distances don't lead to a bound orbit, it will
  // reduce them and try again.
  while((simpchi[0] == LARGERR2 || simpchi[1] == LARGERR2 || simpchi[2] == LARGERR2) && simplex[0][0]>MINHERGETDIST) {
    // Calculate chi-square values for each point in the initial simplex
    // Note that the output vectors fitRA, fitDec, and resid are null-wiped
    // internally, so it isn't necessary to wipe them here.
    for(i=0;i<3;i++) {
      if(verbose>=2) cout << "Calling Hergetchi01 with distances " << simplex[i][0] << " " << simplex[i][1] << " : ";
      simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
      if(verbose>=2) cout << "reduced chi-square value is " << simpchi[i]/obsMJD.size() << "\n";
      simp_eval_ct++;
      simp_total_ct++;
    }
    if(simpchi[0] == LARGERR2 || simpchi[1] == LARGERR2 || simpchi[2] == LARGERR2) {
      if(verbose>=2) cout << "Hergetchi01 returned failure code with simplex:\n";
      for(i=0;i<3;i++) {
	if(verbose>=2) cout << simplex[i][0] << " " << simplex[i][1] << "\n";
	simplex[i][0]*=HERGET_DOWNSCALE;
	simplex[i][1]*=HERGET_DOWNSCALE;
      }
    }
  }
  if(simplex[0][0]<=MINHERGETDIST) {
    if(verbose>=2) cerr << "ERROR: no acceptable solutions found for the Kepler two-point boundary value problem:\n";
    if(verbose>=2) cerr << "Method of Herget cannot proceed with these data\n";
    for(i=0;i<10;i++) orbit.push_back(-1.0l);
    return(LARGERR2);
  }
  if(verbose>=2) cout << "Reduced chi-square value for input distances is " << simpchi[0]/obsMJD.size() << "\n";
  
  // Find best and worst points
  worstpoint=bestpoint=0;
  bestchi = worstchi = simpchi[0];
  for(i=1;i<3;i++) {
    if(simpchi[i]<bestchi) {
      bestchi = simpchi[i];
      bestpoint=i;
    }
    if(simpchi[i]>worstchi) {
      worstchi = simpchi[i];
      worstpoint=i;
    }
  }
  simprange = (worstchi-bestchi)/bestchi;

  // LAUNCH DOWNHILL SIMPLEX SEARCH
  while(simprange>ftol && simp_total_ct <= SIMP_MAXCT_TOTAL) {
    if(verbose>=2) cout << fixed << setprecision(6) << "Eval " << simp_total_ct << ": Best reduced chi-square value is " << bestchi/obsMJD.size() << ", range is " << simprange << ", vector is " << simplex[bestpoint][0] << " "  << simplex[bestpoint][1] << "\n";
    
    // Try to reflect away from worst point
    // Find mean over all the points except the worst one
    refdist[0] = refdist[1] = 0.0L;
    for(i=0;i<3;i++) {
      if(i!=worstpoint) {
	refdist[0] += simplex[i][0]/2.0L;
	refdist[1] += simplex[i][1]/2.0L;
      }
    }
    // Calculate new trial point
    trialdist[0] = refdist[0] - (simplex[worstpoint][0] - refdist[0]);
    trialdist[1] = refdist[1] - (simplex[worstpoint][1] - refdist[1]);
    // Calculate chi-square value at this new point
    chisq = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
    simp_eval_ct++;
    simp_total_ct++;
    if(chisq<bestchi) {
      // Very good result. Let this point replace worstpoint in the simplex
      for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
      simpchi[worstpoint]=chisq;
     // Extrapolate further in this direction: maybe we can do even better
      trialdist[0] = refdist[0] - 2.0L*(simplex[worstpoint][0] - refdist[0]);
      trialdist[1] = refdist[1] - 2.0L*(simplex[worstpoint][1] - refdist[1]);
      newchi = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
      simp_eval_ct++;
      simp_total_ct++;
      if(newchi<chisq) {
	// Let this even better point replace worstpoint in the simplex
	for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
	simpchi[worstpoint]=newchi;
      }
      // This closes the case where reflecting away from the
      // worst point was a big success.
    } else {
      // Reflecting away from the worst point wasn't great, but
      // we'll see what we can manage.
      if(chisq<worstchi) { 
	// The new point was at least better than the previous worst.
	// Add it to the simplex in place of the worst point
	for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
	simpchi[worstpoint]=chisq;
      } else {
	// The new point was really no good.
	// This is the part of the story where we give up on
	// reflecting away from the worst point, and we try
	// something else.
	// First, try contracting away from the bad point,
	// instead of reflecting away from it.
	trialdist[0] = 0.5L*(simplex[worstpoint][0] + refdist[0]);
	trialdist[1] = 0.5L*(simplex[worstpoint][1] + refdist[1]);
	// Calculate chi-square value at this new point
	chisq = Hergetchi01(trialdist[0], trialdist[1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
	simp_eval_ct++;
	simp_total_ct++;
	if(chisq<worstchi) {
	  // The new point is better than the previous worst point
	  // Add it to the simplex in place of the worst point
	  for(j=0;j<2;j++) simplex[worstpoint][j] = trialdist[j];
	  simpchi[worstpoint]=chisq;
	} else {
	  // Even contracting away from the bad point didn't help.
	  // Only one thing left to try: contract toward the best point.
	  // This means each point will become an average of the best
	  // point and its former self.
	  for(i=0;i<3;i++) {
	    if(i!=bestpoint) {
	      simplex[i][0] = 0.5L*(simplex[i][0] + simplex[bestpoint][0]);
	      simplex[i][1] = 0.5L*(simplex[i][1] + simplex[bestpoint][1]);
	      simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
	      simp_eval_ct++;
	      simp_total_ct++;
	    }
	  }
	  // Close case where nothing worked but contracting around the best point.
	}
	// Close case where reflecting away from the worst point did not work. 
      }
      // Close case where reflecting away from the worst point was not a big success.
    }
    // Expand the simplex if we've been running for a long time
    if(simp_eval_ct>SIMP_EXPAND_NUM && simp_total_ct <= SIMP_MAXCT_EXPAND) {
      // Zero the counter
      simp_eval_ct=0;
      // Find center of the simplex
      refdist[0] = (simplex[0][0] + simplex[1][0] + simplex[2][0])/3.0L;
      refdist[1] = (simplex[0][1] + simplex[1][1] + simplex[2][1])/3.0L;
      // Expand the simplex
      for(i=0;i<3;i++) {
	simplex[i][0] = refdist[0] + (simplex[i][0]-refdist[0])*SIMP_EXPAND_FAC;
	simplex[i][1] = refdist[1] + (simplex[i][1]-refdist[1])*SIMP_EXPAND_FAC;
      }
      // Re-evaluate the chi-square values
      for(i=0;i<3;i++) {
	simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
      }
    }

    // Identify best and worst points for next iteration.
    worstpoint=bestpoint=0;
    bestchi = worstchi = simpchi[0];
    for(i=1;i<3;i++) {
      if(simpchi[i]<bestchi) {
	bestchi = simpchi[i];
	bestpoint=i;
	} 
      if(simpchi[i]>worstchi) {
	worstchi = simpchi[i];
	worstpoint=i;
      }
    }
    if(bestchi<global_bestchi) {
      global_bestchi = bestchi;
      global_bestd1 = simplex[bestpoint][0];
      global_bestd2 = simplex[bestpoint][1];
    }

    if(bestchi<LARGERR2) simprange = (worstchi-bestchi)/bestchi;
    else {
      if(verbose>=2) cout << "WARNING: probing a simplex with no valid points!\n";
      // We have problems: expanding the simplex resulted in no
      // acceptable points at all.
      simprange = LARGERR2;
      // Try to find something reasonable.
      while((simpchi[0] == LARGERR2 || simpchi[1] == LARGERR2 || simpchi[2] == LARGERR2) && simplex[0][0]>MINHERGETDIST) {
	// Calculate chi-square values for each point in the initial simplex
	// Note that the output vectors fitRA, fitDec, and resid are null-wiped
	// internally, so it isn't necessary to wipe them here.
	for(i=0;i<3;i++) {
	  if(verbose>=2) cout << "Calling Hergetchi01 with distances " << simplex[i][0] << " " << simplex[i][1] << "\n";
	  simpchi[i] = Hergetchi01(simplex[i][0], simplex[i][1], Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
	  simp_eval_ct++;
	  simp_total_ct++;
	}
	if(simpchi[0] == LARGERR2 || simpchi[1] == LARGERR2 || simpchi[2] == LARGERR2) {
	  if(verbose>=2) cout << "Hergetchi01 returned failure code with simplex:\n";
	  for(i=0;i<3;i++) {
	    if(verbose>=2) cout << simplex[i][0] << " " << simplex[i][1] << "\n";
	    simplex[i][0]*=HERGET_DOWNSCALE;
	    simplex[i][1]*=HERGET_DOWNSCALE;
	  }
	}
      }
      if(simplex[0][0]<=MINHERGETDIST) {
	if(verbose>=2) cerr << "ERROR: no acceptable solutions found for the Kepler two-point boundary value problem:\n";
	if(verbose>=2) cerr << "Method of Herget cannot proceed with these data\n";
	for(i=0;i<10;i++) orbit.push_back(-1.0l);
	return(LARGERR2);
      } else {
	// We did eventually find an acceptable simplex
	// Find best and worst points
	worstpoint=bestpoint=0;
	bestchi = worstchi = simpchi[0];
	for(i=1;i<3;i++) {
	  if(simpchi[i]<bestchi) {
	    bestchi = simpchi[i];
	    bestpoint=i;
	  }
	  if(simpchi[i]>worstchi) {
	    worstchi = simpchi[i];
	    worstpoint=i;
	  }
	}
	simprange = (worstchi-bestchi)/bestchi;
	if(bestchi<global_bestchi) {
	  global_bestchi = bestchi;
	  global_bestd1 = simplex[bestpoint][0];
	  global_bestd2 = simplex[bestpoint][1];
	}
	// Close case where we eventually found a viable simplex
      }
      // Close case where we had an unviable simplex and had to try to fix it.
    }
    // Close main optimization loop.
  }
  
  
  // Perform fit with final best parameters
  chisq = Hergetchi01(global_bestd1, global_bestd2, Hergetpoint1, Hergetpoint2, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, verbose);
  
  cout << fixed << setprecision(6) << "Best chi/N " << global_bestchi/obsMJD.size() << ", geodists " << global_bestd1 << " and " << global_bestd2 << " orbit a, e = " << orbit[0]/AU_KM << ", " << orbit[1] << "\n";

  orbit.push_back(double(simp_total_ct));
  return(chisq);
}


#undef DEBUG_2PTBVP
#undef SIMP_EXPAND_NUM
#undef SIMP_EXPAND_FAC
#undef SIMP_MAXCT_EXPAND
#undef SIMP_MAXCT_TOTAL

// medind_3d_index: January 31, 2023: Calculate the median of
// a vector of points of type point3d_index in x (dim=1), y (dim=2), or z (dim=3)
long medind_3d_index(const vector <point3d_index> &pointvec, int dim)
{
  vector <point3d_index> pvec = pointvec; //Mutable copy of immutable input vector
  for(long i=0; i<long(pvec.size()); i++) pvec[i].index=i; //Redefine indices
  long medpt = pvec.size()/2; // Central point of vector (it will be off by one half
                              // for a vector with even length, but we don't care).
  if(dim%3 == 1) sort(pvec.begin(), pvec.end(), lower_point3d_index_x()); // Sort vector by x
  else if(dim%3 == 2) sort(pvec.begin(), pvec.end(), lower_point3d_index_y()); // Sort vector by y
  else if(dim%3 == 0) sort(pvec.begin(), pvec.end(), lower_point3d_index_z()); // Sort vector by z
  else {
    cerr << "ERROR: medind_3d_index recieved invalid dimension " << dim << "\n";
    return(-1);
  }
  return(pvec[medpt].index); // Output the index of the median point in
                             // the original, unsorted input vector.
}

// split3d_index: January 31, 2023:
// Given a vector of type point3d_index, split it into two halves,
// a left half with all the points lower than or equal to a specified
// split point along the chosen dimension (use dim = 1, 2, or 3 to
// split along x, y, z, respectively).
int split3d_index(const vector <point3d_index> &pointvec, int dim, long splitpoint, vector <point3d_index> &left, vector <point3d_index> &right)
{
  long i=0;
  long double splitval = 0.0L;

  if(dim%3==1) {
    // split on x
    splitval = pointvec[splitpoint].x;
    for(i=0 ; i<long(pointvec.size()); i++) {
      if(i!=splitpoint && pointvec[i].x<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%3==2) {
    // split on y
    splitval = pointvec[splitpoint].y;
    for(i=0 ; i<long(pointvec.size()); i++) {
      if(i!=splitpoint && pointvec[i].y<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%3==0) {
    // split on z
    splitval = pointvec[splitpoint].z;
    for(i=0 ; i<long(pointvec.size()); i++) {
      if(i!=splitpoint && pointvec[i].z<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else {
      cerr << "ERROR: split3d_index asked to split on undefined dimension " << dim << "\n";
      return(1);
  } 
  return(0);
}

// kdtree_3d_index: January 31, 2023
// Given an input root point, presumed to have null
// right and left branches, load the branches and then
// call kdtree_3d_index on them recursively.
int kdtree_3d_index(const vector <point3d_index> &invec, int dim, long splitpoint, long kdroot, vector <KD_point3d_index> &kdvec)
{
  int lmed=0;
  int rmed=0;
  int kdct = kdvec.size()-1;
  long leftrootkd=-1;
  long rightrootkd=-1;
  point3d_index point0 = point3d_index(0.0l,0.0l,0.0l,0);
  KD_point3d_index lp = KD_point3d_index(point0,-1,-1,0,0);
  KD_point3d_index rp = KD_point3d_index(point0,-1,-1,0,0);
  vector <point3d_index> leftvec = {};
  vector <point3d_index> rightvec = {};

  // Basic outline: split the input vector into a left and a right
  // half, where the left half is below (or level with) splitpoint
  // in the dimension specified by dim, and the right half is above
  // splitpoint. Find the median of the left and right vectors,
  // and make the left and right branches from kdroot point to
  // these medians. Then call kdtree_6D01 itself recursively on
  // each of these median points, to peform a new split along a
  // different dimension.
  split3d_index(invec,dim,splitpoint,leftvec,rightvec);

  dim+=1;
  while(dim>3) dim-=3;

  if(leftvec.size()==1) {
    // Left branch is just a single leaf
    lp = KD_point3d_index(leftvec[0],-1,-1,dim,0); // Define new point as a leaf: branches point nowhere
    kdvec.push_back(lp); // Add this new point to the KD tree.
    kdct++; // Keep track of how many point are in the tree
    kdvec[kdroot].left = kdct; // Stick the new point on the left branch of the input root.
  } else if(leftvec.size()<=0) {
    // There is no left branch
    kdvec[kdroot].left = -1;
  }
  if(rightvec.size()==1) {
    // Right branch is just a single leaf
    rp = KD_point3d_index(rightvec[0],-1,-1,dim,0);
    kdvec.push_back(rp);
    kdct++;
    kdvec[kdroot].right = kdct;
  } else if(rightvec.size()<=0) {
    // There is no right branch
    kdvec[kdroot].right = -1;
  }
   
 if(leftvec.size()>1) {
    lmed = medind_3d_index(leftvec,dim);
    lp = KD_point3d_index(leftvec[lmed],-1,-1,dim,0);
    kdvec.push_back(lp);
    kdct++;
    kdvec[kdroot].left = kdct;
    leftrootkd = kdct;
 }
 
  if(rightvec.size()>1) {
    rmed = medind_3d_index(rightvec,dim);
    rp = KD_point3d_index(rightvec[rmed],-1,-1,dim,0);
    kdvec.push_back(rp);
    kdct++;
    kdvec[kdroot].right = kdct;
    rightrootkd = kdct;
  }
  // I moved these down out of the above loops, because I thought
  // that otherwise, a bunch of stuff might get pushed down by the
  // left loop that the right loop didn't know about.
  if(leftvec.size()>1 && leftrootkd>=0) kdtree_3d_index(leftvec,dim,lmed,leftrootkd,kdvec);
  else if(leftvec.size()>1 && leftrootkd<0)
    {
      cerr << "Error, kdtree_3d_index finds leftroot less than zero with leftvec.size() = " << leftvec.size() << "\n";
    }
  if(rightvec.size()>1 && rightrootkd>=0) kdtree_3d_index(rightvec,dim,rmed,rightrootkd,kdvec);
  else if(rightvec.size()>1 && rightrootkd<0)
    {
      cerr << "Error, kdtree_3d_index finds rightroot less than zero with rightvec.size() = " << rightvec.size() << "\n";
    }

  return(0);
}

// point3d_index_dist2: January 31, 2023:
// Calculate the squared distance in 3-dimensional parameter space
// between two points of class point3d_index.
double point3d_index_dist2(const point3d_index &p1, const point3d_index &p2)
{
  return(DSQUARE(p1.x - p2.x) + DSQUARE(p1.y - p2.y) + DSQUARE(p1.z - p2.z));
}	 

// kdrange_3d_index: January 31, 2023:
// Given a k-d tree vector kdvec created by kdtree_3d_index,
// perform a range-query about the specified point. Returns
// a vector indexing all of the points in the input k-d tree
// that lie within the specified range of the input coordinates.
// Assumes that kdvec[0] is the root of the k-d tree.
int kdrange_3d_index(const vector <KD_point3d_index> &kdvec, const point3d_index &querypoint, double range, vector <long> &indexvec)
{
  double rng2 = range*range;
  int notdone=1;
  int dim=1;
  int currentpoint=0;
  int leftpoint=0;
  int rightpoint=0;
  int goleft=0;
  int goright=0;
  double pointdiff = 0.0l;
  double pdist2 = 0.0l;
  vector <long> checkit={};
  int checknum=0;

  while(notdone>0) {
    // Climb to the top of the k-d tree, keeping track
    // of potentially interesting unexplored branches
    // in the vector checkit.
    while(leftpoint>=0 || rightpoint>=0) {
      // Previous step did not end on a leaf.
      leftpoint = kdvec[currentpoint].left;
      rightpoint = kdvec[currentpoint].right;
      dim = kdvec[currentpoint].dim;
      if(dim%3==1) pointdiff = kdvec[currentpoint].point.x - querypoint.x;
      else if(dim%3==2) pointdiff = kdvec[currentpoint].point.y - querypoint.y;
      else if(dim%3==0) pointdiff = kdvec[currentpoint].point.z - querypoint.z;

      goright = (pointdiff <= range); // possible hits lie to the left;
      goleft = (pointdiff >= -range); // possible hits lie to the right;
      if(goleft && goright) {
	// Current point might be within range.
	pdist2 = point3d_index_dist2(querypoint,kdvec[currentpoint].point);
	if(pdist2 <= rng2) {
	  // Current point is within range. Add it to the output vector
	  indexvec.push_back(currentpoint);
	}
	if(leftpoint>=0) {
	  //Explore leftward first.
	  currentpoint = leftpoint;
	  if(rightpoint>=0) {
	    // Rightward branch will also be explored later
	    checknum++;
	    if(checknum>long(checkit.size())) {
	      checkit.push_back(rightpoint);
	    }
	    else {
	      checkit[checknum-1] = rightpoint;
	    }
	  }
	}
	else if(rightpoint>=0) {
	  // Leftward branch is a dead end: explore rightward branch
	  currentpoint = rightpoint;
	}
      }
      else if(goleft) {
	// Current point cannot be in range, but points that
	// are in range may lie along the left branch.
	if(leftpoint>=0) {
	  currentpoint = leftpoint;
	} else rightpoint=-1; // Dead end, make sure while-loop exits.
      } else if(goright) {
	// Current point cannot be in range, but points that
	// are in range may lie along the right branch.
	if(rightpoint>=0) {
	  currentpoint = rightpoint;
	} else leftpoint=-1;  // Dead end, make sure while-loop exits.
      } else {
	// Program concluded it should go neither left nor right.
	// The likely cause is that it encountered a NAN. Give up on this point.
	leftpoint=rightpoint=-1;
	cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
	cerr << "Query point:\n";
	cerr << querypoint.x << ", " << querypoint.y << ", " << querypoint.z << "\n";
	cerr << "Target point:\n";
 	cerr << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ", " << kdvec[currentpoint].point.z << "\n";
     }
      // Close while-loop checking if we've hit a leaf.
    }
    // We have climbed up the tree to a leaf. Go backwards through
    // the checkit vector and see if there is anything to check.
    checknum=checkit.size();
    while(checknum>=1 && checkit[checknum-1]<0) checknum--;
    if(checknum<=0) {
      //There were no valid entries to check: we're done.
      notdone=0;
    } else {
      //Set currentpoint to the last valid entry in checkit
      currentpoint = checkit[checknum-1];
      //Mark this point as used.
      checkit[checknum-1]=-1;
      leftpoint=rightpoint=0;
    }
  }
  return(0);
}


// medind_4d_index: January 31, 2023: Calculate the median of
// a vector of points of type point4d_index in t (dim=1), x (dim=2), y (dim=3), or z (dim=3)
long medind_4d_index(const vector <point4d_index> &pointvec, int dim)
{
  vector <point4d_index> pvec = pointvec; //Mutable copy of immutable input vector
  for(long i=0; i<long(pvec.size()); i++) pvec[i].index=i; //Redefine indices
  long medpt = pvec.size()/2; // Central point of vector (it will be off by one half
                              // for a vector with even length, but we don't care).
  if(dim%4 == 1) sort(pvec.begin(), pvec.end(), lower_point4d_index_t()); // Sort vector by t
  else if(dim%4 == 2) sort(pvec.begin(), pvec.end(), lower_point4d_index_x()); // Sort vector by x
  else if(dim%4 == 3) sort(pvec.begin(), pvec.end(), lower_point4d_index_y()); // Sort vector by y
  else if(dim%4 == 0) sort(pvec.begin(), pvec.end(), lower_point4d_index_z()); // Sort vector by z
  else {
    cerr << "ERROR: medind_4d_index recieved invalid dimension " << dim << "\n";
    return(-1);
  }
  return(pvec[medpt].index); // Output the index of the median point in
                             // the original, unsorted input vector.
}

// split4d_index: January 31, 2023:
// Given a vector of type point4d_index, split it into two halves,
// a left half with all the points lower than or equal to a specified
// split point along the chosen dimension (use dim = 1, 2, 3, or 4 to
// split along t, x, y, or z, respectively).
int split4d_index(const vector <point4d_index> &pointvec, int dim, long splitpoint, vector <point4d_index> &left, vector <point4d_index> &right)
{
  long i=0;
  long double splitval = 0.0L;

  if(dim%4==1) {
    // split on t
    splitval = pointvec[splitpoint].t;
    for(i=0 ; i<long(pointvec.size()); i++) {
      if(i!=splitpoint && pointvec[i].t<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%4==2) {
    // split on x
    splitval = pointvec[splitpoint].x;
    for(i=0 ; i<long(pointvec.size()); i++) {
      if(i!=splitpoint && pointvec[i].x<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%4==3) {
    // split on y
    splitval = pointvec[splitpoint].y;
    for(i=0 ; i<long(pointvec.size()); i++) {
      if(i!=splitpoint && pointvec[i].y<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else if(dim%4==0) {
    // split on z
    splitval = pointvec[splitpoint].z;
    for(i=0 ; i<long(pointvec.size()); i++) {
      if(i!=splitpoint && pointvec[i].z<=splitval) {
	left.push_back(pointvec[i]);
      } else if (i!=splitpoint) {
	right.push_back(pointvec[i]);
      }
    }
  } else {
      cerr << "ERROR: split4d_index asked to split on undefined dimension " << dim << "\n";
      return(1);
  } 
  return(0);
}

// kdtree_4d_index: January 31, 2023
// Given an input root point, presumed to have null
// right and left branches, load the branches and then
// call kdtree_4d_index on them recursively.
int kdtree_4d_index(const vector <point4d_index> &invec, int dim, long splitpoint, long kdroot, vector <KD_point4d_index> &kdvec)
{
  int lmed=0;
  int rmed=0;
  long kdct = kdvec.size()-1;
  long leftrootkd=-1;
  long rightrootkd=-1;
  point4d_index point0 = point4d_index(0.0l,0.0l,0.0l,0.0l,0);
  KD_point4d_index lp = KD_point4d_index(point0,-1,-1,0,0);
  KD_point4d_index rp = KD_point4d_index(point0,-1,-1,0,0);
  vector <point4d_index> leftvec = {};
  vector <point4d_index> rightvec = {};

  // Basic outline: split the input vector into a left and a right
  // half, where the left half is below (or level with) splitpoint
  // in the dimension specified by dim, and the right half is above
  // splitpoint. Find the median of the left and right vectors,
  // and make the left and right branches from kdroot point to
  // these medians. Then call kdtree_6D01 itself recursively on
  // each of these median points, to peform a new split along a
  // different dimension.
  split4d_index(invec,dim,splitpoint,leftvec,rightvec);

  dim+=1;
  while(dim>4) dim-=4;

  if(leftvec.size()==1) {
    // Left branch is just a single leaf
    lp = KD_point4d_index(leftvec[0],-1,-1,dim,0); // Define new point as a leaf: branches point nowhere
    kdvec.push_back(lp); // Add this new point to the KD tree.
    kdct++; // Keep track of how many point are in the tree
    kdvec[kdroot].left = kdct; // Stick the new point on the left branch of the input root.
  } else if(leftvec.size()<=0) {
    // There is no left branch
    kdvec[kdroot].left = -1;
  }
  if(rightvec.size()==1) {
    // Right branch is just a single leaf
    rp = KD_point4d_index(rightvec[0],-1,-1,dim,0);
    kdvec.push_back(rp);
    kdct++;
    kdvec[kdroot].right = kdct;
  } else if(rightvec.size()<=0) {
    // There is no right branch
    kdvec[kdroot].right = -1;
  }
   
 if(leftvec.size()>1) {
    lmed = medind_4d_index(leftvec,dim);
    lp = KD_point4d_index(leftvec[lmed],-1,-1,dim,0);
    kdvec.push_back(lp);
    kdct++;
    kdvec[kdroot].left = kdct;
    leftrootkd = kdct;
 }
 
  if(rightvec.size()>1) {
    rmed = medind_4d_index(rightvec,dim);
    rp = KD_point4d_index(rightvec[rmed],-1,-1,dim,0);
    kdvec.push_back(rp);
    kdct++;
    kdvec[kdroot].right = kdct;
    rightrootkd = kdct;
  }
  // I moved these down out of the above loops, because I thought
  // that otherwise, a bunch of stuff might get pushed down by the
  // left loop that the right loop didn't know about.
  if(leftvec.size()>1 && leftrootkd>=0) kdtree_4d_index(leftvec,dim,lmed,leftrootkd,kdvec);
  else if(leftvec.size()>1 && leftrootkd<0)
    {
      cerr << "Error, kdtree_4d_index finds leftroot less than zero with leftvec.size() = " << leftvec.size() << "\n";
    }
  if(rightvec.size()>1 && rightrootkd>=0) kdtree_4d_index(rightvec,dim,rmed,rightrootkd,kdvec);
  else if(rightvec.size()>1 && rightrootkd<0)
    {
      cerr << "Error, kdtree_4d_index finds rightroot less than zero with rightvec.size() = " << rightvec.size() << "\n";
    }

  return(0);
}

// point4d_index_dist2: January 31, 2023:
// Calculate the squared distance in 4-dimensional parameter space
// between two points of class point4d_index.
double point4d_index_dist2(const point4d_index &p1, const point4d_index &p2)
{
  return(DSQUARE(p1.t - p2.t) + DSQUARE(p1.x - p2.x) + DSQUARE(p1.y - p2.y) + DSQUARE(p1.z - p2.z));
}	 

// kdrange_4d_index: January 31, 2023:
// Given a k-d tree vector kdvec created by kdtree_4d_index,
// perform a range-query about the specified point. Returns
// a vector indexing all of the points in the input k-d tree
// that lie within the specified range of the input coordinates.
// Assumes that kdvec[0] is the root of the k-d tree.
int kdrange_4d_index(const vector <KD_point4d_index> &kdvec, const point4d_index &querypoint, double range, vector <long> &indexvec)
{
  long double rng2 = range*range;
  int notdone=1;
  int dim=1;
  int currentpoint=0;
  int leftpoint=0;
  int rightpoint=0;
  int goleft=0;
  int goright=0;
  long double pointdiff = 0.0L;
  long double pdist2 = 0.0L;
  vector <long> checkit={};
  int checknum=0;

  while(notdone>0) {
    // Climb to the top of the k-d tree, keeping track
    // of potentially interesting unexplored branches
    // in the vector checkit.
    while(leftpoint>=0 || rightpoint>=0) {
      // Previous step did not end on a leaf.
      leftpoint = kdvec[currentpoint].left;
      rightpoint = kdvec[currentpoint].right;
      dim = kdvec[currentpoint].dim;
      if(dim%4==1) pointdiff = kdvec[currentpoint].point.t - querypoint.t;
      else if(dim%4==2) pointdiff = kdvec[currentpoint].point.x - querypoint.x;
      else if(dim%4==3) pointdiff = kdvec[currentpoint].point.y - querypoint.y;
      else if(dim%4==0) pointdiff = kdvec[currentpoint].point.z - querypoint.z;

      goright = (pointdiff <= range); // possible hits lie to the left;
      goleft = (pointdiff >= -range); // possible hits lie to the right;
      if(goleft && goright) {
	// Current point might be within range.
	pdist2 = point4d_index_dist2(querypoint,kdvec[currentpoint].point);
	if(pdist2 <= rng2) {
	  // Current point is within range. Add it to the output vector
	  indexvec.push_back(currentpoint);
	}
	if(leftpoint>=0) {
	  //Explore leftward first.
	  currentpoint = leftpoint;
	  if(rightpoint>=0) {
	    // Rightward branch will also be explored later
	    checknum++;
	    if(checknum>long(checkit.size())) {
	      checkit.push_back(rightpoint);
	    }
	    else {
	      checkit[checknum-1] = rightpoint;
	    }
	  }
	}
	else if(rightpoint>=0) {
	  // Leftward branch is a dead end: explore rightward branch
	  currentpoint = rightpoint;
	}
      }
      else if(goleft) {
	// Current point cannot be in range, but points that
	// are in range may lie along the left branch.
	if(leftpoint>=0) {
	  currentpoint = leftpoint;
	} else rightpoint=-1; // Dead end, make sure while-loop exits.
      } else if(goright) {
	// Current point cannot be in range, but points that
	// are in range may lie along the right branch.
	if(rightpoint>=0) {
	  currentpoint = rightpoint;
	} else leftpoint=-1;  // Dead end, make sure while-loop exits.
      } else {
	// Program concluded it should go neither left nor right.
	// The likely cause is that it encountered a NAN. Give up on this point.
	leftpoint=rightpoint=-1;
	cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
	cerr << "Query point:\n";
	cerr << querypoint.t << ", " << querypoint.x << ", " << querypoint.y << ", " << querypoint.z << "\n";
	cerr << "Target point:\n";
 	cerr << kdvec[currentpoint].point.t << ", " << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ", " << kdvec[currentpoint].point.z << "\n";
     }
      // Close while-loop checking if we've hit a leaf.
    }
    // We have climbed up the tree to a leaf. Go backwards through
    // the checkit vector and see if there is anything to check.
    checknum=checkit.size();
    while(checknum>=1 && checkit[checknum-1]<0) checknum--;
    if(checknum<=0) {
      //There were no valid entries to check: we're done.
      notdone=0;
    } else {
      //Set currentpoint to the last valid entry in checkit
      currentpoint = checkit[checknum-1];
      //Mark this point as used.
      checkit[checknum-1]=-1;
      leftpoint=rightpoint=0;
    }
  }
  return(0);
}

// MPCcal2MJD: Febuary 01, 2022
// Given a calendar date in the MPC format with integer year,
// integer month, and decimal day, calculate the Modified Julian Day
// (MJD). Works on any date after Jan 01, 1900.
double MPCcal2MJD(int year, int month, double day)
{
  int daystojan;
  double totaldays = 15020.0l; // MJD on UT 1900 January 1.0
  int i=1900;
  int isleap=0;
  vector <int> days_per_month;
  
  if(year<1900 || month<1 || month>12 || day<1.0 || day>32.0) {
    cerr << "ERROR: MPCcal2MJD has bad date: " << year << " " << month << " " << day << "\n";
    return(-1.0l);
  }

  // Load days_per_month
  days_per_month={};
  days_per_month.push_back(0); // Null, to get 1-index
  days_per_month.push_back(31); // January
  days_per_month.push_back(28); // February
  days_per_month.push_back(31); // March
  days_per_month.push_back(30); // April
  days_per_month.push_back(31); // May
  days_per_month.push_back(30); // June
  days_per_month.push_back(31); // July
  days_per_month.push_back(31); // August
  days_per_month.push_back(30); // September
  days_per_month.push_back(31); // October
  days_per_month.push_back(30); // November
  days_per_month.push_back(31); // December
  
  // Calculate the number of days up to 00:00 UT,
  // Jan 1, of the given year.
  daystojan=0;
  for(i=1900;i<year;i++) {
    if((i%4 == 0 && i%100 != 0) || i%400 == 0) {
      isleap=1;
    } else isleap=0;
    daystojan += 365 + isleap;
  }
  
  // Find out if the current year is a leapyear
  i=year;
  if((i%4 == 0 && i%100 != 0) || i%400 == 0) {
    isleap=1;
  } else isleap=0;
  
  totaldays += double(daystojan);
  for(i=1; i<month; i++) { // Remember, days_per_month vector is 1-indexed!
    totaldays += double(days_per_month[i]);
  }
  if(month>2 && isleap==1) totaldays += 1.0l; // Add leap day of current year.

  totaldays += day-1.0l; // Subtract 1 because there is no 0th day of the month.
  return(totaldays);
}

// mpc80_parseline: February 01, 2023:
// Read one line from an MPC 80-column formatted file,
// and load the date (MJD), RA (decimal degrees), Dec (decimal degrees),
// and magnitude as double-precision variables, and the object ID,
// band, and observatory codes as strings.
int mpc80_parseline(const string &lnfromfile, string &object, double *MJD, double *RA, double *Dec, double *mag, string &band, string &obscode)
{
  int i;
  int year=1900;
  int month=1;
  double day=1.0l;
  double hour,min,sec,deg;
  string sdat;
  char c='0';
  char decsign='+';
  double raread;
  double decread;
  
  // Read object name
  object = {};
  for(i=0;i<12;i++) {
    c = lnfromfile[i];
    if(c!=' ' && c!='*') object.push_back(c);
  }
  // Read the year
  sdat = {};
  for(i=15;i<19;i++) {
    sdat.push_back(lnfromfile[i]);
  }
  year = stoi(sdat);
  // Read the month
  sdat = {};
  for(i=20;i<22;i++) {
    sdat.push_back(lnfromfile[i]);
  }
  month = stoi(sdat);
  // Read the day
  sdat = {};
  for(i=23;i<32;i++) {
    c=lnfromfile[i];
    if(c!=' ') sdat.push_back(c);
  }
  day = stod(sdat);
  if(year<1900 || month<1 || month>12 || day < 1.0l || day > 32.0l) {
    cerr << "mpc80_readline cannot read a valid date from the line:\n";
    cerr << lnfromfile;
    return(1);
  }
  // Successfully read the date: convert it to MJD.
  *MJD = MPCcal2MJD(year,month,day);
  // Read the Right Ascension
  // hours
  sdat = {};
  for(i=32;i<34;i++) {
    sdat.push_back(lnfromfile[i]);
  }
  hour = stod(sdat);
  // minutes
  sdat = {};
  for(i=35;i<37;i++) {
    sdat.push_back(lnfromfile[i]);
  }
  if(sdat.size()>0) min = stod(sdat);
  else min = 0.0l;
  // seconds
  sdat = {};
  for(i=38;i<44;i++) {
    c = lnfromfile[i];
    if(c != ' ') sdat.push_back(c);
  }
  if(sdat.size()>0) sec = stod(sdat);
  else sec = 0.0l;
  // Convert the Right Ascension to decimal degrees
  raread = 15.0l*hour + min/4.0l + sec/240.0l;
  if(raread<0.0l || raread>360.0l) {
    cerr << "ERROR: mpc80_readline cannot read a valid Right Ascension. Line:\n";
    cerr << lnfromfile;
    return(1);
  }
  *RA = raread;

  // Read the Declination
  // sign
  decsign = lnfromfile[44];
  if(decsign != '+' && decsign != '-') {
    cerr << "ERROR: mpc80_readline cannot read a valid sign for the Declination. Line:\n";
    cerr << lnfromfile;
    return(1);
  }
  // degrees
  sdat = {};
  for(i=45;i<47;i++) {
    sdat.push_back(lnfromfile[i]);
  }
  deg = stod(sdat);
  // minutes
  sdat = {};
  for(i=48;i<50;i++) {
    sdat.push_back(lnfromfile[i]);
  }
  if(sdat.size()>0) min = stod(sdat);
  else min = 0.0l;
  // seconds
  sdat = {};
  for(i=51;i<56;i++) {
    c = lnfromfile[i];
    if(c != ' ') sdat.push_back(c);
  }
  if(sdat.size()>0) sec = stod(sdat);
  else sec = 0.0l;
  // Convert the Declination to decimal degrees
  decread = deg + min/60.0l + sec/3600.0l;
  if(decsign == '-') decread *= -1.0l;
  if(decread<-90.0l || decread>90.0l) {
    cerr << "ERROR: mpc80_readline cannot read a valid Declination. Line:\n";
    cerr << lnfromfile;
    return(1);
  }
  *Dec = decread;
  // Read the magnitude
  sdat = {};
  for(i=65;i<70;i++) {
    c = lnfromfile[i];
    if(c != ' ') sdat.push_back(c);
  }
  if(sdat.size()>0) *mag = stod(sdat);
  else *mag = 0.0l;

  // Read the band
  band = {};
  band.push_back(lnfromfile[70]);

  // Read the obscode
  obscode={};
  for(i=77;i<80;i++) {
    obscode.push_back(lnfromfile[i]);
  }
  return(0);
}
  
// mpc80_mjd: February 01, 2023:
// Extract the date alone from a line from an MPC 80-column formatted file,
// convert it to MJD, and return it at double precision.
double mpc80_mjd(const string &lnfromfile)
{
  int i;
  int year=1900;
  int month=1;
  double day=1.0l;
  double MJD=0.0l;
  string sdat;
  char c='0';
  
  // Read the year
  sdat = {};
  for(i=15;i<19;i++) {
    sdat.push_back(lnfromfile[i]);
  }
  year = stoi(sdat);
  // Read the month
  sdat = {};
  for(i=20;i<22;i++) {
    sdat.push_back(lnfromfile[i]);
  }
  month = stoi(sdat);
  // Read the day
  sdat = {};
  for(i=23;i<32;i++) {
    c=lnfromfile[i];
    if(c!=' ') sdat.push_back(c);
  }
  day = stod(sdat);
  if(year<1900 || month<1 || month>12 || day < 1.0l || day > 32.0l) {
    cerr << "mpc80_mjd cannot read a valid date from the line:\n";
    cerr << lnfromfile;
    return(-1.0l);
  }
  // Successfully read the date: convert it to MJD.
  MJD = MPCcal2MJD(year,month,day);
  return(MJD);
}

// read_obscode_file: March 14, 2023
// Read a uniformly-formatted observatory code file into a vector
// of type 'observatory'. This will only work properly in the input
// observatory code file has a one-line header, and adheres to the
// format of the ones downloadable from https://minorplanetcenter.net/iau/lists/ObsCodesF.html
// with the lines corresponding to space-based observatories and
// roving observers removed, as they do not have the required
// columns for longitude, parallax cosine, and parallax sine (obslon,plxcos,plxsin).
int read_obscode_file(string obscodefile,  vector <observatory> &observatory_list)
{
  string lnfromfile;
  string stest;
  char obscode[MINSTRINGLEN];
  double obslon,plxcos,plxsin;
  obslon = plxcos = plxsin = 0.0l;
  observatory obs1 = observatory("X05",0l,0l,0l);
  ifstream instream1;
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
  return(0);  
}

// read_obscode_file2: March 14, 2023
// Read an observatory code file into a vector of type 'observatory'.
// The input observatory code file must have a one-line header, and
// must adheres to the format of the ones downloadable from
// https://minorplanetcenter.net/iau/lists/ObsCodesF.html
// This version 2 is more flexible than the original read_obscode_file,
// in that it handles the lines corresponding to space-based observatories,
// and roving observers, which do not have the required columns for longitude, 
// parallax cosine, and parallax sine (obslon,plxcos,plxsin).
// Note that although the function READS these lines without generating
// an error, downstream clients making use of the ouput observatory_list
// may not be able to handle observations from space-based observatories
// or roving observers. The expectation is that we are analyzing data from
// ground-based observatories, for which obslon, plxcos, and plxsin are
// all well-defined.
int read_obscode_file2(string obscodefile,  vector <observatory> &observatory_list, int verbose)
{
  string lnfromfile;
  string stest;
  char obscode[MINSTRINGLEN];
  double obslon,plxcos,plxsin;
  obslon = plxcos = plxsin = 0.0l;
  observatory obs1 = observatory("X05",0l,0l,0l);
  ifstream instream1;
  int i=0;
  int badread=0;
  long linenum=0;
  instream1.open(obscodefile);

  if(!instream1) {
    cerr << "can't open input file " << obscodefile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  linenum++;
  while (!instream1.eof() && !instream1.fail() && !instream1.bad())
    {
      // Wipe results of previous read.
      obslon = plxcos = plxsin = 0.0l;
      badread=0;
      // Read the next line
      getline(instream1,lnfromfile);
      linenum++;
      if(lnfromfile.size()>=30) {
	// Parse the line
	stest={};
	for(i=0; i<3; i++) stest.push_back(lnfromfile[i]);
	stringncopy01(obscode,stest,MINSTRINGLEN);
	stest={};
	for(i=4; i<13; i++) stest.push_back(lnfromfile[i]);
	try { obslon = stod(stest); }
	catch(...) { badread=1; }
	stest={};
	for(i=13; i<21; i++) stest.push_back(lnfromfile[i]);
	try { plxcos = stod(stest); }
	catch(...) { badread=1; }
	stest={};
	for(i=21; i<30; i++) stest.push_back(lnfromfile[i]);
	try { plxsin = stod(stest); }
	catch(...) { badread=1; }
	obs1 = observatory(obscode,obslon,plxcos,plxsin);
	observatory_list.push_back(obs1);
	if(badread>=1 && verbose>=1 && !instream1.eof() && !instream1.fail() && !instream1.bad()) {
	  cerr << "WARNING: could not read full data from line " << linenum << " of ObsCode file " << obscodefile << "\n";
	  cerr << "Here is the problem line: " << lnfromfile << "\n";
	}
      } else if(verbose>=1 && !instream1.eof() && !instream1.fail() && !instream1.bad()) {
	cerr << "WARNING: Line " << linenum << " of ObsCode file " << obscodefile << " is too short to hold valid data\n";
	cerr << "Offending line: " << lnfromfile << "\n";
      }
    }
  instream1.close();
  return(0);  
}


// read_detection_filemt: March 14, 2023:
// Read an input detection file for make_tracklets. This
// function is quite specified to the exact needs of
// make_tracklets, and not likely to be very generally
// useful.
int read_detection_filemt(string indetfile, int idcol, int mjdcol, int racol, int deccol, int magcol,int bandcol, int obscodecol, vector <det_obsmag_indvec> &detvec, int forcerun)
{
  det_obsmag_indvec o1 = det_obsmag_indvec(0l,0l,0l,0l,0l,0l,"null",0l,"V","I11",0,{});
  ifstream instream1;
  string lnfromfile,stest;
  int lct,j,reachedeof,idread,mjdread,raread,decread,magread,bandread,obscoderead;
  lct = j = reachedeof = idread = mjdread = raread = decread = magread = bandread = obscoderead = 0;
  long i = 0;
  char c='0';
  double MJD,RA,Dec,mag;
  MJD = RA = Dec = mag = 0.0l;  
  char idstring[SHORTSTRINGLEN];
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  
  instream1.open(indetfile);
  if(!instream1) {
    cerr << "can't open input file " << indetfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  lct++;
  //cout << lnfromfile << "\n";
  reachedeof = 0;
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
    idread = mjdread = raread = decread = magread = bandread = obscoderead = 0;
    while(i<long(lnfromfile.size()) && lnfromfile.size()>=30 && reachedeof == 0) {
      // Note check on line length: it is completely impossible for a
      // line containing all the required quantities at minimum plausible
      // precision to be less than 30 characters long.
      c='0';
      stest="";
      while(i<long(lnfromfile.size()) && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==idcol) {
	stringncopy01(idstring,stest,SHORTSTRINGLEN);
	idread=1;
      } else if(j==mjdcol) {
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
      } else if(j==bandcol) {
	stringncopy01(band,stest,MINSTRINGLEN);
	bandread=1;
      } else if(j==obscodecol) {
	stringncopy01(obscode,stest,MINSTRINGLEN);
	obscoderead=1;
      }
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0 && lnfromfile.size()>=30) {
      if(!mjdread) {
	cerr << "ERROR: MJD not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	return(2);
      }
      if(!raread) {
	cerr << "ERROR: RA not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	return(2);
      }
      if(!decread) {
	cerr << "ERROR: Dec not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	return(2);
      }
      if(!idread) {
	if(forcerun) {
	  stringncopy01(idstring,"null",SHORTSTRINGLEN);
	  cout << "WARNING: ID not read from line " << detvec.size()+1 << " of input detection file " << indetfile << ".\n";
	  cout << "String ID will be set to null.\n";
	} else {
	  cerr << "ERROR: String ID not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	  return(2);
	}
      }
      if(!magread) {
	if(forcerun) {
	  mag = 99.999;
	  cout << "WARNING: magnitude not read from line " << detvec.size()+1 << " of input detection file " << indetfile << ".\n";
	  cout << "magnitude will be set to 99.999\n";
	} else {
	  cerr << "ERROR: magnitude not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	  return(2);
	}
      }
      if(!bandread) {
	if(forcerun) {
	  stringncopy01(band,"V",MINSTRINGLEN);
	  cout << "WARNING: photometric band not read from line " << detvec.size()+1 << " of input detection file " << indetfile << ".\n";
	  cout << "band will be set to V\n";
	} else {
	  cerr << "ERROR: photometric band not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	  return(2);
	}
      }
      if(!obscoderead) {
	if(forcerun) {
	  stringncopy01(obscode,"500",MINSTRINGLEN);
	  cout << "WARNING: observatory code not read from line " << detvec.size()+1 << " of input detection file " << indetfile << ".\n";
	  cout << "observatory code will be set to 500 (Geocentric)\n";
	} else {
	  cerr << "ERROR: observatory code not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	  return(2);
	}
      }
      o1=det_obsmag_indvec(MJD,RA,Dec,0l,0l,0l,idstring,mag,band,obscode,-lct,{});
      detvec.push_back(o1);
    }
  }
  instream1.close();
  
  if(reachedeof==1) { 
    cout << "Input file " << indetfile << " read successfully to the end.\n";
    return(0);
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << indetfile << " before the end\n";
    return(1);
  } else if(reachedeof==-1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else return(reachedeof);
}

// read_detection_filemt2: April 18, 2023:
// Read an input detection file for make_tracklets. This
// function is quite specified to the exact needs of
// make_tracklets, and not likely to be very generally
// useful.
int read_detection_filemt2(string indetfile, int mjdcol, int racol, int deccol, int magcol, int idcol, int bandcol, int obscodecol, int trail_len_col, int trail_PA_col, int sigmag_col, int sig_across_col, int sig_along_col, int known_obj_col, int det_qual_col, vector <hldet> &detvec, int verbose, int forcerun)
{
  double MJD,RA,Dec;
  MJD = RA = Dec = 0.0l;
  float mag, trail_len, trail_PA, sigmag, sig_across, sig_along;
  mag = -99.99;
  trail_len = 0.0;
  trail_PA = 90.0;
  sigmag = 9.999;
  sig_across = sig_along = 1.0;
  char idstring[SHORTSTRINGLEN];
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  stringncopy01(idstring,"",SHORTSTRINGLEN);
  stringncopy01(obscode,"500",MINSTRINGLEN);
  stringncopy01(band,"V",MINSTRINGLEN);
  int image = -1;
  long known_obj, det_qual, index;
  known_obj = det_qual = index = -1;
  hldet o1 = hldet(MJD, RA, Dec, mag, trail_len, trail_PA, sigmag, sig_across, sig_along, image, idstring, band, obscode, known_obj, det_qual, index);
  ifstream instream1;
  string lnfromfile,stest;
  int lct,j,reachedeof,mjdread,raread,decread,magread,idread,bandread,obscoderead;
  lct = j = reachedeof = mjdread = raread = decread = magread = idread = bandread = obscoderead = 0;
  int trail_len_read, trail_PA_read, sigmag_read, sig_across_read, sig_along_read, known_obj_read, det_qual_read;
  trail_len_read = trail_PA_read = sigmag_read = sig_across_read = sig_along_read = known_obj_read = det_qual_read = 0;
  long i = 0;
  char c='0';
  
  
  instream1.open(indetfile);
  if(!instream1) {
    cerr << "can't open input file " << indetfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  lct++;
  //cout << lnfromfile << "\n";
  reachedeof = 0;
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
    mjdread = raread = decread = magread = idread = bandread = obscoderead = 0;
    trail_len_read = trail_PA_read = sigmag_read = sig_across_read = known_obj_read = det_qual_read = 0;
    MJD = RA = Dec = 0.0l;
    mag = -99.99;
    trail_len = 0.0;
    trail_PA = 90.0;
    sigmag = 9.999;
    sig_across = sig_along = 1.0;
    stringncopy01(idstring,"",SHORTSTRINGLEN);
    stringncopy01(obscode,"500",MINSTRINGLEN);
    stringncopy01(band,"V",MINSTRINGLEN);
    image = -1;
    known_obj = det_qual = index = -1;

    while(i<long(lnfromfile.size()) && lnfromfile.size()>=30 && reachedeof == 0) {
      // Note check on line length: it is completely impossible for a
      // line containing all the required quantities at minimum plausible
      // precision to be less than 30 characters long.
      c='0';
      stest="";
      while(i<long(lnfromfile.size()) && c!=',' && c!='\n' && c!=EOF) {
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
      } else if(j==trail_len_col) {
	trail_len=stod(stest);
	trail_len_read=1;
      } else if(j==trail_PA_col) {
	trail_PA=stod(stest);
	trail_PA_read=1;
      } else if(j==sigmag_col) {
	sigmag=stod(stest);
	sigmag_read=1;
      } else if(j==sig_across_col) {
	sig_across=stod(stest);
	sig_across_read=1;
      } else if(j==sig_along_col) {
	sig_along=stod(stest);
	sig_along_read=1;
      } else if(j==idcol) {
	stringncopy01(idstring,stest,SHORTSTRINGLEN);
	idread=1;
      } else if(j==bandcol) {
	stringncopy01(band,stest,MINSTRINGLEN);
	bandread=1;
      } else if(j==obscodecol) {
	stringncopy01(obscode,stest,MINSTRINGLEN);
	obscoderead=1;
      } else if(j==known_obj_col) {
	known_obj=stol(stest);
	known_obj_read=1;
      } else if(j==det_qual_col) {
	det_qual=stol(stest);
	det_qual_read=1;
      } 
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0 && lnfromfile.size()>=30) {
      if(!mjdread) {
	cerr << "ERROR: MJD not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	return(2);
      }
      if(!raread) {
	cerr << "ERROR: RA not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	return(2);
      }
      if(!decread) {
	cerr << "ERROR: Dec not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	return(2);
      }
      if(!magread) {
	if(forcerun) {
	  mag = 99.999;
	  if(verbose>=2) {
	    cerr << "WARNING: magnitude not read from line " << detvec.size()+1 << " of input detection file " << indetfile << ".\n";
	    cerr << "magnitude will be set to 99.999\n";
	  }
	} else {
	  cerr << "ERROR: magnitude not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	  return(2);
	}
      }
      if(!bandread) {
	if(forcerun) {
	  stringncopy01(band,"V",MINSTRINGLEN);
	  if(verbose>=2) {
	    cerr << "WARNING: photometric band not read from line " << detvec.size()+1 << " of input detection file " << indetfile << ".\n";
	    cerr << "band will be set to V\n";
	  }
	} else {
	  cerr << "ERROR: photometric band not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	  cerr << "bandcol = " << bandcol << " line: " << lnfromfile << "\n";
	  return(2);
	}
      }
      if(!obscoderead) {
	if(forcerun) {
	  stringncopy01(obscode,"500",MINSTRINGLEN);
	  if(verbose>=1) {
	    cerr << "WARNING: observatory code not read from line " << detvec.size()+1 << " of input detection file " << indetfile << ".\n";
	    cerr << "observatory code will be set to 500 (Geocentric)\n";
	  }
	} else {
	  cerr << "ERROR: observatory code not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	  return(2);
	}
      }
      if(!idread) {
	if(forcerun) {
	  stringncopy01(idstring,"null",SHORTSTRINGLEN);
	  if(verbose>=2) {
	    cerr << "WARNING: ID not read from line " << detvec.size()+1 << " of input detection file " << indetfile << ".\n";
	    cerr << "String ID will be set to null.\n";
	  }
	} else {
	  cerr << "ERROR: String ID not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
	  return(2);
	}
      }
      if(verbose>=2 && !trail_len_read) {
	cerr << "Warning: trail length not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
      }
      if(verbose>=2 && !trail_PA_read) {
	cerr << "Warning: trail PA not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
      }
      if(verbose>=2 && !sigmag_read) {
	cerr << "Warning: magnitude uncertainty sigmag not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
      }
      if(verbose>=2 && !sig_across_read) {
	cerr << "Warning: cross-trail astrometric uncertainty sig_across not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
      }
      if(verbose>=2 && !known_obj_read) {
	cerr << "Warning: known object specifier not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
      }
      if(verbose>=2 && !det_qual_read) {
	cerr << "Warning: detection quality specifier not read from line " << detvec.size()+1 << " of input detection file " << indetfile << "!\n";
      }
      o1=hldet(MJD, RA, Dec, mag, trail_len, trail_PA, sigmag, sig_across, sig_along, image, idstring, band, obscode, known_obj, det_qual, -lct);
      detvec.push_back(o1);
    }
  }
  instream1.close();
  
  if(reachedeof==1) { 
    cout << "Input file " << indetfile << " read successfully to the end.\n";
    return(0);
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << indetfile << " before the end\n";
    return(1);
  } else if(reachedeof==-1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else return(reachedeof);
}

// read_pairdet_file: April 20, 2023:
// Read a paired detection file produced by make_tracklets_new.
int read_pairdet_file(string pairdetfile, vector <hldet> &detvec, int verbose)
{
  double MJD, RA, Dec;
  MJD = RA = Dec = 0.0l;
  float mag, trail_len, trail_PA, sigmag, sig_across, sig_along;
  mag =  trail_len = trail_PA = sigmag = sig_across = sig_along = 0.0;
  int image = 0;
  char idstring[SHORTSTRINGLEN];
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  stringncopy01(idstring,"",SHORTSTRINGLEN);
  stringncopy01(obscode,"500",MINSTRINGLEN);
  stringncopy01(band,"V",MINSTRINGLEN);
  long known_obj, det_qual, index;
  known_obj = det_qual = index = 0;
  hldet o1 = hldet(MJD, RA, Dec, mag, trail_len, trail_PA, sigmag, sig_across, sig_along, image, idstring, band, obscode, known_obj, det_qual, index);
  ifstream instream1;
  string lnfromfile,stest;
  int badread=0;
  int reachedeof=0;
  int startpoint=0;
  int endpoint=0;

  detvec={};
  
  instream1.open(pairdetfile);
  if(!instream1) {
    cerr << "can't open input file " << pairdetfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  //cout << lnfromfile << "\n";
  reachedeof = 0;
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      // Read on.
      // Read the MJD
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { MJD = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read MJD string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      else badread=1;
      // Read the RA
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { RA = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read RA string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the Dec
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { Dec = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read Dec string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the mag
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { mag = stof(stest); }
	catch(...) { cerr << "ERROR: cannot read mag string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the trail_len
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { trail_len = stof(stest); }
      catch(...) { cerr << "ERROR: cannot read trail_len string " << stest << " from line " << lnfromfile << "\n";
	badread = 1; }
      } else badread=1;
      // Read the trail_PA
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { trail_PA = stof(stest); }
	catch(...) { cerr << "ERROR: cannot read trail_PA string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      // Read the sigmag
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { sigmag = stof(stest); }
	catch(...) { cerr << "ERROR: cannot read sigmag string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the sig_across
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { sig_across = stof(stest); }
	catch(...) { cerr << "ERROR: cannot read sig_across string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      // Read the sig_along
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { sig_along = stof(stest); }
	catch(...) { cerr << "ERROR: cannot read sig_along string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the image
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { image = stoi(stest); }
	catch(...) { cerr << "ERROR: cannot read image string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the idstring
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(idstring,stest,SHORTSTRINGLEN);
      else badread=1;
      // Read the band
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(band,stest,MINSTRINGLEN);
      else badread=1;
      // Read the obscode
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) 	stringncopy01(obscode,stest,MINSTRINGLEN);
      else badread=1;
      // Read the known_obj
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { known_obj = stol(stest); }
	catch(...) { cerr << "ERROR: cannot read known_obj string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the det_qual
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { det_qual = stol(stest); }
	catch(...) { cerr << "ERROR: cannot read det_qual string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the index
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { index = stol(stest); }
	catch(...) { cerr << "ERROR: cannot read origindex string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      if(badread==0) {
	o1 = hldet(MJD, RA, Dec, mag, trail_len, trail_PA, sigmag, sig_across, sig_along, image, idstring, band, obscode, known_obj, det_qual, index);
	detvec.push_back(o1);
      }
    } else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(badread!=0) {
      cerr << "ERROR reading paired detection file " << pairdetfile << "\n";
      cerr << "Last point was " << detvec.size() << "; last file line was " << lnfromfile << "\n";
      return(badread);
    }
  }
  instream1.close();

  if(badread!=0) {
    cerr << "ERROR reading paired detection file " << pairdetfile << "\n";
    return(badread);
  } 
  if(reachedeof==1) { 
    if(verbose>=1) cout << "Input file " << pairdetfile << " read successfully to the end.\n";
    return(0);
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << pairdetfile << " before the end\n";
    return(1);
  } else if(reachedeof==-1) {
    cerr << "ERROR: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else return(reachedeof);
}

// read_tracklet_file: April 20, 2023:
// Read a tracklet file produced by make_tracklets_new.
int read_tracklet_file(string trackletfile, vector <tracklet> &tracklets, int verbose)
{
  long Img1 = 0;
  double RA1 = 0.0l;
  double Dec1 = 0.0l;
  long Img2 = 0;
  double RA2 = 0.0l;
  double Dec2 = 0.0l;
  int npts = 0;
  long trk_ID = 0;
  tracklet one_tracklet = tracklet(Img1, RA1, Dec1, Img2, RA2, Dec2, npts, trk_ID);
  ifstream instream1;
  string stest,lnfromfile;
  int badread=0;
  int reachedeof=0;
  int startpoint=0;
  int endpoint=0;
  
  tracklets={};
  
  instream1.open(trackletfile);
  if(!instream1) {
    cerr << "can't open input file " << trackletfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  //cout << lnfromfile << "\n";
  reachedeof = 0;
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn

    if(reachedeof == 0) {
      // Read Img1
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { Img1 = stol(stest); }
	catch(...) { cerr << "ERROR: cannot read Img1 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      else badread=1;
      // Read RA1
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { RA1 = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read RA1 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read Dec1
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { Dec1 = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read Dec1 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read Img2
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { Img2 = stol(stest); }
	catch(...) { cerr << "ERROR: cannot read Img2 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      else badread=1;
      // Read RA2
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { RA2 = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read RA2 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read Dec2
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { Dec2 = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read Dec2 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read npts
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { npts = stoi(stest); }
	catch(...) { cerr << "ERROR: cannot read npts string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read trk_ID
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { trk_ID = stol(stest); }
	catch(...) { cerr << "ERROR: cannot read trk_ID string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      if(badread==0) {
	one_tracklet = tracklet(Img1, RA1, Dec1, Img2, RA2, Dec2, npts, trk_ID);
	tracklets.push_back(one_tracklet);
      }
      if(!instream1.eof() && !instream1.fail() && !instream1.bad() && badread!=0) {
	cerr << "ERROR reading tracklet file " << trackletfile << "\n";
	return(badread);
      }
    }
  }
  instream1.close();

  if(badread!=0) {
    cerr << "ERROR reading tracklet file " << trackletfile << "\n";
    return(badread);
  } 
  if(reachedeof==1) { 
    if(verbose>=1) cout << "Input file " << trackletfile << " read successfully to the end.\n";
    return(0);
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << trackletfile << " before the end\n";
    return(1);
  } else if(reachedeof==-1) {
    cerr << "ERROR: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else return(reachedeof);
}

// read_longpair_file: April 20, 2023:
// Read a longpair file: e.g. trk2det or clust2det.
int read_longpair_file(string pairfile, vector <longpair> &pairvec, int verbose)
{
  long i1 = 0;
  long i2 = 0;
  longpair onepair = longpair(i1,i2);
  ifstream instream1;
  string stest,lnfromfile;
  int badread=0;
  int reachedeof=0;
  int startpoint=0;
  int endpoint=0;

  pairvec={};
  instream1.open(pairfile);
  if(!instream1) {
    cerr << "can't open input file " << pairfile << "\n";
    return(1);
  }
  reachedeof = 0;
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(!isdigit(lnfromfile[0])) {
      // Non-numerical: cannot be part of long pair.
      // Skip this possible header or comment line.
      continue;
    }
    if(reachedeof == 0) {
      // Read i1
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { i1 = stol(stest); }
	catch(...) { cerr << "ERROR: cannot read i1 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      else badread=1;
      // Read i2
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { i2 = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read i2 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      if(badread==0) {
	onepair = longpair(i1,i2);
	pairvec.push_back(onepair);
      }
      if(!instream1.eof() && !instream1.fail() && !instream1.bad() && badread!=0) {
	cerr << "ERROR reading long pair file " << pairfile << "\n";
	return(badread);
      }
    }
  }
  instream1.close();

  if(badread!=0) {
    cerr << "ERROR reading long pair file " << pairfile << "\n";
    return(badread);
  } 
  if(reachedeof==1) { 
    if(verbose>=1) cout << "Input file " << pairfile << " read successfully to the end.\n";
    return(0);
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << pairfile << " before the end\n";
    return(1);
  } else if(reachedeof==-1) {
    cerr << "ERROR: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else return(reachedeof);
}

// read_radhyp_file: April 20, 2023:
// Read a file containing heliolinc radial motion hypotheses
int read_radhyp_file(string hypfile, vector <hlradhyp> &accelmat, int verbose)
{
  double HelioRad = 0.0l;
  double R_dot = 0.0l;
  double R_dubdot = 0.0l;
  hlradhyp onehyp = hlradhyp(HelioRad,R_dot,R_dubdot);
  ifstream instream1;
  string stest,lnfromfile;
  int badread=0;
  int reachedeof=0;
  int startpoint=0;
  int endpoint=0;
  
  accelmat={};
    
  instream1.open(hypfile);
  if(!instream1) {
    cerr << "can't open input file " << hypfile << "\n";
    return(1);
  }
  reachedeof = 0;
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(!isdigit(lnfromfile[0])) {
      // Non-numerical: cannot be part of heliocentric motion hypothesis
      // Skip this possible header or comment line.
      continue;
    }
    if(reachedeof == 0) {
      // Read HelioRad
      startpoint=0;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { HelioRad = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read HelioRad string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read R_dot
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { R_dot = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read R_dot string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read R_dubdot
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_sv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { R_dubdot = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read R_dubdot string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      if(badread==0) {
	onehyp = hlradhyp(HelioRad,R_dot,R_dubdot);
	accelmat.push_back(onehyp);
      }
      if(!instream1.eof() && !instream1.fail() && !instream1.bad() && badread!=0) {
	cerr << "ERROR reading heliolcentric motion hypothesis file " << hypfile << "\n";
	cerr << "Last line was " << lnfromfile << " , last point was " << accelmat.size() << "\n";
	return(badread);
      }
    }
  }
  instream1.close();

  if(badread!=0) {
    cerr << "ERROR reading long pair file " << hypfile << "\n";
    return(badread);
  } 
  if(reachedeof==1) { 
    if(verbose>=1) cout << "Input file " << hypfile << " read successfully to the end.\n";
    return(0);
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << hypfile << " before the end\n";
    return(1);
  } else if(reachedeof==-1) {
    cerr << "ERROR: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else return(reachedeof);
}

// read_clustersum_file: April 21, 2023:
// Read a cluster summary file produced by heliolinc_new or link_refine_Herget_new.
int read_clustersum_file(string sumfile, vector <hlclust> &clustvec, int verbose)
{
  long clusternum=0;
  double posRMS,velRMS,totRMS,astromRMS;
  posRMS = velRMS = totRMS = astromRMS = 0.0l;
  int pairnum=0;
  double timespan=0.0l;
  int uniquepoints=0;
  int obsnights=0;
  double metric=0.0l;
  char rating[SHORTSTRINGLEN];
  stringncopy01(rating,"NULL",SHORTSTRINGLEN);
  double heliohyp0,heliohyp1,heliohyp2;
  heliohyp0 = heliohyp1 = heliohyp2 = 0.0l;
  double posX,posY,posZ,velX,velY,velZ;
  posX = posY = posZ = velX = velY = velZ = 0.0l;
  double orbit_a,orbit_e,orbit_MJD;
  orbit_a = orbit_e = orbit_MJD = 0.0l;
  double orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ;
  orbitX = orbitY = orbitZ = orbitVX = orbitVY = orbitVZ = 0.0l;
  long orbit_eval_count=0;
  hlclust onecluster = hlclust(clusternum, posRMS, velRMS, totRMS, astromRMS, pairnum, timespan, uniquepoints, obsnights, metric, rating, heliohyp0, heliohyp1, heliohyp2, posX, posY, posZ, velX, velY, velZ, orbit_a, orbit_e, orbit_MJD, orbitX, orbitY, orbitZ, orbitVX, orbitVY,  orbitVZ, orbit_eval_count);
  ifstream instream1;
  string lnfromfile,stest;
  int badread=0;
  int reachedeof=0;
  int startpoint=0;
  int endpoint=0;
  
  clustvec = {};
  
  instream1.open(sumfile);
  if(!instream1) {
    cerr << "can't open input file " << sumfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  //cout << lnfromfile << "\n";
  reachedeof = 0;
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      // Read on.
      // Read the clusternum
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { clusternum = stol(stest); }
	catch(...) { cerr << "ERROR: cannot read clusternum string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      else badread=1;
      // Read the posRMS
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { posRMS = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read posRMS string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the velRMS
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { velRMS = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read velRMS string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
       // Read the totRMS
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { totRMS = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read totRMS string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the astromRMS
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { astromRMS = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read astromRMS string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
     // Read the pairnum
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { pairnum = stoi(stest); }
	catch(...) { cerr << "ERROR: cannot read pairnum string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the timespan
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { timespan = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read timespan string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the number of unique points
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { uniquepoints = stoi(stest); }
      catch(...) { cerr << "ERROR: cannot read uniquepoints " << stest << " from line " << lnfromfile << "\n";
	badread = 1; }
      } else badread=1;
      // Read the obsnights
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { obsnights = stoi(stest); }
	catch(...) { cerr << "ERROR: cannot read obsnights string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      }
      // Read the metric
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { metric = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read metric string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the rating
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(rating,stest,SHORTSTRINGLEN);
      else badread=1;
      // Read heliohyp0
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { heliohyp0 = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read heliohyp0 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read heliohyp1
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { heliohyp1 = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read heliohyp1 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read heliohyp2
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { heliohyp2 = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read heliohyp2 string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read posX
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { posX = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read posX string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read posY
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { posY = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read posY string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read posZ
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { posZ = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read posZ string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read velX
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { velX = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read velX string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read velY
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { velY = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read velY string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read velZ
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { velZ = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read velZ string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbit_a
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbit_a = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read orbit_a string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbit_e
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbit_e = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read orbit_e string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbit_MJD
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbit_MJD = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read orbit_MJD string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbitX
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbitX = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read orbitX string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbitY
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbitY = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read orbitY string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbitZ
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbitZ = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read orbitZ string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbitVX
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbitVX = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read orbitVX string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbitVY
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbitVY = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read orbitVY string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbitVZ
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbitVZ = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read orbitVZ string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      // Read the orbit_eval_count
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) {
	try { orbit_eval_count = stol(stest); }
	catch(...) { cerr << "ERROR: cannot read orbit_eval_Count string " << stest << " from line " << lnfromfile << "\n";
	  badread = 1; }
      } else badread=1;
      if(badread==0) {
	onecluster = hlclust(clusternum, posRMS, velRMS, totRMS, astromRMS, pairnum, timespan, uniquepoints, obsnights, metric, rating, heliohyp0, heliohyp1, heliohyp2, posX, posY, posZ, velX, velY, velZ, orbit_a, orbit_e, orbit_MJD, orbitX, orbitY, orbitZ, orbitVX, orbitVY,  orbitVZ, orbit_eval_count);
	clustvec.push_back(onecluster);
      }
    } else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(badread!=0) {
      cerr << "ERROR reading cluster summary file " << sumfile << "\n";
      cerr << "Last point was " << clustvec.size() << "; last file line was " << lnfromfile << "\n";
      return(badread);
    }
  }
  instream1.close();

  if(badread!=0) {
    cerr << "ERROR reading cluster summary file " << sumfile << "\n";
    return(badread);
  } 
  if(reachedeof==1) { 
    if(verbose>=1) cout << "Input file " << sumfile << " read successfully to the end.\n";
    return(0);
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << sumfile << " before the end\n";
    return(1);
  } else if(reachedeof==-1) {
    cerr << "ERROR: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else return(reachedeof);
}

// avg_extrema: March 14, 2023: Given an input vector x, find and return
// the average of the extrema: that is, (xmin + xmax)/2.
double avg_extrema(const vector <double> &x) {
  long i=0;
  double xmin=x[0];
  double xmax=x[0];

  for(i=1;i<long(x.size());i++) {
    if(x[i]>xmax) xmax=x[i];
    if(x[i]<xmin) xmin=x[i];
  }
  return(xmin*0.5l + xmax*0.5l);
}

// read_image_file: March 14, 2023: Read an input file
// containing MJD, RA, Dec, obscode for a set of images,
// and partially load a vector of type img_log03.
int read_image_file(string inimfile, vector <img_log03> &img_log)
{
  img_log03 imlog = img_log03(0.0,0.0,0.0,"I11",0,0);
  ifstream instream1;
  int reachedeof,i,j;
  reachedeof = i = j = 0;
  char c = '0';
  double MJD, RA, Dec;
  MJD = RA = Dec = 0.0l;
  string lnfromfile,stest;
  char obscode[MINSTRINGLEN];

  img_log={};
  
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
    MJD=0.0l;
    while(i<long(lnfromfile.size()) && reachedeof == 0) {
      stest="";
      c='0';
      while(i<long(lnfromfile.size()) && c!=',' && c!=' ' && c!='\n' && c!=EOF) {
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
    if((reachedeof == 0 || reachedeof == 1) && MJD>0.0l) {
      // Requirement of MJD>0.0 tests that we read a plausibly
      // valid line.
      imlog=img_log03(MJD,RA,Dec,obscode,0,0);
      img_log.push_back(imlog);
    }
  }
  instream1.close();
  if(reachedeof==1) {
    cout << "Input file " << inimfile << " read successfully to the end.\n";
    return(0);
  }
  else if(reachedeof==-1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else {
    cerr << "Warning: unknown file read problem\n";
    return(3);
  }
}

// read_image_file: April 18, 2023: Read an input file
// containing MJD, RA, Dec, obscode for a set of images,
// and partially load a vector of type hlimage.
int read_image_file(string inimfile, vector <hlimage> &img_log)
{
  hlimage imlog = hlimage(0.0l, 0.0l, 0.0l, "500", 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0, 0);
  ifstream instream1;
  int reachedeof,i,j;
  reachedeof = i = j = 0;
  char c = '0';
  double MJD, RA, Dec;
  MJD = RA = Dec = 0.0l;
  string lnfromfile,stest;
  char obscode[MINSTRINGLEN];

  img_log={};
  
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
    if(!isdigit(lnfromfile[0])) {
      // This cannot be a valid MJD -- maybe it's a header line.
      continue;
    }
    i=0;
    j = 0;
    c='0';
    MJD=0.0l;
    while(i<long(lnfromfile.size()) && reachedeof == 0) {
      stest="";
      c='0';
      while(i<long(lnfromfile.size()) && c!=',' && c!=' ' && c!='\n' && c!=EOF) {
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
    if((reachedeof == 0 || reachedeof == 1) && MJD>0.0l) {
      // Requirement of MJD>0.0 tests that we read a plausibly
      // valid line.
      imlog=hlimage(MJD,RA,Dec,obscode, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0, 0);
      img_log.push_back(imlog);
    }
  }
  instream1.close();
  if(reachedeof==1) {
    cout << "Input file " << inimfile << " read successfully to the end.\n";
    return(0);
  }
  else if(reachedeof==-1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else {
    cerr << "Warning: unknown file read problem\n";
    return(3);
  }
}

// read_image_file2: April 20, 2023: Read an input file
// containing MJD, RA, Dec, obscode for a set of images,
// and fully load a vector of type hlimage.
int read_image_file2(string inimfile, vector <hlimage> &img_log)
{
  hlimage imlog = hlimage(0.0l, 0.0l, 0.0l, "500", 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0, 0);
  ifstream instream1;
  int reachedeof,i,j;
  reachedeof = i = j = 0;
  char c = '0';
  double MJD, RA, Dec;
  MJD = RA = Dec = 0.0l;
  string lnfromfile,stest;
  char obscode[MINSTRINGLEN];
  double X, Y, Z, VX, VY, VZ;
  X = Y = Z = VX = VY = VZ = 0.0l;
  long startind=0;
  long endind=0;
  
  img_log={};
  
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
    if(!isdigit(lnfromfile[0])) {
      // This cannot be a valid MJD -- maybe it's a header line.
      continue;
    }
    i=0;
    j = 0;
    c='0';
    MJD=0.0l;
    while(i<long(lnfromfile.size()) && reachedeof == 0) {
      stest="";
      c='0';
      while(i<long(lnfromfile.size()) && c!=',' && c!=' ' && c!='\n' && c!=EOF) {
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
      else if(j==5) X=stod(stest);
      else if(j==6) Y=stod(stest);
      else if(j==7) Z=stod(stest);
      else if(j==8) VX=stod(stest);
      else if(j==9) VY=stod(stest);
      else if(j==10) VZ=stod(stest);
      else if(j==11) startind=stol(stest);
      else if(j==12) endind=stol(stest);
    }
    if((reachedeof == 0 || reachedeof == 1) && MJD>0.0l) {
      // Requirement of MJD>0.0 tests that we read a plausibly
      // valid line.
      imlog=hlimage(MJD,RA,Dec,obscode, X, Y, Z, VX, VY, VZ, startind, endind);
      img_log.push_back(imlog);
    }
  }
  instream1.close();
  if(reachedeof==1) {
    cout << "Input file " << inimfile << " read successfully to the end.\n";
    return(0);
  }
  else if(reachedeof==-1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else {
    cerr << "Warning: unknown file read problem\n";
    return(3);
  }
}

#define DEBUGB 0
#define DEBUG 0
#define INTEGERIZING_SCALEFAC 100.0l // We divide state vectors by this value to integerize
                                     // them. Given a standard integer with a range
                                     // of +/- 2^31 = 2.15e9, this gives the state vectors
                                     // a range of +/- 2.15e11 km = 1400 AU, which is
                                     // comfortably larger than the solar system.
#define BINSEARCHMAX 100
#define REF_GEODIST 1.0 // Value of geocentric distance to which the user-defined
                             // clustering radius is normalized (AU). In general, the
                             // clustering radius is scaled linearly with geocentric distance.
#define NIGHTSTEP 0.3 // Minimum interval in days between successive points
                           // in a tracklet, to enable them to be counted as being
                           // on separate nights.
#define EPH_INTERP_POLYORDER 5 // Order of polynomial for interpolating JPL ephemerides.
#define TIMECONVSCALE 4.0l // The characteristic timescale used to convert velocities
                           // to distance units is equal to the full temporal span
                           // divided by TIMECONVSCALE.
#define FTOL_HERGET_SIMPLEX 1e-5l
#define ESCAPE_SCALE 0.99l // If the input velocity is above escape velocity, we
                           // scale it down by this amount.
#define MAX_SHUTTER_CORR 10.0 // Implied shutter corrections larger than this value,
                              // in seconds, are implausible and will cause link_refine_Herget
                              // to exit with an error.

// load_image_table: March 14, 2023: Construct an
// image table in the form of a vector of type img_log03.
// If the input image log is non-empty, assume it contains
// the correct MJD, RA, and Dec, and augment it with
// index information based on the input detvec. If the input
// image log vector is empty, infer the number of images,
// MJD, and approximate boresight RA, Dec from the entries
// in the detection vector.
int load_image_table(vector <img_log03> &img_log, const vector <det_obsmag_indvec> &detvec)
{
  img_log03 imlog = img_log03(0.0,0.0,0.0,"I11",0,0);
  vector <img_log03> img_log_tmp = img_log;
  img_log = {};
  // We make a copy of the input image log and then wipe the original,
  // because we are going to reload the original only with images that
  // match detections in the detection catalog: we won't track images
  // that had no detections.

  point3d p3 = point3d(0,0,0);
  point3d p3avg = point3d(0,0,0);
  int imct,detct,startind,endind,i;
  imct = detct = startind = endind = i = 0;
  double mjdnorm,mjdmean,tdelt;
  mjdnorm = mjdmean = tdelt = 0.0l;
  vector <double> x;
  vector <double> y;
  vector <double> z;

  if(DEBUGB==1) cout << "Inside load_image_table\n";
  
  if(img_log_tmp.size() > 0) {
    // We received an input image table, and all we have to do is
    // add the detection information to it.

    // Find the indices in the time-sorted detection file
    // that correspond to the earliest and latest detections
    // on each image, and load these values into imglog02.
    detct=0;
    for(imct=0;imct<long(img_log_tmp.size());imct++) {
      while(detct<long(detvec.size()) && detvec[detct].MJD < img_log_tmp[imct].MJD-IMAGETIMETOL/SOLARDAY) detct++; //Not on any image
      if(detct<long(detvec.size()) && fabs(detvec[detct].MJD-img_log_tmp[imct].MJD)<=IMAGETIMETOL/SOLARDAY && stringnmatch01(detvec[detct].obscode,img_log_tmp[imct].obscode,3)==0) {
	// This should be the first detection on image imct.
	img_log_tmp[imct].startind = detct;
	while(detct<long(detvec.size()) && fabs(detvec[detct].MJD-img_log_tmp[imct].MJD)<=IMAGETIMETOL/SOLARDAY && stringnmatch01(detvec[detct].obscode,img_log_tmp[imct].obscode,3)==0) detct++; //Still on this same image
	// This should be the first detection on the next image
	img_log_tmp[imct].endind = detct;
      }
      if(img_log_tmp[imct].startind >= 0 && img_log_tmp[imct].endind > 0) {
	img_log.push_back(img_log_tmp[imct]);
      }
    }
  } else {
    if(DEBUGB==1) cout << "Creating new image table\n";

    // No input image file was supplied: we have to create one from
    // the sorted detection file.
    mjdnorm = 1.0;
    mjdmean = detvec[0].MJD;
    startind=0;
    for(i=1;i<long(detvec.size());i++) {
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
	//Load it into the vector with mean MJD for all images.
	if(DEBUGB==1) cout << "Working on image " << img_log.size() << ", detections from " << startind << " to " << endind << "\n";
	imlog = img_log03(mjdmean,0.0,0.0,detvec[endind-1].obscode,startind,endind);
	img_log.push_back(imlog);
	// Set up for the next image, starting with detvec[i].MJD;
	mjdmean = detvec[i].MJD;
	mjdnorm = 1.0;
	startind=i;
      }
    }
    // Account for the final image.
    if(isnormal(mjdnorm)) {
      endind=detvec.size(); // Used to be endind=i, change eliminated an OS-dependent segfault.
      if(DEBUGB==1) cout << "Working on final image, " << img_log.size() << ", detections from " << startind << " to " << endind << "\n";
      mjdmean /= mjdnorm;
      //Load it into the vector with mean MJD for all images,
      // and increment image count.
      imlog = img_log03(mjdmean,0.0,0.0,detvec[endind-1].obscode,startind,endind);
      img_log.push_back(imlog);
    }

    //We've now loaded the mean MJDs and the starting and ending
    //detection table indices for each image; it still remains to
    //get the mean RA and Dec.
   
    long detnum = detvec.size();
    long imnum = img_log.size();
    cout << img_log.size() << " unique images were identified.\n";
    cout << "Given our total of " << detvec.size() << " detections,\n";
    cout << "we have " << double(detvec.size())/double(img_log.size()) << " detections per image, on average\n";

    // Find the number of detections and the average RA, Dec on each image.
    // We perform the average after projection onto the unit circle, to
    // avoid wrapping issues.
    detct=imct=0;
    while( imct<imnum && detct<detnum ) {
      int num_dets=0;
      p3avg = point3d(0,0,0);
      x = y = z ={};
      while(detct<detnum && detvec[detct].MJD < img_log[imct].MJD + IMAGETIMETOL/SOLARDAY && stringnmatch01(detvec[detct].obscode,img_log[imct].obscode,3)==0) {
	num_dets++; //Keep count of detections on this image
	p3 =  celeproj01(detvec[detct].RA,detvec[detct].Dec); // Project current detection
	x.push_back(p3.x); // Note that the projection from spherical to Cartesian
	y.push_back(p3.y); // coordinates avoids angle-wrapping issues.
	z.push_back(p3.z);
	detct++;
      }
      // If we got here, we must just have finished with an image.
      // Calculate the averages of the extrema:
      if(num_dets>0) {
	p3avg.x = avg_extrema(x); // Because the sources on an image could be distributed
	p3avg.y = avg_extrema(y); // very non-uniformly, the average of the extrema is
	p3avg.z = avg_extrema(z); // a better indicator for the image center than either
	                          // the mean or the median.
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
      imct++;
    }
  }
  return(0);
}

// load_image_table: April 18, 2023: Construct an
// image table in the form of a vector of type hlimage.
// If the input image log is non-empty, assume it contains
// the correct MJD, RA, and Dec, and augment it with
// index information based on the input detvec. If the input
// image log vector is empty, infer the number of images,
// MJD, and approximate boresight RA, Dec from the entries
// in the detection vector.
// Also (THIS DIFFERS FROM EARLIER OVERLOADED FUNCTION),
// calculate the observer's position and velocity at the instant
// of each image. 
int load_image_table(vector <hlimage> &img_log, const vector <hldet> &detvec, const vector <observatory> &observatory_list, const vector <double> &EarthMJD, const vector <point3d> &Earthpos, const vector <point3d> &Earthvel)
{
  hlimage imlog = hlimage(0.0l, 0.0l, 0.0l, "500", 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0, 0);
  vector <hlimage> img_log_tmp = img_log;
  img_log = {};
  // We make a copy of the input image log and then wipe the original,
  // because we are going to reload the original only with images that
  // match detections in the detection catalog: we won't track images
  // that had no detections.
  
  point3d p3 = point3d(0,0,0);
  point3d p3avg = point3d(0,0,0);
  int imct,detct,startind,endind,i;
  imct = detct = startind = endind = i = 0;
  double mjdnorm,mjdmean,tdelt;
  mjdnorm = mjdmean = tdelt = 0.0l;
  vector <double> x;
  vector <double> y;
  vector <double> z;
  double obslon = 0.0l;
  double plxcos = 0.0l;
  double plxsin = 0.0l;
  point3d obspos = point3d(0,0,0);
  point3d obsvel = point3d(0,0,0);
  int status=0;

  if(DEBUGB==1) cout << "Inside load_image_table\n";
  
  if(img_log_tmp.size() > 0) {
    // We received an input image table, and all we have to do is
    // add the detection information to it.

    // Find the indices in the time-sorted detection file
    // that correspond to the earliest and latest detections
    // on each image, and load these values into imglog02.
    detct=0;
    for(imct=0;imct<long(img_log_tmp.size());imct++) {
      while(detct<long(detvec.size()) && detvec[detct].MJD < img_log_tmp[imct].MJD-IMAGETIMETOL/SOLARDAY) detct++; //Not on any image
      if(detct<long(detvec.size()) && fabs(detvec[detct].MJD-img_log_tmp[imct].MJD)<=IMAGETIMETOL/SOLARDAY && stringnmatch01(detvec[detct].obscode,img_log_tmp[imct].obscode,3)==0) {
	// This should be the first detection on image imct.
	img_log_tmp[imct].startind = detct;
	while(detct<long(detvec.size()) && fabs(detvec[detct].MJD-img_log_tmp[imct].MJD)<=IMAGETIMETOL/SOLARDAY && stringnmatch01(detvec[detct].obscode,img_log_tmp[imct].obscode,3)==0) detct++; //Still on this same image
	// This should be the first detection on the next image
	img_log_tmp[imct].endind = detct;
      }
      if(img_log_tmp[imct].startind >= 0 && img_log_tmp[imct].endind > 0) {
	// This image is good: calculate the observer's position and velocity
	// Look up observatory coordinates for this image.
	status = obscode_lookup(observatory_list,img_log_tmp[imct].obscode,obslon,plxcos,plxsin);
	if(status>0) {
	  cerr << "ERROR: obscode_lookup failed for observatory code " << img_log[imct].obscode << "\n";
	  return(3);
	}
	// Calculate observer's exact heliocentric position and velocity.
	observer_baryvel01(img_log_tmp[imct].MJD, 5, obslon, plxcos, plxsin, EarthMJD, Earthpos, Earthvel, obspos, obsvel);
	img_log_tmp[imct].X = obspos.x;
	img_log_tmp[imct].Y = obspos.y;
	img_log_tmp[imct].Z = obspos.z;
	img_log_tmp[imct].VX = obsvel.x;
	img_log_tmp[imct].VY = obsvel.y;
	img_log_tmp[imct].VZ = obsvel.z;
      
	// Load the image to the output vector
	img_log.push_back(img_log_tmp[imct]);
      }
    }
  } else {
    if(DEBUGB==1) cout << "Creating new image table\n";

    // No input image file was supplied: we have to create one from
    // the sorted detection file.
    mjdnorm = 1.0;
    mjdmean = detvec[0].MJD;
    startind=0;
    for(i=1;i<long(detvec.size());i++) {
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
	// Calculate the mean MJD
	if(isnormal(mjdnorm)) mjdmean /= mjdnorm;
	else mjdmean = 0.0;
	if(!isnormal(mjdmean) || mjdmean<=0.0l) {
	  cerr << "ERROR: invalid mean MJD for image " << img_log.size() << "\n";
	  return(4);
	}
	// Look up the observatory code, so we can get the observer's
	// exact position and velocity.
	status = obscode_lookup(observatory_list,detvec[endind-1].obscode,obslon,plxcos,plxsin);
	if(status>0) {
	  cerr << "ERROR: obscode_lookup failed for observatory code " << img_log[imct].obscode << "\n";
	  return(3);
	}
	// Calculate observer's exact heliocentric position and velocity.
	observer_baryvel01(mjdmean, 5, obslon, plxcos, plxsin, EarthMJD, Earthpos, Earthvel, obspos, obsvel);
	imlog = hlimage(mjdmean,0.0,0.0,detvec[endind-1].obscode,obspos.x,obspos.y,obspos.z,obsvel.x,obsvel.y,obsvel.z,startind,endind);
	//Load it into the vector with mean MJD for all images.
	if(DEBUGB==1) cout << "Working on image " << img_log.size() << ", detections from " << startind << " to " << endind << "\n";
	img_log.push_back(imlog);
	// Set up for the next image, starting with detvec[i].MJD;
	mjdmean = detvec[i].MJD;
	mjdnorm = 1.0;
	startind=i;
      }
    }
    // Account for the final image.
    if(isnormal(mjdnorm)) {
      endind=detvec.size(); // Used to be endind=i, change eliminated an OS-dependent segfault.
      if(DEBUGB==1) cout << "Working on final image, " << img_log.size() << ", detections from " << startind << " to " << endind << "\n";
      mjdmean /= mjdnorm;
      // Look up the observatory code, so we can get the observer's
      // exact position and velocity.
      status = obscode_lookup(observatory_list,detvec[endind-1].obscode,obslon,plxcos,plxsin);
      if(status>0) {
	cerr << "ERROR: obscode_lookup failed for observatory code " << img_log[imct].obscode << "\n";
	return(3);
      }
      // Calculate observer's exact heliocentric position and velocity.
      observer_baryvel01(mjdmean, 5, obslon, plxcos, plxsin, EarthMJD, Earthpos, Earthvel, obspos, obsvel);
      imlog = hlimage(mjdmean,0.0,0.0,detvec[endind-1].obscode,obspos.x,obspos.y,obspos.z,obsvel.x,obsvel.y,obsvel.z,startind,endind);
      
      //Load it into the vector with mean MJD for all images,
      // and increment image count.
      imlog = hlimage(mjdmean,0.0,0.0,detvec[endind-1].obscode,obspos.x,obspos.y,obspos.z,obsvel.x,obsvel.y,obsvel.z,startind,endind);
      img_log.push_back(imlog);
    }

    //We've now loaded the mean MJDs and the starting and ending
    //detection table indices for each image; it still remains to
    //get the mean RA and Dec.
   
    long detnum = detvec.size();
    long imnum = img_log.size();
    cout << img_log.size() << " unique images were identified.\n";
    cout << "Given our total of " << detvec.size() << " detections,\n";
    cout << "we have " << double(detvec.size())/double(img_log.size()) << " detections per image, on average\n";

    // Find the number of detections and the average RA, Dec on each image.
    // We perform the average after projection onto the unit circle, to
    // avoid wrapping issues.
    detct=imct=0;
    while( imct<imnum && detct<detnum ) {
      int num_dets=0;
      p3avg = point3d(0,0,0);
      x = y = z ={};
      while(detct<detnum && detvec[detct].MJD < img_log[imct].MJD + IMAGETIMETOL/SOLARDAY && stringnmatch01(detvec[detct].obscode,img_log[imct].obscode,3)==0) {
	num_dets++; //Keep count of detections on this image
	p3 =  celeproj01(detvec[detct].RA,detvec[detct].Dec); // Project current detection
	x.push_back(p3.x); // Note that the projection from spherical to Cartesian
	y.push_back(p3.y); // coordinates avoids angle-wrapping issues.
	z.push_back(p3.z);
	detct++;
      }
      // If we got here, we must just have finished with an image.
      // Calculate the averages of the extrema:
      if(num_dets>0) {
	p3avg.x = avg_extrema(x); // Because the sources on an image could be distributed
	p3avg.y = avg_extrema(y); // very non-uniformly, the average of the extrema is
	p3avg.z = avg_extrema(z); // a better indicator for the image center than either
	                          // the mean or the median.
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
      imct++;
    }
  }
  return(0);
}


#undef DEBUGB

// load_image_indices: March 23, 2023: Load the starting and
// ending indices in an image table of the form used in the
// python-wrapped version of make_tracklets. 
int load_image_indices(vector <hlimage> &img_log, vector <hldet> &detvec, double imagetimetol, int forcerun)
{
  long imnum = img_log.size();
  long detnum = detvec.size();
  long imct,detct,startind,endind,i;
  imct = detct = startind = endind = i = 0;

  // Load indices in detvec. This is incidental to the main
  // purpose of load_image_indices(), but it's necessary and
  // this is a convenient place to do it. 
  for(detct=0;detct<detnum;detct++) detvec[detct].index = -detct;
  
  detct=0;
  for(imct=0;imct<imnum;imct++) {
    if(detct<detnum && fabs(detvec[detct].MJD-img_log[imct].MJD)<=imagetimetol && stringnmatch01(detvec[detct].obscode,img_log[imct].obscode,3)==0) {
      // This should be the first detection on image imct.
      img_log[imct].startind = detct;
      detvec[detct].image = imct;
      detct++;
      while(detct<detnum && fabs(detvec[detct].MJD-img_log[imct].MJD)<=imagetimetol && stringnmatch01(detvec[detct].obscode,img_log[imct].obscode,3)==0) {
	// Still on this same image
	detvec[detct].image = imct;
	detct++;
      }
      // This should be the first detection on the next image
      img_log[imct].endind = detct;
    } else if(detct<detnum && detvec[detct].MJD > img_log[imct].MJD+imagetimetol) {
      // The next detection is after this image in the ordered time sequence.
      // Therefore, no detections were found on this image
      img_log[imct].startind = img_log[imct].endind = 0;
    } else if(detct<detnum && fabs(detvec[detct].MJD-img_log[imct].MJD)<=imagetimetol && stringnmatch01(detvec[detct].obscode,img_log[imct].obscode,3)>0) {
      // The next detection overlaps this image in the ordered time sequence,
      // but comes after it in the alphabetical listing of obscodes.
      // Therefore, no detections were found on this image.
      img_log[imct].startind = img_log[imct].endind = 0;
    } else if(detct<detnum && (detvec[detct].MJD < img_log[imct].MJD-imagetimetol ||
			       (fabs(detvec[detct].MJD-img_log[imct].MJD)<=imagetimetol && stringnmatch01(detvec[detct].obscode,img_log[imct].obscode,3)<0))) {
      // The next detection is before this image in the ordered sequence, either
      // before it in time OR overlapping in time but before it in the alphabetical listing of obscodes.
      // Therefore, this detection must not appear on any image in the sequence.
      if(forcerun) {
	// With forcerun, we allow detections that aren't on any image,
	// even though the caller really should have made sure this couldn't happen.
	while(detct<detnum && (detvec[detct].MJD < img_log[imct].MJD-imagetimetol ||
			       (fabs(detvec[detct].MJD-img_log[imct].MJD)<=imagetimetol && stringnmatch01(detvec[detct].obscode,img_log[imct].obscode,3)<0))) {
	  detvec[detct].image = -1;
	  detct++;
	}
	// Now we must consider the possibility that we've advanced to
	// a detection that is on the current image.
	// This is the case where we advanced through a series of
	// bad detections until we arrived at a good detection.
	if(detct<detnum && fabs(detvec[detct].MJD-img_log[imct].MJD)<=imagetimetol && stringnmatch01(detvec[detct].obscode,img_log[imct].obscode,3)==0) {
	  // This should be the first detection on image imct.
	  img_log[imct].startind = detct;
	  detvec[detct].image = imct;
	  detct++;
	  while(detct<detnum && fabs(detvec[detct].MJD-img_log[imct].MJD)<=imagetimetol && stringnmatch01(detvec[detct].obscode,img_log[imct].obscode,3)==0) {
	    // Still on this same image
	    detvec[detct].image = imct;
	    detct++;
	  }
	  // This should be the first detection on the next image
	  img_log[imct].endind = detct;
	  // End special case that we advanced through a series of bad detections
	  // to arrive at a good detection.
	}
      } else {
	// forcerun is not on, meaning that it's not acceptable
	// for the image catalog not to span all the detections.
	cerr << "ERROR in load_image_indices: detection " << detct << " not on any image!\n";
	return(1);
      }
    } else if(detct<detnum) {
      // Logically excluded case
      cerr << "ERROR: logically excluded case 1 in load_image_indices\n";
      cerr << "Detection " << detct << " time, obscode: " << detvec[detct].MJD << " " << detvec[detct].obscode << "\n";
      cerr << "Image " << imct << " time, obscode: " << img_log[imct].MJD << " " << img_log[imct].obscode << "\n";
      if(!forcerun) return(1);
    } else if(detct>=detnum) {
      // We ran past the end of the detection catalog, apparently without running out of images
      // All remaining images will have no detections.
      cerr << "WARNING: The image log continues past the end of the detection catalog\n";
      img_log[imct].startind = img_log[imct].endind = 0;
    } else {
      cerr << "ERROR: logically excluded case 2 in load_image_indices\n";
      cerr << "Detection " << detct << " time, obscode: " << detvec[detct].MJD << " " << detvec[detct].obscode << "\n";
      cerr << "Image " << imct << " time, obscode: " << img_log[imct].MJD << " " << img_log[imct].obscode << "\n";
      if(!forcerun) return(1);
    }
  }
  // Deal with any left-over detections after the last image.
  if(detct<detnum) {
    cerr << "ERROR: Ran out of images at detection " << detct << ", short of the last detection at " << detnum << "\n";
    if(!forcerun) return(1);
    else {
      while(detct<detnum) {
	detvec[detct].image = -1;
	detct++;
      }
    }
  }
  return(0);
}


//find_pairs: March 24, 2023:  Create pairs, output a vector pairdets of type hldet;
// a vector indvecs of type vector <long>, with the same length as pairdets,
// giving the indices of all the detections paired with a given detection;
// and the vector pairvec of type longpair, giving all the pairs of detections.
int find_pairs(vector <hldet> &detvec, const vector <hlimage> &img_log, vector <hldet> &pairdets, vector <vector <long>> &indvecs, vector <longpair> &pairvec, double mintime, double maxtime, double imrad, double maxvel, int verbose)
{
  int imnum = img_log.size();
  int imct=0;
  long detct=0;
  long pdct=0; // count of detections that have been paired
  long pairct=0; // count of actual pairs
  xy_index xyind=xy_index(0.0, 0.0, 0);
  vector <xy_index> axyvec = {};
  double dist,pa;
  dist = pa = 0.0l;
  long dettarg=0;
  longpair onepair = longpair(0,0);
  vector <long> ivec1;

  pairvec={};
  pairdets={};
  indvecs = {};
  ivec1={};
  
  // Loop over images for image A
  for(imct=0;imct<imnum;imct++) {
    if(img_log[imct].endind<=0 || img_log[imct].endind<=img_log[imct].startind) continue; // No detections on this image.
    int apct=0;
    int adetct=0;
    // See if there are any images that might match
    vector <int> imagematches = {};
    int imatchcount = 0;
    int imtarg=imct+1;
    while(imtarg<imnum && img_log[imtarg].MJD < img_log[imct].MJD + maxtime) {
      double timediff = img_log[imtarg].MJD-img_log[imct].MJD;
      if(!isnormal(timediff) || timediff<0.0) {
	cerr << "WARNING: Negative time difference " << timediff << " encountered between images " << imct << " and " << imtarg << "\n";
      }
      // See if the images are close enough on the sky.
      double imcendist = distradec01(img_log[imct].RA, img_log[imct].Dec, img_log[imtarg].RA, img_log[imtarg].Dec);
      if(imcendist<2.0*imrad+maxvel*timediff && timediff>=mintime && img_log[imtarg].endind>0 && img_log[imtarg].endind>img_log[imtarg].startind) {
	if(DEBUG>=1) cout << "  pairs may exist between images " << imct << " and " << imtarg << ": dist = " << imcendist << ", timediff = " << timediff << "\n";
	imagematches.push_back(imtarg);
      }
      imtarg++;
    }
    if(verbose>=1) cout << "Looking for pairs for image " << imct << ": " << imagematches.size() << " later images are worth searching\n";
    int imatchnum = imagematches.size();
    if(imatchnum>0) {
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
      for(imatchcount=0;imatchcount<imatchnum;imatchcount++)
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
	for(detct=0 ; detct<long(axyvec.size()) ; detct++) {
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
	      ivec1={};
	      indvecs.push_back(ivec1);  // Load empty index vector
	      detvec[axyvec[detct].index].index = pdct; // Re-assign index to apply to paired detection vector
	      pdct++; // Increment count of paired detections
	      adetct++;
	      if(pdct!=long(pairdets.size()) || pdct!=long(indvecs.size())) {
		cerr << "\nERROR: PAIRED DETECTION MISMATCH: " << pdct << " vs " << pairdets.size() << " vs " << indvecs.size() << "\n";
		return(1);
	      }
	    }
	    // Record image B detections
	    for(matchct=0;matchct<matchnum;matchct++) {
	      matchpt = indexvec[matchct];
	      if(detvec[kdvec[matchpt].point.index].index<0) {
		//This detection has not yet been paired with any other.
		detvec[kdvec[matchpt].point.index].index *= -1; // Mark as paired by changing to positive sign
		pairdets.push_back(detvec[kdvec[matchpt].point.index]); // Load into paired detection vector
		ivec1={};
		indvecs.push_back(ivec1); // Load empty index vector
		detvec[kdvec[matchpt].point.index].index = pdct; // Re-assign index to apply to paired detection vector
		pdct++; // Increment count of paired detections
		if(pdct!=long(pairdets.size()) || pdct!=long(indvecs.size())) {
		  cerr << "\nERROR: PAIRED DETECTION MISMATCH: " << pdct << " vs " << pairdets.size() << " vs " << indvecs.size() << "\n";
		  return(1);
		}
	      }
	      // Write index values for both components of the
	      // new pair to the pair vector, regardless of whether
	      // the index values are pre-existing or newly assigned.
	      onepair = longpair(detvec[axyvec[detct].index].index,detvec[kdvec[matchpt].point.index].index);
	      pairvec.push_back(onepair);
	      pairct++;
	      apct++;
	      // Load index of each detection into the paired index vector of the other
	      if(kdvec[matchpt].point.index>=0 && kdvec[matchpt].point.index < long(detvec.size()) && axyvec[detct].index >=0 && axyvec[detct].index < long(detvec.size())) {
		if(detvec[kdvec[matchpt].point.index].index >= 0 && detvec[kdvec[matchpt].point.index].index < long(detvec.size()) && detvec[axyvec[detct].index].index >= 0 && detvec[axyvec[detct].index].index < long(detvec.size())) {
		  indvecs[detvec[axyvec[detct].index].index].push_back(detvec[kdvec[matchpt].point.index].index);
		  indvecs[detvec[kdvec[matchpt].point.index].index].push_back(detvec[axyvec[detct].index].index);
		} else {
		  cerr << "ERROR: trying to load out-of-range points to indvecs\n";
		  cerr << "Points are " <<  detvec[kdvec[matchpt].point.index].index << " and " << detvec[axyvec[detct].index].index  << "\n";
		  cerr << "Permitted range is 0 to " << detvec.size() << "\n";
		  return(8);
		}
	      } else {
		cerr << "ERROR: attempting to access out-of-range values in detvec\n";
		cerr << "Indices are " << kdvec[matchpt].point.index << " and " << axyvec[detct].index << "\n";
		cerr << "Permitted range is 0 to " << detvec.size() << "\n";
		return(9);
	      }				  
	    }
	    // Close if-statement checking if image A detection was matched to anything.
	  }
	  // Close loop over detections on source image (image A)
	}
	// Close loop over image B candidates
      }
      // Close if-statement checking if any images could match image A      
    }
    if(verbose>=1) cout << "Image " << imct << ": found " << adetct << " newly paired detections and a total of " << apct << " pairs.\n";
    // Close loop over images for image A
  }
  if(verbose>=1) cout << "Test count of paired detections: " << pdct << " " << pairdets.size() << "\n";
  if(verbose>=1) cout << "Test count of pairs: " << pairct << " " << pairvec.size() << "\n";
  return(0);
}

//merge_pairs: March 24, 2023: Given the output from find_pairs,
//merge pairs into tracklets with more than two points, if possible
int merge_pairs(const vector <hldet> &pairdets, vector <vector <long>> &indvecs, const vector <longpair> &pairvec, vector <tracklet> &tracklets, vector <longpair> &trk2det, int mintrkpts, double maxgcr, double minarc, double minvel, double maxvel, int verbose)
{
  long detnum = pairdets.size();
  long detct=0;
  long i = 0;
  long_index ppn = long_index(0,0);
  vector <long_index> pair_partner_num;
  long pdct=0;
  int istracklet=0;
  vector <hldet> ppset = {};
  vector <vector <long>> ppind = {};
  xy_index xyind=xy_index(0.0, 0.0, 0);
  vector <xy_index> axyvec = {};
  double dist = 0.0l;
  double pa = 0.0l;
  tracklet track1 = tracklet(0,0.0l,0.0l,0,0.0l,0.0l,0,0);
  longpair onepair = longpair(0,0);
  long j=0;
  long k=0;
  int biggest_tracklet = -1;
  int tracklet_size = 0;
  point3d_index p3di = point3d_index(0.0l,0.0l,0.0l,0);
  vector <point3d_index>   track_mrdi_vec;
  int trkptnum=0;
  int istimedup=1;
  vector <double> timevec;
  vector <double> xvec;
  vector <double> yvec;
  vector <long> detindexvec;
  double slopex,slopey,interceptx,intercepty,worsterr;
  slopex = slopey = interceptx = intercepty = worsterr = 0.0l;
  vector <double> fiterr = {};
  vector <double> fiterr2 = {};
  int worstpoint=-1;
  double dtref,dt,dx,dy,angvel;
  dtref = dt = dx = dy = angvel = 0.0l;
  double outra1,outdec1,outra2,outdec2;
  outra1 = outdec1 = outra2 = outdec2 = 0.0l;
  int rp1,rp2,instep;
  rp1=rp2=instep=0;

  if(detnum != long(indvecs.size())) {
    cerr << "ERROR: merge_pairs received input vectors pairdets and indvecs\n";
    cerr << "with different lengths (" << detnum << " and " << indvecs.size() << "\n";
    return(4);
  }
  if(verbose) {
    cout << "Input vector lengths: pairdets: " << detnum << ", indvecs: " << indvecs.size() << ", pairvec: " << pairvec.size() << "\n";
  }

  // Sanity-check indvecs
  cout << "merge_pairs is sanity-checking indvecs\n";
  for(detct=0; detct<detnum; detct++) {
    for(i=0; i<long(indvecs[detct].size()); i++) {
      if(indvecs[detct][i]<0 || indvecs[detct][i]>=detnum) {
	cerr << "ERROR: indvecs[" << detct << "][" << i << "] out of range: " << indvecs[detct][i] << "\n";
	cerr << "Acceptable range is 0 to " << detnum << "\n";
	return(9);
      }
    }
  }
  cout << "Sanity-check finished\n";
  
  tracklets={};
  trk2det={}; // Wipe output vectors.
  
  // Load a vector storing the number of pair-partners found for each detection.
  for(i=0;i<detnum;i++) {
    ppn = long_index(indvecs[i].size(),i);
    pair_partner_num.push_back(ppn);
  }
  if(detnum != long(pair_partner_num.size())) {
    cerr << "ERROR: newly constructed vector pair_partner_num in  merge_pairs\n";
    cerr << "is not the same length as pairdets vector (" << detnum << " vs " << pair_partner_num.size() << "\n";
    return(4);
  }
  // Sort the new vector by number of pair-partners
  sort(pair_partner_num.begin(), pair_partner_num.end(), lower_long_index());
  
  // Analyze paired detections in order of decreasing number of partners.
  // At the same time, load the tracklet file and the trk2det file
  cout << "Constructing tracklets, and loading output vectors\n";
  
  for(i=detnum-1; i>=0 ;i--) {
    pdct=pair_partner_num[i].index; // Decode from pair_partner_num sorted list to actual pairdets index.
    istracklet=0; // Assume there is no tracklet unless one is confirmed to exist.
    if(long(indvecs[pdct].size()) > mintrkpts-1) { // Use mintrkpts-1 because the root detection pdct
                                                     // is itself is a potential point in the tracklet
      if(verbose>=1) {
	cout << "Working on detection " << i << " = " << pdct << " of " << detnum << ", with " << pair_partner_num[i].lelem << " = " << indvecs[pdct].size() << " pair partners";
	if(verbose>=3) {
	  cout << ":\n";
	  for(j=0; j<long(indvecs[pdct].size()); j++) {
	    cout << indvecs[pdct][j] << ", ";
	  }
	  cout << "\n";
	} else cout << "\n";
      }
      // Detection number pdct is paired with more than one
      // other detection.
      // Project all of these pairs relative to detection pdct,
      // storing x,y projected coordinates in axyvec.
      axyvec={};
      ppset={};
      ppind={};
      for(j=0; j<long(indvecs[pdct].size()); j++) { // Loop over the pair-partners of detection pdct.
	detct = indvecs[pdct][j]; // detct is the pairdets index of a pair-partner to detection pdct.
	if(detct<0 || detct>=long(indvecs.size())) {
	  cerr << "Error: merge_pairs attempting to query detct=" << detct << ", out of range 0-" << long(indvecs.size()) << " = " << detnum << "\n";
	  return(8);
	}
	if(indvecs[detct].size()>0) {
	  // Detection detct hasn't already been allocated to a tracklet,
	  // and hence is available for inclusion in a new tracklet anchored by pdct.
	  // Project detct into an arc-WCS style x,y coords centered on pdct.
	  distradec02(pairdets[pdct].RA, pairdets[pdct].Dec, pairdets[detct].RA, pairdets[detct].Dec, &dist, &pa);
	  dist *= 3600.0L; // Convert distance from degrees to arcsec.
	  xyind = xy_index(dist*sin(pa/DEGPRAD),dist*cos(pa/DEGPRAD),detct);
	  axyvec.push_back(xyind);
	  ppset.push_back(pairdets[detct]);
	  ppind.push_back({}); // We need these vectors mainly just to have some way to store the
	                       // indices of mutually consistent pair partners on the next step
	}
      }
      if(verbose>=2) cout << "Loaded axyvec and ppset vectors OK, with sizes " << axyvec.size() << " and " << ppset.size() << "\n";
      if(axyvec.size() != ppset.size() || axyvec.size() != ppind.size()) {
	cerr << "ERROR: vectors of projected and original\n";
	cerr << "pair partner candidates do not have the same length!\n";
	cerr << axyvec.size() << ", " << ppset.size() << ", and " << ppind.size() << " must all be the same, and are not!\n";
	return(3);
      }
      // Perform n^2 search on the projected points stored in axyvec
      // to find the largest subset that lie along a consistent line.
      for(j=0; j<long(axyvec.size()); j++) {
	dtref = ppset[j].MJD - pairdets[pdct].MJD; // Time from anchor detection pdct to pair-partner j.
	if(dtref == 0) {
	  cerr << "ERROR: paired detections with no time separation!\n";
	  return(4);
	}
	// Make sure corresponding index vector in ppset is empty
	ppind[j] = {};
	// Count addition pair partners (besides ppset[j]) that plausibly
	// lie along the line defined by pdct and ppset[j].
	if(DEBUG>=2) cout << "Counting consistent pair partners\n";
	for(k=0; k<long(axyvec.size()); k++) {
	  if(j!=k) {
	    dt = ppset[k].MJD - pairdets[pdct].MJD; // Time from anchor detection pdct to pair-partner j.
	    // Find out if the projected x,y coords scale with time from pdct
	    // in a consistent way for detections j and k.
	    dx = axyvec[k].x - axyvec[j].x*(dt/dtref);
	    dy = axyvec[k].y - axyvec[j].y*(dt/dtref);
	    dist = sqrt(dx*dx + dy*dy); 
	    if(verbose>=3) cout << "Detection " << axyvec[j].index << ":" << axyvec[k].index << " dist = " << dist << "\n";
	    if(dist < 2.0*maxgcr) {
	      // Detections j, k, and pdct all lie along a plausibly
	      // linear, constant-velocity trajectory on the sky.
	      ppind[j].push_back(k); // Store detection k as possible tracklet partner for pdct and j.
	    }
	  }
	}
      }
      // Now, ppset stores all the possible pair-partners of detection pdct,
      // and ppind stores, for each one of these, the ppset indices of ADDITIONAL ones
      // lie on a potentially consistent trajectory with it: that is, are tracklet partners.
      // Find which detection in ppset has the largest number of possible tracklet partners.
      biggest_tracklet=-1;
      tracklet_size=0;
      for(j=0; j<long(ppset.size()); j++) {
	if(long(ppind[j].size())+2 > tracklet_size) {
	  tracklet_size = ppind[j].size()+2; //We add one for pdct, one for j, to get actual tracklet size
	  biggest_tracklet = j;
	  if(DEBUG>=2) cout << "bt = " << biggest_tracklet << ", size = " << tracklet_size << "\n";
	} else if(DEBUG>=2) cout << "not the biggest\n";
      }
      if(verbose>=2 && biggest_tracklet>=0) cout << "Biggest tracklet is " << biggest_tracklet << ", which corresponds to " << axyvec[biggest_tracklet].index << ", with size " << tracklet_size << "\n";
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
 	for(j=0; j<long(ppind[biggest_tracklet].size()); j++) {
	  p3di = point3d_index(ppset[ppind[biggest_tracklet][j]].MJD - pairdets[pdct].MJD,axyvec[ppind[biggest_tracklet][j]].x,axyvec[ppind[biggest_tracklet][j]].y,axyvec[ppind[biggest_tracklet][j]].index);
	  track_mrdi_vec.push_back(p3di);
	}
	// Sort track_mrdi_vec by time.
	sort(track_mrdi_vec.begin(), track_mrdi_vec.end(), lower_point3d_index_x());
	// Load time, x, y, and index vectors from sorted track_mrdi_vec.
	timevec=xvec=yvec={};
	detindexvec={};
	for(j=0;j<long(track_mrdi_vec.size());j++)
	  {
	    timevec.push_back(track_mrdi_vec[j].x);
	    xvec.push_back(track_mrdi_vec[j].y);
	    yvec.push_back(track_mrdi_vec[j].z);
	    detindexvec.push_back(track_mrdi_vec[j].index);
	  }
	if(track_mrdi_vec.size() != timevec.size() || track_mrdi_vec.size() != xvec.size() || track_mrdi_vec.size() != yvec.size() || track_mrdi_vec.size() != detindexvec.size()) {
	  cerr << "ERROR: vector length mismatch in vectors for tracklet-fitting!\n";
	  cerr << "Lengths of track_mrdi_vec, timevec, xvec, yvec, and detindexvec:\n";
	  cerr << track_mrdi_vec.size() << ", " << timevec.size()  << ", " << xvec.size()  << ", " << yvec.size()  << ", " << detindexvec.size() << "\n";
	  return(6);
	}
 	if(DEBUG>=2) {
	  cout << "First iteration linear fit vectors:\n";
	  for(j=0; j<long(timevec.size()); j++) {
	    cout << detindexvec[j] << " " << timevec[j] << " " << xvec[j] << " " << yvec[j] << "\n";
	  }
	}

	// Perform fit to projected x coordinate as a function of time
	linfituw01(timevec, xvec, slopex, interceptx);
 	// Perform fit to projected y coordinate as a function of time
	linfituw01(timevec, yvec, slopey, intercepty);
	// Load vector of residuals
	fiterr = {};
	for(j=0; j<long(timevec.size()); j++) {
	  fiterr.push_back(sqrt(DSQUARE(timevec[j]*slopex+interceptx-xvec[j]) + DSQUARE(timevec[j]*slopey+intercepty-yvec[j])));
	}
	// Ditch duplicate times, if there are any
	istimedup=1; // Guilty until proven innocent
	while(istimedup==1 && long(timevec.size())>=mintrkpts+1) {
	  istimedup=0;
	  j=1;
	  while(j<long(timevec.size()) && istimedup==0) {
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
	    if(timevec.size() != xvec.size() || timevec.size() != yvec.size() || timevec.size() != detindexvec.size()) {
	      cerr << "ERROR: vector length mismatch in vectors for tracklet-fitting!\n";
	      cerr << "Lengths of timevec, xvec, yvec, and detindexvec:\n";
	      cerr  << timevec.size()  << ", " << xvec.size()  << ", " << yvec.size()  << ", " << detindexvec.size() << "\n";
	      return(6);
	    }
	    // Re-do linear fit
	    // Perform fit to projected x coordinate as a function of time
	    linfituw01(timevec, xvec, slopex, interceptx);
	    // Perform fit to projected y coordinate as a function of time
	    linfituw01(timevec, yvec, slopey, intercepty);
	    // Load vector of residuals
	    fiterr = {};
	    for(j=0; j<long(timevec.size()); j++) {
	      fiterr.push_back(sqrt(DSQUARE(timevec[j]*slopex+interceptx-xvec[j]) + DSQUARE(timevec[j]*slopey+intercepty-yvec[j])));
	    }
	  }
	}
	// Find worst error.  
	worsterr = 0.0l;
	for(j=0; j<long(timevec.size()); j++) {
	  if(fiterr[j]>worsterr) {
	    worsterr = fiterr[j];
	    worstpoint = j;
	  }
	}
	// Reject successive points until either there are only three left
	// or the worst error drops below maxgcr.
	while(worsterr>maxgcr && timevec.size()>3 && long(timevec.size())>=mintrkpts) {
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
	  if(timevec.size() != xvec.size() || timevec.size() != yvec.size() || timevec.size() != detindexvec.size()) {
	    cerr << "ERROR: vector length mismatch in vectors for tracklet-fitting!\n";
	    cerr << "Lengths of timevec, xvec, yvec, and detindexvec:\n";
	    cerr  << timevec.size()  << ", " << xvec.size()  << ", " << yvec.size()  << ", " << detindexvec.size() << "\n";
	    return(6);
	  }
	  // Perform fit to projected x coordinate as a function of time
	  linfituw01(timevec, xvec, slopex, interceptx);
	  // Perform fit to projected y coordinate as a function of time
	  linfituw01(timevec, yvec, slopey, intercepty);
	  // Load vector of residuals
	  fiterr = {};
	  for(j=0; j<long(timevec.size()); j++) {
	    fiterr.push_back(sqrt(DSQUARE(timevec[j]*slopex+interceptx-xvec[j]) + DSQUARE(timevec[j]*slopey+intercepty-yvec[j])));
	  }
	  // Find worst error.  
	  worsterr = 0.0l;
	  if(fiterr.size() != timevec.size()) {
	    cerr << "Error: fiterr and timevec have different sizes: " << fiterr.size() << "vs. " << timevec.size() << "\n";
	    return(7);
	  }
	  for(j=0; j<long(timevec.size()); j++) {
	    if(fiterr[j]>worsterr) {
	      worsterr = fiterr[j];
	      worstpoint = j;
	    }
	  }
	}
	if(worsterr<=maxgcr && timevec.size()>=3 && long(timevec.size())>=mintrkpts) {
	  // We succeeded in finding a tracklet with no time-duplicates, and
	  // no outliers beyond maxgcr. Prepare to write it to the pair file.
	  // Select points that will represent this tracklet.
	  instep = (timevec.size()-1)/4;
	  rp1 = instep;
	  rp2 = timevec.size()-1-instep;
	  if(rp1==rp2) {
	    cerr << "ERROR: both representative points for a tracklet are the same!\n";
	    cerr << "size, instep, rp1, rp2: " << timevec.size() << " " << instep << " " << rp1 << " " << rp2 << "\n";
	    return(5);
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
	    track1 = tracklet(pairdets[detindexvec[rp1]].image,outra1,outdec1,pairdets[detindexvec[rp2]].image,outra2,outdec2,detindexvec.size(),tracklets.size());
	    tracklets.push_back(track1);
	    for(j=0; j<long(detindexvec.size()); j++) {
	      onepair = longpair(tracklets[tracklets.size()-1].trk_ID,detindexvec[j]);
	      indvecs[detindexvec[j]] = {};
	      trk2det.push_back(onepair);
	    }
	    istracklet=1;
	    // Close if-statement confirming that a bona fide,
	    // aligned tracklet was found and written to the output file.
	  } else {
	    istracklet=0;
	    if(verbose>=1) cout << "A tracklet was rejected: arc = " << setprecision(3) << fixed << dist << " < " << minarc << " or angvel = " << setprecision(5) << fixed << angvel << " not in range " << setprecision(3) << fixed << minvel << "-" << maxvel << "\n";
	  }
	} else istracklet=0;
	// Close else-statement confirming there was a candidate for
	// being an aligned tracklet.
      }
      // Close if-statement checking that detection i has more than
      // one pair-partner, and hence COULD be part of a tracklet
    } else istracklet=0;
    if((istracklet==0 || long(indvecs[pdct].size())>0) && mintrkpts==2) {
      // Either there was no tracklet (istracklet==0) or there was a tracklet,
      // but the original root point pdct got rejected from it.
      // In either case, it's necessary to write out the (surviving) pairs
      // associated with detection pdct.
      for(j=0; j<long(indvecs[pdct].size()); j++) {
	k=indvecs[pdct][j];
	// Calculate angular arc and angular velocity
	distradec02(pairdets[pdct].RA,pairdets[pdct].Dec,pairdets[k].RA,pairdets[k].Dec, &dist, &pa);
	angvel = dist/fabs(pairdets[pdct].MJD-pairdets[k].MJD); // Degrees per day
	dist *= 3600.0l; // Arcseconds
	if(indvecs[k].size()>0 && pairdets[k].MJD>pairdets[pdct].MJD && angvel>=minvel && dist>=minarc && angvel<=maxvel) {
	  track1 = tracklet(pairdets[pdct].image,pairdets[pdct].RA,pairdets[pdct].Dec,pairdets[k].image,pairdets[k].RA,pairdets[k].Dec,2,tracklets.size());
	  tracklets.push_back(track1);
	  onepair = longpair(tracklets.size()-1,pdct);
	  trk2det.push_back(onepair);
	  onepair = longpair(tracklets.size()-1,k);
	  trk2det.push_back(onepair);
	} else if(angvel<minvel || dist<minarc) {
	  if(verbose>=1) cout << "A pair was rejected: arc = " << setprecision(3) << fixed << dist << " < " << minarc << " or angvel = " << setprecision(5) << fixed << angvel << " not in range " << setprecision(3) << fixed << minvel << "-" << maxvel << "\n";
	}
      }
    }
    // Close loop over all detections
  }
  return(0);
}

// make_tracklets: April 07, 2023: dummy wrapper for make_tracklets

int make_tracklets(vector <hldet> &detvec, vector <hlimage> &image_log, MakeTrackletsConfig config, vector <hldet> &pairdets,vector <tracklet> &tracklets, vector <longpair> &trk2det)
{
 
  long unsigned int i=0;
  std::vector <longpair> pairvec;
  std::vector <vector <long>> indvecs;
  
  // Echo config struct
  cout << "Configuration parameters:\n";
  cout << "Min. number of tracklet points: " << config.mintrkpts << "\n";
  cout << "Time-tolerance for matching detections on the same image: " << config.imagetimetol << " days (" << config.imagetimetol*SOLARDAY << " seconds)\n";
  cout << "Maximum angular velocity: " << config.maxvel << " deg/day\n";
  cout << "Minimum angular velocity: " << config.minvel << " deg/day\n";
  cout << "Minimum angular arc: " << config.minarc << " arcsec\n";
  cout << "Maximum inter-image time interval: " << config.maxtime << " days (" << config.maxtime*1440.0 << " minutes)\n";
  cout << "Minimum inter-image time interval: " << config.mintime << " days (" << config.mintime*1440.0 << " minutes)\n";
  cout << "Image radius: " << config.imagerad << " degrees\n";
  cout << "Maximum Great Circle Residual for tracklets with more than two points: " << config.maxgcr << " arcsec\n";
  if(config.forcerun) {
    cout << "forcerun has been invoked: execution will attempt to push through\n";
    cout << "any errors that are not immediately fatal, including those that\n";
    cout << "could produce inaccurate final results.\n";
  }
  if(config.verbose) cout << "Verbose output has been requested\n";
  
  int status = load_image_indices(image_log, detvec, config.imagetimetol, config.forcerun);
  if(status!=0) {
    cerr << "ERROR: failed to load_image_indices from detection vector\n";
    return(status);
  }
  
  // Echo detection vector
  //for(i=0;i<detvec.size();i++) {
  //  cout << "det " << i << " " << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << " " << detvec[i].mag  << " " << detvec[i].obscode << " " << detvec[i].image << "\n";
  //}
  // Echo image log
  for(i=0;i<image_log.size();i++) {
    cout << "image " << i << " " << image_log[i].MJD << " " << image_log[i].RA << " " << image_log[i].Dec << " " << image_log[i].X << " " << image_log[i].obscode  << " " << image_log[i].startind  << " " << image_log[i].endind << "\n";
  }

  // Create pairs, output a vector pairdets of type hldet; a vector indvecs of type vector <long>,
  // with the same length as pairdets, giving the indices of all the detections paired with a given detection;
  // and the vector pairvec of type longpair, giving all the pairs of detections.
  status = find_pairs(detvec, image_log, pairdets, indvecs, pairvec, config.mintime, config.maxtime, config.imagerad, config.maxvel, config.verbose);
  
  if(status!=0) {
    cerr << "ERROR: find_pairs reports failure status " << status << "\n";
    return(status);
  }
  status = merge_pairs(pairdets, indvecs, pairvec, tracklets, trk2det, config.mintrkpts, config.maxgcr, config.minarc, config.minvel, config.maxvel, config.verbose);
  if(status!=0) {
    cerr << "ERROR: merge_pairs reports failure status " << status << "\n";
    return(status);
  } else cout << "merge_pairs finished OK\n";
  return(0);
}

int trk2statevec(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double heliodist, double heliovel, double helioacc, double chartimescale, vector <point6ix2> &allstatevecs, double mjdref, double mingeoobs, double minimpactpar)
{
  allstatevecs={};
  long imnum = image_log.size();
  long imct=0;
  long pairnum = tracklets.size();
  long pairct=0;
  int badpoint=0;
  int status1=0;
  int status2=0;
  int num_dist_solutions=0;
  int solnct=0;
  double mjdavg=0l;
  vector <double> heliodistvec;
  double delta1 = 0.0l;
  double RA,Dec;
  long i1,i2;
  i1=i2=0;
  point6dx2 statevec1 = point6dx2(0l,0l,0l,0l,0l,0l,0,0);
  point6ix2 stateveci = point6ix2(0,0,0,0,0,0,0,0);
  point3d observerpos1 = point3d(0l,0l,0l);
  point3d observerpos2 = point3d(0l,0l,0l);
  point3d targpos1 = point3d(0l,0l,0l);
  point3d targpos2 = point3d(0l,0l,0l);
  point3d targvel1 = point3d(0l,0l,0l);
  point3d targvel2 = point3d(0l,0l,0l);
  point3d unitbary = point3d(0l,0l,0l);
  vector <point3d> targposvec1;
  vector <point3d> targposvec2;
  int glob_warning=0;
  vector <double> deltavec1;
  vector <double> deltavec2;
  double absvelocity=0l;
  double impactpar=0l;
  double timediff=0l;
 
  // Calculate approximate heliocentric distances from the
  // input quadratic approximation.
  heliodistvec={};
  for(imct=0;imct<imnum;imct++) {
    delta1 = image_log[imct].MJD - mjdref;
      heliodistvec.push_back(heliodist + heliovel*delta1 + 0.5*helioacc*delta1*delta1);
      if(heliodistvec[imct]<=0.0l) {
	badpoint=1;
	return(1);
      }
  }
  if(badpoint==0 && long(heliodistvec.size())!=imnum) {
    cerr << "ERROR: number of heliocentric distance values does\n";
    cerr << "not match the number of input images!\n";
    return(2);
  }
  for(pairct=0; pairct<pairnum; pairct++) {
    badpoint=0;
    // Obtain indices to the image_log and heliocentric distance vectors.
    i1=tracklets[pairct].Img1;
    i2=tracklets[pairct].Img2;
    // Project the first point
    RA = tracklets[pairct].RA1;
    Dec = tracklets[pairct].Dec1;
    celestial_to_stateunit(RA,Dec,unitbary);
    observerpos1 = point3d(image_log[i1].X,image_log[i1].Y,image_log[i1].Z);
    targposvec1={};
    deltavec1={};
    status1 = helioproj02(unitbary,observerpos1, heliodistvec[i1], deltavec1, targposvec1);
    RA = tracklets[pairct].RA2;
    Dec = tracklets[pairct].Dec2;
    celestial_to_stateunit(RA,Dec,unitbary);
    observerpos2 = point3d(image_log[i2].X,image_log[i2].Y,image_log[i2].Z);
    targposvec2={};
    deltavec2={};
    status2 = helioproj02(unitbary, observerpos2, heliodistvec[i2], deltavec2, targposvec2);
    if(status1 > 0 && status2 > 0 && badpoint==0) {
      // Calculate time difference between the observations
      timediff = (image_log[i2].MJD - image_log[i1].MJD)*SOLARDAY;
      // Did helioproj find two solutions in both cases, or only one?
      num_dist_solutions = status1;
      if(num_dist_solutions > status2) num_dist_solutions = status2;
      // Loop over solutions (num_dist_solutions can only be 1 or 2).
      for(solnct=0; solnct<num_dist_solutions; solnct++) {
	// Begin new stuff added to eliminate 'globs'
	// These are spurious linkages of unreasonably large numbers (typically tens of thousands)
	// of detections that arise when the hypothetical heliocentric distance at a time when
	// many observations are acquired is extremely close to, but slightly greater than,
	// the heliocentric distance of the observer. Then detections over a large area of sky
	// wind up with projected 3-D positions in an extremely small volume -- and furthermore,
	// they all have similar velocities because the very small geocentric distance causes
	// the inferred velocities to be dominated by the observer's motion and the heliocentric
	// hypothesis, with only a negligible contribution from the on-sky angular velocity.
	glob_warning=0;
	if(deltavec1[solnct]<mingeoobs*AU_KM && deltavec2[solnct]<mingeoobs*AU_KM) {
	  // New-start
	  // Load target positions
	  targpos1 = targposvec1[solnct];
	  targpos2 = targposvec2[solnct];
	  // Calculate positions relative to observer
	  targpos1.x -= observerpos1.x;
	  targpos1.y -= observerpos1.y;
	  targpos1.z -= observerpos1.z;
	    
	  targpos2.x -= observerpos2.x;
	  targpos2.y -= observerpos2.y;
	  targpos2.z -= observerpos2.z;
	    
	  // Calculate velocity relative to observer
	  targvel1.x = (targpos2.x - targpos1.x)/timediff;
	  targvel1.y = (targpos2.y - targpos1.y)/timediff;
	  targvel1.z = (targpos2.z - targpos1.z)/timediff;
   
	  // Calculate impact parameter (past or future).
	  absvelocity = vecabs3d(targvel1);
	  impactpar = dotprod3d(targpos1,targvel1)/absvelocity;
	  // Effectively, we've projected targpos1 onto the velocity
	  // vector, and impactpar temporarily holds the magnitude of this projection.
	  // Subtract off the projection of the distance onto the velocity unit vector
	  targpos1.x -= impactpar*targvel1.x/absvelocity;
	  targpos1.y -= impactpar*targvel1.y/absvelocity;
	  targpos1.z -= impactpar*targvel1.z/absvelocity;
	  // Now targpos1 is the impact parameter vector at projected closest approach.
	  impactpar  = vecabs3d(targpos1); // Now impactpar is really the impact parameter
	  if(impactpar<=minimpactpar) {
	    // The hypothesis implies the object already passed with minimpactpar km of the Earth
	    // in the likely case that minimpactpar has been set to imply an actual impact,
	    // it's not our problem anymore.
	    glob_warning=1;
	  }
	}
	if(!glob_warning) {
	  targpos1 = targposvec1[solnct];
	  targpos2 = targposvec2[solnct];
	  
	  targvel1.x = (targpos2.x - targpos1.x)/timediff;
	  targvel1.y = (targpos2.y - targpos1.y)/timediff;
	  targvel1.z = (targpos2.z - targpos1.z)/timediff;

	  targpos1.x = 0.5L*targpos2.x + 0.5L*targpos1.x;
	  targpos1.y = 0.5L*targpos2.y + 0.5L*targpos1.y;
	  targpos1.z = 0.5L*targpos2.z + 0.5L*targpos1.z;
      
	  // Integrate orbit to the reference time.
	  mjdavg = 0.5l*image_log[i1].MJD + 0.5l*image_log[i2].MJD;
	  status1 = Keplerint(GMSUN_KM3_SEC2,mjdavg,targpos1,targvel1,mjdref,targpos2,targvel2);
	  if(status1 == 0 && badpoint==0) {
	    statevec1 = point6dx2(targpos2.x,targpos2.y,targpos2.z,chartimescale*targvel2.x,chartimescale*targvel2.y,chartimescale*targvel2.z,pairct,0);
	    // Note that the multiplication by chartimescale converts velocities in km/sec
	    // to units of km, for apples-to-apples comparison with the positions.
	    stateveci = conv_6d_to_6i(statevec1,INTEGERIZING_SCALEFAC);
	    allstatevecs.push_back(stateveci);
	  } else {
	    // Kepler integration encountered unphysical situation.
	    continue;
	  }
	}
      }
    } else {
      badpoint=1;
      // Heliocentric projection found no physical solution.
      continue;
    }
  }
  return(0);
}

// tracklet_lookup: Given a vector of type longpair that is a catalog
// of the form trk2det, find and return all of the entries corresponding
// to tracklet number trknum. The form of the input catalog is that it
// is monotonically sorted by the first index (e.g. tracklet count), but
// with the number of entries of the same tracklet count -- and hence,
// the vector index where a given tracklet trknum will start -- not known ahead of time. 
vector <long> tracklet_lookup(const vector <longpair> &trk2det, long trknum)
{
  vector <long> outvec;
  long catnum = trk2det.size();
  long i=catnum/2;
  long ilo=0;
  long ihi=catnum-1;
  int itnum=0;

  if(DEBUG>=2) cout << "Looking up tracklet " << trknum << "\n";
		
  while(itnum<BINSEARCHMAX && trk2det[i].i1 != trknum) {
    if(trk2det[i].i1 < trknum) {
      if(DEBUG>=2) cout << "Guess = " << i << " trknum = " << trk2det[i].i1 << ": too low\n";
      // Guess is too low. Make it the new lower bound
      ilo = i;
      // Reset to midway between the current low and high bounds
      i = (ilo+ihi)/2;
      itnum++;
      if(i<0) i=0;
      else if(i>=catnum) i=catnum-1;
    } else if(trk2det[i].i1 > trknum) {
      if(DEBUG>=2) cout << "Guess = " << i << " trknum = " << trk2det[i].i1 << ": too high\n";
      // Guess is too high. Make it the new upper bound
      ihi = i;
      // Reset to midway between the current low and high bounds
      i = (ilo+ihi)/2;
      itnum++;
      if(i<0) i=0;
      else if(i>=catnum) i=catnum-1;
    } 
  }
  if(trk2det[i].i1 != trknum) {
    cerr << "ERROR: lookup failed for tracklet number " << trknum <<"\n";
    outvec={};
    return(outvec);
  }
  // If we get here, we must have found the tracklet. Move upward to find where it begins.
  while(i>=0 && trk2det[i].i1 == trknum) i--;
  // Now point i+1 must be the start of the tracklet
  i+=1;
  if(DEBUG>0) cout << "Tracklet " << trknum << " begins at line number " << i << "\n";
  outvec={};
  while(i<catnum && trk2det[i].i1 == trknum) {
    if(DEBUG>0) cout << "Loading point " << outvec.size() << " of tracklet " << trknum << ", which corresponds to detection " << trk2det[i].i2 << "\n";
    outvec.push_back(trk2det[i].i2);
    i++;
  }
  return(outvec);
}

// earthpos01: March 28, 2023: wrapper to get an old-style 3D
// position for the Earth from a vector of the new EarthState struct.
point3d earthpos01(const vector <EarthState> &earthpos, double mjd)
{
  point3d earthnow = point3d(0,0,0);
  int polyorder=EPH_INTERP_POLYORDER;
  int status = planetpos01(mjd, polyorder, earthpos, earthnow); 
  if(status==0) return(earthnow);
  else {
    cerr << "ERROR: ephemeris interpolation code planetpos01,\n";
    cerr << "called by earthpos01, returned bad output\n";
     return(earthnow);
  }
}


int form_clusters(const vector <point6ix2> &allstatevecs, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const point3d &Earthrefpos, double heliodist, double heliovel, double helioacc, double chartimescale, vector <hlclust> &outclust, vector <longpair> &clust2det, long &realclusternum, double cluster_radius, double dbscan_npt, double mingeodist, double geologstep, double maxgeodist, int mintimespan, int minobsnights, int verbose)
{
  int geobinct = 0;
  long i=0;
  long statevecnum = allstatevecs.size();
  double georadcen = mingeodist*intpowD(geologstep,geobinct);
  double georadmin=0l;
  double georadmax=0l;
  vector <point6ix2> binstatevecs;
  point6dx2 statevec1 = point6dx2(0l,0l,0l,0l,0l,0l,0,0);
  point6ix2 stateveci = point6ix2(0,0,0,0,0,0,0,0);
  double geodist=0l;
  long kdroot=0;
  long splitpoint=0;
  int gridpoint_clusternum=0;
  KD_point6ix2 kdpoint = KD_point6ix2(stateveci,-1,-1,1,0);
  vector <KD_point6ix2> kdvec;
  vector <KD6i_clust> kdclust;
  vector <long> pointind;
  vector <long> pointjunk;
  vector <double> clustmjd;
  vector <double> mjdstep;
  double timespan=0;
  int daysteps = 0;
  int obsnights = 0;
  longpair c2d = longpair(0,0);
  hlclust onecluster = hlclust(0, 0.0l, 0.0l, 0.0l, 0.0l, 0, 0.0l, 0, 0, 0.0l, "NULL", 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0);
  double posRMS = 0.0l;
  double velRMS = 0.0l;
  double totRMS = 0.0l;
  double astromRMS = 0.0l;
  int pairnum = 0;
  int uniquepoints = 0;
  double posX, posY, posZ, velX, velY, velZ;
  posX = posY = posZ = velX = velY = velZ = 0.0l;
  double orbit_a, orbit_e, orbit_MJD, orbitX, orbitY, orbitZ, orbitVX, orbitVY, orbitVZ;
  orbit_a = orbit_e = orbit_MJD = orbitX = orbitY = orbitZ = orbitVX = orbitVY = orbitVZ = 0.0l;
  long orbit_eval_count = 0;
  int clusterct=0;
  double clustmetric=0.0l;
  string rating;
  int pairct=0;
  int j=0;

  // Loop over geocentric bins, selecting the subset of state-vectors
  // in each bin, and running DBSCAN only on those, with clustering radius
  // adjusted accordingly.
  while(georadcen<=maxgeodist) {
    georadcen = mingeodist*intpowD(geologstep,geobinct);
    georadmin = georadcen/geologstep;
    georadmax = georadcen*geologstep;
    // Load new array of state vectors, limited to those in the current geocentric bin
    binstatevecs={};
    for(i=0; i<statevecnum; i++) {
      // Reverse integerization of the state vector.
      // This is only possible to a crude approximation, of course.
      statevec1 = conv_6i_to_6d(allstatevecs[i],INTEGERIZING_SCALEFAC);
      // Calculate geocentric distance in AU
      geodist = sqrt(DSQUARE(statevec1.x-Earthrefpos.x) + DSQUARE(statevec1.y-Earthrefpos.y) + DSQUARE(statevec1.z-Earthrefpos.z))/AU_KM;
      if(geodist >= georadmin && geodist <= georadmax) {
	// This state vector is in the geocentric radius bin we are currently considering.
	// Add it to binstatevecs.
	binstatevecs.push_back(allstatevecs[i]);
      }
    }
    if(verbose>=1) cout << "Found " << binstatevecs.size() << " state vectors in geocentric bin from " << georadmin << " to " << georadmax << " AU\n";
    if(binstatevecs.size()<=1) {
      geobinct++;
      continue; // No clusters possible, skip to the next step.
    }
      
    kdvec={};
    kdroot = splitpoint = 0;
    splitpoint=medind_6ix2(binstatevecs,1);
    kdpoint = KD_point6ix2(binstatevecs[splitpoint],-1,-1,1,0);
    kdvec.push_back(kdpoint);
    kdtree_6i01(binstatevecs,1,splitpoint,kdroot,kdvec);
    
    if(verbose>=1) cout << "Created a KD tree with " << kdvec.size() << " branches\n";

    kdclust={};
    int clusternum = DBSCAN_6i01(kdvec, cluster_radius*(georadcen/REF_GEODIST)/INTEGERIZING_SCALEFAC, dbscan_npt, INTEGERIZING_SCALEFAC, kdclust, verbose);
    if(verbose>=1) cout << "DBSCAN_6i01 finished, with " << clusternum << " = " << kdclust.size() << " clusters found\n";
    for(clusterct=0; clusterct<long(kdclust.size()); clusterct++) {
      // Scale cluster RMS down to reference geocentric distance
      if(DEBUG >= 2) cout << "scaling kdclust rms for cluster " << clusterct << " out of " << kdclust.size() << "\n";
      fflush(stdout);
      for(i=0;i<9;i++) {
	if(DEBUG >= 2) cout << "scaling rmsvec point " << i << " out of " << kdclust[clusterct].rmsvec.size() << "\n";
	if(DEBUG >= 2) cout << "RMS = " << kdclust[clusterct].rmsvec[i];
	kdclust[clusterct].rmsvec[i] *= REF_GEODIST/georadcen;
	if(DEBUG >= 2) cout << ", scales to " << kdclust[clusterct].rmsvec[i] << "\n";
      }
      // Note that RMS is scaled down for more distant clusters, to
      // avoid bias against them in post-processing.
	
      // Map cluster to individual detections.
      // create vector of unique detection indices.
      if(DEBUG >= 1) cout << "Loading pointind for " << kdclust[clusterct].numpoints << " of cluster #" << clusterct <<  " out of " << kdclust.size() << "\n";
      fflush(stdout);
      pointind = {};
      for(i=0;i<kdclust[clusterct].numpoints;i++) {
	pairct=kdvec[kdclust[clusterct].clustind[i]].point.i1;
	if(DEBUG >= 2) cout << "Looking up tracklet " << pairct << " out of " << tracklets.size() << "\n";
	pointjunk = {};
	pointjunk = tracklet_lookup(trk2det, pairct);
	if(DEBUG >= 2) cout << "Found " << pointjunk.size() << " detections for tracklet " << pairct << "\n";
	for(j=0; j<long(pointjunk.size()); j++) {
	  pointind.push_back(pointjunk[j]);
	}
      }
      // Sort vector of detection indices
      sort(pointind.begin(), pointind.end());
      // Cull out duplicate entries
      pointjunk = pointind;
      pointind={};
      pointind.push_back(pointjunk[0]);
      for(i=1; i<long(pointjunk.size()); i++) {
	if(pointjunk[i]!=pointjunk[i-1]) pointind.push_back(pointjunk[i]);
      }
      uniquepoints = pointind.size();

      // Load vector of detection MJD's
      clustmjd = {};
      for(i=0; i<long(pointind.size()); i++) {
	clustmjd.push_back(detvec[pointind[i]].MJD);
      }

      // Sort vector of MJD's
      sort(clustmjd.begin(), clustmjd.end());
      timespan = clustmjd[clustmjd.size()-1] - clustmjd[0];
      // Load vector of MJD steps
      mjdstep={};
      for(i=1; i<long(clustmjd.size()); i++) {
	mjdstep.push_back(clustmjd[i]-clustmjd[i-1]);
      }
      // Count steps large enough to suggest a daytime period between nights.
      daysteps=0;	
      for(i=0; i<long(mjdstep.size()); i++) {
	if(mjdstep[i]>NIGHTSTEP) daysteps++;
      }
      obsnights = daysteps+1;
      // Does cluster pass the criteria for a linked detection?
      if(timespan >= mintimespan && obsnights >= minobsnights) {
	if(verbose >= 1) cout << "Loading good cluster " << realclusternum << " with timespan " << timespan << " and obsnights " << obsnights << "\n";
	fflush(stdout);
	if(verbose >= 1) cout << "Cluster passes discovery criteria: will be designated as cluster " << realclusternum << "\n";
	// Check whether cluster is composed purely of detections from
	// a single simulated object (i.e., would be a real discovery) or is a mixture
	// of detections from two or more different simulated objects (i.e., spurious).
	rating="PURE";
	for(i=0; i<long(pointind.size()); i++) {
	  if(i>0 && stringnmatch01(detvec[pointind[i]].idstring,detvec[pointind[i-1]].idstring,SHORTSTRINGLEN)!=0) rating="MIXED";
	}
	if(DEBUG >= 1) cout << "Rating is found to be " << rating << "\n";
	fflush(stdout);
	// Write all individual detections in this cluster to the clust2det array
	for(i=0; i<long(pointind.size()); i++) {
	  c2d = longpair(realclusternum,pointind[i]);
	  clust2det.push_back(c2d);
	}

	// Calculate values for the statistics in the output array (class hlclust) that have
	// not been caculated already.
	clustmetric = double(pointind.size())*double(obsnights)*timespan/kdclust[clusterct].rmsvec[8];
	// Note contents of rmsvec: [0] xrms, [1] yrms, [2] zrms, [3] vxrms, [4] vyrms, [5] vzrms,
	// [6] overall position rms, [7] overall velocity rms, [8] overall rms
	posRMS = kdclust[clusterct].rmsvec[6];
	velRMS = kdclust[clusterct].rmsvec[7];
	totRMS = kdclust[clusterct].rmsvec[8];
	pairnum = kdclust[clusterct].numpoints; // This is the original total number of pairs/tracklets assigned
	                                        // to this cluster, which might have a lot of overlap in terms of 
	                                        // individual detections (of which the non-overlapping count has
	                                        // already been saved in 'uniquepoints').
	// Now save the values of the mean state vectors at the reference time.
	posX = kdclust[clusterct].meanvec[0];
	posY = kdclust[clusterct].meanvec[1];
	posZ = kdclust[clusterct].meanvec[2];
	velX = kdclust[clusterct].meanvec[3]/chartimescale;
	velY = kdclust[clusterct].meanvec[4]/chartimescale;
	velZ = kdclust[clusterct].meanvec[5]/chartimescale;
	// Some of the statistics in the hlclust class relate to orbit-fitting,
	// and are meant for later use. For now, set them all to zero.
	astromRMS = orbit_a = orbit_e = 0.0l;
	orbit_MJD = orbitX = orbitY = orbitZ = orbitVX = orbitVY = orbitVZ = 0.0l;
	orbit_eval_count = 0;
	// Write overall cluster statistics to the outclust array.	
	onecluster = hlclust(realclusternum, posRMS, velRMS, totRMS, astromRMS, pairnum, timespan, uniquepoints, obsnights, clustmetric, rating, heliodist/AU_KM, heliovel/SOLARDAY, helioacc*1000.0/SOLARDAY/SOLARDAY, posX, posY, posZ, velX, velY, velZ, orbit_a, orbit_e, orbit_MJD, orbitX, orbitY, orbitZ, orbitVX, orbitVY, orbitVZ, orbit_eval_count);
	// cout << "kdload velrms: " << velRMS << " " << kdclust[clusterct].rmsvec[7] << " " << onecluster.velRMS << "\n";
	outclust.push_back(onecluster);
	realclusternum++;
	gridpoint_clusternum++;
      }
    }
    // Move on to the next bin in geocentric distance
    geobinct++;
  }
  if(verbose>=0) cout << "Identified " << gridpoint_clusternum << " candidate linkages\n";
  return(0);
}

// heliolinc_alg: April 11, 2023: dummy wrapper for heliolinc,
// calls all important algorithmic stuff.

int heliolinc_alg(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det)
{
  outclust = {};
  clust2det = {};
   
  point3d Earthrefpos = point3d(0l,0l,0l);
  long imnum = image_log.size();
  long pairnum = tracklets.size();
  long trk2detnum = trk2det.size();
  long accelnum = radhyp.size();
  long accelct=0;

  vector <double> heliodist;
  vector <double> heliovel;
  vector <double> helioacc;
  long realclusternum, gridpoint_clusternum, status;
  realclusternum = gridpoint_clusternum = status = 0;
  vector <point6ix2> allstatevecs;

  // Echo config struct
  cout << "Configuration parameters:\n";
  cout << "MJD of reference time: " << config.MJDref << "\n";
  cout << "DBSCAN clustering radius: " << config.clustrad << " km\n";
  cout << "DBSCAN npt: " << config.dbscan_npt << "\n";
  cout << "Min number of distinct observing nights for a valid linkage: " << config.minobsnights << "\n";
  cout << "Min time span for a valid linkage: " << config.mintimespan << " days\n";
  cout << "Min geocentric distance (center of innermost bin): " << config.mingeodist << " AU\n";
  cout << "Max geocentric distance (will be exceeded by center only of the outermost bin): " << config.maxgeodist << " AU\n";
  cout << "Logarthmic step size (and bin width) for geocentric distance bins: " << config.geologstep << "\n";
  cout << "Minimum inferred geocentric distance for a valid tracklet: " << config.mingeoobs << " AU\n";
  cout << "Minimum inferred impact parameter (w.r.t. Earth) for a valid tracklet: " << config.minimpactpar << " Earth radii\n";
  if(config.verbose) cout << "Verbose output selected\n";
  
  if(imnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty image catalog\n";
    return(1);
  } else if(pairnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty tracklet array\n";
    return(1);
  } else if(trk2detnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty trk2det array\n";
    return(1);
  } else if(accelnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty heliocentric hypothesis array\n";
    return(1);
  }
  
  double MJDmin = image_log[0].MJD;
  double MJDmax = image_log[imnum-1].MJD;
  if(config.MJDref<MJDmin || config.MJDref>MJDmax) {
    // Input reference MJD is invalid. Suggest a better value before exiting.
    cerr << "\nERROR: reference MJD, supplied as " << config.MJDref << ",\n";
    cerr << "must fall in the time interval spanned by the data (" << MJDmin << " to " << MJDmax << "\n";
    cerr << fixed << setprecision(2) << "Suggested value is " << MJDmin*0.5l + MJDmax*0.5l << "\n";
    cout << "based on your input image catalog\n";
    return(2);
  }

  double chartimescale = (MJDmax - MJDmin)*SOLARDAY/TIMECONVSCALE; // Note that the units are seconds.
  Earthrefpos = earthpos01(earthpos, config.MJDref);

  // Convert heliocentric radial motion hypothesis matrix
  // from units of AU, AU/day, and GMSun/R^2
  // to units of km, km/day, and km/day^2.
  heliodist = heliovel = helioacc = {};
  for(accelct=0; accelct<accelnum; accelct++) {
    heliodist.push_back(radhyp[accelct].HelioRad * AU_KM);
    heliovel.push_back(radhyp[accelct].R_dot * AU_KM);
    helioacc.push_back(radhyp[accelct].R_dubdot * (-GMSUN_KM3_SEC2*SOLARDAY*SOLARDAY/heliodist[accelct]/heliodist[accelct]));
  }

  // Begin master loop over heliocentric hypotheses
  outclust={};
  clust2det={};
  realclusternum=0;  
  for(accelct=0;accelct<accelnum;accelct++) {
    gridpoint_clusternum=0;
    
    // Covert all tracklets into state vectors at the reference time, under
    // the assumption that the heliocentric distance hypothesis is correct.
    status = trk2statevec(image_log, tracklets, heliodist[accelct], heliovel[accelct], helioacc[accelct], chartimescale, allstatevecs, config.MJDref, config.mingeoobs, config.minimpactpar);
    
    if(status==1) {
      cerr << "WARNING: hypothesis " << accelct << ": " << radhyp[accelct].HelioRad << " " << radhyp[accelct].R_dot << " " << radhyp[accelct].R_dubdot << " led to\nnegative heliocentric distance or other invalid result: SKIPPING\n";
      continue;
    } else if(status==2) {
      // This is a weirder error case and is fatal.
      cerr << "Fatal error case from trk2statevec.\n";
      return(3);
    }
    // If we get here, trk2statevec probably ran OK.
    if(allstatevecs.size()<=1) continue; // No clusters possible, skip to the next step.
    if(config.verbose>=0) cout << pairnum << " input pairs/tracklets led to " << allstatevecs.size() << " physically reasonable state vectors\n";

    status = form_clusters(allstatevecs, detvec, tracklets, trk2det, Earthrefpos, heliodist[accelct], heliovel[accelct], helioacc[accelct], chartimescale, outclust, clust2det, realclusternum, config.clustrad, config.dbscan_npt, config.mingeodist, config.geologstep, config.maxgeodist, config.mintimespan, config.minobsnights, config.verbose);
  }
  return(0);    
}

// link_refine_Herget: April 11, 2023:
// Algorithmic portion to be called by wrappers.
int link_refine_Herget(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <hlclust> &inclust, const vector  <longpair> &inclust2det, LinkRefineConfig config, vector <hlclust> &outclust, vector <longpair> &outclust2det)
{
  long i=0;
  long imnum = image_log.size();
  long imct=0;
  long detnum = detvec.size();
  long inclustnum = inclust.size();
  long inclustct=0;
  long clusternum=0;
  double clustmetric = 0.0l;
  hlclust onecluster = hlclust(0, 0.0l, 0.0l, 0.0l, 0.0l, 0, 0.0l, 0, 0, 0.0l, "NULL", 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0);
  vector <long> clustind;
  vector <hldet> clusterdets;
  vector <hlclust> holdclust;
  int ptnum,ptct,istimedup,detsused;
  ptnum=ptct=istimedup=detsused=0;
  point3d onepoint = point3d(0.0L,0.0L,0.0L);
  vector <point3d> observerpos;
  vector <double> obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit;
  double geodist1,geodist2, v_escape, v_helio, astromrms, chisq;
  double ftol = FTOL_HERGET_SIMPLEX;
  double simplex_scale = SIMPLEX_SCALEFAC;
  double X, Y, Z, dt;
  X = Y = Z = dt = 0l;
  point3d startpos = point3d(0.0l,0.0l,0.0l);
  point3d startvel = point3d(0.0l,0.0l,0.0l);
  point3d endpos = point3d(0.0l,0.0l,0.0l);
  point3d endvel = point3d(0.0l,0.0l,0.0l);
  vector <int> detusedvec = {};
  vector <vector <long>> clustindmat = {};
  vector <double_index> metric_index;
  long clusterct,clusterct2,goodclusternum;
  clusterct = clusterct2 = goodclusternum = 0;
  double_index dindex = double_index(0l,0);
  vector <double_index> sortclust;
  longpair c2d = longpair(0,0);
  char rating[SHORTSTRINGLEN];
 
  make_ivec(detnum, detusedvec); // All the entries are guaranteed to be 0.

  // Wipe output arrays
  holdclust={};
  clustindmat={};
  outclust={};
  outclust2det={};
 
  cout << "Reference MJD: " << config.MJDref << "\n";
  cout << "Maximum RMS in km: " << config.maxrms << "\n";
  cout << "In calculating the cluster quality metric, the number of\nunique points will be raised to the power of " << config.ptpow << ";\n";
  cout << "the number of unique nights will be raised to the power of " << config.nightpow << ";\n";
  cout << "the total timespan will be raised to the power of " << config.timepow << ";\n";
  cout << "and the astrometric RMS will be raised to the power of (negative) " << config.rmspow << "\n";
  if(config.verbose>=1) cout << "verbose output has been selected\n";

  // Launch master loop over all the input clusters.
  for(inclustct=0; inclustct<inclustnum; inclustct++) {
    onecluster = inclust[inclustct];
    if(inclustct!=onecluster.clusternum) {
      cerr << "ERROR: cluster index mismatch " << inclustct << " != " << onecluster.clusternum << " at input cluster " << inclustct << "\n";
      return(5);
    }
    if(onecluster.totRMS<=config.maxrms) {
      // This cluster passes the initial cut. Analyze it.
      // Load a vector with the indices to detvec
      clustind = {};
      clustind = tracklet_lookup(inclust2det, inclustct);
      ptnum = clustind.size();
      if(ptnum!=onecluster.uniquepoints) {
	cerr << "ERROR: point number mismatch " << ptnum << " != " << onecluster.uniquepoints << " at input cluster " << inclustct << "\n";
	return(6);
      }
      // Load vector of detections for this cluster
      clusterdets={};
      for(i=0; i<ptnum; i++) {
	clusterdets.push_back(detvec[clustind[i]]);
      }
      sort(clusterdets.begin(), clusterdets.end(), early_hldet());
      istimedup=0;
      for(ptct=1; ptct<ptnum; ptct++) {
	if(clusterdets[ptct-1].MJD == clusterdets[ptct].MJD && stringnmatch01(clusterdets[ptct-1].obscode,clusterdets[ptct].obscode,3)==0) istimedup=1;
      }
      if(istimedup==0) {
	// The cluster is good so far.
	// Recalculate clustermetric
	clustmetric = intpowD(double(onecluster.uniquepoints),config.ptpow)*intpowD(double(onecluster.obsnights),config.nightpow)*intpowD(onecluster.timespan,config.timepow);
	// Note that the value of clustermetric just calculated
	// will later be divided by the reduced chi-square value of the
	// astrometric fit, before it is ultimately used as a selection criterion.
	    
	// Perform orbit fitting using the method of Herget, to get astrometric residuals
	// Load observational vectors
	observerpos = {};
	obsMJD = obsRA = obsDec = sigastrom = fitRA = fitDec = resid = orbit = {};
	for(ptct=0; ptct<ptnum; ptct++) {
	  obsMJD.push_back(clusterdets[ptct].MJD);
	  obsRA.push_back(clusterdets[ptct].RA);
	  obsDec.push_back(clusterdets[ptct].Dec);
	  sigastrom.push_back(1.0L); // WARNING, THIS IS CRUDE AND NEEDS FIXING
	  imct = clusterdets[ptct].image;
	  if(imct>=imnum) {
	    cerr << "ERROR: attempting to access image " << imct << " of only " << imnum << " available\n";
	    return(8);
	  }
	  X = image_log[imct].X;
	  Y = image_log[imct].Y;
	  Z = image_log[imct].Z;
	  if(image_log[imct].MJD!=clusterdets[ptct].MJD) {
	    // A shutter-travel correction must have been applied to the
	    // detection time relative to the image time. Use the stored
	    // velocity info from the image log to apply a correction to
	    // the observer position.
	    dt = clusterdets[ptct].MJD - image_log[imct].MJD;
	    if(dt*SOLARDAY > MAX_SHUTTER_CORR) {
	      cerr << "ERROR: detection vs. image time mismatch of " << dt*SOLARDAY << " seconds.\n";
	      cerr << "Something has gone wrong: no shutter is that slow\n";
	      return(4);
	    }
	    X += image_log[imct].VX*dt;
	    Y += image_log[imct].VY*dt;
	    Z += image_log[imct].VZ*dt;
	  }
	  onepoint = point3d(X,Y,Z);
	  observerpos.push_back(onepoint);
	}
	// Use mean state vectors to estimate positions
	startpos.x = onecluster.posX;
	startpos.y = onecluster.posY;
	startpos.z = onecluster.posZ;
	startvel.x = onecluster.velX;
	startvel.y = onecluster.velY;
	startvel.z = onecluster.velZ;
	    
	// Check if the velocity is above escape
	v_escape = sqrt(2.0L*GMSUN_KM3_SEC2/vecabs3d(startpos));
	v_helio = vecabs3d(startvel);
	if(v_helio>=v_escape) {
	  cerr << "WARNING: mean state vector velocity was " << v_helio/v_escape << " times higher than solar escape\n";
	  startvel.x *= ESCAPE_SCALE*v_escape/v_helio;
	  startvel.y *= ESCAPE_SCALE*v_escape/v_helio;
	  startvel.z *= ESCAPE_SCALE*v_escape/v_helio;
	}
	// Calculate position at first observation
	Keplerint(GMSUN_KM3_SEC2, config.MJDref, startpos, startvel, obsMJD[0], endpos, endvel);
	// Find vector relative to the observer by subtracting off the observer's position.
	endpos.x -= observerpos[0].x;
	endpos.y -= observerpos[0].y;
	endpos.z -= observerpos[0].z;
	geodist1 = vecabs3d(endpos)/AU_KM;
	// Calculate position at last observation
	Keplerint(GMSUN_KM3_SEC2, config.MJDref, startpos, startvel, obsMJD[ptnum-1], endpos, endvel);
	endpos.x -= observerpos[ptnum-1].x;
	endpos.y -= observerpos[ptnum-1].y;
	endpos.z -= observerpos[ptnum-1].z;
	geodist2 = vecabs3d(endpos)/AU_KM;
	simplex_scale = SIMPLEX_SCALEFAC;
	if(config.verbose>=2) {
	  cout << "Calling Hergetfit01 with dists " << geodist1 << " and " << geodist2 << "\n";
	}
	if(config.verbose>=2) {
	  cout << "Launching Hergetfit01 for cluster " << inclustct << ":\n";
	  for(i=0;i<=ptnum;i++) {
	    cout << "Point " << i << ": " << obsMJD[i] << " " << obsRA[i] << " "  << obsDec[i] << " " << sigastrom[i] << "\n";
	  }
	}
	if(config.verbose>=1) cout << "Cluster " << inclustct << " of " << inclustnum << " is good: ";
	if(config.verbose>=2) cout << "\n";
	chisq = Hergetfit01(geodist1, geodist2, simplex_scale, config.simptype, ftol, 1, ptnum, observerpos, obsMJD, obsRA, obsDec, sigastrom, fitRA, fitDec, resid, orbit, config.verbose);
	// orbit vector contains: semimajor axis [0], eccentricity [1],
	// mjd at epoch [2], the state vectors [3-8], and the number of
	// orbit evaluations (~iterations) required to reach convergence [9].
      
	chisq /= double(ptnum); // Now it's the reduced chi square value
	astromrms = sqrt(chisq); // This gives the actual astrometric RMS in arcseconds if all the
	// entries in sigastrom are 1.0. Otherwise it's a measure of the
	// RMS in units of the typical uncertainty.
	// Include this astrometric RMS value in the cluster metric and the RMS vector
	onecluster.astromRMS = astromrms; // rmsvec[3]: astrometric rms in arcsec.
	onecluster.metric = clustmetric/intpowD(astromrms,config.rmspow); // Under the default value rmspow=2, this is equivalent
	                                                                  // to dividing by the chi-square value rather than just
	                                                                  // the astrometric RMS, which has the desireable effect of
	                                                                  // prioritizing low astrometric error even more.
	onecluster.orbit_a = orbit[0]/AU_KM;
	onecluster.orbit_e = orbit[1];
	onecluster.orbit_MJD = orbit[2];
	onecluster.orbitX = orbit[3];
	onecluster.orbitY = orbit[4];
	onecluster.orbitZ = orbit[5];
	onecluster.orbitVX = orbit[6];
	onecluster.orbitVY = orbit[7];
	onecluster.orbitVZ = orbit[8];
	onecluster.orbit_eval_count = long(round(orbit[9]));
	// Push new cluster on to holding vector holdclust
	holdclust.push_back(onecluster);
	clustindmat.push_back(clustind);
	// Close if-statement checking for duplicate MJDs
      }
      // Close if-statement checking the RMS was low enough.
    }	
    // Close loop on input cluster arrays.
  }
  clusternum = holdclust.size();
  cout << "Finished loading input clusters: " << clusternum << " out of " << inclustnum << " passed initial screening.\n";

  // Load just clustermetric values and indices from
  // clustanvec into metric_index
  metric_index = {};
  // Record indices so information won't be lost on sort
  for(clusterct=0; clusterct<clusternum; clusterct++) {
    dindex = double_index(holdclust[clusterct].metric,clusterct);
    metric_index.push_back(dindex);
  }
  // Sort metric_index
  sort(metric_index.begin(), metric_index.end(), lower_double_index());
  
  // Loop on clusters, starting with the best (highest metric),
  // and eliminating duplicates
  goodclusternum=0;
  for(clusterct2=clusternum-1; clusterct2>=0; clusterct2--) {
    // Look up the correct cluster from metric_index
    clusterct = metric_index[clusterct2].index;
    inclustct = holdclust[clusterct].clusternum;
    onecluster = holdclust[clusterct];
    // Pull out the vector of detection indices from clustindmat
    clustind = clustindmat[clusterct];
    ptnum = clustind.size();
    // Sanity-check the count of unique detections in this cluster
    if(onecluster.uniquepoints>0 && ptnum!=onecluster.uniquepoints) {
      cerr << "ERROR: 2nd-stage point number mismatch " << ptnum << " != " << onecluster.uniquepoints << " at input cluster " << inclustct << " (" << clusterct << ")\n";
      return(7);
    }
    // See if all of them are still unused
    detsused = 0;
    for(ptct=1; ptct<ptnum; ptct++) {
      if(detusedvec[clustind[ptct]]!=0) detsused+=1;
    }
    if(onecluster.uniquepoints>0 && onecluster.totRMS<=config.maxrms && detsused==0) {
      // This is a good cluster not already marked as used.
      goodclusternum++;
      cout << "Accepted good cluster " << goodclusternum << " with metric " << onecluster.metric << "\n";
      // See whether cluster is pure or mixed.
      stringncopy01(rating,"PURE",SHORTSTRINGLEN);
      for(ptct=1; ptct<ptnum; ptct++) {
	if(stringnmatch01(detvec[clustind[ptct]].idstring,detvec[clustind[ptct-1]].idstring,SHORTSTRINGLEN) != 0) {
	  stringncopy01(rating,"MIXED",SHORTSTRINGLEN);
	}
      }
      // Figure out the time order of cluster points, so we can write them out in order.
      sortclust = {};
      for(ptct=0; ptct<ptnum; ptct++) {
	dindex = double_index(detvec[clustind[ptct]].MJD,clustind[ptct]);
	sortclust.push_back(dindex);
	// Also mark each detection as used
	detusedvec[clustind[ptct]]=1;
      }
      sort(sortclust.begin(), sortclust.end(), lower_double_index());
      // Load new clusternumber for onecluster
      onecluster.clusternum = outclust.size();
      // Push back the new, vetted cluster on to outclust
      outclust.push_back(onecluster);
      // Push back the detection indices 
      for(ptct=0; ptct<ptnum; ptct++) {
	c2d = longpair(onecluster.clusternum,sortclust[ptct].index);
	outclust2det.push_back(c2d);
      }
      // Close if-statement checking if this is a good, non-duplicated cluster
    }
    // Close loop on all clusters passing the initial cuts
  }
  return(0);
}


#undef DEBUGB
#undef DEBUG
#undef INTEGERIZING_SCALEFAC
#undef BINSEARCHMAX
#undef REF_GEODIST
#undef NIGHTSTEP
#undef EPH_INTERP_POLYORDER
#undef TIMECONVSCALE
#undef TIMECONVSCALE
#undef FTOL_HERGET_SIMPLEX
#undef ESCAPE_SCALE

// parse_clust2det: April 26, 2023:
// Use the clust2det mapping vector produced by heliolinc or
// link_refine to parse out clusters into a cluster detection
// vector -- that is, a vector that includes detection records
// only for the valid clusters. In the output hldet vector,
// 'index' gets replaced with the cluster number.
int parse_clust2det(const vector <hldet> &detvec, const vector <longpair> &inclust2det, vector <hldet> &clustdet)
{
  long detnum = detvec.size();
  long c2dnum = inclust2det.size();
  long detct=0;
  long c2dct=0;
  long clustct=0;
  clustdet = {};
  hldet o1 = hldet(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, "null", "V", "500", 0, 0, 0);
  
  for(c2dct=0; c2dct<c2dnum; c2dct++) {
    clustct = inclust2det[c2dct].i1;
    detct = inclust2det[c2dct].i2;
    if(detct>=0 && detct<detnum) {
      o1 = detvec[detct];
      o1.index = clustct;
      clustdet.push_back(o1);
    } else {
      cerr << "Warning: parse_clust2det tried to access out-of-range detection\n";
      cerr << "number " << detct << ", while valid range is 0-" << detnum << "\n";
    }
  }
  return(0);
}
