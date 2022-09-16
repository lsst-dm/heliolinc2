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

void make_ivec(int nx, vector <int> &ivec)
{
  int i=0;
  ivec={};
  for(i=0;i<=nx;i++) ivec.push_back(0);
}

void make_imat(int nx, int ny, vector <vector <int>> &imat)
{
  int i=0;
  int j=0;
  vector <int> tvec;
  imat = {};
  
  for(i=0;i<=nx;i++) {
    tvec={};
    for(j=0;j<=ny;j++) tvec.push_back(0);
    imat.push_back(tvec);
  }
}

void make_lvec(int nx, vector <long> &lvec)
{
  int i=0;
  lvec={};
  for(i=0;i<=nx;i++) lvec.push_back(0);
}

void make_lmat(int nx, int ny, vector <vector <long>> &lmat)
{
  int i=0;
  int j=0;
  vector <long> tvec;
  lmat = {};
  
  for(i=0;i<=nx;i++) {
    tvec={};
    for(j=0;j<=ny;j++) tvec.push_back(0);
    lmat.push_back(tvec);
  }
}

void make_cvec(int nx, vector <char> &cvec)
{
  int i=0;
  cvec={};
  for(i=0;i<=nx;i++) cvec.push_back('\0');
}

void make_cmat(int nx, int ny, vector <vector <char>> &cmat)
{
  int i=0;
  int j=0;
  vector <char> tvec;
  cmat = {};
  
  for(i=0;i<=nx;i++) {
    tvec={};
    for(j=0;j<=ny;j++) tvec.push_back('\0');
    cmat.push_back(tvec);
  }
}

void make_dvec(int nx, vector <double> &dvec)
{
  int i=0;
  dvec={};
  for(i=0;i<=nx;i++) dvec.push_back(0.0);
}

void make_dmat(int nx, int ny, vector <vector <double>> &dmat)
{
  int i=0;
  int j=0;
  vector <double> tvec;
  dmat = {};
  
  for(i=0;i<=nx;i++) {
    tvec={};
    for(j=0;j<=ny;j++) tvec.push_back(0.0);
    dmat.push_back(tvec);
  }
}

void make_LDvec(int nx, vector <long double> &ldvec)
{
  int i=0;
  ldvec={};
  for(i=0;i<=nx;i++) ldvec.push_back(0.0);
}

void make_LDmat(int nx, int ny, vector <vector <long double>> &ldmat)
{
  int i=0;
  int j=0;
  vector <long double> tvec;
  ldmat = {};
  
  for(i=0;i<=nx;i++) {
    tvec={};
    for(j=0;j<=ny;j++) tvec.push_back(0.0);
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
  for(int i=0; i<xyv.size(); i++) xyv[i].index=i; //Redefine indices
  long medpt = xyv.size()/2;
  if(dim%2==1) sort(xyv.begin(), xyv.end(), xyind_lower_x());
  else sort(xyv.begin(), xyv.end(), xyind_lower_y());
  return(xyv[medpt].index);
}

int splitxy(const vector <xy_index> &xyvec, int dim, long splitpoint, vector <xy_index> &left, vector <xy_index> &right)
{
  long i=0;
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
int kdtree01(const vector <xy_index> &xyvec, int dim, long rootptxy, long rootptkd, vector <kdpoint> &kdvec)
{
  int status=0;
  int lmed=0;
  int rmed=0;
  int kdct = kdvec.size()-1;
  int i=0;
  long leftrootkd=-1;
  long rightrootkd=-1;
  xy_index xyi = xy_index(0.0,0.0,0);
  kdpoint root = kdvec[kdct];
  kdpoint lp = kdpoint(xyi,-1,-1,0);
  kdpoint rp = kdpoint(xyi,-1,-1,0);
  kdpoint kdtest = kdpoint(xyi,-1,-1,0);
  vector <xy_index> leftvec = {};
  vector <xy_index> rightvec = {};

  status = splitxy(xyvec,dim,rootptxy,leftvec,rightvec);
  
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
  int branchct=0;
  double rng2 = range*range;
  int notdone=1;
  int kdveclen = kdvec.size();
  int dim=1;
  int currentpoint=0;
  int leftpoint=0;
  int rightpoint=0;
  int goleft=0;
  int goright=0;
  double xdiff=0.0;
  double ydiff=0.0;
  vector <long> checkit={};
  int i=0;
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
  for(long i=0; i<pvec.size(); i++) pvec[i].i1=i; //Redefine indices
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
int splitLDx2(const vector <point6LDx2> &pointvec, int dim, long splitpoint, vector <point6LDx2> &left, vector <point6LDx2> &right)
{
  long i=0;
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
int kdtree_6D01(const vector <point6LDx2> &invec, int dim, long splitpoint, long kdroot, vector <KD_point6LDx2> &kdvec)
{
  int status=0;
  int lmed=0;
  int rmed=0;
  int kdct = kdvec.size()-1;
  int i=0;
  long leftrootkd=-1;
  long rightrootkd=-1;
  point6LDx2 point0 = point6LDx2(0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0,0);
  KD_point6LDx2 root = kdvec[kdct];
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
  status = splitLDx2(invec,dim,splitpoint,leftvec,rightvec);

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
  int branchct=0;
  long double rng2 = range*range;
  int notdone=1;
  int kdveclen = kdvec.size();
  int dim=1;
  int currentpoint=0;
  int leftpoint=0;
  int rightpoint=0;
  int goleft=0;
  int goright=0;
  long double pointdiff = 0.0L;
  long double pdist2 = 0.0L;
  vector <long> checkit={};
  int i=0;
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
  
  for(int i=0; i<cluster.size(); i++) {
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
  
  for(int i=0; i<cluster.size(); i++) {
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
  int clustptnum=0;
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
  int i=0;
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
      else if(queryout.size() >= npt) {
	// This is a core point of a new cluster.
	cout << "Point " << kdct << ": cluster core with " << queryout.size() << " neighbors.\n";
	fakeclusternum++;
	kdtree[kdct].flag = fakeclusternum;
	// Begin loading cluster
	clusterind.push_back(kdct);
	cluster.push_back(kdtree[kdct]);
	// Loop on points in cluster.
	clustptct=0;
	while(clustptct<queryout.size()) {
	  if(kdtree[queryout[clustptct]].flag != 0) {
	    // Current point has already been considered: skip to the next.
	    clustptct++;
	  } else {
	    // Range-query current cluster point.
	    querypoint = kdtree[queryout[clustptct]].point;
	    subquery={};
	    kdrange_6D01(kdtree, querypoint, clustrad, subquery);
	    if(subquery.size()>=npt) {
	      // This point is a core point.
	      kdtree[queryout[clustptct]].flag = fakeclusternum;
	      clusterind.push_back(queryout[clustptct]);
	      cluster.push_back(kdtree[queryout[clustptct]]);
	      // Add additional points to queryout as appropriate
	      for(i=0;i<subquery.size();i++) {
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
	for(i=0;i<clusterind.size();i++) {
	  pointind.push_back(kdtree[clusterind[i]].point.i1);
	  pointind.push_back(kdtree[clusterind[i]].point.i2);
	}
	// Sort vector of detection indices
	sort(pointind.begin(), pointind.end());
	// Cull out duplicate entries
	pointjunk = pointind;
	pointind={};
	pointind.push_back(pointjunk[0]);
	for(i=1; i<pointjunk.size(); i++) {
	  if(pointjunk[i]!=pointjunk[i-1]) pointind.push_back(pointjunk[i]);
	}
	// Load vector of detection MJD's
	clustmjd = {};
	for(i=0; i<pointind.size(); i++) {
	  clustmjd.push_back(detvec[pointind[i]].MJD);
	}
	// Sort vector of MJD's
	sort(clustmjd.begin(), clustmjd.end());
	timespan = clustmjd[clustmjd.size()-1] - clustmjd[0];
	// Load vector of MJD steps
	mjdstep={};
	for(i=1; i<clustmjd.size(); i++) {
	  mjdstep.push_back(clustmjd[i]-clustmjd[i-1]);
	}
	// Count steps large enough to suggest a daytime period between nights.
	numdaysteps=0;	
	for(i=0; i<mjdstep.size(); i++) {
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
	  for(i=0; i<pointind.size(); i++) {
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
  int clustptnum=0;
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
  long double trms = 0.0L;
  int i=0;

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
      else if(queryout.size() >= npt) {
	// This is a core point of a new cluster.
	cout << "Point " << kdct << ": cluster core with " << queryout.size() << " neighbors.\n";
	clusternum++;
	kdtree[kdct].flag = clusternum;
	// Begin loading cluster
	clusterind.push_back(kdct);
	cluster.push_back(kdtree[kdct]);
	// Loop on points in cluster.
	clustptct=0;
	while(clustptct<queryout.size()) {
	  if(kdtree[queryout[clustptct]].flag != 0) {
	    // Current point has already been considered: skip to the next.
	    clustptct++;
	  } else {
	    // Range-query current cluster point.
	    querypoint = kdtree[queryout[clustptct]].point;
	    subquery={};
	    kdrange_6D01(kdtree, querypoint, clustrad, subquery);
	    if(subquery.size()>=npt) {
	      // This point is a core point.
	      kdtree[queryout[clustptct]].flag = clusternum;
	      clusterind.push_back(queryout[clustptct]);
	      cluster.push_back(kdtree[queryout[clustptct]]);
	      // Add additional points to queryout as appropriate
	      for(i=0;i<subquery.size();i++) {
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
	trms = cluster_stats6D01(cluster, meanvec, rmsvec);
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

long medind_6ix2(const vector <point6ix2> &pointvec, int dim)
{
  vector <point6ix2> pvec = pointvec; //Mutable copy of immutable input vector
  for(long i=0; i<pvec.size(); i++) pvec[i].i1=i; //Redefine indices
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
int splitix2(const vector <point6ix2> &pointvec, int dim, long splitpoint, vector <point6ix2> &left, vector <point6ix2> &right)
{
  long i=0;
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
int kdtree_6i01(const vector <point6ix2> &invec, int dim, long splitpoint, long kdroot, vector <KD_point6ix2> &kdvec)
{
  int status=0;
  int lmed=0;
  int rmed=0;
  int kdct = kdvec.size()-1;
  int i=0;
  long leftrootkd=-1;
  long rightrootkd=-1;
  point6ix2 point0 = point6ix2(0,0,0,0,0,0,0,0);
  KD_point6ix2 root = kdvec[kdct];
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
  status = splitix2(invec,dim,splitpoint,leftvec,rightvec);

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
  int branchct=0;
  long rng2 = range*range;
  int notdone=1;
  int kdveclen = kdvec.size();
  int dim=1;
  int currentpoint=0;
  int leftpoint=0;
  int rightpoint=0;
  int goleft=0;
  int goright=0;
  long pointdiff = 0;
  long pdist2 = 0;
  vector <long> checkit={};
  int i=0;
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
  
  for(int i=0; i<cluster.size(); i++) {
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
  
  for(int i=0; i<cluster.size(); i++) {
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
int DBSCAN_6i01(vector <KD_point6ix2> &kdtree, double clustrad, int npt, double intconvscale, vector <KD6i_clust> &outclusters)
{
  long kdnum = kdtree.size();
  long kdct=0;
  int clustptnum=0;
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
  double trms = 0.0l;
  int i=0;

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
	// cout << "Point " << kdct << ": noise\n";
      }
      else if(queryout.size() >= npt) {
	// This is a core point of a new cluster.
	cout << "Point " << kdct << ": cluster core with " << queryout.size() << " neighbors.\n";
	clusternum++;
	kdtree[kdct].flag = clusternum;
	// Begin loading cluster
	clusterind.push_back(kdct);
	cluster.push_back(kdtree[kdct]);
	// Loop on points in cluster.
	clustptct=0;
	while(clustptct<queryout.size()) {
	  if(kdtree[queryout[clustptct]].flag != 0) {
	    // Current point has already been considered: skip to the next.
	    clustptct++;
	  } else {
	    // Range-query current cluster point.
	    querypoint = kdtree[queryout[clustptct]].point;
	    subquery={};
	    kdrange_6i01(kdtree, querypoint, clustrad, subquery);
	    if(subquery.size()>=npt) {
	      // This point is a core point.
	      kdtree[queryout[clustptct]].flag = clusternum;
	      clusterind.push_back(queryout[clustptct]);
	      cluster.push_back(kdtree[queryout[clustptct]]);
	      // Add additional points to queryout as appropriate
	      for(i=0;i<subquery.size();i++) {
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
	if(cluster.size()>=npt) {	  
	  // This cluster has enough points to be considered.
	  // Calculate some cluster statistics.
	  meanvec = rmsvec = {};
	  trms = cluster_stats6i01(cluster, intconvscale, meanvec, rmsvec);
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
  // -x and y are switched above becuase we are rotating by 90 degrees
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
  // -x and y are switched above becuase we are rotating by 90 degrees
  // after the pole-switch, to get the old North Celestial Pole
  // on the +y axis where it should be.
  return(0);
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
  int i=0;
  char c = '0';
  int reachedend=0;
  string teststring, lnfromfile;
  double x,y,z,vx,vy,vz,MJD;
  MJD=x=y=z=vz=vy=vz=0.0;
  
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
  int i=0;
  char c = '0';
  int reachedend=0;
  string teststring, lnfromfile;
  long double x,y,z,vx,vy,vz,MJD;
  MJD=x=y=z=vz=vy=vz=0.0;
  
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


/*solvematrix01: November 23, 2021
Given a matrix with dimensions eqnum,eqnum+1, interpret
it as a system of eqnum linear equations in eqnum unknowns,
with the first term in each equation being the constant term
and the others being the coefficients of x1,x2,x3,etc;
solve for the vector of x values or report the matrix to
be singular.*/
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
  int i=0;
  int j=0;
  int k=0;
  int status=0;
  int npoints = x.size();
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
  int i=0;
  int j=0;
  int k=0;
  int status=0;
  int npoints = x.size();
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
  int fitnum = polyorder+1;
  int pointsbefore = fitnum - fitnum/2;
  int pbf=0;
  vector <double> xvec;
  vector <double> yvec;
  vector <double> fitvec;
  double tdelt=0;
  double sumvar=0;
  int i=0;
  int j=0;
  int k=0;
  make_dvec(fitnum,fitvec);

  // Convert input time from UT1 standard to the dynamical time (TT) used
  // by JPL Horizons for the state vectors.
  detmjd += TTDELTAT/SOLARDAY;
  // TT is ahead of UT1 because the Earth's rotation is slowing down.
  
  //Interpolate to find the planet's exact position at the time
  //of the detection.
  pbf=0;
  i=posmjd.size();
  if(planetpos.size()!=i) {
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
  int fitnum = polyorder+1;
  int pointsbefore = fitnum - fitnum/2;
  int pbf=0;
  vector <long double> xvec;
  vector <long double> yvec;
  vector <long double> fitvec;
  long double tdelt=0;
  long double sumvar=0;
  int i=0;
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
  i=posmjd.size();
  if(planetpos.size()!=i) {
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
  int fitnum = polyorder+1;
  int pointsbefore = fitnum - fitnum/2;
  int pbf=0;
  vector <long double> xvec;
  vector <long double> yvec;
  vector <long double> fitvec;
  long double tdelt=0;
  long double sumvar=0;
  int i=0;
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
  i=posmjd.size();
  if(planetpos.size()!=i) {
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
  double alphapos,alphaneg;
  double opdistcos,sunelong,barydist;
  double obsdot,opelong;
  
  //cout << fixed << setprecision(9) << "helioproj01 input observer pos: " << obsbary.x << " " << obsbary.y << " " << obsbary.z << "\n";
  //cout << "barycentric unit vector: " << unitbary.x << " " << unitbary.y << " " << unitbary.z << "\n";

  barydist = sqrt(dotprod3d(obsbary,obsbary));
  obsdot = dotprod3d(unitbary,obsbary);
  opdistcos = obsdot/barydist;
  opelong = acos(opdistcos)*DEGPRAD;
  if(opelong < 0.0) opelong = 90.0 - opelong;
  sunelong = 180.0 - opelong;

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
  long double alphapos,alphaneg,obsdot,barydist;
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
  int endnear=0;
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
  if(fabs(planetmjd[lastpoint]-mjdend) < fabs(planetmjd[lastpoint+1]-mjdend)) endnear = lastpoint;
  else endnear = lastpoint+1;

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
    
// Keplerint: December 08, 2021:
// Integrate an orbit assuming we have a Keplerian 2-body problem
// with all the mass in the Sun, and the input position and velocity
// are relative to the Sun.
int Keplerint(const long double MGsun, const long double mjdstart, const point3LD &startpos, const point3LD &startvel, const long double mjdend, point3LD &endpos, point3LD &endvel)
{
  long double e,E,a,lscalar,r0,v0,r1,v1;
  point3LD lvec = point3LD(0L,0L,0L);
  point3LD lunit = point3LD(0L,0L,0L);
  point3LD r0unit = point3LD(0L,0L,0L);
  point3LD r1unit = point3LD(0L,0L,0L);
  point3LD v1unit = point3LD(0L,0L,0L);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0L;
  long double costheta,theta0,theta1,radvel,cospsi,psi;
  costheta = theta0 = theta1 = radvel = cospsi = psi = 0L;
  long double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0L;
  long double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0L;
  int status = 0;
  long double junkra,junkdec,sinev,thetav,v1ra,v1dec;
  junkra = junkdec = sinev = thetav = v1ra = v1dec = 0L;
    
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
    cerr << "WARNING: Keplerint finds costheta+1.0 = " << costheta+1.0 << "\n";
    costheta = -1.0L;
    theta0 = M_PI;
  } else if(costheta>1.0) {
    cerr << "WARNING: Keplerint finds costheta-1.0 = " << costheta-1.0 << "\n";
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
    cerr << "WARNING: Keplerint finds cospsi+1.0 = " << cospsi+1.0 << "\n";
    cospsi = -1.0L;
    psi = M_PI;
  } else if(cospsi>1.0) {
    cerr << "WARNING: Keplerint finds cospsi-1.0 = " << cospsi-1.0 << "\n";
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
  int verbose=1;
  long double vtheta1,vtheta2,posRA,posDec,velRA,velDec;
  long double GMsun=0.0L;
  long double omega=0.0L;
  long double E=0.0L;
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
  point3LD lunit = point3LD(0L,0L,0L);
  point3LD r0unit = point3LD(0L,0L,0L);
  point3LD r1unit = point3LD(0L,0L,0L);
  point3LD v1unit = point3LD(0L,0L,0L);
  e = E = a = lscalar = r0 = v0 = r1 = v1 = 0L;
  long double coshH,H0,theta0,theta1,radvel,H;
  coshH = H0 = theta0 = theta1 = radvel = H = 0L;
  long double omega, t0omega, t0, t1omega, t1;
  omega = t0omega = t0 = t1omega = t1 = 0L;
  long double lra,ldec,r0ra,r0dec,r1ra,r1dec,newra,newdec;
  lra = ldec = r0ra = r0dec = r1ra = r1dec = newra = newdec = 0L;
  int status = 0;
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
  int planetpointct = 0;
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
  int i=0;
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
  int i=0;
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
  int i=0;
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
  int i=0;
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
  int i=0;
  int lnct=0;
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
  int i=0;
  int lnct=0;
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
  int pnum = poscluster.size();
  int i=0;
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
  double delta,xal,yal,xty,xsq,nsum,rms,err,errmax;
  double siga,sigb;

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
  for(j=0;j<str.size();j++) outstring.push_back(str[j]);
  return(outstring);
}

// get_csv_string01: Given a line read from a csv file, and an
// starting point along that line, read the next comma-separated value,
// and put it into the output string. If the read was successful, return
// the line index of the comma, newline, or EOF at the end of the value read.
// Otherwise, return -1 as an error code.
int get_csv_string01(const string &lnfromfile, string &outstring, int startpoint)
{
  int i=startpoint;
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
  int i=0;
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
  int endhere=-1;
  int planetpointnum = planetmjd.size();
  vector <long double> forwardmjd;
  vector <long double> backwardmjd;
   int planetpointct = 0;
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
      for(i=0;i<forwardmjd.size();i++) {
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
    while(latestpoint<forwardmjd.size()-1) {
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
    for(j=0;j<temptime.size();j++) {
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
      for(i=0;i<backwardmjd.size();i++) {
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
    while(latestpoint<backwardmjd.size()-1) {
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
  int k=0;
  int outnum = endpoint-startpoint+1;
  int outct=0;
  int endhere=-1;
  int planetpointnum = planetmjd.size();
   int planetpointct = 0;
  int pointafter=0;
  int pointbefore=0;
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
  } else if(startpoint<0 || endpoint>=planetmjd.size()) {
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
  double x;
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
  long double barydist,normdot1,normaldist;
  point3LD plane_normvec = point3LD(0L,0L,0L);
  point3LD plane_to_obs = point3LD(0L,0L,0L);
  point3LD periobs = point3LD(0L,0L,0L);
  
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
  
  // 3d. Subtract the resulting physical vector from the instantaneous position of the observer,
  //     yielding the position of the point on the plane closest to the observer.
  //     Call this point 'periobs'.
  periobs.x = obsbary.x - plane_to_obs.x;
  periobs.y = obsbary.y - plane_to_obs.y;
  periobs.z = obsbary.z - plane_to_obs.z;
  
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
