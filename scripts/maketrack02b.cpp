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

class det_indx{ // Detection of some astronomical source with
                  // cross-reference index.
public:
  double MJD;
  double RA;
  double Dec;
  long indx;
  det_indx(double mjd, double ra, double dec, long indx) :MJD(mjd), RA(ra), Dec(dec), indx(indx) { }
};

class xy_indx{ // Double-precision x,y point plus long index
public:
  double x;
  double y;
  long indx;
  xy_indx(double x, double y, long indx) :x(x), y(y), indx(indx) { }
};

class kdpoint { // Point in an xy_indx k-d tree designed as a vector,
                // with left and right giving the index of the left and
                // right branches, and the dimension of splitting
                // specified at each node.
public:
  xy_indx point;
  long left;
  long right;
  int dim;
  kdpoint(xy_indx point, long left, long right, int dim) :point(point), left(left), right(right), dim(dim) {}
};
  
class xyind_lower_x{ // Function for sorting vectors of type xy_indx by x
public:
  inline bool operator() (const xy_indx& p1, const xy_indx& p2) {
    return(p1.x < p2.x);
  }
};

class xyind_lower_y{ // Function for sorting vectors of type xy_indx by y
public:
  inline bool operator() (const xy_indx& p1, const xy_indx& p2) {
    return(p1.y < p2.y);
  }
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
  inline bool operator() (const det_indx& o1, const det_indx& o2) {
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
  double asinepa=0.0;
  
  // Calculate the distance d, in radians, between
  // the two points.
  x1=cos(dec1/DEGPRAD)*cos(ra1/DEGPRAD);
  y1=cos(dec1/DEGPRAD)*sin(ra1/DEGPRAD);
  z1=sin(dec1/DEGPRAD);
  x2=cos(dec2/DEGPRAD)*cos(ra2/DEGPRAD);
  y2=cos(dec2/DEGPRAD)*sin(ra2/DEGPRAD);
  z2=sin(dec2/DEGPRAD);
  h=sqrt(DSQUARE(x1-x2)+DSQUARE(y1-y2)+DSQUARE(z1-z2));
  if(fabs(h)<=2.0) d=2.0*asin(h/2.0);
  else cerr << "WARNING: argument of asin is " << h/2.0 << "\n";
  
  *dist = d*DEGPRAD;
  if(!isnormal(*dist)) cerr << "WARNING: found arc distance " << *dist << "\n";

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
    if(fabs(sinepa)<1.0) asinepa=asin(sinepa);
    else if(sinepa >= 1.0) asinepa = M_PI/2.0;
    else if(sinepa <= -1.0) asinepa = -M_PI/2.0;
    if(cosinepa>=0.0) celpa = asinepa;
    else celpa = M_PI - asinepa;
    *pa = celpa*DEGPRAD;
    if(!isnormal(*pa)) {
      cerr << "WARNING: found position angle " << *pa << "\n";
      cerr << "ra1,2: " << ra1 << " " << ra2 << "dec1,2: " << dec1 << " " << dec2 << " poleangle: " << poleangle << " colats: " << colat1 << " " << colat2 << " sinepa, cosinepa: " << sinepa << " " << cosinepa << "\n";
    }
  }
  else if(ra1<ra2) {
    // Wrapped case with point 2 west of point 1
    // across zero RA: the PA should be greater
    // than 180 degrees.
    poleangle = (ra1+(double)360.0-ra2)/DEGPRAD;
    sinepa = sin(colat2)*sin(poleangle)/sin(d);
    if(fabs(sinepa)<1.0) asinepa=asin(sinepa);
    else if(sinepa >= 1.0) asinepa = M_PI/2.0;
    else if(sinepa <= -1.0) asinepa = -M_PI/2.0;
    if(cosinepa>=0.0) celpa = asinepa;
    else celpa = M_PI - asinepa;
    *pa = (double)360.0 - celpa*DEGPRAD;
    if(!isnormal(*pa)) {
      cerr << "WARNING: found position angle " << *pa << "\n";
      cerr << "ra1,2: " << ra1 << " " << ra2 << "dec1,2: " << dec1 << " " << dec2 << " poleangle: " << poleangle << " colats: " << colat1 << " " << colat2 << " sinepa, cosinepa: " << sinepa << " " << cosinepa << "\n";
    }
  }
  else if(ra1>ra2 && (ra1-ra2)<=180.0) {
    // Simple case, point 2 is west of point 1,
    // so PA should be greater than 180 degrees.
    poleangle = (ra1-ra2)/DEGPRAD;
    sinepa = sin(colat2)*sin(poleangle)/sin(d);
    if(fabs(sinepa)<1.0) asinepa=asin(sinepa);
    else if(sinepa >= 1.0) asinepa = M_PI/2.0;
    else if(sinepa <= -1.0) asinepa = -M_PI/2.0;
    if(cosinepa>=0.0) celpa = asinepa;
    else celpa = M_PI - asinepa;
    *pa = (double)360.0 - celpa*DEGPRAD;
    if(!isnormal(*pa)) {
      cerr << "WARNING: found position angle " << *pa << "\n";
      cerr << "ra1,2: " << ra1 << " " << ra2 << "dec1,2: " << dec1 << " " << dec2 << " poleangle: " << poleangle << " colats: " << colat1 << " " << colat2 << " sinepa, cosinepa: " << sinepa << " " << cosinepa << "\n";
    }
  }
  else if(ra1>ra2) {
    // Wrapped case with point 2 east of point 1
    // across zero RA: the PA should be less
    // than 180.0 degrees.
    poleangle = (ra2+(double)360.0-ra1)/DEGPRAD;
    sinepa = sin(colat2)*sin(poleangle)/sin(d);
    if(fabs(sinepa)<1.0) asinepa=asin(sinepa);
    else if(sinepa >= 1.0) asinepa = M_PI/2.0;
    else if(sinepa <= -1.0) asinepa = -M_PI/2.0;
    if(cosinepa>=0.0) celpa = asinepa;
    else celpa = M_PI - asinepa;
    *pa = celpa*DEGPRAD;
    if(!isnormal(*pa)) {
      cerr << "WARNING: found position angle " << *pa << "\n";
      cerr << "ra1,2: " << ra1 << " " << ra2 << "dec1,2: " << dec1 << " " << dec2 << " poleangle: " << poleangle << " colats: " << colat1 << " " << colat2 << " sinepa, cosinepa: " << sinepa << " " << cosinepa << "\n";
    }
  }
  return(0);
}

long medindex(const vector <xy_indx> &xyvec, int dim)
{
  vector <xy_indx> xyv = xyvec; //Mutable copy of immutable input vector
  for(int i=0; i<xyv.size(); i++) xyv[i].indx=i; //Redefine indices
  long medpt = xyv.size()/2;
  if(dim%2==1) sort(xyv.begin(), xyv.end(), xyind_lower_x());
  else sort(xyv.begin(), xyv.end(), xyind_lower_y());
  return(xyv[medpt].indx);
}

int splitxy(const vector <xy_indx> &xyvec, int dim, long splitpoint, vector <xy_indx> &left, vector <xy_indx> &right)
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
int kdtree01(const vector <xy_indx> &xyvec, int dim, long rootptxy, long rootptkd, vector <kdpoint> &kdvec)
{
  int status=0;
  int lmed=0;
  int rmed=0;
  int kdct = kdvec.size()-1;
  int i=0;
  long leftrootkd=-1;
  long rightrootkd=-1;
  xy_indx xyi = xy_indx(0.0,0.0,0);
  kdpoint root = kdvec[kdct];
  kdpoint lp = kdpoint(xyi,-1,-1,0);
  kdpoint rp = kdpoint(xyi,-1,-1,0);
  kdpoint kdtest = kdpoint(xyi,-1,-1,0);
  vector <xy_indx> leftvec = {};
  vector <xy_indx> rightvec = {};

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
int kdrange01(const vector <kdpoint> &kdvec,double x,double y,double range,vector <long> &indxvec)
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
	      indxvec.push_back(currentpoint);
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
	      indxvec.push_back(currentpoint);
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
      
      
static void show_usage()
{
  cerr << "Usage: maketrack02a -dets detfile -imgs imfile -outimgs output image file/ \n";
  cerr << "-pairs pairfile -pairdets paired detection file/ \n";
  cerr << "-imrad image radius(deg) -maxtime max inter-image time interval (hr)/ \n";
  cerr << "-maxvel maximum angular velocity (deg/day)\n";
}
    
int main(int argc, char *argv[])
{
  det_indx o1 = det_indx(0,0,0,0);
  vector <det_indx> detvec = {};
  vector <det_indx> pairdets = {};
  img_log02 imlog = img_log02(0.0,0.0,0.0,0,0);
  vector <img_log02> img_log = {};
  longpair onepair = longpair(0,0);
  vector <longpair> pairvec ={};
  point3d p3 = point3d(0,0,0);
  point3d p3avg = point3d(0,0,0);
  double tdelt = 0;
  double mjdmean = 0;
  double mjdnorm = 0;
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
  double MJD,RA,Dec;
  double maxvel = MAXVEL; // Max angular velocity in deg/day
  double maxtime = MAXTIME; // Max time interval a tracklet could span,
                            // in days.
  double maxdist = MAXVEL*MAXTIME/24.0; // Max angular distance a tracklet
                                   // could span, in degrees.
  double imrad = IMAGERAD; // radius from image center to most distant corner (deg).
  string indetfile;
  string inimfile;
  string outimfile;
  string outpairfile="outpairfile01.txt";
  string pairdetfile="pairdetfile01.txt";
  

  if(argc<3)
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
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }
  cout << "indet file " << indetfile << "\n";
  cout << "inimage file " << inimfile << "\n";
  cout << "output image file " << outimfile << "\n";
  cout << "pairfile file " << outpairfile << "\n";
  cout << "paired detection file " << pairdetfile << "\n";
  cout << "image radius " << imrad << "\n";
  cout << "max time interval " << maxtime << "\n";
  cout << "maxvel " << maxvel << "\n";
  maxtime/=24.0; /*Unit conversion from hours to days*/
  maxdist = maxtime*maxvel;
  
  ifstream instream1 {indetfile};
  if(!instream1) error("can't open input file ",indetfile);
  // Skip one-line header
  getline(instream1,lnfromfile);
  cout << lnfromfile << "\n";
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
      while(i<lnfromfile.size() && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==MJDCOL) MJD=stod(stest);
      else if(j==RACOL) RA=stod(stest);
      else if(j==DECCOL) Dec=stod(stest);
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0) {
      // cout<<"MJD, RA, Dec: " << MJD-floor(MJD) << " " << RA << " " << Dec << "\n";
      o1=det_indx(MJD,RA,Dec,-1);
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
  sort(detvec.begin(), detvec.end(), early_detindx());

  // Get image information.
  if(inimfile.size()>0)
    {
      // Read input image file: MJD, RA, Dec:
        ifstream instream1 {inimfile};
	if(!instream1) error("can't open input file ",inimfile);
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
      vector <det_indx> imobs = {};
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
      xy_indx xyind=xy_indx(0.0, 0.0, 0);
      vector <xy_indx> axyvec = {};
      double dist,pa;
      int dettarg=0;
      for(detct=img_log[imct].startind ; detct<img_log[imct].endind ; detct++) {
	distradec02(img_log[imct].RA, img_log[imct].Dec,detvec[detct].RA,detvec[detct].Dec,&dist,&pa);
	xyind = xy_indx(dist*sin(pa/DEGPRAD),dist*cos(pa/DEGPRAD),detct);
	axyvec.push_back(xyind);
      }
      // Loop over images with potential matches (image B's)
      for(imatchcount=0;imatchcount<imagematches.size();imatchcount++)
      {
	imtarg = imagematches[imatchcount];
	double range = (img_log[imtarg].MJD-img_log[imct].MJD)*maxvel;
	vector <xy_indx> bxyvec = {};
	// Project all detections on image B
	for(dettarg=img_log[imtarg].startind ; dettarg<img_log[imtarg].endind ; dettarg++) {
	  distradec02(img_log[imct].RA, img_log[imct].Dec,detvec[dettarg].RA,detvec[dettarg].Dec,&dist,&pa);
	  xyind = xy_indx(dist*sin(pa/DEGPRAD),dist*cos(pa/DEGPRAD),dettarg);
	  bxyvec.push_back(xyind);
	}
	// Create k-d tree of detections on image B (imtarg).
	int dim=1;
	xy_indx xyi = bxyvec[0];
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
	  vector <long> indxvec = {};
	  if(isnormal(axyvec[detct].x) &&isnormal(axyvec[detct].y)) {
	     kdrange01(kdvec,axyvec[detct].x,axyvec[detct].y,range,indxvec);
	  } else {
	    cout << "WARNING: detection " << detct << " on image " << imct << " not normal: " << axyvec[detct].x << " " << axyvec[detct].y << "\n";
	  }
	  int matchnum=indxvec.size();
	  long matchpt=0;
	  int matchct=0;
	  if(matchnum>0) {
	    // Record image A detection as paired, if not already recorded.
	    if(detvec[axyvec[detct].indx].indx<0) {
	      //This detection has not yet been paired with any other.
	      detvec[axyvec[detct].indx].indx = pdct; // Mark as paired
	      pairdets.push_back(detvec[axyvec[detct].indx]); // Load into paired detection vector
	      pairdets[pdct].indx = axyvec[detct].indx; // Record index from original detection vector
	      pdct++; // Increment count of paired detections
	      adetct++;
	      if(pdct!=pairdets.size()) cerr << "\nWARNING: PAIRED DETECTION MISMATCH: " << pdct << " vs " << pairdets.size() << "\n\n";
	    }
	    // Record image B detections
	    for(matchct=0;matchct<matchnum;matchct++) {
	      matchpt = indxvec[matchct];
	      if(detvec[kdvec[matchpt].point.indx].indx<0) {
		//This detection has not yet been paired with any other.
		detvec[kdvec[matchpt].point.indx].indx = pdct; // Mark as paired
		pairdets.push_back(detvec[kdvec[matchpt].point.indx]); // Load into paired detection vector
		pairdets[pdct].indx = kdvec[matchpt].point.indx; // Record index from original detection vector
		pdct++; // Increment count of paired detections
		if(pdct!=pairdets.size()) cerr << "\nWARNING: PAIRED DETECTION MISMATCH: " << pdct << " vs " << pairdets.size() << "\n\n";
	      }
	      // Write indx values for both components of the
	      // new pair to the pair vector, regardless of whether
	      // the indx values are pre-existing or newly assigned.
	      onepair = longpair(detvec[axyvec[detct].indx].indx,detvec[kdvec[matchpt].point.indx].indx);
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
  for(i=0;i<pairdets.size();i++) {
    outstream3 << fixed << setprecision(6) << pairdets[i].MJD << " " << pairdets[i].RA << " " << pairdets[i].Dec << " " << pairdets[i].indx << "\n"; 
  }
  // Write pair vector to file
  cout << "Writing pairs to output file\n";
  ofstream outstream4 {outpairfile};
  for(i=0;i<pairvec.size();i++) {
    outstream4 <<  pairvec[i].i1 << " " <<  pairvec[i].i2 << "\n";
  }
  
  // Write indx'ed input detection vector to file as a test.
  cout << "Writing indexed input detection vector as a test\n";
  ofstream outstream5 {"indextest01.txt"};
  for(i=0;i<detvec.size();i++) {
    outstream5 << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << " " << detvec[i].indx << "\n"; 
  }

  return(0);
}
