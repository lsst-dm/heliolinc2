// December 07, 2021: solarsyst_dyn_geo01.cpp
// Library of functions useful for Solar System dynamics and
// geometry. Includes functions originally developed in orbint02a.cpp,
// maketrack02b.cpp, projtest01b.cpp, and other places.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

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

double dotprod3LD(point3LD p1, point3LD p2)
{
  return(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z);
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

int celestial_to_statevec(double RA, double Dec,double delta,point3d &baryvec)
{
  double x,y,z,theta,phi,thetapole,phipole;
  x = y = z = theta = phi = thetapole = phipole = 0.0;
  theta = Dec/DEGPRAD;
  phi = RA/DEGPRAD;
  thetapole = NEPDEC/DEGPRAD;
      
  z = sin(theta);
  x = -cos(theta)*sin(phi);
  y = cos(theta)*cos(phi);
  baryvec.z = delta*(z*sin(thetapole) + x*cos(thetapole));
  baryvec.y = delta*(z*cos(thetapole) - x*sin(thetapole));
  baryvec.x = delta*y;
  return(0);
}

int celestial_to_statevecLD(long double RA, long double Dec,long double delta,point3LD &baryvec)
{
  long double x,y,z,theta,phi,thetapole,phipole;
  x = y = z = theta = phi = thetapole = phipole = 0.0;
  theta = Dec/DEGPRAD;
  phi = RA/DEGPRAD;
  thetapole = NEPDEC/DEGPRAD;
      
  z = sin(theta);
  x = -cos(theta)*sin(phi);
  y = cos(theta)*cos(phi);
  baryvec.z = delta*(z*sin(thetapole) + x*cos(thetapole));
  baryvec.y = delta*(z*cos(thetapole) - x*sin(thetapole));
  baryvec.x = delta*y;
  return(0);
}

int celestial_to_stateunit(double RA, double Dec,point3d &baryvec)
{
  double x,y,z,theta,phi,thetapole,phipole;
  x = y = z = theta = phi = thetapole = phipole = 0.0;
  theta = Dec/DEGPRAD;
  phi = RA/DEGPRAD;
  thetapole = NEPDEC/DEGPRAD;
      
  z = sin(theta);
  x = -cos(theta)*sin(phi);
  y = cos(theta)*cos(phi);
  baryvec.z = z*sin(thetapole) + x*cos(thetapole);
  baryvec.y = z*cos(thetapole) - x*sin(thetapole);
  baryvec.x = y;
  return(0);
}

int celestial_to_stateunitLD(long double RA, long double Dec, point3LD &baryvec)
{
  long double x,y,z,theta,phi,thetapole,phipole;
  x = y = z = theta = phi = thetapole = phipole = 0.0;
  theta = Dec/DEGPRAD;
  phi = RA/DEGPRAD;
  thetapole = NEPDEC/DEGPRAD;
      
  z = sin(theta);
  x = -cos(theta)*sin(phi);
  y = cos(theta)*cos(phi);
  baryvec.z = z*sin(thetapole) + x*cos(thetapole);
  baryvec.y = z*cos(thetapole) - x*sin(thetapole);
  baryvec.x = y;
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
int poleswitch01(const double &inRA, const double &inDec, const double &poleRA, const double &poleDec, const double &oldpoleRA, double *newRA, double *newDec)
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
  *newRA = phip;
  *newDec = thetap;
  return(0);
}


// poleswitch01LD: December 9, 2021: given a double precision input
// celestial position IN RADIANS, and the celestial position of the pole
// of a new coordinate system, calculates the position
// of the point in the new coordinate system.  The desired
// RA for the old pole in the new coordinates is also required.
int poleswitch01LD(const long double &inRA, const long double &inDec, const long double &poleRA, const long double &poleDec, const long double &oldpoleRA, long double *newRA, long double *newDec)
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
  *newRA = phip;
  *newDec = thetap;
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
  outpos.x = planetpos[pbf].x + fitvec[0]*tdelt;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.x += sumvar;
    sumvar *= tdelt; // One more power of tdelt, for the position.
    outpos.x += sumvar;
  }
  // Load vector to fit y-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetvel[i].y);
  // Solve for polynomial interpolation
  perfectpoly01LD(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outvel.y = fitvec[0];
  outpos.y = planetpos[pbf].y + fitvec[0]*tdelt;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.y += sumvar;
    sumvar *= tdelt; // One more power of tdelt, for the position.
    outpos.y += sumvar;
  }
  // Load vector to fit z-coordinate
  yvec={};
  for(i=pbf;i<pbf+fitnum;i++) yvec.push_back(planetvel[i].z);
  // Solve for polynomial interpolation
  perfectpoly01LD(xvec,yvec,fitvec);
  // Calculate interpolated position.
  outvel.z = fitvec[0];
  outpos.z = planetpos[pbf].z + fitvec[0]*tdelt;
  for(j=1;j<fitnum;j++) {
    sumvar = fitvec[j]*tdelt;
    for(k=2;k<=j;k++) sumvar*=tdelt;
    outvel.z += sumvar;
    sumvar *= tdelt; // One more power of tdelt, for the position.
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
  cout << "obs_from_geocen: " << obs_from_geocen.x << " " << obs_from_geocen.y << " " << obs_from_geocen.z << " \n";
  cout << "geocen_from_barycen: " << geocen_from_barycen.x << " " << geocen_from_barycen.y << " " << geocen_from_barycen.z << "\n";
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
  
  cout << fixed << setprecision(9) << "helioproj01 input observer pos: " << obsbary.x << " " << obsbary.y << " " << obsbary.z << "\n";
  cout << "barycentric unit vector: " << unitbary.x << " " << unitbary.y << " " << unitbary.z << "\n";

  barydist = sqrt(dotprod3d(obsbary,obsbary));
  obsdot = dotprod3d(unitbary,obsbary);
  opdistcos = obsdot/barydist;
  opelong = acos(opdistcos)*DEGPRAD;
  if(opelong < 0.0) opelong = 90.0 - opelong;
  sunelong = 180.0 - opelong;

  cout << "barydist = " << barydist << ", obsdot = " << obsdot << ", opdistcos = " << opdistcos << ", opelong = " << opelong << ", sunelong = " << sunelong << "\n";
  
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
      cout << "a = " << a << ", b = " << b << ", c = " << c << "\n";
      alphapos = (-b + sqrt(b*b-4.0*a*c))/2.0/a;
      if(alphapos>0.0) {
	geodist=alphapos;
	projbary.x = obsbary.x + alphapos*unitbary.x;
	projbary.y = obsbary.y + alphapos*unitbary.y;
	projbary.z = obsbary.z + alphapos*unitbary.z;
	cout << "Distance from Earth: " << alphapos << "\n";
	cout << "Barycentric coordinates: " << projbary.x << " " << projbary.y << " " << projbary.z << "\n";
	return(0);
      }
      else return(-1);
    }
  return(-1);
}


int helioproj01LD(point3LD unitbary, point3LD obsbary, long double heliodist, long double &geodist, point3LD &projbary)
{
  long double a,b,c;
  long double alphapos,alphaneg;
  long double opdistcos,sunelong,barydist;
  long double obsdot,opelong;
  
  cout << fixed << setprecision(9) << "helioproj01 input observer pos: " << obsbary.x << " " << obsbary.y << " " << obsbary.z << "\n";
  cout << "barycentric unit vector: " << unitbary.x << " " << unitbary.y << " " << unitbary.z << "\n";

  barydist = sqrt(dotprod3LD(obsbary,obsbary));
  obsdot = dotprod3LD(unitbary,obsbary);
  opdistcos = obsdot/barydist;
  opelong = acos(opdistcos)*DEGPRAD;
  if(opelong < 0.0L) opelong = 90.0L - opelong;
  sunelong = 180.0L - opelong;

  cout << "barydist = " << barydist << ", obsdot = " << obsdot << ", opdistcos = " << opdistcos << ", opelong = " << opelong << ", sunelong = " << sunelong << "\n";
  
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
      cout << "a = " << a << ", b = " << b << ", c = " << c << "\n";
      alphapos = (-b + sqrt(b*b-4.0L*a*c))/2.0L/a;
      if(alphapos>0.0L) {
	geodist=alphapos;
	projbary.x = obsbary.x + alphapos*unitbary.x;
	projbary.y = obsbary.y + alphapos*unitbary.y;
	projbary.z = obsbary.z + alphapos*unitbary.z;
	cout << "Distance from Earth: " << alphapos << "\n";
	cout << "Barycentric coordinates: " << projbary.x << " " << projbary.y << " " << projbary.z << "\n";
	return(0);
      }
      else return(-1);
    }
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
    cout << "Warning: kep_trancendental used " << itct << " iterations and still had an error of " << fpsi << ", greater than tolerance of " << tol << "\n";
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
      cerr << "ERROR: Keplerint finds positive total energy: " << E << "\n";
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
 
  // The new time t1 for which we want to re-evaluate psi is
  // given by t0 + mjdend-mjdstart.
  t1omega = t0omega + (mjdend-mjdstart)*SOLARDAY*omega;
  //cout << " t1omega = " << t1omega;
  while(t1omega > 2.0L*M_PI) t1omega -= 2.0L*M_PI;
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
  poleswitch01LD(r0ra/DEGPRAD,r0dec/DEGPRAD,lra/DEGPRAD,ldec/DEGPRAD,0.0L,&newra,&newdec); // Output is radians
  //cout << "Orbital plane coords of original position: " << newra*DEGPRAD << " " << newdec*DEGPRAD << "\n";
  // Rotate starting unit vector around the angular momentum axis by
  // the calculated angle.
  newra += theta1-theta0;
  // cout << "Orbital plane coords of new position " << newra*DEGPRAD << " " << newdec*DEGPRAD << "\n";
  // The unit vector for the new position r1 is on the equator at this RA,
  // in the coordinate system that has the angular momentum vector at the pole.
  // Convert back to the original coordinate system.
  poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,&r1ra,&r1dec); // Output is radians
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
  poleswitch01LD(newra,0.0L,0.0L,ldec/DEGPRAD,lra/DEGPRAD,&v1ra,&v1dec); // Output is radians

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
	// Weird case where the requesed mjdend falls exactly on a timestep
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



