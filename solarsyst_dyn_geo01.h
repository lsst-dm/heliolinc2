// December 07, 2021: solarsyst_dyn_geo01.h:
// Library of functions useful for Solar System dynamics and
// geometry. Includes functions originally developed in orbint02a.cpp,
// maketrack02b.cpp, projtest01b.cpp, and other places.

using namespace std;

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <list>
#include <forward_list>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <array>
#include <regex>
#include <random>
#include <stdexcept>
#include <cstdio>

#define DEGPRAD (180.0L/M_PI) /*Degrees per radian*/
#define LSQUARE(x) long(x)*long(x)
#define DSQUARE(x) double(x)*double(x)
#define LDSQUARE(x) ((long double)x)*((long double)x)

#define SOLARDAY 86400.0L
#define NEPRA 270.0L //Right ascension of the North Ecliptic Pole.
#define NEPDEC 66.560708333333L // Declination of the North Ecliptic Pole
#define MJDOFF 2400000.5L // Offset from Julian Days to Modified Julian days.
#define EARTHEQUATRAD 6378.140L // Earth's equatorial radius.
#define AU_KM 1.495978700e8L /*One AU in km*/
#define AU 1.49597870700e11L /*One AU in meters*/

#define ZET0 2.5976176L /*precession formula constants in units of arcsec*/
#define ZET1 2306.0809506L
#define ZET2 0.3019015L
#define ZET3 0.0179663L
#define ZET4 -32.7e-6L
#define ZET5 -0.2e-6L
#define Z0 -2.5976176L
#define Z1 2306.0803226L
#define Z2 1.0947790L
#define Z3 0.0182273L
#define Z4 47.0e-6L
#define Z5 -0.3e-6L
#define THET1 2004.1917476L
#define THET2 -0.4269353L
#define THET3 -0.0418251L
#define THET4 -60.1e-6L
#define THET5 -0.1e-6L
#define SMALLANG 1.0e-7L /*angle in radians within which declinations are 
                           collapsed to the pole in precess01 routine.*/

#define TTDELTAT 70.0L // Time in seconds by which TT differs from UT. Thus,
                      // a light-travel-time corrected dynamical position at
                      // TT =  T + TTDELTAT corresponds to the celestial
                      // coordinates observed at UT = T.*/
#define CLIGHT 2.99792458e8L // Speed of light in m/sec
#define CLIGHT_AUDAY CLIGHT*SOLARDAY/AU // Speed of light in AU/day
#define KEPTRANSITMAX 50 // Maximum number of iterations to use in Newton's
                         // method solution of the trancendental Kepler Equation.
#define KEPTRANSTOL 1e-15L // Maximum error for an acceptable solution of the
                           // trancendental Kepler Equation.
#define LARGERR 1e30L // Large number supposed to be a safe initialization
                      // for most minimum-finding problems.
#define GMSUN_KM3_SEC2 132712440041.279419L // GM for the Sun: that is, the Universal
                                           // Gravitational Constant times the solar mass,
                                           // in units of km^3/sec^2. 


class det_bsc{ // Basic, minimal detection of some astronomical source 
public:
  double MJD;
  double RA;
  double Dec;
  det_bsc(double mjd, double ra, double dec) :MJD(mjd), RA(ra), Dec(dec) { }
};

class det_index{ // Detection of some astronomical source with
                  // cross-reference index.
public:
  double MJD;
  double RA;
  double Dec;
  long index;
  det_index(double mjd, double ra, double dec, long index) :MJD(mjd), RA(ra), Dec(dec), index(index) { }
};

class det_OC_index{ // Detection of some astronomical source with
                    // the observer's X, Y, Z coordinates, an index
                    // for cross-referencing, and a string identifier.
public:
  long double MJD;
  long double RA;
  long double Dec;
  long double x;
  long double y;
  long double z;
  string idstring;
  long index;
  det_OC_index(long double mjd, long double ra, long double dec, long double x, long double y, long double z, string idstring, long index) :MJD(mjd), RA(ra), Dec(dec), x(x), y(y), z(z), idstring(idstring), index(index) { }
};

class xy_index{ // Double-precision x,y point plus long index
public:
  double x;
  double y;
  long index;
  xy_index(double x, double y, long index) :x(x), y(y), index(index) { }
};

class kdpoint { // Point in an xy_index k-d tree designed as a vector,
                // with left and right giving the index of the left and
                // right branches, and the dimension of splitting
                // specified at each node.
public:
  xy_index point;
  long left;
  long right;
  int dim;
  kdpoint(xy_index point, long left, long right, int dim) :point(point), left(left), right(right), dim(dim) {}
};
  
class xyind_lower_x{ // Function for sorting vectors of type xy_index by x
public:
  inline bool operator() (const xy_index& p1, const xy_index& p2) {
    return(p1.x < p2.x);
  }
};

class xyind_lower_y{ // Function for sorting vectors of type xy_index by y
public:
  inline bool operator() (const xy_index& p1, const xy_index& p2) {
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

class early_det_index{
public:
  inline bool operator() (const det_index& o1, const det_index& o2) {
    return(o1.MJD < o2.MJD);
  }
};

class early_det_OC_index{
public:
  inline bool operator() (const det_OC_index& o1, const det_OC_index& o2) {
    return(o1.MJD < o2.MJD);
  }
};

class stringsort_det_OC_index{
public:
  inline bool operator() (const det_OC_index& o1, const det_OC_index& o2) {
    return(o1.idstring < o2.idstring);
  }
};

class early_imlg2{
public:
  inline bool operator() (const img_log02& i1, const img_log02& i2) {
    return(i1.MJD < i2.MJD);
  }
};

class longpair{ // Pair of long integers
public:
  long i1;
  long i2;
  longpair(long i1, long i2) :i1(i1), i2(i2) { }
};

class point3d{ // Double-precision 3-D point
public:
  double x;
  double y;
  double z;
  point3d(double x, double y, double z) :x(x), y(y), z(z) { }
};

class point3LD{ // Long double-precision 3-D point
public:
  long double x;
  long double y;
  long double z;
  point3LD(long double x, long double y, long double z) :x(x), y(y), z(z) { }
};

class point6LDx2{ // Long double-precision 6-D point plus 2 long integer indices
                  // The purpose of the indices is to specify the pair of detections
                  // from which the position and velocity state vectors originated.
public:
  long double x;
  long double y;
  long double z;
  long double vx;
  long double vy;
  long double vz;
  long i1;
  long i2;
  point6LDx2(long double x, long double y, long double z, long double vx, long double vy, long double vz, long i1, long i2) :x(x), y(y), z(z), vx(vx), vy(vy), vz(vz), i1(i1), i2(i2) {}
};

class lower_point6LDx2_x{ // Sort point6LDx2's by x
public:
  inline bool operator() (const point6LDx2& p1, const point6LDx2& p2) {
    return(p1.x < p2.x);
  }
};

class lower_point6LDx2_y{ // Sort point6LDx2's by y
public:
  inline bool operator() (const point6LDx2& p1, const point6LDx2& p2) {
    return(p1.y < p2.y);
  }
};

class lower_point6LDx2_z{ // Sort point6LDx2's by z
public:
  inline bool operator() (const point6LDx2& p1, const point6LDx2& p2) {
    return(p1.z < p2.z);
  }
};

class lower_point6LDx2_vx{ // Sort point6LDx2's by vx
public:
  inline bool operator() (const point6LDx2& p1, const point6LDx2& p2) {
    return(p1.vx < p2.vx);
  }
};

class lower_point6LDx2_vy{ // Sort point6LDx2's by vy
public:
  inline bool operator() (const point6LDx2& p1, const point6LDx2& p2) {
    return(p1.vy < p2.vy);
  }
};

class lower_point6LDx2_vz{ // Sort point6LDx2's by vz
public:
  inline bool operator() (const point6LDx2& p1, const point6LDx2& p2) {
    return(p1.vz < p2.vz);
  }
};


class KD_point6LDx2{ // 6-dimensional KD tree based on points of type point6LDx2.
                     // These points carry two long integer indices not used by
                     // the KD tree, which itself uses two additional long integer
                     // indices to specify the left and right branches, and one
                     // regular integer to specify the dimension of the current
                     // branching. We also provide a final integer, flag, to mark
                     // points as already-checked for DBSCAN or other algorithms.
public:
  point6LDx2 point;
  long left;
  long right;
  int dim;
  int flag;
  KD_point6LDx2(point6LDx2 point, long left, long right, int dim, int flag) :point(point), left(left), right(right), dim(dim), flag(flag) {}
};

class KD6_clust{ // Cluster of points produced by DBSCAN_6D01.
public:
  int numpoints;
  vector <long> clustind;
  vector <long double> meanvec;
  vector <long double> rmsvec;
  KD6_clust(int numpoints, vector <long> clustind, vector <long double> meanvec, vector <long double> rmsvec) :numpoints(numpoints), clustind(clustind), meanvec(meanvec), rmsvec(rmsvec) {}
};

class KD6i_clust{ // Cluster of points produced by DBSCAN_6i01.
public:
  int numpoints;
  vector <long> clustind;
  vector <double> meanvec;
  vector <double> rmsvec;
  KD6i_clust(int numpoints, vector <long> clustind, vector <double> meanvec, vector <double> rmsvec) :numpoints(numpoints), clustind(clustind), meanvec(meanvec), rmsvec(rmsvec) {}
};

class point6ix2{  // integer 6-D point plus 2 long integer indices
                  // The integer 6-D point is supposed to hold integerized
                  // versions of dynamical state vectors (xyz position and velocity).
                  // the purpose of the indices is to specify the pair of detections
                  // from which the position and velocity state vectors originated.
public:
  int x;
  int y;
  int z;
  int vx;
  int vy;
  int vz;
  long i1;
  long i2;
  point6ix2(int x, int y, int z, int vx, int vy, int vz, long i1, long i2) :x(x), y(y), z(z), vx(vx), vy(vy), vz(vz), i1(i1), i2(i2) {}
};

class lower_point6ix2_x{ // Sort point6ix2's by x
public:
  inline bool operator() (const point6ix2& p1, const point6ix2& p2) {
    return(p1.x < p2.x);
  }
};

class lower_point6ix2_y{ // Sort point6ix2's by y
public:
  inline bool operator() (const point6ix2& p1, const point6ix2& p2) {
    return(p1.y < p2.y);
  }
};

class lower_point6ix2_z{ // Sort point6ix2's by z
public:
  inline bool operator() (const point6ix2& p1, const point6ix2& p2) {
    return(p1.z < p2.z);
  }
};

class lower_point6ix2_vx{ // Sort point6ix2's by vx
public:
  inline bool operator() (const point6ix2& p1, const point6ix2& p2) {
    return(p1.vx < p2.vx);
  }
};

class lower_point6ix2_vy{ // Sort point6ix2's by vy
public:
  inline bool operator() (const point6ix2& p1, const point6ix2& p2) {
    return(p1.vy < p2.vy);
  }
};

class lower_point6ix2_vz{ // Sort point6ix2's by vz
public:
  inline bool operator() (const point6ix2& p1, const point6ix2& p2) {
    return(p1.vz < p2.vz);
  }
};


class KD_point6ix2{ // 6-dimensional KD tree based on points of type point6ix2.
                     // These points carry two long integer indices not used by
                     // the KD tree, which itself uses two additional long integer
                     // indices to specify the left and right branches, and one
                     // regular integer to specify the dimension of the current
                     // branching. We also provide a final integer, flag, to mark
                     // points as already-checked for DBSCAN or other algorithms.
public:
  point6ix2 point;
  long left;
  long right;
  int dim;
  int flag;
  KD_point6ix2(point6ix2 point, long left, long right, int dim, int flag) :point(point), left(left), right(right), dim(dim), flag(flag) {}
};


void make_dvec(int nx, vector <double> &dvec);
void make_dmat(int nx, int ny, vector <vector <double>> &dmat);
void make_LDvec(int nx, vector <long double> &ldvec);
void make_LDmat(int nx, int ny, vector <vector <long double>> &ldmat);
double dotprod3d(point3d p1, point3d p2);
double dotprod3LD(point3LD p1, point3LD p2);
point3d crossprod3d(point3d p1, point3d p2);
point3LD crossprod3LD(point3LD p1, point3LD p2);
long double intpowLD(long double x, int p);
long double factorialLD(int p);
point3d celeproj01(double RA, double Dec);
int celedeproj01(point3d p3, double *RA, double *Dec);
point3LD celeproj01LD(long double RA, long double Dec);
int celedeproj01LD(point3LD p3, long double *RA, long double *Dec);
double angspeed01(det_bsc o1, det_bsc o2);
double distradec01(double RA1, double Dec1, double RA2, double Dec2);
int distradec02(double ra1,double dec1,double ra2,double dec2,double *dist,double *pa);
long medindex(const vector <xy_index> &xyvec, int dim);
int splitxy(const vector <xy_index> &xyvec, int dim, long splitpoint, vector <xy_index> &left, vector <xy_index> &right);
int kdtree01(const vector <xy_index> &xyvec, int dim, long rootptxy, long rootptkd, vector <kdpoint> &kdvec);
int kdrange01(const vector <kdpoint> &kdvec,double x,double y,double range,vector <long> &indexvec);
long medind_6LDx2(const vector <point6LDx2> &pointvec, int dim);
int splitLDx2(const vector <point6LDx2> &pointvec, int dim, long splitpoint, vector <point6LDx2> &left, vector <point6LDx2> &right);
int kdtree_6D01(const vector <point6LDx2> &invec, int dim, long splitpoint, long kdroot, vector <KD_point6LDx2> &kdvec);
long double point6LDx2_dist(const point6LDx2 &p1, const point6LDx2 &p2);
long double point6LDx2_dist2(const point6LDx2 &p1, const point6LDx2 &p2);
int kdrange_6D01(const vector <KD_point6LDx2> &kdvec, const point6LDx2 &querypoint, long double range, vector <long> &indexvec);
long double cluster_stats6D01(const vector <KD_point6LDx2> &cluster, vector <long double> &meanvals, vector <long double> &rmsvals);
int DBSCAN_6D01(vector <KD_point6LDx2> &kdtree, long double clustrad, int npt, const vector <det_bsc> &detvec, const vector <string> &det_id_vec, vector <KD6_clust> &outclusters, string rmsfile);
int DBSCAN_6D02(vector <KD_point6LDx2> &kdtree, long double clustrad, int npt, vector <KD6_clust> &outclusters);
point6ix2 conv_6LD_to_6i(point6LDx2 p1, long double scale);
long medind_6ix2(const vector <point6ix2> &pointvec, int dim);
int splitix2(const vector <point6ix2> &pointvec, int dim, long splitpoint, vector <point6ix2> &left, vector <point6ix2> &right);
int kdtree_6i01(const vector <point6ix2> &invec, int dim, long splitpoint, long kdroot, vector <KD_point6ix2> &kdvec);
long point6ix2_dist2(const point6ix2 &p1, const point6ix2 &p2);
int kdrange_6i01(const vector <KD_point6ix2> &kdvec, const point6ix2 &querypoint, long range, vector <long> &indexvec);
double cluster_stats6i01(const vector <KD_point6ix2> &cluster, double intconvscale, vector <double> &meanvals, vector <double> &rmsvals);
int DBSCAN_6i01(vector <KD_point6ix2> &kdtree, double clustrad, int npt, double intconvscale, vector <KD6i_clust> &outclusters);
int celestial_to_statevec(double RA, double Dec,double delta,point3d &baryvec);
int celestial_to_statevecLD(long double RA, long double Dec,long double delta,point3LD &baryvec);
int celestial_to_stateunit(double RA, double Dec,point3d &baryvec);
int celestial_to_stateunitLD(long double RA, long double Dec, point3LD &baryvec);
int read_horizons_file(string infile, vector <double> &mjdvec, vector <point3d> &pos, vector <point3d> &vel);
int read_horizons_fileLD(string infile, vector <long double> &mjdvec, vector <point3LD> &pos, vector <point3LD> &vel);
int poleswitch01(const double &inRA, const double &inDec, const double &poleRA, const double &poleDec, const double &oldpoleRA, double *newRA, double *newDec);
int poleswitch01LD(const long double &inRA, const long double &inDec, const long double &poleRA, const long double &poleDec, const long double &oldpoleRA, long double *newRA, long double *newDec);
int precess01a(double ra1,double dec1,double mjd,double *ra2,double *dec2,int precesscon);
int precess01aLD(long double ra1,long double dec1,long double mjd,long double *ra2,long double *dec2,int precesscon);
int solvematrix01(const vector <vector <double>> &inmat, int eqnum, vector <double> &outvec, int verbose);
int solvematrix01LD(const vector <vector <long double>> &inmat, int eqnum, vector <long double> &outvec, int verbose);
int perfectpoly01(const vector <double> &x, const vector <double> &y, vector <double> &fitvec);
int perfectpoly01LD(const vector <long double> &x, const vector <long double> &y, vector <long double> &fitvec);
int planetpos01(double detmjd, int polyorder, const vector <double> &posmjd, const vector <point3d> &planetpos, point3d &outpos);
int planetpos01LD(long double detmjd, int polyorder, const vector <long double> &posmjd, const vector <point3LD> &planetpos, point3LD &outpos);
int planetposvel01LD(long double detmjd, int polyorder, const vector <long double> &posmjd, const vector <point3LD> &planetpos, const vector <point3LD> &planetvel, point3LD &outpos, point3LD &outvel);
int nplanetpos01LD(long double detmjd, int planetnum, int polyorder, const vector <long double> &posmjd, const vector <point3LD> &planetpos, vector <point3LD> &outpos);
int nplanetgrab01LD(int pointrequest, int planetnum, const vector <long double> &posmjd, const vector <point3LD> &planetpos, vector <point3LD> &outpos);
int observer_barycoords01(double detmjd, int polyorder, double lon, double obscos, double obssine, const vector <double> &posmjd, const vector <point3d> &planetpos, point3d &outpos);
int observer_barycoords01LD(long double detmjd, int polyorder, long double lon, long double obscos, long double obssine, const vector <long double> &posmjd, const vector <point3LD> &planetpos, point3LD &outpos);
int helioproj01(point3d unitbary, point3d obsbary,double heliodist,double &geodist, point3d &projbary);
int helioproj01LD(point3LD unitbary, point3LD obsbary, long double heliodist, long double &geodist, point3LD &projbary);
int accelcalc01LD(int planetnum, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const point3LD &targpos, point3LD &accel);
int integrate_orbit_constac(int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, long double mjdstart, point3LD startpos, point3LD startvel, long double mjdend, point3LD &endpos, point3LD &endvel);
long double kep_transcendental(long double q, long double e, long double tol);
int Keplerint(const long double MGsun, const long double mjdstart, const point3LD &startpos, const point3LD &startvel, const long double mjdend, point3LD &endpos, point3LD &endvel);
int integrate_orbit01LD(int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, long double mjdstart, point3LD startpos, point3LD startvel, long double mjdend, point3LD &endpos, point3LD &endvel);
int integrate_orbit02LD(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, long double mjdstart, point3LD startpos, point3LD startvel, long double mjdend, point3LD &endpos, point3LD &endvel);
int iswhitespace(int c);
int readconfigLD(ifstream &instream1, long double *ldval);
int readconfigd(ifstream &instream1, double *dval);
int readconfigint(ifstream &instream1, int *ival);
int readconfigstring(ifstream &instream1, string &sval);
int read_accel_fileLD(string accelfile, vector <long double> &heliodist, vector <long double> &heliovel, vector <long double> &helioacc);
long double weight_posvel_rms(const vector <point3LD> &poscluster,const vector <point3LD> &velcluster,const long double dtime, vector <long double> &rmsvec);


