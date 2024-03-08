// December 07, 2021: solarsyst_dyn_geo01.h:
// Library of functions useful for Solar System dynamics and
// geometry. Includes functions originally developed in orbint02a.cpp,
// maketrack02b.cpp, projtest01b.cpp, and other places.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cstring>
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
#include <cassert>

#include <omp.h>

using namespace std;

#define DEGPRAD (180.0L/M_PI) /*Degrees per radian*/
#define ASECPRAD (648000.0L/M_PI) /*Arcseconds per radian*/
#define LSQUARE(x) (long(x)*long(x))
#define DSQUARE(x) (double(x)*double(x))
#define LDSQUARE(x) (((long double)x)*((long double)x))

#define SOLARDAY 86400.0L
#define SIDEREALDAY 86164.09053L
#define NEPRA 270.0L //Right ascension of the North Ecliptic Pole.
#define NEPDEC 66.560708333333L // Declination of the North Ecliptic Pole
#define NGPRA 192.728L //Right ascension of the North Galactic Pole (Karim and Mamajek 2017)
#define NGPDEC 26.863L // Declination of the North Galactic Pole (Karim and Mamajek 2017)
#define NCPGAL_LON 122.928L // Galactic longitude of the North Celestial Pole (Karim and Mamajek 2017)
#define MJDOFF 2400000.5L // Offset from Julian Days to Modified Julian days.
#define EARTHEQUATRAD 6378.140L // Earth's equatorial radius.
#define AU_KM 1.495978700e8L /*One AU in km*/
#define AU 1.49597870700e11L /*One AU in meters*/
#define EFOLDS_PER_MAG 0.921034037197618 // E-foldings per astronomical magnitude
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
#define VSMALLANG 1.0e-9L /*angle difference below which angles are collapsed
                            to zero in distradec routines.*/

#define TTDELTAT 70.0L // Time in seconds by which TT differs from UT. Thus,
                      // a light-travel-time corrected dynamical position at
                      // TT =  T + TTDELTAT corresponds to the celestial
                      // coordinates observed at UT = T.*/
#define CLIGHT 2.99792458e8L // Speed of light in m/sec
#define CLIGHT_AUDAY (CLIGHT*SOLARDAY/AU) // Speed of light in AU/day
#define KEPTRANSITMAX 50 // Maximum number of iterations to use in Newton's
                         // method solution of the trancendental Kepler Equation.
#define KVSTAR_ITMAX 10 // Maximum number of iterations for solving the
                        // Kepler two-point boundary value problem with
                        // universal variables.
#define KEPTRANSTOL 1e-15L // Maximum error for an acceptable solution of the
                           // trancendental Kepler Equation.
#define KEPTRANSTOL2 (double(1e-10)) // Maximum error for an acceptable solution of the
                            // trancendental Kepler Equation.
#define HYPTRANSTOL 1e-13l // Maximum error for an acceptable solution of the
                           // Kepler Equation with universal variables
#define KEP2PBVPTOL 1e-12L // Maximum error for an acceptable solution of the
                           // Kepler two-point boundary value problem.
#define KEP2PBVPTOL2 (double(1e-8l)) // Maximum error for an acceptable solution of the
                           // Kepler two-point boundary value problem.
#define QSTAR_XTOL 1e-4l  // Transition to use series rather than analytical forms
                         // for the Q(x*) function from Danby 6.11
#define QSTAR_N  10 // Number of terms to calculate for the series expansion 
                   // of the Q(x*) function from Danby 6.11
#define LARGERR 1e30L // Large number supposed to be a safe initialization
                      // for most minimum-finding problems.
#define LARGERR2 (double(1e30l)) // Large number supposed to be a safe initialization
                                // for most minimum-finding problems.
#define LARGERR3 1e60L // Large number supposed to be a safe initialization
                       // for most minimum-finding problems.
#define GMSUN_KM3_SEC2 132712440041.279419L // GM for the Sun: that is, the Universal
                                            // Gravitational Constant times the solar mass,
                                            // in units of km^3/sec^2. 
#define KCONST 0.0172020989484485L // Mean daily motion in radians for an orbit around
                                   // the sun with a semimajor axis of 1AU.
#define MINSTRINGLEN 5 // Minimum size of character array we use: e.g., for filter bandpass or obscode.
#define SHORTSTRINGLEN 20 // Standard size for a short-ish string, used, e.g. for detection idstring
#define MEDSTRINGLEN 80 // Medium string length, should hold most file paths
#define LONGSTRINGLEN 256 // Should hold any reasonable file path, use if not pressed for memory.
#define RAND_MAX_64 18446744073709551616.0L

#define PHASECONST_A1 3.33L // Constants from H-G system phase equation
#define PHASECONST_A2 1.87L // for calculating asteroid apparent magnitudes.
#define PHASECONST_B1 0.63L
#define PHASECONST_B2 1.22L
#define WARN_INVERSE_TRIG 0 // If 0, don't print warning messages about taking arcsin or arccos of values
                            // or arccos of values infinitesimally outside the range [-1,1]. Just collapse
                            // them to exactly 1.0 or -1.0 like a good robot.
#define HERGET_DOWNSCALE 0.9L // If your Herget orbit fit is leading to hyperbolic orbits,
                              // shrink the distances by this factor to see if you can bring them in range.
#define MINHERGETDIST 0.0001L // Make sure distances and distance separations in Herget
                              // orbit fitting don't get smaller than this value (in AU).
#define SIMPLEX_SCALEFAC 0.2L // Initial scaling for Herget simplex. For input guess x, y points will be
                              // x ( 1 + SIMPLEX_SCALEFAC - SIMPLEX_SCALEFACE^2),
                              // y ( 1 + SIMPLEX_SCALEFAC + SIMPLEX_SCALEFACE^2) and
                              // x ( 1 - SIMPLEX_SCALEFAC + SIMPLEX_SCALEFACE^2),
                              // y ( 1 - SIMPLEX_SCALEFAC - SIMPLEX_SCALEFACE^2),
                              // The point here is to create a guaranteed linearly independent
                              // simplex where the distances nevertheless tend to vary together,
                              // reducing the chance of hyperbolic orbits.
#define SIMPLEX_SCALE_LIMIT ((sqrt(5.0L)-1.0L)/2.0L) // Can get negative distances if simplex_scale
                                                     // exceeds this value.
#define LOG10_e 0.434294481903252L
#define LN10 2.30258509299405L
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds

// Parameters related to heliolinc clustering
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
#define MAX_SHUTTER_CORR 10.0 // Implied shutter corrections larger than this value,
                              // in seconds, are implausible and will cause link_refine_Herget
                              // to exit with an error.
// End parameters related to heliolinc clustering


// String-handling stuff that has to be declared early because other things depend on it.
void stringncopy01(char *dest, const string &source, int n);
int stringnmatch01(const char *string1, const char *string2, int n);

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
  double RA;
  double Dec;
  long double x;
  long double y;
  long double z;
  string idstring;
  long index;
  det_OC_index(long double mjd, double ra, double dec, long double x, long double y, long double z, string idstring, long index) :MJD(mjd), RA(ra), Dec(dec), x(x), y(y), z(z), idstring(idstring), index(index) { }
};

class det_OC_indvec{ // Detection of some astronomical source with
                    // the observer's X, Y, Z coordinates, a string
                    // identifier, an index for cross-referencing,
                    // and an integer vector intended to hold the
                    // indices of other detections with which the
                    // current one has been paired.
public:
  long double MJD;
  double RA;
  double Dec;
  long double x;
  long double y;
  long double z;
  string idstring;
  long index;
  vector <int> indvec;
  det_OC_indvec(long double mjd, double ra, double dec, long double x, long double y, long double z, string idstring, long index, vector <int> indvec) :MJD(mjd), RA(ra), Dec(dec), x(x), y(y), z(z), idstring(idstring), index(index) , indvec(indvec) { }
};


class det_obsmag_indvec{ // Detection of some astronomical source with
                         // the observer's X, Y, Z coordinates, a string
                         // identifier, the magnitude and photometric band,
                         // the observatory code, an index for cross-referencing,
                         // and an integer vector intended to hold the
                         // indices of other detections with which the
                         // current one has been paired.
                         // Note well that the calling function should ensure
                         // the input idstring, band, and obcode strings are
                         // not too long.
public:
  double MJD;
  double RA;
  double Dec;
  double x;
  double y;
  double z;
  char idstring[SHORTSTRINGLEN];
  double mag;
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  long index;
  vector <int> indvec;
  det_obsmag_indvec(long double mjd, double ra, double dec, double x, double y, double z, const string &idstring, double mag, const string &band, const string &obscode, long index, vector <int> indvec) :MJD(mjd), RA(ra), Dec(dec), x(x), y(y), z(z), mag(mag), index(index) , indvec(indvec) {
    // Copy input value of idstring, making sure it's not too long
    assert(idstring.size() < sizeof(this->idstring));
    std::strncpy(this->idstring, idstring.c_str(), sizeof(this->idstring));
    this->idstring[sizeof(this->idstring)-1] = 0;

    // Copy input value for photometric band, making sure it's not too long
    assert(band.size() < sizeof(this->band));
    std::strncpy(this->band, band.c_str(), sizeof(this->band));
    this->band[sizeof(this->band)-1] = 0;
    
    // Copy input value for obscode, making sure it's not too long
    assert(obscode.size() < sizeof(this->obscode));
    std::strncpy(this->obscode, obscode.c_str(), sizeof(this->obscode));
    this->obscode[sizeof(this->obscode)-1] = 0;
  }
  det_obsmag_indvec() = default;
};


class hldet{ // Detection of some astronomical source with
             // info required for new version of heliolinc:
public:
  double MJD;
  double RA;
  double Dec;
  float mag;
  float trail_len;  // In arcseconds. Use zero unless object is certainly trailed.
  float trail_PA;   // In degrees. Default to 90 deg for un-trailed sources.
  float sigmag;
  float sig_across; // Astrometric uncertainty cross-trail (or in RA for untrailed sources), in arcsec
  float sig_along;  // Astrometric uncertainty cross-trail (or in Dec for untrailed sources), in arcsec
  int image;        // Vitally important index to image catalog.
  char idstring[SHORTSTRINGLEN];  // E.g., diaSourceID
  char band[MINSTRINGLEN];        // photometric band
  char obscode[MINSTRINGLEN];     // Observatory code
  long known_obj;                 // Identifier for a known object, or prob. of being a known object (e.g. ATLAS Pkn)
  long det_qual;                  // Detection quality -- e.g., a metric for the likelihood that this detection is real.
  long index;                     // Line count in original input structure
  hldet(double mjd, double ra, double dec, float mag, float trail_len, float trail_PA, float sigmag, float sig_across, float sig_along, int image, const string &idstring, const string &band, const string &obscode, long known_obj, long det_qual, long index) :MJD(mjd), RA(ra), Dec(dec), mag(mag), trail_len(trail_len), trail_PA(trail_PA), sigmag(sigmag), sig_across(sig_across), sig_along(sig_along), image(image), known_obj(known_obj), det_qual(det_qual), index(index) {
    // Copy input value of idstring, making sure it's not too long
    assert(idstring.size() < sizeof(this->idstring));
    std::strncpy(this->idstring, idstring.c_str(), sizeof(this->idstring));
    this->idstring[sizeof(this->idstring)-1] = 0;

    // Copy input value for photometric band, making sure it's not too long
    assert(band.size() < sizeof(this->band));
    std::strncpy(this->band, band.c_str(), sizeof(this->band));
    this->band[sizeof(this->band)-1] = 0;
    
    // Copy input value for obscode, making sure it's not too long
    assert(obscode.size() < sizeof(this->obscode));
    std::strncpy(this->obscode, obscode.c_str(), sizeof(this->obscode));
    this->obscode[sizeof(this->obscode)-1] = 0;
  }
  hldet() = default;
};

class hlimage{ // Astronomical image with MJD of mid-exposure, boresight RA and Dec,
               // observatory code, and observer X, Y, Z, VX, VY, VZ, and exposure time in seconds
public:
  double MJD;
  double RA;
  double Dec;
  char obscode[MINSTRINGLEN];     // Observatory code
  double X;
  double Y;
  double Z;
  double VX;
  double VY;
  double VZ;
  long startind;
  long endind;
  double exptime; // Exposure time in seconds
  hlimage(double mjd, double ra, double dec, const string &obscode, double x, double y, double z, double vx, double vy, double vz, long startind, long endind, double exptime) :MJD(mjd), RA(ra), Dec(dec), X(x), Y(y), Z(z), VX(vx), VY(vy), VZ(vz), startind(startind), endind(endind), exptime(exptime) {
    // Copy input value for obscode, making sure it's not too long
    assert(obscode.size() < sizeof(this->obscode));
    std::strncpy(this->obscode, obscode.c_str(), sizeof(this->obscode));
    this->obscode[sizeof(this->obscode)-1] = 0;
  }
  hlimage() = default;
};

struct MakeTrackletsConfig {
  int mintrkpts = 2;
  double imagetimetol = IMAGETIMETOL/SOLARDAY;    // Tolerance for matching image time, in days
  double maxvel = 1.5l;         // Default max angular velocity in deg/day.
  double minvel = 0.0l;         // Min angular velocity in deg/day
  double minarc = 0.0l;         // Min total angular arc in arcseconds
  double maxtime = 1.5 / 24.0;  // Default max inter-image time interval
                                // for tracklets, in days.
  double mintime = 1.0 / SOLARDAY; // Minimum inter-image time interval, in days.
  double imagerad = 2.0;        // radius from image center to most distant corner (deg)
  double maxgcr = 0.5;          // Default maximum Great Circle Residual allowed for a valid tracklet
  double exptime = 30.0;        // Defulat exposure time in seconds.
  double siglenscale = 0.5;     // Default scaling from trail length to its uncertainty
  double sigpascale = 1.0;      // Default size in arcseconds used to estimate trail
                                // PA uncertainty, using sigPA = DEGPRAD*sigpascale/trail_len
  int forcerun = 0; // Pushes through all but the immediately fatal errors.
  int verbose = 0;  // Prints monitoring output
};

struct HeliolincConfig {
  double MJDref = 0.0l;        // MJD of reference time. Must be specified at runtime:
                               //   no sensible default is possible, because the reference
                               //   time MUST be near the center of the time spanned by
                               //   the input detection catalog.
  double clustrad = 1.0e5l;    // Clustering radius (km) at the reference geocentric 
                               // distance, REF_GEODIST, which is declared as 1.0 AU in
                               // solarsyst_dyn_geo01.cpp, as of 2024 March 08.
                               // The clustering radius scales linearly with geocentric
                               // distance, except inward of clustchangerad (see below).
  double clustchangerad = 0.5l; // Geocentric distance (AU), within which the clustering
                                // radius remains fixed at a minimum value and does not
                                // continue to decrease for samller geocentric distances.
  int dbscan_npt = 3;            // Number of points npt for the DBSCAN algorithm
  int minobsnights = 3;        // Minimum number of distinct observing nights for a valid linkage
  double mintimespan = 1.0l;   // Minimum timespan for a valid linkage
  double mingeodist = 0.10l;   // Geocentric distance (AU) at the center of the innermost
                               //   geocentric distance bin.
  double maxgeodist = 100.0l;    // Minimum value in AU for the center of the outermost distance bin
  double geologstep = 1.5l;      // Factor by which the center of a geocentric distance bin
                               //   moves outward from one bin to the next: e.g. 0.1, 0.15, 0.225, 0.338.
                               //   This is also the logarithmic width of each bin: e.g., the bin centered
                               //   at 0.15 AU ranges from 0.1 to 0.225 AU.
  double mingeoobs = 0.0l;     // Minimum inferred geocentric distance for a valid tracklet.
  double minimpactpar = 0.0l;  // Minimum inferred impact parameter for a valid tracklet.
                               // These last two parameters can be raised above the default
                               //   value of zero to limit the formation of 'globs' -- enormous
                               //   spurious linkages -- when probing heliocentric hypothesis
                               //   that have the potential to imply near-zero distance from Earth
                               //   for some observations.
  int use_univar = 0;          // Use the universal variable formulation of the Kepler equations,
                               // rather than the default fg function formulation.
  double max_v_inf = 0.0l;     // Maximum value of v_infinity relative to the sun. Setting this
                               // to zero restricts us to bound orbits. To find interstellar objects,
                               // we set it to positive values. In this case it is also necessary
                               // to set use_univar=1, since only the universal variable formulation
                               // can handle unbound orbits.
  int verbose=0;
};

struct HeliovaneConfig {
  double MJDref = 0.0l;        // MJD of reference time. Must be specified at runtime:
                               //   no sensible default is possible, because the reference
                               //   time MUST be near the center of the time spanned by
                               //   the input detection catalog.
  double clustrad = 1.0e5l;    // Clustering radius for the DBSCAN algorithm, in km.
  int dbscan_npt = 3;            // Number of points npt for the DBSCAN algorithm
  int minobsnights = 3;        // Minimum number of distinct observing nights for a valid linkage
  double mintimespan = 1.0l;   // Minimum timespan for a valid linkage
  double mingeodist = 0.05l;   // Geocentric distance (AU) at the center of the innermost
                               //   geocentric distance bin.
  double maxgeodist = 2.0l;    // Minimum value in AU for the center of the outermost distance bin
  double geologstep = 1.5l;    // Factor by which the center of a geocentric distance bin
                               //   moves outward from one bin to the next: e.g. 0.1, 0.15, 0.225, 0.338.
                               //   This is also the logarithmic width of each bin: e.g., the bin centered
                               //   at 0.15 AU ranges from 0.1 to 0.225 AU.
  double mingeoobs = 0.0l;     // Minimum inferred geocentric distance for a valid tracklet.
  double minimpactpar = 0.0l;  // Minimum inferred impact parameter for a valid tracklet.
                               // These last two parameters can be raised above the default
                               //   value of zero to limit the formation of 'globs' -- enormous
                               //   spurious linkages -- when probing heliocentric hypothesis
                               //   that have the potential to imply near-zero distance from Earth
                               //   for some observations.
  int use_univar = 0;          // Use the universal variable formulation of the Kepler equations,
                               // rather than the default fg function formulation.
  double minsunelong = 0.0l;  // Minimum solar elongation in degrees
  double maxsunelong = 180.0l;  // Maximum solar elongation in degrees
  double min_incid_angle = 40.0l; // Minimum value for the angle at which the observer-to-target unitvector
                                 // intersects the heliocentric vane, in degrees.
  double maxheliodist = 2.0l; // Maximum heliocentric distance in AU. Partially redundant with
                               // min_proj_sine.
  double max_v_inf = 0.0l;     // Maximum value of v_infinity relative to the sun. Setting this
                               // to zero restricts us to bound orbits. To find interstellar objects,
                               // we set it to positive values. In this case it is also necessary
                               // to set use_univar=1, since only the universal variable formulation
                               // can handle unbound orbits.
  int verbose=0;
};

struct LinkRefineConfig {
  double MJDref = 0.0l;        // MJD of reference time. Must be specified at runtime:
                               //   no sensible default is possible, because the reference
                               //   time MUST be near the center of the time spanned by
                               //   the input detection catalog.
  int simptype = 0;            // Defines how Herget_simplex_int constructs the starting simplex
                               //   in the 2-D parameter space of geocentric distance at first detection
                               //   (geodist1) and geocentric distance at last detection (geodist2)
                               //   simptype=0 uses multiplicative scaling to create an approximately equilateral triangle:
                               //   simptype=1 creates a simplex elongated along the direction defined by geodist1=geodist2
                               //   simptype=2 uses subtraction to create a precisely equilateral triangle
  int ptpow = 1;               // Power to which we raise the number of unique detections, when calculating the cluster quality metric.
  int nightpow = 1;            // Power to which we raise the number of distinct observing nights, when calculating the cluster quality metric.
  int timepow = 0;             // Power to which we raise the total temporal span of the linkage, when calculating the cluster quality metric.
  int rmspow = 2;              // Power to which we raise the RMS astrometric residual when calculating the cluster quality metric.
  double maxrms = 1.0e5;       // Maximum scaled RMS in km for a viable cluster
  int verbose=0;
};

struct LinkPurifyConfig {
  double MJDref = 0.0l;        // MJD of reference time. Must be specified at runtime:
                               //   no sensible default is possible, because the reference
                               //   time MUST be near the center of the time spanned by
                               //   the input detection catalog.
  int simptype = 0;            // Defines how Herget_simplex_int constructs the starting simplex
                               //   in the 2-D parameter space of geocentric distance at first detection
                               //   (geodist1) and geocentric distance at last detection (geodist2)
                               //   simptype=0 uses multiplicative scaling to create an approximately equilateral triangle:
                               //   simptype=1 creates a simplex elongated along the direction defined by geodist1=geodist2
                               //   simptype=2 uses subtraction to create a precisely equilateral triangle
  int ptpow = 1;               // Power to which we raise the number of unique detections, when calculating the cluster quality metric.
  int nightpow = 1;            // Power to which we raise the number of distinct observing nights, when calculating the cluster quality metric.
  int timepow = 0;             // Power to which we raise the total temporal span of the linkage, when calculating the cluster quality metric.
  int rmspow = 2;              // Power to which we raise the RMS astrometric residual when calculating the cluster quality metric.
  double maxrms = 1.0e5;       // Maximum scaled RMS in km for a viable cluster
  double max_oop = 1000.0;     // Maximum scaled out-of-plane RMS in km for a viable cluster
  double rejfrac = 0.5;        // Maximum fraction of points that can be rejected.
  double max_astrom_rms = 1.0; // Maximum RMS astrometric residual, in arcsec.
  int minobsnights = 3;        // Minimum number of distinct observing nights for a valid linkage
  int minpointnum = 6;         // Minimum number of individual detections for a valid linkage
  int use_heliovane = 0;       // Are we analyzing data from heliovane rather than heliolinc?
  int verbose=0;
};

class observatory{ // observatory code and location values.
public:
  char obscode[MINSTRINGLEN];
  double obslon;
  double plxcos;
  double plxsin;
  observatory(const string &obscode, double obslon, double plxcos, double plxsin) :obslon(obslon), plxcos(plxcos), plxsin(plxsin) {
    // Copy input value of obscode, making sure it's not too long
    assert(obscode.size() < sizeof(this->obscode));
    std::strncpy(this->obscode, obscode.c_str(), sizeof(this->obscode));
    this->obscode[sizeof(this->obscode)-1] = 0;
  }
  observatory() = default;
};

class tracklet{ // Pair or tracklet for heliolinc 
public:
  long Img1;
  double RA1;
  double Dec1;
  long Img2;
  double RA2;
  double Dec2;
  int npts;
  long trk_ID;
  tracklet(long img1, double ra1, double dec1, long img2, double ra2, double dec2, int npts, long trk_id) :Img1(img1), RA1(ra1), Dec1(dec1), Img2(img2), RA2(ra2), Dec2(dec2), npts(npts), trk_ID(trk_id) { }
  tracklet() = default;
};

class hlradhyp{ // Heliolinc heliocentric radial motion hypothesis
public:
  double HelioRad;
  double R_dot;
  double R_dubdot;
  hlradhyp(double heliorad, double rdot, double rdubdot) :HelioRad(heliorad), R_dot(rdot), R_dubdot(rdubdot) { }
  hlradhyp() = default;
};

class hlclust{ // cluster (linkage) type used as output by heliolinc and I/O by link_refine
public:
  long clusternum;
  double posRMS;
  double velRMS;
  double totRMS;
  double astromRMS;
  int pairnum;
  double timespan;
  int uniquepoints;
  int obsnights;
  double metric;
  char rating[SHORTSTRINGLEN];
  double heliohyp0;
  double heliohyp1;
  double heliohyp2;
  double posX;
  double posY;
  double posZ;
  double velX;
  double velY;
  double velZ;
  double orbit_a;
  double orbit_e;
  double orbit_MJD;
  double orbitX;
  double orbitY;
  double orbitZ;
  double orbitVX;
  double orbitVY;
  double orbitVZ;
  long orbit_eval_count;
  hlclust(long clusternum, double posRMS, double velRMS, double totRMS, double astromRMS, int pairnum, double timespan, int uniquepoints, int obsnights, double metric, const string &rating, double heliohyp0, double heliohyp1, double heliohyp2, double posX, double posY, double posZ, double velX, double velY, double velZ, double orbit_a, double orbit_e, double orbit_MJD, double orbitX, double orbitY, double orbitZ, double orbitVX, double orbitVY, double orbitVZ, long orbit_eval_count) :clusternum(clusternum), posRMS(posRMS), velRMS(velRMS), totRMS(totRMS), astromRMS(astromRMS), pairnum(pairnum), timespan(timespan), uniquepoints(uniquepoints), obsnights(obsnights), metric(metric), heliohyp0(heliohyp0), heliohyp1(heliohyp1), heliohyp2(heliohyp2), posX(posX), posY(posY), posZ(posZ), velX(velX), velY(velY), velZ(velZ), orbit_a(orbit_a), orbit_e(orbit_e), orbit_MJD(orbit_MJD), orbitX(orbitX), orbitY(orbitY), orbitZ(orbitZ), orbitVX(orbitVX), orbitVY(orbitVY), orbitVZ(orbitVZ), orbit_eval_count(orbit_eval_count) { 
    // Copy input value of rating, making sure it's not too long
    assert(rating.size() < sizeof(this->rating));
    std::strncpy(this->rating, rating.c_str(), sizeof(this->rating));
    this->rating[sizeof(this->rating)-1] = 0;
  }
  hlclust() = default;
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

class img_log03{ // Image log that indexes a time-sorted detection vector including observatory code:
                 // MJD, RA, Dec, obscode, starting index, ending index
public:
  double MJD;
  double RA;
  double Dec;
  char obscode[MINSTRINGLEN];
  long startind;
  long endind;
  img_log03(double mjd, double ra, double dec, const string &obscode, long startind, long endind) :MJD(mjd), RA(ra), Dec(dec), startind(startind), endind(endind) {
    // Copy input value of obscode, making sure it's not too long
    assert(obscode.size() < sizeof(this->obscode));
    std::strncpy(this->obscode, obscode.c_str(), sizeof(this->obscode));
    this->obscode[sizeof(this->obscode)-1] = 0;
  }
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

class early_det_OC_indvec{
public:
  inline bool operator() (const det_OC_indvec& o1, const det_OC_indvec& o2) {
    return(o1.MJD < o2.MJD);
  }
};

class stringsort_det_OC_index{
public:
  inline bool operator() (const det_OC_index& o1, const det_OC_index& o2) {
    return(o1.idstring < o2.idstring);
  }
};

class early_det_obsmag_indvec{
public:
  inline bool operator() (const det_obsmag_indvec& o1, const det_obsmag_indvec& o2) {
    return(o1.MJD < o2.MJD || (o1.MJD==o2.MJD && stringnmatch01(o1.obscode,o2.obscode,3)==-1) || (o1.MJD==o2.MJD && stringnmatch01(o1.obscode,o2.obscode,3)==0 && o1.RA<o2.RA));
  }
};

class early_hldet{
public:
  inline bool operator() (const hldet& o1, const hldet& o2) {
    return(o1.MJD < o2.MJD || (o1.MJD==o2.MJD && stringnmatch01(o1.obscode,o2.obscode,3)==-1) || (o1.MJD==o2.MJD && stringnmatch01(o1.obscode,o2.obscode,3)==0 && o1.RA<o2.RA));
  }
};

class early_imlg2{
public:
  inline bool operator() (const img_log02& i1, const img_log02& i2) {
    return(i1.MJD < i2.MJD);
  }
};

class early_imlg3{
public:
  inline bool operator() (const img_log03& i1, const img_log03& i2) {
    return(i1.MJD < i2.MJD || (i1.MJD == i2.MJD && stringnmatch01(i1.obscode,i2.obscode,3)==-1));
  }
};

class early_hlimage{
public:
  inline bool operator() (const hlimage& i1, const hlimage& i2) {
    return(i1.MJD < i2.MJD || (i1.MJD == i2.MJD && stringnmatch01(i1.obscode,i2.obscode,3)==-1));
  }
};

class longpair{ // Pair of long integers
public:
  long i1;
  long i2;
  longpair(long i1, long i2) :i1(i1), i2(i2) { }
  longpair() = default;
};

class point2d{ // Double-precision 2-D point
public:
  double x;
  double y;
  point2d(double x, double y) :x(x), y(y) { }
};

class point3d{ // Double-precision 3-D point
public:
  double x;
  double y;
  double z;
  point3d(double x, double y, double z) :x(x), y(y), z(z) { }
};

class point3d_index{ // Double-precision 3-D point with long-integer idex
public:
  double x;
  double y;
  double z;
  long index;
  point3d_index(double x, double y, double z, long index) :x(x), y(y), z(z), index(index) { }
};

class lower_point3d_index_x{ // Sort point3d_index by x
public:
  inline bool operator() (const point3d_index& p1, const point3d_index& p2) {
    return(p1.x < p2.x);
  }
};

class lower_point3d_index_y{ // Sort point3d_index by y
public:
  inline bool operator() (const point3d_index& p1, const point3d_index& p2) {
    return(p1.y < p2.y);
  }
};

class lower_point3d_index_z{ // Sort point3d_index by z
public:
  inline bool operator() (const point3d_index& p1, const point3d_index& p2) {
    return(p1.z < p2.z);
  }
};

class point4d_index{ // Double-precision 4-D point with long-integer idex
public:
  double t;
  double x;
  double y;
  double z;
  long index;
  point4d_index(double t, double x, double y, double z, long index) :t(t), x(x), y(y), z(z), index(index) { }
};

class lower_point4d_index_t{ // Sort point4d_index by t
public:
  inline bool operator() (const point4d_index& p1, const point4d_index& p2) {
    return(p1.t < p2.t);
  }
};

class lower_point4d_index_x{ // Sort point4d_index by x
public:
  inline bool operator() (const point4d_index& p1, const point4d_index& p2) {
    return(p1.x < p2.x);
  }
};

class lower_point4d_index_y{ // Sort point4d_index by y
public:
  inline bool operator() (const point4d_index& p1, const point4d_index& p2) {
    return(p1.y < p2.y);
  }
};

class lower_point4d_index_z{ // Sort point4d_index by z
public:
  inline bool operator() (const point4d_index& p1, const point4d_index& p2) {
    return(p1.z < p2.z);
  }
};

class point3LD{ // Long double-precision 3-D point
public:
  long double x;
  long double y;
  long double z;
  point3LD(long double x, long double y, long double z) :x(x), y(y), z(z) { }
};

class point6dx2{ // Double-precision 6-D point plus 2 long integer indices
                  // The purpose of the indices is to specify the pair of detections
                  // from which the position and velocity state vectors originated.
public:
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  long i1;
  long i2;
  point6dx2(double x, double y, double z, double vx, double vy, double vz, long i1, long i2) :x(x), y(y), z(z), vx(vx), vy(vy), vz(vz), i1(i1), i2(i2) {}
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

class KD_point3d_index{ // 3-dimensional KD tree based on points of type point point3d_index.
                     // These points carry a long integer index not used by
                     // the KD tree, which itself uses two additional long integer
                     // indices to specify the left and right branches, and one
                     // regular integer to specify the dimension of the current
                     // branching. We also provide a final integer, flag, to mark
                     // points as already-checked for DBSCAN or other algorithms.
public:
  point3d_index point;
  long left;
  long right;
  int dim;
  int flag;
  KD_point3d_index(point3d_index point, long left, long right, int dim, int flag) :point(point), left(left), right(right), dim(dim), flag(flag) {}
};

class KD_point4d_index{ // 4-dimensional KD tree based on points of type point point4d_index.
                     // These points carry a long integer index not used by
                     // the KD tree, which itself uses two additional long integer
                     // indices to specify the left and right branches, and one
                     // regular integer to specify the dimension of the current
                     // branching. We also provide a final integer, flag, to mark
                     // points as already-checked for DBSCAN or other algorithms.
public:
  point4d_index point;
  long left;
  long right;
  int dim;
  int flag;
  KD_point4d_index(point4d_index point, long left, long right, int dim, int flag) :point(point), left(left), right(right), dim(dim), flag(flag) {}
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

class dumbdet { // Dumb representation of a detection, needed only
                // for ancluster01a.cpp to fix idiocy of projectpairs04c.cpp
                // not writing out detection indices -- to be superseded as
                // soon as projectpairs04c.cpp is fixed.
public:
  double mjd;
  string detid;
  dumbdet(double mjd, string detid) :mjd(mjd), detid(detid) { }
};
  
class early_dumbdet{
public:
  inline bool operator() (const dumbdet& o1, const dumbdet& o2) {
    return(o1.mjd < o2.mjd || (o1.mjd == o2.mjd && o1.detid < o2.detid) );
  }
};

class det_svec { // Detection with string ID, integer idex, and vector of
                 // indices. The vector is intended to contain indices for
                 // all the clusters of which this detection is a memeber.
public:
  double mjd;
  double RA;
  double Dec;
  string detid;
  int index;
  vector <int> indvec;
  det_svec(double mjd, double RA, double Dec, string detid, int index, vector <int> indvec) :mjd(mjd), RA(RA), Dec(Dec), detid(detid), index(index), indvec(indvec) { }
};

class early_det_svec{
public:
  inline bool operator() (const det_svec& o1, const det_svec& o2) {
    return(o1.mjd < o2.mjd || (o1.mjd == o2.mjd && o1.detid < o2.detid) );
  }
};

class clusteran01{ // Analysis of a heliolinc cluster (candidate discovery).
public:
  string recfile;
  int accelct;
  int clusternum;
  int numpoints;
  double timespan;
  int daysteps;
  vector <double> heliopar;
  vector <double> rmsvec; // pos(3), vel(3), postotal, veltotal, total 
  vector <int> clustind;
  clusteran01( string recfile, int accelct, int clusternum, int numpoints, double timespan, int daysteps, vector <double> heliopar, vector <double> rmsvec, vector <int> clustind) :recfile(recfile), accelct(accelct), clusternum(clusternum), numpoints(numpoints), timespan(timespan), daysteps(daysteps), heliopar(heliopar), rmsvec(rmsvec), clustind(clustind) {}
};

class clusteran02{ // Analysis of a heliolinc cluster (candidate discovery).
                   // Differs from clusteran01 in that it includes a quality metric
public:
  string recfile;
  int accelct;
  int clusternum;
  int numpoints;
  double timespan;
  int daysteps;
  vector <double> heliopar;
  vector <double> rmsvec; // pos(3), vel(3), postotal, veltotal, total 
  vector <int> clustind;
  double clustermetric;
  clusteran02( string recfile, int accelct, int clusternum, int numpoints, double timespan, int daysteps, vector <double> heliopar, vector <double> rmsvec, vector <int> clustind, double clustermetric) :recfile(recfile), accelct(accelct), clusternum(clusternum), numpoints(numpoints), timespan(timespan), daysteps(daysteps), heliopar(heliopar), rmsvec(rmsvec), clustind(clustind), clustermetric(clustermetric) {}
};

class lowermetric_clusteran02{
public:
  inline bool operator() (const clusteran02& o1, const clusteran02& o2) {
    return(o1.clustermetric < o2.clustermetric);
  }
};

class clusteran03{ // Analysis of a heliolinc cluster (candidate discovery).
                   // Differs from clusteran02 in that it uses floats rather
                   // than double, to preserve space.
public:
  string recfile;
  int accelct;
  int clusternum;
  int numpoints;
  float timespan;
  int daysteps;
  vector <float> heliopar;
  vector <float> rmsvec; // pos(3), vel(3), postotal, veltotal, total 
  vector <int> clustind;
  float clustermetric;
  clusteran03( string recfile, int accelct, int clusternum, int numpoints, float timespan, int daysteps, vector <float> heliopar, vector <float> rmsvec, vector <int> clustind, float clustermetric) :recfile(recfile), accelct(accelct), clusternum(clusternum), numpoints(numpoints), timespan(timespan), daysteps(daysteps), heliopar(heliopar), rmsvec(rmsvec), clustind(clustind), clustermetric(clustermetric) {}
};

class lowermetric_clusteran03{
public:
  inline bool operator() (const clusteran03& o1, const clusteran03& o2) {
    return(o1.clustermetric < o2.clustermetric);
  }
};

class clusteran04{ // New, more complete storage of important data on a heliolinc cluster

public:
  int clusterct;
  vector <float> rmsvec; // postotal, veltotal, total 
  int pairnum;
  float timespan;
  int uniquepoints;
  int obsnights;
  float clustermetric;
  char rating[SHORTSTRINGLEN];
  vector <float> heliopar;
  vector <double> statevecs;
  vector <int> clustind;
  
  clusteran04(int clusterct, vector <float> rmsvec, int pairnum, float timespan, int uniquepoints, int obsnights, float clustermetric, const string &rating, vector <float> heliopar, vector <double> statevecs, vector <int> clustind) :clusterct(clusterct), rmsvec(rmsvec), pairnum(pairnum), timespan(timespan), uniquepoints(uniquepoints), obsnights(obsnights), clustermetric(clustermetric), heliopar(heliopar), statevecs(statevecs), clustind(clustind) {
    // Copy input value for rating, making sure it's not too long
    assert(rating.size() < sizeof(this->rating));
    std::strncpy(this->rating, rating.c_str(), sizeof(this->rating));
    this->rating[sizeof(this->rating)-1] = 0;
  }
  clusteran04() = default;
};

class lowermetric_clusteran04{
public:
  inline bool operator() (const clusteran04& o1, const clusteran04& o2) {
    return(o1.clustermetric < o2.clustermetric);
  }
};


class heliogridpoint{ // 
public:
  long double dist;
  long double vel;
  vector <long double> acc;
  heliogridpoint(long double dist, long double vel, vector <long double> acc) :dist(dist), vel(vel), acc(acc) {}
};

class string_index{ // Pairs a string with an index, intended for use
                    // accessing elements in a vector of a more complex
                    // class in an order sorted by a string element,
                    // without actually re-sorting the vector.
public:
  string selem;
  long index;
  string_index(string selem, long index) :selem(selem), index(index) {}
};

class lower_string_index{
public:
  inline bool operator() (const string_index& i1, const string_index& i2) {
    return(i1.selem < i2.selem);
  }
};

class char1000_index{ // Pairs a string with an index, intended for use
                      // accessing elements in a vector of a more complex
                      // class in an order sorted by a string element,
                      // without actually re-sorting the vector.
public:
  char selem[1000];
  long index;
  char1000_index(const string &selem, long index) :index(index) {
    // Copy input value of selem, making sure it's not too long
    assert(selem.size() < sizeof(this->selem));
    std::strncpy(this->selem, selem.c_str(), sizeof(this->selem));
    this->selem[sizeof(this->selem)-1] = 0;
  }
};

string char1000_getstring(const char1000_index &ci);

class lower_char1000_index{
public:
  inline bool operator() (const char1000_index& i1, const char1000_index& i2) {
    string s1 = char1000_getstring(i1);
    string s2 = char1000_getstring(i2);
    return(s1 < s2);
  }
};

class char500_index{ // Pairs a string with an index, intended for use
                      // accessing elements in a vector of a more complex
                      // class in an order sorted by a string element,
                      // without actually re-sorting the vector.
public:
  char selem[1000];
  long index;
  char500_index(const string &selem, long index) :index(index) {
    // Copy input value of selem, making sure it's not too long
    assert(selem.size() < sizeof(this->selem));
    std::strncpy(this->selem, selem.c_str(), sizeof(this->selem));
    this->selem[sizeof(this->selem)-1] = 0;
  }
};

string char500_getstring(const char500_index &ci);

class lower_char500_index{
public:
  inline bool operator() (const char500_index& i1, const char500_index& i2) {
    string s1 = char500_getstring(i1);
    string s2 = char500_getstring(i2);
    return(s1 < s2);
  }
};

class char256_index{ // Pairs a string with an index, intended for use
                      // accessing elements in a vector of a more complex
                      // class in an order sorted by a string element,
                      // without actually re-sorting the vector.
public:
  char selem[1000];
  long index;
  char256_index(const string &selem, long index) :index(index) {
    // Copy input value of selem, making sure it's not too long
    assert(selem.size() < sizeof(this->selem));
    std::strncpy(this->selem, selem.c_str(), sizeof(this->selem));
    this->selem[sizeof(this->selem)-1] = 0;
  }
};

string char256_getstring(const char256_index &ci);

class lower_char256_index{
public:
  inline bool operator() (const char256_index& i1, const char256_index& i2) {
    string s1 = char256_getstring(i1);
    string s2 = char256_getstring(i2);
    return(s1 < s2);
  }
};

class char128_index{ // Pairs a string with an index, intended for use
                      // accessing elements in a vector of a more complex
                      // class in an order sorted by a string element,
                      // without actually re-sorting the vector.
public:
  char selem[1000];
  long index;
  char128_index(const string &selem, long index) :index(index) {
    // Copy input value of selem, making sure it's not too long
    assert(selem.size() < sizeof(this->selem));
    std::strncpy(this->selem, selem.c_str(), sizeof(this->selem));
    this->selem[sizeof(this->selem)-1] = 0;
  }
};

string char128_getstring(const char128_index &ci);

class lower_char128_index{
public:
  inline bool operator() (const char128_index& i1, const char128_index& i2) {
    string s1 = char128_getstring(i1);
    string s2 = char128_getstring(i2);
    return(s1 < s2);
  }
};

class char64_index{ // Pairs a string with an index, intended for use
                      // accessing elements in a vector of a more complex
                      // class in an order sorted by a string element,
                      // without actually re-sorting the vector.
public:
  char selem[1000];
  long index;
  char64_index(const string &selem, long index) :index(index) {
    // Copy input value of selem, making sure it's not too long
    assert(selem.size() < sizeof(this->selem));
    std::strncpy(this->selem, selem.c_str(), sizeof(this->selem));
    this->selem[sizeof(this->selem)-1] = 0;
  }
};

string char64_getstring(const char64_index &ci);

class lower_char64_index{
public:
  inline bool operator() (const char64_index& i1, const char64_index& i2) {
    string s1 = char64_getstring(i1);
    string s2 = char64_getstring(i2);
    return(s1 < s2);
  }
};

class long_index{   // Pairs a long integer with an index, intended for use
                    // accessing elements in a vector of a more complex
                    // class in an order sorted by a long integer element,
                    // without actually re-sorting the vector.
public:
  long lelem;
  long index;
  long_index(long lelem, long index) :lelem(lelem), index(index) {}
};

class lower_long_index{
public:
  inline bool operator() (const long_index& i1, const long_index& i2) {
    return(i1.lelem < i2.lelem);
  }
};

class float_index{ // Pairs a float with an index, intended for use
                    // accessing elements in a vector of a more complex
                    // class in an order sorted by a double-precision component,
                    // without actually re-sorting the vector.
public:
  float felem;
  long index;
  float_index(float felem, long index) :felem(felem), index(index) {}
};

class lower_float_index{
public:
  inline bool operator() (const float_index& i1, const float_index& i2) {
    return(i1.felem < i2.felem);
  }
};

class double_index{ // Pairs a double with an index, intended for use
                    // accessing elements in a vector of a more complex
                    // class in an order sorted by a double-precision component,
                    // without actually re-sorting the vector.
public:
  double delem;
  long index;
  double_index(double delem, long index) :delem(delem), index(index) {}
};

class lower_double_index{
public:
  inline bool operator() (const double_index& i1, const double_index& i2) {
    return(i1.delem < i2.delem);
  }
};


class ldouble_index{ // Pairs a long double with an index, intended for use
                    // accessing elements in a vector of a more complex
                    // class in an order sorted by a long double component
                    // without actually re-sorting the vector.
public:
  long double ldelem;
  long index;
  ldouble_index(long double ldelem, long index) :ldelem(ldelem), index(index) {}
};

class lower_ldouble_index{
public:
  inline bool operator() (const ldouble_index& i1, const ldouble_index& i2) {
    return(i1.ldelem < i2.ldelem);
  }
};

class longpd_index{   // Pairs a long integer and a double with an index, intended for use
                      // accessing elements in a vector of a more complex
                      // class in an order sorted by a combination of a
                      // long integer element and a double-precision element,
                      // without actually re-sorting the vector.
public:
  long lelem;
  double delem;
  long index;
  longpd_index(long lelem, double delem, long index) :lelem(lelem), delem(delem), index(index) {}
};

class lower_longpd_index{
public:
  inline bool operator() (const longpd_index& i1, const longpd_index& i2) {
    return(i1.lelem < i2.lelem || (i1.lelem == i2.lelem && i1.delem > i2.delem));
    // Note that this will put the lowest values of lelem,
    // but the HIGHEST values of delem, near the beginning of the vector.
    // This is the specific behavior required by form_clusters_kd4().
  }
};


class keplerian_orbit{
public:
  long double semimaj_axis; // in AU
  long double eccentricity; // unitless
  long double inclination;  // in degrees
  long double long_ascend_node; // Longitude of the ascending node, in degrees
  long double arg_perihelion;   // Argument of perihelion, in degrees
  long double mean_anom;        // Mean anomaly at the epoch, in degrees
  long double mjd_epoch;        // Epoch for the orbit in MJD
  long double mean_daily_motion; // in degrees/day
  keplerian_orbit(long double semimaj_axis, long double eccentricity, long double inclination, long double long_ascend_node, long double arg_perihelion, long double mean_anom, long double mjd_epoch, long double mean_daily_motion) :semimaj_axis(semimaj_axis), eccentricity(eccentricity), inclination(inclination), long_ascend_node(long_ascend_node), arg_perihelion(arg_perihelion), mean_anom(mean_anom), mjd_epoch(mjd_epoch), mean_daily_motion(mean_daily_motion) {}
};

// asteroid_orbit is the same as keplerian_orbit, except that
// it uses doubles rather than long doubles, and
// it also holds a string designation, the absolute magnitude H,
// and the phase parameter G.

class asteroid_orbit{
public:
  string desig;
  double semimaj_axis; // in AU
  double eccentricity; // unitless
  double inclination;  // in degrees
  double long_ascend_node; // Longitude of the ascending node, in degrees
  double arg_perihelion;   // Argument of perihelion, in degrees
  double mean_anom;        // Mean anomaly at the epoch, in degrees
  double mjd_epoch;        // Epoch for the orbit in MJD
  double mean_daily_motion; // in degrees/day
  double H;                // absolute magnitude
  double G;                // phase slope parameter
  asteroid_orbit(string desig, double semimaj_axis, double eccentricity, double inclination, double long_ascend_node, double arg_perihelion, double mean_anom, double mjd_epoch, double mean_daily_motion, double H, double G) :desig(desig), semimaj_axis(semimaj_axis), eccentricity(eccentricity), inclination(inclination), long_ascend_node(long_ascend_node), arg_perihelion(arg_perihelion), mean_anom(mean_anom), mjd_epoch(mjd_epoch), mean_daily_motion(mean_daily_motion), H(H), G(G) {}
};

// asteroid_orbitLD is the same as asteroid_orbit, except that it uses
// long doubles
class asteroid_orbitLD{
public:
  string desig;
  long double semimaj_axis; // in AU
  long double eccentricity; // unitless
  long double inclination;  // in degrees
  long double long_ascend_node; // Longitude of the ascending node, in degrees
  long double arg_perihelion;   // Argument of perihelion, in degrees
  long double mean_anom;        // Mean anomaly at the epoch, in degrees
  long double mjd_epoch;        // Epoch for the orbit in MJD
  long double mean_daily_motion; // in degrees/day
  long double H;                // absolute magnitude
  long double G;                // phase slope parameter
  asteroid_orbitLD(string desig, long double semimaj_axis, long double eccentricity, long double inclination, long double long_ascend_node, long double arg_perihelion, long double mean_anom, long double mjd_epoch, long double mean_daily_motion, long double H, long double G) :desig(desig), semimaj_axis(semimaj_axis), eccentricity(eccentricity), inclination(inclination), long_ascend_node(long_ascend_node), arg_perihelion(arg_perihelion), mean_anom(mean_anom), mjd_epoch(mjd_epoch), mean_daily_motion(mean_daily_motion), H(H), G(G) {}
};

struct EarthState {
    double MJD;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;

    EarthState() = default;

    EarthState(double MJD, double x, double y, double z,
               double vx, double vy, double vz)
            : MJD(MJD), x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) {}
};

void make_ivec(long nx, vector <int> &ivec);
void make_imat(int nx, int ny, vector <vector <int>> &imat);
void make_lvec(int nx, vector <long> &lvec);
void make_lmat(int nx, int ny, vector <vector <long>> &lmat);
void make_cvec(int nx, vector <char> &cvec);
void make_cmat(int nx, int ny, vector <vector <char>> &cmat);
void make_dvec(int nx, vector <double> &dvec);
void make_dmat(int nx, int ny, vector <vector <double>> &dmat);
void make_LDvec(int nx, vector <long double> &ldvec);
void make_LDmat(int nx, int ny, vector <vector <long double>> &ldmat);
double dotprod3d(point3d p1, point3d p2);
long double dotprod3LD(point3LD p1, point3LD p2);
double vecabs3d(point3d p1);
long double vecabs3LD(point3LD p1);
int vecnorm3d(point3d &p1);
int vecnorm3LD(point3LD &p1);
point3d crossprod3d(point3d p1, point3d p2);
point3LD crossprod3LD(point3LD p1, point3LD p2);
long double intpowLD(long double x, int p);
long double factorialLD(int p);
double intpowD(double x, int p);
double factorialD(int p);
double dmean01(const vector <double> &invec);
double drms01(const vector <double> &invec);
int dmeanrms01(const vector <double> &invec, double *mean, double *rms);
point3d celeproj01(double RA, double Dec);
int celedeproj01(point3d p3, double *RA, double *Dec);
point3LD celeproj01LD(long double RA, long double Dec);
int celedeproj01LD(point3LD p3, long double *RA, long double *Dec);
double angspeed01(det_bsc o1, det_bsc o2);
double distradec01(double RA1, double Dec1, double RA2, double Dec2);
int distradec02(double ra1,double dec1,double ra2,double dec2,double *dist,double *pa);
long medindex(const vector <xy_index> &xyvec, int dim);
int splitxy(const vector <xy_index> &xyvec, int dim, long unsigned int splitpoint, vector <xy_index> &left, vector <xy_index> &right);
int kdtree01(const vector <xy_index> &xyvec, int dim, long unsigned int rootptxy, long unsigned int rootptkd, vector <kdpoint> &kdvec);
int kdrange01(const vector <kdpoint> &kdvec,double x,double y,double range,vector <long> &indexvec);
long medind_3d_index(const vector <point3d_index> &pointvec, int dim);
int split3d_index(const vector <point3d_index> &pointvec, int dim, long splitpoint, vector <point3d_index> &left, vector <point3d_index> &right);
int kdtree_3d_index(const vector <point3d_index> &invec, int dim, long splitpoint, long kdroot, vector <KD_point3d_index> &kdvec);
double point3d_index_dist2(const point3d_index &p1, const point3d_index &p2);
int kdrange_3d_index(const vector <KD_point3d_index> &kdvec, const point3d_index &querypoint, double range, vector <long> &indexvec);
long medind_4d_index(const vector <point4d_index> &pointvec, int dim);
int split4d_index(const vector <point4d_index> &pointvec, int dim, long splitpoint, vector <point4d_index> &left, vector <point4d_index> &right);
int kdtree_4d_index(const vector <point4d_index> &invec, int dim, long splitpoint, long kdroot, vector <KD_point4d_index> &kdvec);
double point4d_index_dist2(const point4d_index &p1, const point4d_index &p2);
int kdrange_4d_index(const vector <KD_point4d_index> &kdvec, const point4d_index &querypoint, double range, vector <long> &indexvec);
long medind_6LDx2(const vector <point6LDx2> &pointvec, int dim);
int splitLDx2(const vector <point6LDx2> &pointvec, int dim, long unsigned int splitpoint, vector <point6LDx2> &left, vector <point6LDx2> &right);
int kdtree_6D01(const vector <point6LDx2> &invec, int dim, long unsigned int splitpoint, long unsigned int kdroot, vector <KD_point6LDx2> &kdvec);
long double point6LDx2_dist(const point6LDx2 &p1, const point6LDx2 &p2);
long double point6LDx2_dist2(const point6LDx2 &p1, const point6LDx2 &p2);
int kdrange_6D01(const vector <KD_point6LDx2> &kdvec, const point6LDx2 &querypoint, long double range, vector <long> &indexvec);
long double cluster_stats6D01(const vector <KD_point6LDx2> &cluster, vector <long double> &meanvals, vector <long double> &rmsvals);
int DBSCAN_6D01(vector <KD_point6LDx2> &kdtree, long double clustrad, int npt, const vector <det_bsc> &detvec, const vector <string> &det_id_vec, vector <KD6_clust> &outclusters, string rmsfile);
int DBSCAN_6D02(vector <KD_point6LDx2> &kdtree, long double clustrad, int npt, vector <KD6_clust> &outclusters);
point6ix2 conv_6LD_to_6i(point6LDx2 p1, long double scale);
point6LDx2 conv_6i_to_6LD(point6ix2 p1, long double scale);
point6ix2 conv_6d_to_6i(point6dx2 p1, double scale);
point6dx2 conv_6i_to_6d(point6ix2 p1, double scale);
long medind_6ix2(const vector <point6ix2> &pointvec, int dim);
int splitix2(const vector <point6ix2> &pointvec, int dim, long unsigned int splitpoint, vector <point6ix2> &left, vector <point6ix2> &right);
int kdtree_6i01(const vector <point6ix2> &invec, int dim, long unsigned int splitpoint, long unsigned int kdroot, vector <KD_point6ix2> &kdvec);
long point6ix2_dist2(const point6ix2 &p1, const point6ix2 &p2);
int kdrange_6i01(const vector <KD_point6ix2> &kdvec, const point6ix2 &querypoint, long range, vector <long> &indexvec);
double cluster_stats6i01(const vector <KD_point6ix2> &cluster, double intconvscale, vector <double> &meanvals, vector <double> &rmsvals);
int DBSCAN_6i01(vector <KD_point6ix2> &kdtree, double clustrad, int npt, double intconvscale, vector <KD6i_clust> &outclusters, int verbose);
int KDRclust_6i01(vector <KD_point6ix2> &kdtree, double clustrad, int npt, double intconvscale, vector <KD6i_clust> &outclusters, int verbose);
int celestial_to_statevec(double RA, double Dec,double delta,point3d &baryvec);
int celestial_to_statevecLD(long double RA, long double Dec,long double delta,point3LD &baryvec);
int celestial_to_stateunit(double RA, double Dec,point3d &baryvec);
int celestial_to_stateunitLD(long double RA, long double Dec, point3LD &baryvec);
int get_csv_string01(const string &lnfromfile, string &outstring, int startpoint);
int get_sv_string01(const string &lnfromfile, string &outstring, int startpoint);
int read_horizons_file(string infile, vector <double> &mjdvec, vector <point3d> &pos, vector <point3d> &vel);
int read_horizons_fileLD(string infile, vector <long double> &mjdvec, vector <point3LD> &pos, vector <point3LD> &vel);
int read_horizons_csv(string infile, vector <double> &mjdvec, vector <point3d> &pos, vector <point3d> &vel);
int read_horizons_csv(string infile, vector <EarthState> &earthpos);
int poleswitch01(const double &inRA, const double &inDec, const double &poleRA, const double &poleDec, const double &oldpoleRA, double &newRA, double &newDec);
int poleswitch01LD(const long double &inRA, const long double &inDec, const long double &poleRA, const long double &poleDec, const long double &oldpoleRA, long double &newRA, long double &newDec);
int poleswitch02(const double &inRA, const double &inDec, const double &poleRA, const double &poleDec, const double &oldpoleRA, double &newRA, double &newDec);
int poleswitch02LD(const long double &inRA, const long double &inDec, const long double &poleRA, const long double &poleDec, const long double &oldpoleRA, long double &newRA, long double &newDec);
int precess01a(double ra1,double dec1,double mjd,double *ra2,double *dec2,int precesscon);
int precess01aLD(long double ra1,long double dec1,long double mjd,long double *ra2,long double *dec2,int precesscon);
int solvematrix01(const vector <vector <double>> &inmat, int eqnum, vector <double> &outvec, int verbose);
int solvematrix01LD(const vector <vector <long double>> &inmat, int eqnum, vector <long double> &outvec, int verbose);
int perfectpoly01(const vector <double> &x, const vector <double> &y, vector <double> &fitvec);
int perfectpoly01LD(const vector <long double> &x, const vector <long double> &y, vector <long double> &fitvec);
int planetpos01(double detmjd, int polyorder, const vector <double> &posmjd, const vector <point3d> &planetpos, point3d &outpos);
int planetpos01(double detmjd, int polyorder, const vector <EarthState> &planetpos, point3d &outpos);
int planetpos01LD(long double detmjd, int polyorder, const vector <long double> &posmjd, const vector <point3LD> &planetpos, point3LD &outpos);
int planetposvel01(double detmjd, int polyorder, const vector <double> &posmjd, const vector <point3d> &planetpos, const vector <point3d> &planetvel, point3d &outpos, point3d &outvel);
int planetposvel01(double detmjd, int polyorder, const vector <EarthState> &planetpos, point3d &outpos, point3d &outvel);
int planetposvel01LD(long double detmjd, int polyorder, const vector <long double> &posmjd, const vector <point3LD> &planetpos, const vector <point3LD> &planetvel, point3LD &outpos, point3LD &outvel);
int nplanetpos01LD(long double detmjd, int planetnum, int polyorder, const vector <long double> &posmjd, const vector <point3LD> &planetpos, vector <point3LD> &outpos);
int nplanetgrab01LD(int pointrequest, int planetnum, const vector <long double> &posmjd, const vector <point3LD> &planetpos, vector <point3LD> &outpos);
int observer_barycoords01(double detmjd, int polyorder, double lon, double obscos, double obssine, const vector <double> &posmjd, const vector <point3d> &planetpos, point3d &outpos);
int observer_barycoords01LD(long double detmjd, int polyorder, long double lon, long double obscos, long double obssine, const vector <long double> &posmjd, const vector <point3LD> &planetpos, point3LD &outpos);
int observer_baryvel01(double detmjd, int polyorder, double lon, double obscos, double obssine, const vector <double> &posmjd, const vector <point3d> &planetpos, const vector <point3d> &planetvel, point3d &outpos, point3d &outvel);
int observer_baryvel01(double detmjd, int polyorder, double lon, double obscos, double obssine, const vector <EarthState> &earthpos, point3d &outpos, point3d &outvel);
int helioproj01(point3d unitbary, point3d obsbary,double heliodist,double &geodist, point3d &projbary);
int helioproj01LD(point3LD unitbary, point3LD obsbary, long double heliodist, long double &geodist, point3LD &projbary);
int helioproj02LD(point3LD unitbary, point3LD obsbary, long double heliodist, vector <long double> &geodist, vector <point3LD> &projbary);
int helioproj02(point3d unitbary, point3d obsbary, double heliodist, vector <double> &geodist, vector <point3d> &projbary);
int accelcalc01LD(int planetnum, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const point3LD &targpos, point3LD &accel);
int integrate_orbit_constac(int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, long double mjdstart, point3LD startpos, point3LD startvel, long double mjdend, point3LD &endpos, point3LD &endvel);
long double kep_transcendental(long double q, long double e, long double tol);
double kep_transcendental(double q, double e, double tol);
int Stumpff_func(const double xc, double *c0, double *c1, double *c2, double *c3);
int Stumpff_func_cf(const double xc, double *c0, double *c1, double *c2, double *c3);
int Keplerint(const long double MGsun, const long double mjdstart, const point3LD &startpos, const point3LD &startvel, const long double mjdend, point3LD &endpos, point3LD &endvel);
int Keplerint(const double MGsun, const double mjdstart, const point3d &startpos, const point3d &startvel, const double mjdend, point3d &endpos, point3d &endvel);
int Kepler_fg_func_int(const double MGsun, const double mjdstart, const point3d &startpos, const point3d &startvel, const double mjdend, point3d &endpos, point3d &endvel);
int Kepler_univ_int(const double MGsun, const double mjdstart, const point3d &startpos, const point3d &startvel, const double mjdend, point3d &endpos, point3d &endvel);
int Kepler2dyn(const long double mjdnow, const keplerian_orbit &keporb, point3LD &outpos,  point3LD &outvel);
int Kepler2dyn(const double mjdnow, const asteroid_orbit &oneorb, point3d &outpos,  point3d &outvel);
int Kepler2dyn(const long double mjdnow, const asteroid_orbitLD &keporb, point3LD &outpos,  point3LD &outvel);
long double hyp_transcendental(long double q, long double e, long double tol);
double hyp_transcendental(double q, double e, double tol);
int Hyper_Kepint(const long double MGsun, const long double mjdstart, const point3LD &startpos, const point3LD &startvel, const long double mjdend, point3LD &endpos, point3LD &endvel);
int Hyper_Kepint(const double MGsun, const double mjdstart, const point3d &startpos, const point3d &startvel, const double mjdend, point3d &endpos, point3d &endvel);
int integrate_orbit01LD(int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, long double mjdstart, point3LD startpos, point3LD startvel, long double mjdend, point3LD &endpos, point3LD &endvel);
int integrate_orbit02LD(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, long double mjdstart, point3LD startpos, point3LD startvel, long double mjdend, point3LD &endpos, point3LD &endvel);
int iswhitespace(int c);
int readconfigLD(ifstream &instream1, long double *ldval);
int readconfigd(ifstream &instream1, double *dval);
int readconfigint(ifstream &instream1, int *ival);
int readconfigstring(ifstream &instream1, string &sval);
int read_accel_fileLD(string accelfile, vector <long double> &heliodist, vector <long double> &heliovel, vector <long double> &helioacc);
int read_longitude_fileLD(string accelfile, vector <long double> &longitude_vel, vector <long double> &longitude_acc);
long double weight_posvel_rms(const vector <point3LD> &poscluster,const vector <point3LD> &velcluster,const long double dtime, vector <long double> &rmsvec);
int linfituw01(const vector <double> &x, const vector <double> &y, double &slope, double &intercept);
int linfit01(const vector <double> &x, const vector <double> &y, const vector <double> &yerr, double &slope, double &intercept);
int multilinfit01(const vector <double> &yvec, const vector <vector <double>> &xmat, int pnum, int fitnum, vector <double> &avec, int verbose);
int multilinfit02(const vector <double> &yvec, const vector <double> &sigvec, const vector <vector <double>> &xmat, int pnum, int fitnum, vector <double> &avec, int verbose);
int multilinfit02b(const vector <double> &yvec, const vector <double> &varvec, const vector <vector <double>> &xmat, int pnum, int fitnum, vector <double> &avec, int verbose);
int arc2cel01(double racenter,double deccenter,double dist,double pa,double &outra,double &outdec);
int obscode_lookup(const vector <observatory> &observatory_list, const char* obscode, double &obslon, double &plxcos,double &plxsin);
string intzero01i(const int i, const int n);
int get_col_vector01(const string &lnfromfile, vector <string> &outvec);
int mjd2mpcdate(double MJD,int &year,int &month,double &day);
int stringline01(const string &lnfromfile, vector <string> &outstrings);
long double ldmedian(const vector <long double> &invec);
int ldmedian_minmax(const vector <long double> &invec, long double &median, long double &min, long double &max);
double dmedian(const vector <double> &invec);
int dmedian_minmax(const vector <double> &invec, double &median, double &min, long double &max);
float fmedian(const vector <float> &invec);
int fmedian_minmax(const vector <float> &invec, float &median, float &min, float &max);
int stateunit_to_celestial(point3d &baryvec, double &RA, double &Dec);
int stateunitLD_to_celestial(point3LD &baryvec, long double &RA, long double &Dec);
int integrate_orbit03LD(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const vector <long double> &obsMJD, point3LD startpos, point3LD startvel, long double mjdstart, vector <point3LD> &obspos, vector <point3LD> &obsvel);
long double tortoisechi01(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, const vector <long double> &scalestate, long double timescale, long double mjdstart, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid);
int integrate_orbit04LD(int polyorder, int planetnum, const vector <long double> &planetmjd, const vector <long double> &planetmasses, const vector <point3LD> &planetpos, point3LD startpos, point3LD startvel, int startpoint, int endpoint, vector <long double> &outMJD, vector <point3LD> &outpos, vector <point3LD> &outvel);
double gaussian_deviate();
int uvw_to_galcoord(const double &u, const double &v, const double &w, double &RA, double &Dec);
long double unitvar(mt19937_64 &generator);
double gaussian_deviate_mt(mt19937_64 &generator);
int multilinfit01(const vector <double> &yvec, const vector <double> &sigvec, const vector <vector <double>> &xmat, int pnum, int fitnum, vector <double> &avec);
int polyfit01(const vector <double> &yvec, const vector <double> &sigvec, const vector <double> &xvec, int pnum, int polyorder, vector <double> &avec);
int vaneproj01LD(point3LD unitbary, point3LD obsbary, long double ecliplon, long double &geodist, point3LD &projbary);
int vaneproj01d(point3d unitbary, point3d obsbary, double ecliplon, double min_proj_sine, double &geodist, point3d &projbary);
double Twopoint_KepQ(double x);
//int Twopoint_Qf(double z, double MGsun, double lambda1, double lambda2, double X, double Y, double deltat, double *f, double *fprime);
int Twopoint_Kepler_vel(const double MGsun, const point3d startpoint, const point3d endpoint, const double timediff, point3d &startvel, double *semimaj, int itmax);
double Twopoint_KepQstar(double x);
int Twopoint_Kepler_vstar(const double MGsun, const point3d startpoint, const point3d endpoint, const double timediff, point3d &startvel, int itmax);
int Keplerint_multipoint01(const long double MGsun, const long double mjdstart, const vector <long double> &obsMJD, const point3LD &startpos, const point3LD &startvel, vector <point3LD> &obspos, vector <point3LD> &obsvel);
int Keplerint_multipoint02(const long double MGsun, const long double mjdstart, const vector <long double> &obsMJD, const point3LD &startpos, const point3LD &startvel, vector <point3LD> &obspos, vector <point3LD> &obsvel, long double *semimajor_axis, long double *eccen, long double *angperi);
int Keplerint_multipoint02(const double MGsun, const double mjdstart, const vector <double> &obsMJD, const point3d &startpos, const point3d &startvel, vector <point3d> &obspos, vector <point3d> &obsvel, double *semimajor_axis, double *eccen, double *angperi);
int Keplerint_multipoint_fgfunc(const double MGsun, const double mjdstart, const vector <double> &obsMJD, const point3d &startpos, const point3d &startvel, vector <point3d> &obspos, vector <point3d> &obsvel, double *semimajor_axis, double *eccen);
int Keplerint_multipoint_univar(const double MGsun, const double mjdstart, const vector <double> &obsMJD, const point3d &startpos, const point3d &startvel, vector <point3d> &obspos, vector <point3d> &obsvel, double *semimajor_axis, double *eccen);
long double orbitchi01(const point3LD &objectpos, const point3LD &objectvel, const long double mjdstart, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid);
long double orbitchi02(const point3LD &objectpos, const point3LD &objectvel, const long double mjdstart, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid, long double *semimajor_axis, long double *eccen, long double *angperi);
double orbitchi02(const point3d &objectpos, const point3d &objectvel, const double mjdstart, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid, double *semimajor_axis, double *eccen, double *angperi);
double orbitchi_fgfunc(const point3d &objectpos, const point3d &objectvel, const double mjdstart, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid, double *semimajor_axis, double *eccen);
double orbitchi_univar(const point3d &objectpos, const point3d &objectvel, const double mjdstart, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid, double *semimajor_axis, double *eccen);
long double TwopointF(long double a, long double k, long double lambda1, long double lambda2, long double deltat, long double Xsign, long double Ysign);
double TwopointF(double a, double k, double lambda1, double lambda2, double deltat, double Xsign, double Ysign);
long double TwopointFprime(long double a, long double k, long double lambda1, long double lambda2, long double deltat, long double Xsign, long double Ysign);
double TwopointFprime(double a, double k, double lambda1, double lambda2, double deltat, double Xsign, double Ysign);
int eccen_calc_fast(long double a, point3LD rvec1, point3LD rvec2, long double *e, long double *theta, long double Xsign, long double Ysign);
point3LD Twopoint_Kepler_v1(const long double GMsun, const point3LD startpos, const point3LD endpos, const long double timediff, long double Y, long double *a, int itmax, int verbose);
point3d Twopoint_Kepler_v1(const double GMsun, const point3d startpos, const point3d endpos, const double timediff, const double Ysign, double *a, int itmax, int verbose);
point3LD geodist_to_3Dpos01(long double RA, long double Dec, point3LD observerpos, long double geodist);
point3d geodist_to_3dpos01(double RA, double Dec, point3d observerpos, double geodist);
int Herget_unboundcheck01(long double geodist1, long double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec);
long double Hergetchi01(long double geodist1, long double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid, vector <long double> &orbit, int verbose);
double Hergetchi01(double geodist1, double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid, vector <double> &orbit, int verbose);
double Hergetchi_vstar(double geodist1, double geodist2, int Hergetpoint1, int Hergetpoint2, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid, vector <double> &orbit, int verbose);
int Herget_simplex_int(long double geodist1, long double geodist2, long double simpscale, long double simplex[3][2], int simptype);
int Herget_simplex_int(double geodist1, double geodist2, double simpscale, double simplex[3][2], int simptype);
long double Hergetfit01(long double geodist1, long double geodist2, long double simplex_scale, int simptype, long double ftol, int point1, int point2, const vector <point3LD> &observerpos, const vector <long double> &obsMJD, const vector <long double> &obsRA, const vector <long double> &obsDec, const vector <long double> &sigastrom, vector <long double> &fitRA, vector <long double> &fitDec, vector <long double> &resid, vector <long double> &orbit, int verbose);
double Hergetfit_vstar(double geodist1, double geodist2, double simplex_scale, int simptype, double ftol, int point1, int point2, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, vector <double> &fitRA, vector <double> &fitDec, vector <double> &resid, vector <double> &orbit, int verbose);
int wrap_Hergetfit01(double simplex_scale, int simptype, double ftol, int point1, int point2, const vector <point3d> &observerpos, const vector <double> &obsMJD, const vector <double> &obsRA, const vector <double> &obsDec, const vector <double> &sigastrom, double MJDref, int rmspow, int verbose, hlclust &onecluster);
int wrap_Hergetfit02(double simplex_scale, int simptype, double ftol, const vector <vector <point3d>> &observerpos, const vector <vector <double>> &obsMJD, const vector <vector <double>> &obsRA, const vector <vector <double>> &obsDec, const vector <vector <double>> &sigastrom, double MJDref, int rmspow, int verbose, vector <hlclust> &outclust, int threadct);
int wrap_Hergetfit_vstar(double simplex_scale, int simptype, double ftol, const vector <vector <point3d>> &observerpos, const vector <vector <double>> &obsMJD, const vector <vector <double>> &obsRA, const vector <vector <double>> &obsDec, const vector <vector <double>> &sigastrom, double MJDref, int rmspow, int verbose, vector <hlclust> &outclust, int threadct);
double MPCcal2MJD(int year, int month, double day);
int mpc80_parseline(const string &lnfromfile, string &object, double *MJD, double *RA, double *Dec, double *mag, string &band, string &obscode);
double mpc80_mjd(const string &lnfromfile);
int read_obscode_file(string obscodefile,  vector <observatory> &observatory_list);
int read_obscode_file2(string obscodefile,  vector <observatory> &observatory_list, int verbose);
int read_detection_filemt(string indetfile, int idcol, int mjdcol, int racol, int deccol, int magcol,int bandcol, int obscodecol, vector <det_obsmag_indvec> &detvec, int forcerun);
int read_detection_filemt2(string indetfile, int mjdcol, int racol, int deccol, int magcol, int idcol, int bandcol, int obscodecol, int trail_len_col, int trail_PA_col, int sigmag_col, int sig_across_col, int sig_along_col, int known_obj_col, int det_qual_col, vector <hldet> &detvec, int verbose, int forcerun);
int read_detection_file_MPC80(string indetfile, vector <hldet> &detvec);
int read_pairdet_file(string pairdetfile, vector <hldet> &detvec, int verbose);
int read_tracklet_file(string trackletfile, vector <tracklet> &tracklets, int verbose);
int read_longpair_file(string pairfile, vector <longpair> &pairvec, int verbose);
int append_longpair_file(string pairfile, long oldsize, vector <longpair> &pairvec, int verbose);
int read_radhyp_file(string hypfile, vector <hlradhyp> &accelmat, int verbose);
int read_clustersum_file(string sumfile, vector <hlclust> &clustvec, int verbose);
int append_clustersum_file(string sumfile, vector <hlclust> &clustvec, int verbose);
double avg_extrema(const vector <double> &x);
int read_image_file(string inimfile, vector <img_log03> &img_log);
int read_image_file(string inimfile, vector <hlimage> &img_log);
int read_image_file2(string inimfile, vector <hlimage> &img_log);
int load_image_table(vector <img_log03> &img_log, const vector <det_obsmag_indvec> &detvec);
int load_image_table(vector <hlimage> &img_log, const vector <hldet> &detvec, const vector <observatory> &observatory_list, const vector <double> &Earthmjd, const vector <point3d> &Earthpos, const vector <point3d> &Earthvel);
int load_image_indices(vector <hlimage> &img_log, vector <hldet> &detvec, double imagetimetol, int forcerun);
int load_image_indices2(vector <hlimage> &img_log, vector <hldet> &detvec, double imagetimetol, int forcerun);
int find_pairs(vector <hldet> &detvec, const vector <hlimage> &img_log, vector <hldet> &pairdets, vector <vector <long>> &indvecs, vector <longpair> &pairvec, double mintime, double maxtime, double imrad, double maxvel, int verbose);
int find_trailpairs(vector <hldet> &detvec, const vector <hlimage> &img_log, vector <hldet> &pairdets, vector <vector <long>> &indvecs, vector <longpair> &pairvec, double mintime, double maxtime, double imrad, double maxvel, double siglenscale, double sigpascale, int verbose);
int merge_pairs(const vector <hldet> &pairdets, vector <vector <long>> &indvecs, const vector <longpair> &pairvec, vector <tracklet> &tracklets, vector <longpair> &trk2det, int mintrkpts, double maxgcr, double minarc, double minvel, double maxvel, int verbose);
int merge_trailpairs(const vector <hldet> &pairdets, vector <vector <long>> &indvecs, const vector <longpair> &pairvec, vector <tracklet> &tracklets, vector <longpair> &trk2det, int mintrkpts, double maxgcr, double minarc, double minvel, double maxvel, int verbose);
int record_pairs(vector <hldet> &detvec, vector <hldet> &detvec_fixed, vector <tracklet> &tracklets, vector <longpair> &trk2det, int verbose);
int make_tracklets(vector <hldet> &detvec, vector <hlimage> &image_log, MakeTrackletsConfig config, vector <hldet> &pairdets,vector <tracklet> &tracklets, vector <longpair> &trk2det);
int make_trailed_tracklets(vector <hldet> &detvec, vector <hlimage> &image_log, MakeTrackletsConfig config, vector <hldet> &pairdets,vector <tracklet> &tracklets, vector <longpair> &trk2det);
int remake_tracklets(vector <hldet> &detvec, vector <hldet> &detvec_fixed, vector <hlimage> &image_log,vector <tracklet> &tracklets, vector <longpair> &trk2det, int verbose);
int trk2statevec(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double heliodist, double heliovel, double helioacc, double chartimescale, vector <point6ix2> &allstatevecs, double mjdref, double mingeoobs, double minimpactpar);
int trk2statevec_fgfunc(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double heliodist, double heliovel, double helioacc, double chartimescale, vector <point6ix2> &allstatevecs, double mjdref, double mingeoobs, double minimpactpar, double max_v_inf);
int trk2statevec_clusterprobe(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double heliodist, double heliovel, double helioacc, double chartimescale, vector <point6dx2> &allstatevecs, double mjdref);
int trk2statevec_clusterprobe_innea(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double heliodist, double heliovel, double helioacc, double chartimescale, vector <point6dx2> &allstatevecs, double mjdref);
int trk2statevec_univar(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double heliodist, double heliovel, double helioacc, double chartimescale, vector <point6ix2> &allstatevecs, double mjdref, double mingeoobs, double minimpactpar, double max_v_inf);
int trk2statevane_fgfunc(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double lambda0, double lambda_dot, double lambda_ddot, double chartimescale, double minsunelong, double maxsunelong, double min_proj_sine, double max_heliodist, vector <point6ix2> &allstatevecs, double mjdref, double mingeoobs, double minimpactpar, double max_v_inf);
int trk2statevane_univar(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double lambda0, double lambda_dot, double lambda_ddot, double chartimescale, double minsunelong, double maxsunelong, double min_proj_sine, double maxheliodist, vector <point6ix2> &allstatevecs, double mjdref, double mingeoobs, double minimpactpar, double max_v_inf);
void lastroot(const vector <double> &intvec, vector <double> &rootvec, long N);
int trk2statevec_omp(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double heliodist, double heliovel, double helioacc, double chartimescale, vector <point6ix2> &allstatevecs, double mjdref, double mingeoobs, double minimpactpar);
int trk2statevec_omp2(const vector <hlimage> &image_log, const vector <tracklet> &tracklets, double heliodist, double heliovel, double helioacc, double chartimescale, vector <point6ix2> &allstatevecs, double mjdref, double mingeoobs, double minimpactpar);
vector <long> tracklet_lookup(const vector <longpair> &trk2det, long trknum);
point3d earthpos01(const vector <EarthState> &earthpos, double mjd);
int form_clusters(const vector <point6ix2> &allstatevecs, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const point3d &Earthrefpos, double heliodist, double heliovel, double helioacc, double chartimescale, vector <hlclust> &outclust, vector <longpair> &clust2det, long &realclusternum, double cluster_radius, double dbscan_npt, double mingeodist, double geologstep, double maxgeodist, int mintimespan, int minobsnights, int verbose);
int form_clusters_kd(const vector <point6ix2> &allstatevecs, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const point3d &Earthrefpos, double heliodist, double heliovel, double helioacc, double chartimescale, vector <hlclust> &outclust, vector <longpair> &clust2det, long &realclusternum, double cluster_radius, double dbscan_npt, double mingeodist, double geologstep, double maxgeodist, int mintimespan, int minobsnights, int verbose);
long intvec_lower(const vector <long> &ivec1, const vector <long> &ivec2);
long intvec_findplace(const vector <vector <long>> &imat, const vector <long> &ivec, int &dup);
long blend_vector(vector<long> vec);
int form_clusters_kd2(const vector <point6ix2> &allstatevecs, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const point3d &Earthrefpos, double heliodist, double heliovel, double helioacc, double chartimescale, vector <hlclust> &outclust, vector <longpair> &clust2det, long &realclusternum, double cluster_radius, double dbscan_npt, double mingeodist, double geologstep, double maxgeodist, int mintimespan, int minobsnights, int verbose);
int form_clusters_kd3(const vector <point6ix2> &allstatevecs, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const point3d &Earthrefpos, double heliodist, double heliovel, double helioacc, double chartimescale, vector <hlclust> &outclust, vector <longpair> &clust2det, long &realclusternum, double cluster_radius, double dbscan_npt, double mingeodist, double geologstep, double maxgeodist, int mintimespan, int minobsnights, int verbose);
int form_clusters_kd4(const vector <point6ix2> &allstatevecs, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const point3d &Earthrefpos, double heliodist, double heliovel, double helioacc, double chartimescale, vector <hlclust> &outclust, vector <longpair> &clust2det, long &realclusternum, double cluster_radius, double clustchangerad, double dbscan_npt, double mingeodist, double geologstep, double maxgeodist, int mintimespan, int minobsnights, int verbose);
int heliolinc_alg(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_fgfunc(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_univar(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_danby(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliovane_alg_danby(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &lambdahyp, const vector <EarthState> &earthpos, HeliovaneConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_omp(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_omp2(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_omp3(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_ompdanby(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_ompkd(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_ompkd3(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_ompkd4(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliolinc_alg_kd(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &radhyp, const vector <EarthState> &earthpos, HeliolincConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int heliovane_alg_ompdanby(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <tracklet> &tracklets, const vector <longpair> &trk2det, const vector <hlradhyp> &lambdahyp, const vector <EarthState> &earthpos, HeliovaneConfig config, vector <hlclust> &outclust, vector <longpair> &clust2det);
int link_refine_Herget(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <hlclust> &inclust, const vector  <longpair> &inclust2det, LinkRefineConfig config, vector <hlclust> &outclust, vector <longpair> &outclust2det);
int link_refine_Herget_univar(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <hlclust> &inclust, const vector  <longpair> &inclust2det, LinkRefineConfig config, vector <hlclust> &outclust, vector <longpair> &outclust2det);
int link_dedup(const vector <hlclust> &inclust, const vector  <longpair> &inclust2det, vector <hlclust> &outclust, vector  <longpair> &outclust2det);
int planepolefind(const vector <point3d> &invecs, point3d &polevec);
int link_purify(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <hlclust> &inclust, const vector  <longpair> &inclust2det, LinkPurifyConfig config, vector <hlclust> &outclust, vector <longpair> &outclust2det);
int link_planarity(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <hlclust> &inclust1, const vector  <longpair> &inclust2det1, LinkPurifyConfig config, vector <hlclust> &outclust, vector <longpair> &outclust2det);
int link_refine_Herget_omp(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <hlclust> &inclust, const vector  <longpair> &inclust2det, LinkRefineConfig config, vector <hlclust> &outclust, vector <longpair> &outclust2det);
int link_refine_Herget_omp2(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <hlclust> &inclust, const vector  <longpair> &inclust2det, LinkRefineConfig config, vector <hlclust> &outclust, vector <longpair> &outclust2det);
int link_refine_Herget_omp3(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <hlclust> &inclust, const vector  <longpair> &inclust2det, LinkRefineConfig config, vector <hlclust> &outclust, vector <longpair> &outclust2det);
int link_refine_Herget_omp4(const vector <hlimage> &image_log, const vector <hldet> &detvec, const vector <hlclust> &inclust, const vector  <longpair> &inclust2det, LinkRefineConfig config, vector <hlclust> &outclust, vector <longpair> &outclust2det);
int parse_clust2det(const vector <hldet> &detvec, const vector <longpair> &inclust2det, vector <hldet> &clustdet);
int greatcircfit(const vector <hldet> &trackvec, double &poleRA, double &poleDec,double &angvel,double &pa,double &crosstrack,double &alongtrack);
int read_orbline(ifstream &instream1, asteroid_orbit &oneorb);
int read_orbline(ifstream &instream1, asteroid_orbitLD &oneorb);
int read_orbline(string lnfromfile, asteroid_orbit &oneorb);

