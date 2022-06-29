#ifndef HELA_MAKETRACKLETS_H
#define HELA_MAKETRACKELTS_H

#define DEBUG 0

#include "constants.h"
#include "types.h"
#include <tuple>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <cmath>

using namespace std;

struct MakeTrackletsConfig {
    int mintrkpts = 2;
    double imagetimetol = 1.0;    // Tolerance for matching image time, in seconds
    double maxvel = 1.5l;         // Default max angular velocity in deg/day.
    double minvel = 0.0l;         // Min angular velocity in deg/day
    double minarc = 0.0l;         // Min total angular arc in arcseconds
    double maxtime = 1.5 / 24.0;  // Default max inter-image time interval
                                // for tracklets, in days.
    double mintime = 1.0 / SOLARDAY;
    double angvel = 0.0l;
    double maxdist = 1.5 * 1.5 / 24.0;  // Max angular distance a tracklet
                                        // could span, in degrees.
    double imagerad = 2.0;        // radius from image center to most distant corner (deg)
    double maxgcr = 0.5;         // Default maximum Great Circle Residual allowed for a valid tracklet

    // Input files names
    string indetfile;
    string inimfile;
    string earthfile;
    string obscodefile;
};

using Point3LD = Point3<long double>;
using Point3D = Point3<double>;

void stringncopy01(char *dest, string const& source, int n);
int stringnmatch01(char const* string1, char const* string2, int n);
int read_horizons_fileLD(string infile, std::vector<long double> &mjdvec, std::vector<Point3LD> &pos, std::vector<Point3LD> &vel);
int obscode_lookup(std::vector<Observatory> const& observatory_list, char const* obscode, double &obslon, double &plxcos, double &plxsin);
void make_LDvec(int nx, std::vector<long double> &ldvec);
void make_LDmat(int nx, int ny, std::vector<std::vector<long double>> &ldmat);
int solvematrix01LD(std::vector<std::vector<long double>> const& inmat, int eqnum, std::vector<long double> &outvec, int verbose);
int perfectpoly01LD(std::vector<long double> const& x, std::vector<long double> const& y, std::vector<long double> &fitvec);
int planetpos01LD(long double detmjd, int polyorder, std::vector<long double> const& posmjd, std::vector<Point3LD> const& planetpos, Point3LD &outpos);
int celestial_to_statevecLD(long double RA, long double Dec, long double delta, Point3LD &baryvec);
int precess01aLD(long double ra1, long double dec1,long double mjd, long double *ra2, long double *dec2, int precesscon);
int observer_barycoords01LD(long double detmjd, int polyorder, long double lon, long double obscos, long double obssine, std::vector<long double> const& posmjd, std::vector<Point3LD> const& planetpos, Point3LD &outpos);
long medindex(std::vector<XYIndex> const& xyvec, int dim);
Point3D celeproj01(double RA, double Dec);
int celedeproj01(Point3D p3, double *RA, double *Dec);
double distradec01(double RA1, double Dec1, double RA2, double Dec2);
int distradec02(double ra1,double dec1,double ra2,double dec2,double *dist,double *pa);
int splitxy(std::vector<XYIndex> const& xyvec, int dim, long splitpoint, std::vector<XYIndex> &left, std::vector<XYIndex> &right);
int kdtree01(std::vector<XYIndex> const& xyvec, int dim, long rootptxy, long rootptkd, std::vector<KDPoint> &kdvec);
int kdrange01(std::vector<KDPoint> const& kdvec, double x, double y, double range, std::vector<long> &indexvec);
int linfituw01(std::vector<double> const& x, std::vector<double> const& y, double &slope, double &intercept);
int arc2cel01(double racenter, double deccenter, double dist, double pa, double &outra,double &outdec);

inline bool timeCompareDetections(Detection const& d1, Detection const& d2) {
    return (d1.MJD < d2.MJD || (d1.MJD==d1.MJD && stringnmatch01(d1.obscode, d2.obscode, 3) == -1) || (d1.MJD==d2.MJD && stringnmatch01(d1.obscode, d2.obscode, 3) == 0 && d1.RA < d2.RA));
}
inline bool timeCompareImageLog(ImageLog const& i1, ImageLog const& i2) {
    return (i1.MJD < i2.MJD || (i1.MJD == i2.MJD && stringnmatch01(i1.obscode, i2.obscode, 3) == -1));
}
inline bool compareXYIndexX(XYIndex const& p1, XYIndex const& p2) { return (p1.x < p2.x); }
inline bool compareXYIndexY(XYIndex const& p1, XYIndex const& p2) { return (p1.y < p2.y); }
inline bool comparePoint3DIndexX(Point3DIndex const& p1, Point3DIndex const& p2) { return (p1.x < p2.x); }
inline bool compareLongIndex(LongIndex const& i1, LongIndex const& i2) { return (i1.lelem < i2.lelem); }

std::vector<Observatory> readObscodeFile(string obscodefile);

std::vector<Detection> readDetectionsFile(
    string indetfile,
    int idcol,
    int mjdcol,
    int racol,
    int deccol,
    int magcol,
    int bandcol,
    int obscodecol
);

std::vector<ImageLog> readImageFile(
    MakeTrackletsConfig const& config,
    std::vector<Detection> const& detvec,
    string inimfile
);

std::vector<ImageLog> makeImageLogs(
    MakeTrackletsConfig const& config,
    std::vector<Detection> const& detvec
);

std::tuple<std::vector<long double>, std::vector<Point3LD>, std::vector<Point3LD>> readEarthEphemerides(
    string earthfile
);

void computeHelioPositions(
    std::vector<Detection> &detvec,
    std::vector<ImageLog> const& img_log,
    std::vector<Observatory> const& observatory_list,
    std::vector<long double> const& EarthMJD,
    std::vector<Point3LD> const& Earthpos
);

std::tuple<std::vector<Detection>, std::vector<LongPair>> buildTracklets(
    MakeTrackletsConfig const& config,
    std::vector<Detection> &detvec,
    std::vector<ImageLog> &img_log
);

void refineTracklets(
    MakeTrackletsConfig const& config,
    std::vector<Detection> &pairdets,
    string outpairfile
);

#endif // HELA_MAKETRACKLETS_H