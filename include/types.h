// Common types
#ifndef HELA_TYPES_H
#define HELA_TYPES_H

#include <cassert>
#include <string>
#include <cstring>
#include <vector>

using namespace std;

class Detection {
public:
    long double MJD;
    double RA;
    double Dec;
    char idstring[SHORTSTRINGLEN];
    double mag;
    char band[MINSTRINGLEN];
    char obscode[MINSTRINGLEN];
    long double x;
    long double y;
    long double z;
    long index;

    Detection() = default;

    Detection(long double mjd, double ra, double dec, long double x, long double y, long double z,
              string const& idstring, double mag, string const& band, string const& obscode, long index)
            : MJD(mjd), RA(ra), Dec(dec), mag(mag), x(x), y(y), z(z), index(index) {
        // Copy input value of idstring, making sure it's not too long
        assert(idstring.size() < sizeof(this->idstring));
        std::strncpy(this->idstring, idstring.c_str(), sizeof(this->idstring));
        this->idstring[sizeof(this->idstring) - 1] = 0;

        // Copy input value for photometric band, making sure it's not too long
        assert(band.size() < sizeof(this->band));
        std::strncpy(this->band, band.c_str(), sizeof(this->band));
        this->band[sizeof(this->band) - 1] = 0;

        // Copy input value for obscode, making sure it's not too long
        assert(obscode.size() < sizeof(this->obscode));
        std::strncpy(this->obscode, obscode.c_str(), sizeof(this->obscode));
        this->obscode[sizeof(this->obscode) - 1] = 0;
    }
};

class Observatory {
public:
    char obscode[MINSTRINGLEN];
    double obslon;
    double plxcos;
    double plxsin;

    Observatory() = default;

    Observatory(string const& obscode, double obslon, double plxcos, double plxsin)
            : obslon(obslon), plxcos(plxcos), plxsin(plxsin) {
        // Copy input value of obscode, making sure it's not too long
        assert(obscode.size() < sizeof(this->obscode));
        std::strncpy(this->obscode, obscode.c_str(), sizeof(this->obscode));
        this->obscode[sizeof(this->obscode) - 1] = 0;
    }
};

class ImageLog {
public:
    double MJD;
    double RA;
    double Dec;
    char obscode[MINSTRINGLEN];
    long startind;
    long endind;

    ImageLog() = default;

    ImageLog(double mjd, double ra, double dec, const string& obscode, long startind, long endind)
            : MJD(mjd), RA(ra), Dec(dec), startind(startind), endind(endind) {
        // Copy input value of obscode, making sure it's not too long
        assert(obscode.size() < sizeof(this->obscode));
        std::strncpy(this->obscode, obscode.c_str(), sizeof(this->obscode));
        this->obscode[sizeof(this->obscode) - 1] = 0;
    }
};

struct EarthState {
    long double MJD;
    long double x;
    long double y;
    long double z;
    long double vx;
    long double vy;
    long double vz;

    EarthState() = default;

    EarthState(long double MJD, long double x, long double y, long double z,
               long double vx, long double vy, long double vz)
            : MJD(MJD), x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) {}
};

struct LongPair {  // Pair of long integers
    long i1;
    long i2;
    LongPair(long i1, long i2) : i1(i1), i2(i2) {}
};

template <typename T>
struct Point3 {  // 3-D point in arbitrary precision
    T x;
    T y;
    T z;
    Point3<T>(T x, T y, T z) : x(x), y(y), z(z) {}
};

struct Point3DIndex {
    double x;
    double y;
    double z;
    long index;
    Point3DIndex(double x, double y, double z, long index) : x(x), y(y), z(z), index(index) {}
};

struct LongIndex {
    long lelem;
    long index;
    LongIndex(long lelem, long index) : lelem(lelem), index(index) {}
};

struct XYIndex {
    double x;
    double y;
    long index;
    XYIndex(double x, double y, long index) : x(x), y(y), index(index) {}
};

struct KDPoint {
    XYIndex point;
    long left;
    long right;
    int dim;
    KDPoint(XYIndex point, long left, long right, int dim)
            : point(point), left(left), right(right), dim(dim) {}
};

#endif // HELA_TYPES_H