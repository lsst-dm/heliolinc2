#include "solarsyst_dyn_geo01.h"
#include <tuple>

#define DEBUG 0

// make some type aliases for now
typedef det_obsmag_indvec Detection;
typedef img_log03 ImageLog;

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

struct StateVector {
    long double mjd = 0.0;
    point3LD pos = point3LD(0, 0, 0);
    point3LD vel = point3LD(0, 0, 0);
};

std::vector<observatory> readObscodeFile(string obscodefile);

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

std::tuple<std::vector<long double>, std::vector<point3LD>, std::vector<point3LD>> readEarthEphemerides(
    string earthfile
);

void computeHelioPositions(
    std::vector<Detection> &detvec,
    std::vector<ImageLog> const& img_log,
    std::vector<observatory> const& observatory_list,
    std::vector<long double> const& EarthMJD,
    std::vector<point3LD> const& Earthpos
);

std::tuple<std::vector<Detection>, std::vector<longpair>> buildTracklets(
    MakeTrackletsConfig const& config,
    std::vector<Detection> &detvec,
    std::vector<ImageLog> &img_log
);

void refineTracklets(
    MakeTrackletsConfig const& config,
    std::vector<Detection> &pairdets,
    string outpairfile
);

/* These are used types/functions

TYPES
[] det_obsmag_indvec
[] observatory
[] img_log03
[] longpair
[] point3d
[] point3LD
[] long_index
[] xy_index
[] point3d_index
[] kdpoint
[] medindex

METHODS
[] distradec01
[] distradec02
[] kdtree01
    [] splitxy
[] kdrange01
[] linfituw01
[] arc2cel01

# Methods for fast sorting
[] early_det_obsmag_indvec
[] lower_long_index
[] lower_point3d_index_x
[] early_imlg3
[] lower_point6LDx2_x,y,z,vx,vy,vz

# Methods used in IO
[] read_horizons_fileLD
[] stringncopy01
[] stringnmatch01
[] celeproj01
[] celedeproj01
[] obscode_lookup
[] observer_barycoords01LD
    [] precess01a
    [] celestial_to_statevec
    [] planetpos01
        [] make_dvec
        [] make_dmat
        [] perfectpoly01
            [] solvematrix01


*/