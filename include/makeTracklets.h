#include "solarsyst_dyn_geo01.h"

#define DEBUG 0

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
};

void buildTracklets(
    MakeTrackletsConfig config,
    std::vector<det_obsmag_indvec> detvec,
    std::vector<img_log03> img_log,
    std::vector<longpair> &pairvec,
    std::vector<det_obsmag_indvec> &pairdets
);

void refineTracklets(
    MakeTrackletsConfig config,
    std::vector<det_obsmag_indvec> &pairdets,
    string outpairfile
);

static void show_usage() {
    cerr << "Usage: maketrack04b -dets detfile -imgs imfile -outimgs output image file/ \n";
    cerr << "-pairs pairfile -pairdets paired detection file -colformat column format file/ \n";
    cerr << "-imrad image radius(deg) -maxtime max inter-image time interval (hr)/ \n";
    cerr << "-mintime min inter-image time interval (hr) -maxGCR maximum GRC -mintrkpts min. num. of "
            "tracklet points/\n";
    cerr << "-minvel minimum angular velocity (deg/day) -maxvel maximum angular velocity (deg/day) \n";
    cerr << "-minarc minimum total angular arc (arcsec) -earth earthfile -obscode obscodefile\n";
    cerr << "\nor, at minimum\n\n";
    cerr << "maketrack04b -dets detfile -earth earthfile -obscode obscodefile\n";
    cerr << "Note well that the minimum invocation will leave a bunch of things\n";
    cerr << "set to defaults that may not be what you want.\n";
};

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