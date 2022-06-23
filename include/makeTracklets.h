#include "solarsyst_dyn_geo01.h"

int const NUMPOS = 3;
int const IDCOL = 1;
int const MJDCOL = 3;
int const RACOL = 6;
int const DECCOL = 8;
int const MAGCOL = 32;
int const BANDCOL = 26;
int const OBSCODECOL = 38;
int const COLS_TO_READ = 7;
double const IMAGETIMETOL = 1.0;    // Tolerance for matching image time, in seconds
double const MAXVEL = 1.5;          // Default max angular velocity in deg/day.
double const MAXTIME = 1.5 / 24.0;  // Default max inter-image time interval
                                    // for tracklets, in days.
double const IMAGERAD = 2.0;        // radius from image center to most distant corner (deg)
double const MAX_GCR = 0.5;         // Default maximum Great Circle Residual allowed for a valid tracklet

#define DEBUG 0

void makeTracklets(
    std::vector<det_obsmag_indvec> detvec,
    std::vector<img_log03> img_log,
    std::vector<longpair> &pairvec,
    std::vector<det_obsmag_indvec> &pairdet,
    string outpairfile,
    string pairdetfile
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