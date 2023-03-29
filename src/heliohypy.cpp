// March 21, 2023: heliohypy (heliocentric hypothesis code for python)
// Implementation of make_tracklets for python.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void fill_struct(hldet & out, hldet const& in) {
    out.MJD = in.MJD;
    out.RA = in.RA;
    out.Dec = in.Dec;
    out.mag = in.mag;
    out.trail_len = in.trail_len;
    out.trail_PA = in.trail_PA;
    out.sigmag = in.sigmag;
    out.sig_across = in.sig_across;
    out.sig_along = in.sig_along;
    out.image = in.image;
    memcpy(out.idstring, in.idstring, sizeof(in.idstring));
    memcpy(out.band, in.band, sizeof(in.band));
    memcpy(out.obscode, in.obscode, sizeof(in.obscode));
    out.index = in.index;
}

void fill_struct(EarthState & out, EarthState const& in) {
    out.MJD = in.MJD;
    out.x = in.x;
    out.y = in.y;
    out.z = in.z;
    out.vx = in.vx;
    out.vy = in.vy;
    out.vz = in.vz;
}

void fill_struct(hlimage & out, hlimage const& in) {
    out.MJD = in.MJD;
    out.RA = in.RA;
    out.Dec = in.Dec;
    memcpy(out.obscode, in.obscode, sizeof(in.obscode));
    out.X = in.X;
    out.Y = in.Y;
    out.Z = in.Z;
    out.VX = in.VX;
    out.VY = in.VY;
    out.VZ = in.VZ;
}

void fill_struct(tracklet & out, tracklet const& in) {
    out.Img1 = in.Img1;
    out.RA1 = in.RA1;
    out.Dec1 = in.Dec1;
    out.Img2 = in.Img2;
    out.RA2 = in.RA2;
    out.Dec2 = in.Dec2;
    out.npts = in.npts;
    out.trk_ID = in.trk_ID;
}

void fill_struct(longpair & out, longpair const& in) {
    out.i1 = in.i1;
    out.i2 = in.i2;
}

void fill_struct(hlradhyp & out, hlradhyp const& in) {
    out.HelioRad = in.HelioRad;
    out.R_dot = in.R_dot;
    out.R_dubdot = in.R_dubdot;
}

void fill_struct(hlclust & out, hlclust const& in) {
    out.clusternum = in.clusternum;
    out.posRMS = in.posRMS;
    out.astromRMS = in.astromRMS;
    out.totRMS = in.totRMS;
    out.pairnum = in.pairnum;
    out.timespan = in.timespan;
    out.uniquepoints = in.uniquepoints;
    out.obsnights = in.obsnights;
    out.metric = in.metric;
    memcpy(out.rating, in.rating, sizeof(in.rating));
    out.heliohyp0 = in.heliohyp0;
    out.heliohyp1 = in.heliohyp1;
    out.heliohyp2 = in.heliohyp2;
    out.posX = in.posX;
    out.posY = in.posY;
    out.posZ = in.posZ;
    out.velX = in.velX;
    out.velY = in.velY;
    out.velZ = in.velZ;
    out.orbit_a = in.orbit_a;
    out.orbit_e = in.orbit_e;
    out.orbit_MJD = in.orbit_MJD;
    out.orbitX = in.orbitX;
    out.orbitY = in.orbitY;
    out.orbitZ = in.orbitZ;
    out.orbitVX = in.orbitVX;
    out.orbitVY = in.orbitVY;
    out.orbitVZ = in.orbitVZ;
    out.orbit_eval_count = in.orbit_eval_count;
}

template<typename T>
std::vector<T> ndarray_to_vec(py::array_t<T> py_vec) {
    std::vector<T> vec = {};

    // Get a reference to the py_array data
    auto data_ref = py_vec.unchecked();

    // Place the numpy data into the c++ type
    for (long int i = 0; i < data_ref.size(); i++) {
        T data_out;
        auto &data_in = data_ref[i];

        fill_struct(data_out, data_in);

        vec.push_back(data_out);
    }

    return vec;
}

template<typename T>
py::array vec_to_ndarray(std::vector<T> const& vec) {
    // Allocate a structured numpy array of type T
    auto py_vec = py::array_t<T>(vec.size());

    // Get a mutable reference to the ndarray data
    auto data_ref = py_vec.mutable_unchecked();

    // Place vector data into numpy array
    for (long int i = 0; i < data_ref.size(); i++) {
        auto data_in = vec[i];
        auto &data_out = data_ref[i];

        fill_struct(data_out, data_in);
    }
    return py_vec;
}

py::array iotest02(py::array_t<hldet> py_ioin)
{
  std::vector <hldet> iovec = ndarray_to_vec(py_ioin);
  long unsigned int i=0;
  for(i=0; i<iovec.size(); i++) {
    std::cout << iovec[i].MJD << " " << iovec[i].RA << " " << iovec[i].Dec << " " << iovec[i].mag << " " << iovec[i].trail_len << " " << iovec[i].trail_PA << " " << iovec[i].sigmag  << " " << iovec[i].sig_across << " " << iovec[i].sig_along << " " << iovec[i].image << " " << iovec[i].idstring << " " << iovec[i].band << " " << iovec[i].obscode << " " << iovec[i].index << "\n";
    iovec[i].MJD += MJDOFF;
  }

  auto py_ioout = vec_to_ndarray<hldet>(iovec);
  return(py_ioout);
}

py::array observer_coords(double detmjd, double lon, double obscos, double obssine, py::array_t<EarthState> earthin)
{
  std::vector <EarthState> earthpos = ndarray_to_vec(earthin);
  int polyorder=5;
  std::vector <double> posmjd;
  std::vector <point3d> planetpos;
  int i=0;
  int earthnum = earthpos.size();
  point3d earthnow = point3d(0,0,0);
  std::vector <double> outvec;
  
  posmjd={};
  planetpos={};
  for(i=0; i<earthnum; i++) {
    posmjd.push_back(earthpos[i].MJD);
    earthnow = point3d(earthpos[i].x, earthpos[i].y, earthpos[i].z);
    planetpos.push_back(earthnow);
  }  
  observer_barycoords01(detmjd, polyorder, lon, obscos, obssine, posmjd, planetpos, earthnow);
  outvec={};
  outvec.push_back(earthnow.x);
  outvec.push_back(earthnow.y);
  outvec.push_back(earthnow.z);
  
  auto py_outvec = py::array_t<double>(outvec.size());
  auto data_ref = py_outvec.mutable_unchecked();
  // Place vector data into numpy array
  for (long int i = 0; i < data_ref.size(); i++) {
    auto data_in = outvec[i];
    auto &data_out = data_ref[i];
    data_out = data_in;
  }

  return(py_outvec);
}

py::array observer_vel(double detmjd, double lon, double obscos, double obssine, py::array_t<EarthState> earthin)
{
  std::vector <EarthState> earthpos = ndarray_to_vec(earthin);
  int polyorder=5;
  std::vector <double> posmjd;
  std::vector <point3d> planetpos;
  std::vector <point3d> planetvel;
  int i=0;
  int earthnum = earthpos.size();
  point3d earthnow = point3d(0,0,0);
  point3d earthvel = point3d(0,0,0);
  std::vector <double> outvec;

  
  posmjd={};
  planetpos={};
  planetvel={};
  for(i=0; i<earthnum; i++) {
    posmjd.push_back(earthpos[i].MJD);
    earthnow = point3d(earthpos[i].x, earthpos[i].y, earthpos[i].z);
    planetpos.push_back(earthnow);
    earthvel = point3d(earthpos[i].vx, earthpos[i].vy, earthpos[i].vz);
    planetvel.push_back(earthvel);
  }  
  observer_baryvel01(detmjd, polyorder, lon, obscos, obssine, posmjd, planetpos, planetvel, earthnow, earthvel);

  outvec={};
  outvec.push_back(earthnow.x);
  outvec.push_back(earthnow.y);
  outvec.push_back(earthnow.z);
  outvec.push_back(earthvel.x);
  outvec.push_back(earthvel.y);
  outvec.push_back(earthvel.z);
 
  auto py_outvec = py::array_t<double>(outvec.size());
  auto data_ref = py_outvec.mutable_unchecked();
  // Place vector data into numpy array
  for (long int i = 0; i < data_ref.size(); i++) {
    auto data_in = outvec[i];
    auto &data_out = data_ref[i];
    data_out = data_in;
  }

  return(py_outvec);
}

std::tuple<py::array, py::array, py::array> makeTracklets(
    MakeTrackletsConfig config,
    py::array_t<hldet> py_detvec,
    py::array_t<hlimage> py_imglog
  ) {
  cout << "C++ wrapper for make_tracklets, now fully functional\n";
  
  std::vector <hldet> detvec = ndarray_to_vec(py_detvec);
  std::vector <hlimage> image_log = ndarray_to_vec(py_imglog);
  long unsigned int i=0;
  std::vector <longpair> pairvec;
  std::vector <vector <long>> indvecs;
  std::vector <hldet> pairdets;
  std::vector <tracklet> tracklets;
  std::vector <longpair> trk2det;  
   
  // Echo config struct
  cout << "Configuration parameters:\n";
  cout << "Min. number of tracklet points: " << config.mintrkpts << "\n";
  cout << "Time-tolerance for matching detections on the same image: " << config.imagetimetol << " days (" << config.imagetimetol*SOLARDAY << " seconds)\n";
  cout << "Maximum angular velocity: " << config.maxvel << " deg/day\n";
  cout << "Minimum angular velocity: " << config.minvel << " deg/day\n";
  cout << "Minimum angular arc: " << config.minarc << " arcsec\n";
  cout << "Maximum inter-image time interval: " << config.maxtime << " days (" << config.maxtime*1440.0 << " minutes)\n";
  cout << "Minimum inter-image time interval: " << config.mintime << " days (" << config.mintime*1440.0 << " minutes)\n";
  cout << "Image radius: " << config.imagerad << " degrees\n";
  cout << "Maximum Great Circle Residual for tracklets with more than two points: " << config.maxgcr << " arcsec\n";
  if(config.forcerun) {
    cout << "forcerun has been invoked: execution will attempt to push through\n";
    cout << "any errors that are not immediately fatal, including those that\n";
    cout << "could produce inaccurate final results.\n";
  }
  
  int status = load_image_indices(image_log, detvec, config.imagetimetol, config.forcerun);
  if(status!=0) {
    cerr << "ERROR: failed to load_image_indices from detection vector\n";
    auto py_detout = vec_to_ndarray<hldet>({});
    return(std::make_tuple(py_detout, py_detout, py_detout));
  }
  
  // Echo detection vector
  for(i=0;i<detvec.size();i++) {
    cout << "det " << i << " " << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << " " << detvec[i].mag  << " " << detvec[i].obscode << " " << detvec[i].image << "\n";
  }
  // Echo image log
  for(i=0;i<image_log.size();i++) {
    cout << "image " << i << " " << image_log[i].MJD << " " << image_log[i].RA << " " << image_log[i].Dec << " " << image_log[i].X << " " << image_log[i].obscode  << " " << image_log[i].startind  << " " << image_log[i].endind << "\n";
  }

  // Create pairs, output a vector pairdets of type hldet; a vector indvecs of type vector <long>,
  // with the same length as pairdets, giving the indices of all the detections paired with a given detection;
  // and the vector pairvec of type longpair, giving all the pairs of detections.
  status = find_pairs(detvec, image_log, pairdets, indvecs, pairvec, config.mintime, config.maxtime, config.imagerad, config.maxvel);
  
  if(status!=0) {
    cerr << "ERROR: find_pairs reports failure status " << status << "\n";
    auto py_detout = vec_to_ndarray<hldet>({});
    return(std::make_tuple(py_detout, py_detout, py_detout));
  }
  status = merge_pairs(pairdets, indvecs, pairvec, tracklets, trk2det, config.mintrkpts, config.maxgcr, config.minarc, config.minvel, config.maxvel);

  auto py_detout1 = vec_to_ndarray<hldet>(pairdets);
  auto py_detout2 = vec_to_ndarray<tracklet>(tracklets);
  auto py_detout3 = vec_to_ndarray<longpair>(trk2det);

  return(std::make_tuple(py_detout1, py_detout2, py_detout3));
}


#define TIMECONVSCALE 4.0L // The characteristic timescale used to convert velocities
                           // to distance units is equal to the full temporal span
                           // divided by TIMECONVSCALE.
#define INTRANIGHTSTEP 0.3 // Minimum interval in days between successive points
                           // in a tracklet, to enable them to be counted as being
                           // on separate nights.
#define INTEGERIZING_SCALE 100.0L // We divide state vectors by this value to integerize
                                   // them. Given a standard integer with a range
                                   // of +/- 2^31 = 2.15e9, this gives the state vectors
                                   // a range of +/- 2.15e11 km = 1400 AU, which is
                                   // comfortably larger than the solar system.
#define CRAD_REF_GEODIST 1.0 // Value of geocentric distance to which the user-defined
                             // clustering radius is normalized (AU). In general, the
                             // clustering radius is scaled linearly with geocentric distance.

std::tuple<py::array, py::array> heliolinc(
    HeliolincConfig config,
    py::array_t<hlimage> py_imglog,
    py::array_t<tracklet> py_tracklets,
    py::array_t<longpair> py_trk2det,
    py::array_t<hlradhyp> py_radhyp,
    py::array_t<EarthState> earthin
  ) {
  cout << "C++ wrapper for heliolinc\n";
  
  std::vector <hlimage> image_log = ndarray_to_vec(py_imglog);
  std::vector <tracklet> tracklets = ndarray_to_vec(py_tracklets);
  std::vector <longpair> trk2det = ndarray_to_vec(py_trk2det);
  std::vector <hlradhyp> radhyp = ndarray_to_vec(py_radhyp);
  std::vector <EarthState> earthpos = ndarray_to_vec(earthin);
  int polyorder=5;
  point3d Earthrefpos = point3d(0l,0l,0l);

  long imnum = image_log.size();
  long pairnum = tracklets.size();
  long trk2detnum = trk2det.size();
  long accelnum = radhyp.size();
  long accelct=0;

  vector <double> heliodist;
  vector <double> heliovel;
  vector <double> helioacc;
  vector <hlclust> outclust;
  vector <longpair> clust2det;
  long realclusternum, gridpoint_clusternum, status;
  realclusternum = gridpoint_clusternum = status = 0;
  vector <point6ix2> allstatevecs;

  // Echo config struct
  cout << "Configuration parameters:\n";
  cout << "MJD of reference time: " << config.MJDref << "\n";
  cout << "DBSCAN clustering radius: " << config.clustrad << " km\n";
  cout << "DBSCAN npt: " << config.dbscan_npt << "\n";
  cout << "Min number of distinct observing nights for a valid linkage: " << config.minobsnights << "\n";
  cout << "Min time span for a valid linkage: " << config.mintimespan << " days\n";
  cout << "Min geocentric distance (center of innermost bin): " << config.mingeodist << " AU\n";
  cout << "Max geocentric distance (will be exceeded by center only of the outermost bin): " << config.maxgeodist << " AU\n";
  cout << "Logarthmic step size (and bin width) for geocentric distance bins: " << config.geologstep << "\n";
  cout << "Minimum inferred geocentric distance for a valid tracklet: " << config.mingeoobs << " AU\n";
  cout << "Minimum inferred impact parameter (w.r.t. Earth) for a valid tracklet: " << config.minimpactpar << " Earth radii\n";
  if(config.verbose) cout << "Verbose output selected\n";
  
  if(imnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty image catalog\n";
    auto py_clustout = vec_to_ndarray<hlclust>({});
    return(std::make_tuple(py_clustout, py_clustout));
  } else if(pairnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty tracklet array\n";
    auto py_clustout = vec_to_ndarray<hlclust>({});
    return(std::make_tuple(py_clustout, py_clustout));
  } else if(trk2detnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty trk2det array\n";
    auto py_clustout = vec_to_ndarray<hlclust>({});
    return(std::make_tuple(py_clustout, py_clustout));
  } else if(accelnum<=0) {
    cerr << "ERROR: heliolinc supplied with empty heliocentric hypothesis array\n";
    auto py_clustout = vec_to_ndarray<hlclust>({});
    return(std::make_tuple(py_clustout, py_clustout));
  }
  
  double MJDmin = image_log[0].MJD;
  double MJDmax = image_log[imnum-1].MJD;
  if(config.MJDref<MJDmin || config.MJDref>MJDmax) {
    // Input reference MJD is invalid. Suggest a better value before exiting.
    cerr << "\nERROR: reference MJD must fall in the time interval spanned by the data\n";
    cerr << fixed << setprecision(2) << "Suggested value is " << MJDmin*0.5l + MJDmax*0.5l << "\n";
    cout << "based on your input image catalog\n";
    auto py_clustout = vec_to_ndarray<hlclust>({});
    return(std::make_tuple(py_clustout, py_clustout));
  }

  double chartimescale = (MJDmax - MJDmin)*SOLARDAY/TIMECONVSCALE; // Note that the units are seconds.
  point3d earthpos01(const vector <EarthState> &earthpos, double mjd)
{
  point3d earthnow = point3d(0,0,0);

  // Convert heliocentric radial motion hypothesis matrix
  // from units of AU, AU/day, and GMSun/R^2
  // to units of km, km/day, and km/day^2.
  heliodist = heliovel = helioacc = {};
  for(accelct=0; accelct<accelnum; accelct++) {
    heliodist.push_back(radhyp[accelct].HelioRad * AU_KM);
    heliovel.push_back(radhyp[accelct].R_dot * AU_KM);
    helioacc.push_back(radhyp[accelct].R_dubdot * (-GMSUN_KM3_SEC2*SOLARDAY*SOLARDAY/heliodist[accelct]/heliodist[accelct]));
  }

  // Begin master loop over heliocentric hypotheses
  outclust={};
  clust2det={};
  realclusternum=0;  
  for(accelct=0;accelct<accelnum;accelct++) {
    gridpoint_clusternum=0;
    
    // Covert all tracklets into state vectors at the reference time, under
    // the assumption that the heliocentric distance hypothesis is correct.
    status = trk2statevec(image_log, tracklets, heliodist[accelct], heliovel[accelct], helioacc[accelct], chartimescale, allstatevecs, config.MJDref, config.mingeoobs, config.minimpactpar);
    
    if(status==1) {
      cerr << "WARNING: hypothesis " << accelct << ": " << radhyp[accelct].HelioRad << " " << radhyp[accelct].R_dot << " " << radhyp[accelct].R_dubdot << " led to\nnegative heliocentric distance or other invalid result: SKIPPING\n";
      continue;
    } else if(status==2) {
      // This is a weirder error case and is fatal.
      cerr << "Fatal error case from trk2statevec.\n";
      auto py_clustout = vec_to_ndarray<hlclust>({});
      return(std::make_tuple(py_clustout, py_clustout));
    }
    // If we get here, trk2statevec probably ran OK.
    if(allstatevecs.size()<=1) continue; // No clusters possible, skip to the next step.
    if(config.verbose>=0) cout << pairnum << " input pairs/tracklets led to " << allstatevecs.size() << " physically reasonable state vectors\n";

    status = form_clusters(allstatevecs, tracklets, trk2det, Earthrefpos, heliodist[accelct], heliovel[accelct], helioacc[accelct], outclust, clust2det, config.mingeodist, config.geologstep, config.maxgeodist, config.mintimespan, config.minobsnights, config.verbose);
  }
    
  auto py_detout1 = vec_to_ndarray<hlclust>(outclust);
  auto py_detout2 = vec_to_ndarray<longpair>(clust2det);

  return(std::make_tuple(py_detout1, py_detout2));
}

#undef TIMECONVSCALE
#undef INTRANIGHTSTEP
#undef INTEGERIZING_SCALE
#undef CRAD_REF_GEODIST


PYBIND11_MODULE(heliohypy, m) {
    m.doc() = "pybind11 I/O test"; // optional module docstring
    
    PYBIND11_NUMPY_DTYPE(hldet, MJD, RA, Dec, mag, trail_len, trail_PA, sigmag, sig_across, sig_along, image, idstring, band, obscode, index);
    PYBIND11_NUMPY_DTYPE(EarthState, MJD, x, y, z, vx, vy, vz);
    PYBIND11_NUMPY_DTYPE(hlimage, MJD, RA, Dec, obscode, X, Y, Z, VX, VY, VZ, startind, endind);
    PYBIND11_NUMPY_DTYPE(longpair, i1, i2);
    PYBIND11_NUMPY_DTYPE(tracklet, Img1, RA1, Dec1, Img2, RA2, Dec2, npts, trk_ID);
    PYBIND11_NUMPY_DTYPE(hlradhyp, HelioRad, R_dot, R_dubdot);
    PYBIND11_NUMPY_DTYPE(hlclust, clusternum, posRMS, velRMS, totRMS, astromRMS, pairnum, timespan, uniquepoints, obsnights, metric, rating, heliohyp0, heliohyp1, heliohyp2, posX, posY, posZ, velX, velY, velZ, orbit_a, orbit_e, orbit_MJD, orbitX, orbitY, orbitZ, orbitVX, orbitVY, orbitVZ, orbit_eval_count);

    
    py::class_<hldet>(m, "hldet")
      .def(py::init<double &, double &, double &, float &, float &, float &, float &, float &, float &, int &, std::string &, std::string &, std::string &, long &>());
    
    py::class_<EarthState>(m, "EarthState")
      .def(py::init<double &, double &, double &, double &, double &, double &, double &>());
    
    py::class_<hlimage>(m, "hlimage")
      .def(py::init<double &, double &, double &, std::string &, double &, double &, double &, double &, double &, double &, int &, int &>());

    py::class_<longpair>(m, "longpair")
      .def(py::init<long &, long &>());

    py::class_<tracklet>(m, "tracklet")
      .def(py::init<long &, double &, double &, long &, double &, double &, int &, long &>());

    py::class_<hlradhyp>(m, "hlradhyp")
      .def(py::init<double &, double &, double &>());
    
    py::class_<hlclust>(m, "hlclust")
      .def(py::init<long &, double &, double &, double &, double &, int &, double &, int &, int &, double &, std::string &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, long &>());


    // Config class for MakeTracklets
    py::class_<MakeTrackletsConfig>(m, "MakeTrackletsConfig")
      .def(py::init<>())
      .def_readwrite("mintrkpts", &MakeTrackletsConfig::mintrkpts)
      .def_readwrite("imagetimetol", &MakeTrackletsConfig::imagetimetol)
      .def_readwrite("maxvel", &MakeTrackletsConfig::maxvel)
      .def_readwrite("minvel", &MakeTrackletsConfig::minvel)
      .def_readwrite("minarc", &MakeTrackletsConfig::minarc)
      .def_readwrite("maxtime", &MakeTrackletsConfig::maxtime)
      .def_readwrite("mintime", &MakeTrackletsConfig::mintime)
      .def_readwrite("imagerad", &MakeTrackletsConfig::imagerad)
      .def_readwrite("maxgcr", &MakeTrackletsConfig::maxgcr)
      .def_readwrite("forcerun", &MakeTrackletsConfig::forcerun);

    // Config class for Heliolinc
    py::class_<HeliolincConfig>(m, "HeliolincConfig")
      .def(py::init<>())
      .def_readwrite("MJDref", &HeliolincConfig::MJDref)
      .def_readwrite("clustrad", &HeliolincConfig::clustrad) 
      .def_readwrite("dbscan_npt", &HeliolincConfig::dbscan_npt) 
      .def_readwrite("minobsnights", &HeliolincConfig::minobsnights) 
      .def_readwrite("mintimespan", &HeliolincConfig::mintimespan)
      .def_readwrite("mingeodist", &HeliolincConfig::mingeodist)
      .def_readwrite("maxgeodist", &HeliolincConfig::maxgeodist)
      .def_readwrite("geologstep", &HeliolincConfig::geologstep)
      .def_readwrite("mingeoobs", &HeliolincConfig::mingeoobs)
      .def_readwrite("minimpactpar", &HeliolincConfig::minimpactpar);
      .def_readwrite("verbose", &HeliolincConfig::verbose);

    
    m.def("iotest02", &iotest02, "A function to test python I/O");
    m.def("observer_coords", &observer_coords, "calculate position of an observer on Earth");
    m.def("observer_vel", &observer_vel, "calculate position of an observer on Earth");
    m.def("makeTracklets", &makeTracklets,  "Make tracklets from set of detections.");
 }


