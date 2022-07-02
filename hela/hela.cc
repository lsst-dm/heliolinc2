#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "makeTracklets.h"

namespace py = pybind11;

void fill_struct(Observatory & out, Observatory const& in, long i) {
    memcpy(out.obscode, in.obscode, sizeof(in.obscode));
    out.obslon = in.obslon;
    out.plxcos = in.plxcos;
    out.plxsin = in.plxsin;
}

void fill_struct(Detection & out, Detection const& in, long i) {
    out.MJD = in.MJD;
    out.RA = in.RA;
    out.Dec = in.Dec;
    memcpy(out.idstring, in.idstring, sizeof(in.idstring));
    out.mag = in.mag;
    memcpy(out.band, in.band, sizeof(in.band));
    memcpy(out.obscode, in.obscode, sizeof(in.obscode));
    out.index = -i-2;
}

template<typename T>
std::vector<T> ndarray_to_vec(py::array_t<T> py_vec) {
    std::vector<T> vec = {};

    // Get a reference to the py_array data
    auto data_ref = py_vec.unchecked();

    // Place the numpy data into the c++ type
    for (long i = 0; i < data_ref.size(); i++) {
        T data_out;
        auto &data_in = data_ref[i];

        fill_struct(data_out, data_in, i);

        vec.push_back(data_out);
    }

    return vec;
}

void makeTracklets(
    MakeTrackletsConfig config
) {
    cout << "Testing C++ wrapper" << endl;

    int idcol = 1;
    int mjdcol = 3;
    int racol = 6;
    int deccol = 8;
    int magcol = 32;
    int bandcol = 26;
    int obscodecol = 38;

    cout << "Reading obscode file... ";
    std::vector<Observatory> observatory_list = readObscodeFile(config.obscodefile);
    cout << "Done." << endl;

    cout << "Reading detections file... ";
    std::vector<Detection> detvec = readDetectionsFile(config.indetfile, idcol, mjdcol, racol, deccol, magcol, bandcol, obscodecol);
    cout << "Done." << endl;

    cout << "Read/Make image file... ";
    std::vector<ImageLog> img_log = {};
    if (config.inimfile.size() > 0) {
        cout << "Reading file... ";
        img_log = readImageFile(config, detvec, config.inimfile);
        cout << "Done." << endl;
    } else{
        cout << "Making image info from detections... ";
        img_log = makeImageLogs(config, detvec);
        cout << "Done." << endl;
    }
    cout << "Done." << endl;

    cout << "Read Earth ephemerides... ";
    std::vector<long double> EarthMJD = {};
    std::vector<Point3LD> Earthpos = {};
    std::vector<Point3LD> Earthvel = {};
    std::tie(EarthMJD, Earthpos, Earthvel) = readEarthEphemerides(config.earthfile);
    cout << "Done." << endl;

    cout << "Writing test files...";
    // Output for testing
    ofstream outstream;
    if (DEBUG >= 2) {
        // Test: print out time-sorted detection table.
        outstream.open("testjunk01.txt");
        for (size_t i = 0; i < detvec.size(); i++) {
            outstream << fixed << setprecision(6) << detvec[i].MJD << " " << detvec[i].RA << " "
                       << detvec[i].Dec << " " << detvec[i].mag << " " << detvec[i].band << " "
                       << detvec[i].obscode << "\n";
        }
        outstream.close();
    }

    // Hard-code files for now
    string outimfile = "imfile_month04a.txt";
    if (outimfile.size() > 0) {
        // Write and print image log table
        outstream.open(outimfile);
        for (size_t imct = 0; imct < img_log.size(); imct++) {
            //	  cout << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
            //	  cout << " " << img_log[imct].Dec << " " << img_log[imct].startind << " " <<
            // img_log[imct].endind << " "; 	  cout << img_log[imct].endind-img_log[imct].startind << "\n";
            outstream << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
            outstream << fixed << setprecision(6) << " " << img_log[imct].Dec << " " << img_log[imct].obscode
                       << " ";
            outstream << img_log[imct].startind << " " << img_log[imct].endind << "\n";
        }
        outstream.close();
    }
    cout << "Done." << endl;

    cout << "Compute heliocentric positions... ";
    computeHelioPositions(detvec, img_log, observatory_list, EarthMJD, Earthpos);
    cout << "Done." << endl;

    cout << "Make tracklets... ";
    std::vector<LongPair> pairvec = {};
    std::vector<Detection> pairdets = {};
    std::vector<std::vector<int>> indvecs = {};
    std::tie(pairdets, pairvec, indvecs) = buildTracklets(config, detvec, img_log);
    cout << "Done." << endl;

    cout << "Writing paired detections file... ";
    string pairdetfile = "pairdets_sol_month04a.csv";
    outstream.open(pairdetfile);
    outstream << "#MJD,RA,Dec,observerX,observerY,observerZ,stringID,mag,band,obscode,origindex\n";
    for (size_t i = 0; i < pairdets.size(); i++) {
        outstream << fixed << setprecision(7) << pairdets[i].MJD << "," << pairdets[i].RA << ","
                   << pairdets[i].Dec << ",";
        outstream << fixed << setprecision(3) << pairdets[i].x << "," << pairdets[i].y << ","
                   << pairdets[i].z << ",";
        outstream << fixed << setprecision(3) << pairdets[i].idstring << "," << pairdets[i].mag << ","
                   << pairdets[i].band << ",";
        outstream << pairdets[i].obscode << "," << pairdets[i].index << "\n";
    }
    outstream.close();
    cout << "Done." << endl;

    cout << "Refining tracklets... ";
    string outpairfile = "pair_sol_month04a.csv";
    refineTracklets(config, pairdets, indvecs, outpairfile);
    cout << "Done." << endl;

    return;
}

void makeTracklets(
    MakeTrackletsConfig config,
    py::array_t<Observatory> py_obsv,
    py::array_t<Detection> py_detvec
) {
    cout << "Testing C++ wrapper" << endl;

    //int idcol = 1;
    //int mjdcol = 3;
    //int racol = 6;
    //int deccol = 8;
    //int magcol = 32;
    //int bandcol = 26;
    //int obscodecol = 38;

    cout << "Reading obscode file... ";
    std::vector<Observatory> observatory_list = ndarray_to_vec(py_obsv);
    //std::vector<Observatory> observatory_list = readObscodeFile(config.obscodefile);
    cout << "Done." << endl;
    fflush(stdout);

    cout << "Reading detections file... ";
    std::vector<Detection> detvec = ndarray_to_vec(py_detvec);
    sort(detvec.begin(), detvec.end(), timeCompareDetections);
    //std::vector<Detection> detvec =
    //        readDetectionsFile(config.indetfile, idcol, mjdcol, racol, deccol, magcol, bandcol, obscodecol);
    cout << "Done." << endl;
    fflush(stdout);

    cout << "Read/Make image file... ";
    std::vector<ImageLog> img_log = {};
    if (config.inimfile.size() > 0) {
        cout << "Reading file... ";
        img_log = readImageFile(config, detvec, config.inimfile);
        cout << "Done." << endl;
    } else {
        cout << "Making image info from detections... ";
        img_log = makeImageLogs(config, detvec);
        cout << "Done." << endl;
    }
    fflush(stdout);

    cout << "Read Earth ephemerides... ";
    std::vector<long double> EarthMJD = {};
    std::vector<Point3LD> Earthpos = {};
    std::vector<Point3LD> Earthvel = {};
    std::tie(EarthMJD, Earthpos, Earthvel) = readEarthEphemerides(config.earthfile);
    cout << "Done." << endl;
    fflush(stdout);

    cout << "Writing test files...";
    // Output for testing
    ofstream outstream;
    if (DEBUG >= 2) {
        // Test: print out time-sorted detection table.
        outstream.open("testjunk01.txt");
        for (size_t i = 0; i < detvec.size(); i++) {
            outstream << fixed << setprecision(6) << detvec[i].MJD << " " << detvec[i].RA << " "
                      << detvec[i].Dec << " " << detvec[i].mag << " " << detvec[i].band << " "
                      << detvec[i].obscode << "\n";
        }
        outstream.close();
    }

    // Hard-code files for now
    string outimfile = "imfile_month04a.txt";
    if (outimfile.size() > 0) {
        // Write and print image log table
        outstream.open(outimfile);
        for (size_t imct = 0; imct < img_log.size(); imct++) {
            //	  cout << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
            //	  cout << " " << img_log[imct].Dec << " " << img_log[imct].startind << " " <<
            // img_log[imct].endind << " "; 	  cout << img_log[imct].endind-img_log[imct].startind << "\n";
            outstream << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
            outstream << fixed << setprecision(6) << " " << img_log[imct].Dec << " " << img_log[imct].obscode
                      << " ";
            outstream << img_log[imct].startind << " " << img_log[imct].endind << "\n";
        }
        outstream.close();
    }
    cout << "Done." << endl;
    fflush(stdout);

    cout << "Compute heliocentric positions... ";
    computeHelioPositions(detvec, img_log, observatory_list, EarthMJD, Earthpos);
    cout << "Done." << endl;
    fflush(stdout);

    cout << "Make tracklets... ";
    std::vector<LongPair> pairvec = {};
    std::vector<Detection> pairdets = {};
    std::vector<std::vector<int>> indvecs = {};
    std::tie(pairdets, pairvec, indvecs) = buildTracklets(config, detvec, img_log);
    cout << "Done." << endl;
    fflush(stdout);

    cout << "Writing paired detections file... ";
    string pairdetfile = "pairdets_sol_month04a.csv";
    outstream.open(pairdetfile);
    outstream << "#MJD,RA,Dec,observerX,observerY,observerZ,stringID,mag,band,obscode,origindex\n";
    for (size_t i = 0; i < pairdets.size(); i++) {
        outstream << fixed << setprecision(7) << pairdets[i].MJD << "," << pairdets[i].RA << ","
                  << pairdets[i].Dec << ",";
        outstream << fixed << setprecision(3) << pairdets[i].x << "," << pairdets[i].y << "," << pairdets[i].z
                  << ",";
        outstream << fixed << setprecision(3) << pairdets[i].idstring << "," << pairdets[i].mag << ","
                  << pairdets[i].band << ",";
        outstream << pairdets[i].obscode << "," << pairdets[i].index << "\n";
    }
    outstream.close();
    cout << "Done." << endl;
    fflush(stdout);

    cout << "Refining tracklets... ";
    string outpairfile = "pair_sol_month04a.csv";
    refineTracklets(config, pairdets, indvecs, outpairfile);
    cout << "Done." << endl;
    fflush(stdout);

    return;
}

PYBIND11_MODULE(hela, m) {
    m.doc() = "HelioLINC Advanced (hela)";

    // Structures/Classes
    PYBIND11_NUMPY_DTYPE(Observatory, obscode, obslon, plxcos, plxsin);
    PYBIND11_NUMPY_DTYPE(Detection, MJD, RA, Dec, x, y, z, idstring, mag, band, obscode, index);
    // PYBIND11_NUMPY_DTYPE(ImageLog, MJD, RA, Dec, obscode, startind, endind);

    // Config class
    py::class_<MakeTrackletsConfig>(m, "MakeTrackletsConfig")
        .def(py::init<>())
        .def_readwrite("mintrkpts", &MakeTrackletsConfig::mintrkpts)
        .def_readwrite("imagetimetol", &MakeTrackletsConfig::imagetimetol)
        .def_readwrite("maxvel", &MakeTrackletsConfig::maxvel)
        .def_readwrite("minvel", &MakeTrackletsConfig::minvel)
        .def_readwrite("minarc", &MakeTrackletsConfig::minarc)
        .def_readwrite("maxtime", &MakeTrackletsConfig::maxtime)
        .def_readwrite("mintime", &MakeTrackletsConfig::mintime)
        .def_readwrite("angvel", &MakeTrackletsConfig::angvel)
        .def_readwrite("maxdist", &MakeTrackletsConfig::maxdist)
        .def_readwrite("imagerad", &MakeTrackletsConfig::imagerad)
        .def_readwrite("maxgcr", &MakeTrackletsConfig::maxgcr)
        .def_readwrite("indetfile", &MakeTrackletsConfig::indetfile)
        .def_readwrite("inimfile", &MakeTrackletsConfig::inimfile)
        .def_readwrite("earthfile", &MakeTrackletsConfig::earthfile)
        .def_readwrite("obscodefile", &MakeTrackletsConfig::obscodefile);

    m.def("makeTracklets", py::overload_cast<MakeTrackletsConfig>(&makeTracklets), "Make tracklets from set of detections.");
    m.def(
        "makeTracklets",
        py::overload_cast<MakeTrackletsConfig, py::array_t<Observatory>, py::array_t<Detection>>(&makeTracklets),
        "Make tracklets from set of detections."
    );
}