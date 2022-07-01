#include "makeTracklets.h"
#include "cmath"

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

int main(int argc, char *argv[]) {
    MakeTrackletsConfig config;
    Detection o1 = Detection(0L, 0l, 0l, 0L, 0L, 0L, "null", 0l, "V", "I11", 0);
    vector<Detection> detvec = {};
    vector<Detection> pairdets = {};
    vector<Detection> ppset = {};
    std::vector<std::vector<int>> indvecs = {};
    Observatory obs1 = Observatory("Ill", 0l, 0l, 0l);
    vector<Observatory> observatory_list = {};
    ImageLog imlog = ImageLog(0.0, 0.0, 0.0, "I11", 0, 0);
    vector<ImageLog> img_log_tmp = {};
    vector<ImageLog> img_log = {};
    vector<LongPair> pairvec = {};
    Point3D p3 = Point3D(0.0, 0.0, 0.0);
    Point3D p3avg = Point3D(0.0, 0.0, 0.0);
    vector<Point3LD> Earthpos;
    vector<Point3LD> Earthvel;
    vector<Point3LD> observer_heliopos;
    vector<long double> EarthMJD;
    Point3LD outpos = Point3LD(0.0, 0.0, 0.0);
    double tdelt = 0;
    double mjdmean = 0;
    double mjdnorm = 0;
    char idstring[SHORTSTRINGLEN];
    char band[MINSTRINGLEN];
    char obscode[MINSTRINGLEN];
    string lnfromfile;
    int status = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int imct = 0;
    int imctp = 0;
    int imnum = 0;
    long detnum = 0;
    long num_dets = 0;
    long detct = 0;
    int startind = 0;
    int endind = 0;
    int reachedeof = 0;
    char c = '0';
    long double MJD, RA, Dec;
    MJD = RA = Dec = 0.0L;
    double mag = 0l;
    double maxvel = config.maxvel;  // Max angular velocity in deg/day
    double minvel = config.minvel;    // Min angular velocity in deg/day
    double minarc = config.minarc;    // Min total angular arc in arcseconds
    double angvel = config.angvel;
    double maxtime = config.maxtime;           // Max time interval a tracklet could span,
                                        // in days.
    double maxdist = maxvel * maxtime;  // Max angular distance a tracklet
                                        // could span, in degrees.
    double imrad = config.imagerad;            // radius from image center to most distant corner (deg).
    string indetfile;
    string inimfile;
    string outimfile;
    string earthfile;
    string obscodefile;
    string colformatfile;
    string outpairfile = "outpairfile01.txt";
    string pairdetfile = "pairdetfile01.txt";
    double obslon = 289.26345L;
    double plxcos = 0.865020L;
    double plxsin = -0.500901L;
    long lct = 0;
    double mintime = config.mintime;
    double maxgcr = config.maxgcr;
    int idcol = 1;
    int mjdcol = 3;
    int racol = 6;
    int deccol = 8;
    int magcol = 32;
    int bandcol = 26;
    int obscodecol = 38;
    int cols_to_read = 7;
    int colreadct = 0;
    ifstream instream1;
    ofstream outstream1;
    string stest;
    int mintrkpts = config.mintrkpts;

    if (argc < 7) {
        show_usage();
        return (1);
    }

    i = 1;
    while (i < argc) {
        cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
        if (string(argv[i]) == "-d" || string(argv[i]) == "-dets" || string(argv[i]) == "-det" ||
            string(argv[i]) == "--dets" || string(argv[i]) == "--det" || string(argv[i]) == "--detection" ||
            string(argv[i]) == "--detections") {
            if (i + 1 < argc) {
                // There is still something to read;
                indetfile = argv[++i];
                i++;
            } else {
                cerr << "Detection file keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" ||
                   string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" ||
                   string(argv[i]) == "--image" || string(argv[i]) == "--images") {
            if (i + 1 < argc) {
                // There is still something to read;
                inimfile = argv[++i];
                i++;
            } else {
                cerr << "Image file keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-outim" || string(argv[i]) == "-outimg" ||
                   string(argv[i]) == "-outimgs" || string(argv[i]) == "--outimages" ||
                   string(argv[i]) == "--outimage" || string(argv[i]) == "--outimgs" ||
                   string(argv[i]) == "--outimg") {
            if (i + 1 < argc) {
                // There is still something to read;
                outimfile = argv[++i];
                i++;
            } else {
                cerr << "Output image file keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-pairs" || string(argv[i]) == "-pairfile") {
            if (i + 1 < argc) {
                // There is still something to read;
                outpairfile = argv[++i];
                i++;
            } else {
                cerr << "Output pair file keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" ||
                   string(argv[i]) == "-detpairs") {
            if (i + 1 < argc) {
                // There is still something to read;
                pairdetfile = argv[++i];
                i++;
            } else {
                cerr << "Output paired detection file keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-imrad") {
            if (i + 1 < argc) {
                // There is still something to read;
                imrad = stod(argv[++i]);
                config.imagerad = imrad;
                i++;
                if (!isnormal(imrad) || imrad <= 0.0) {
                    cerr << "Error: invalid image radius (" << imrad << " deg) supplied.\n";
                    cerr << "Image radius must be strictly positive!\n";
                    return (2);
                }
            } else {
                cerr << "Output image radius keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-maxtime") {
            if (i + 1 < argc) {
                // There is still something to read;
                maxtime = stod(argv[++i]);
                config.maxtime = maxtime;
                i++;
                if (isnormal(maxtime) && maxtime > 0.0) {
                    maxtime /= 24.0;  // Convert from hours to days.
                    config.maxtime /= 24.0;
                } else {
                    cerr << "Error: invalid maximum inter-image time interval\n";
                    cerr << "(" << maxtime << " hr) supplied: must be strictly positive.\n";
                    return (2);
                }
            } else {
                cerr << "Maximum inter-image time interval\nkeyword supplied with no corresponding "
                        "argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-mintime") {
            if (i + 1 < argc) {
                // There is still something to read;
                mintime = stod(argv[++i]);
                config.mintime = mintime;
                i++;
                if ((isnormal(mintime) || mintime == 0.0) && mintime >= 0.0) {
                    mintime /= 24.0;  // Convert from hours to days
                } else {
                    cerr << "Error: invalid minimum inter-image time interval\n";
                    cerr << "(" << mintime << " hr) supplied: must be non-negative.\n";
                    return (2);
                }
            } else {
                cerr << "Minimum inter-image time interval\nkeyword supplied with no corresponding "
                        "argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-minvel") {
            if (i + 1 < argc) {
                // There is still something to read;
                minvel = stod(argv[++i]);
                config.minvel = minvel;
                i++;
                if (!isnormal(minvel) && minvel != 0.0l) {
                    cerr << "Error: invalid minimum angular velocity\n";
                    cerr << "(" << minvel << "deg/day) supplied.\n";
                    return (2);
                }
            } else {
                cerr << "Minimum angular velocity\nkeyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-maxvel") {
            if (i + 1 < argc) {
                // There is still something to read;
                maxvel = stod(argv[++i]);
                config.maxvel = maxvel;
                i++;
                if (!isnormal(maxvel) || maxvel <= 0.0) {
                    cerr << "Error: invalid maximum angular velocity\n";
                    cerr << "(" << maxvel << "deg/day) supplied: must be strictly positive.\n";
                    return (2);
                }
            } else {
                cerr << "Maximum angular velocity\nkeyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-maxGCR" || string(argv[i]) == "-maxgcr") {
            if (i + 1 < argc) {
                // There is still something to read;
                maxgcr = stod(argv[++i]);
                config.maxgcr = maxgcr;
                i++;
                if (!isnormal(maxgcr) || maxgcr <= 0.0) {
                    cerr << "Error: invalid maximum Great Circle residual\n";
                    cerr << "(" << maxgcr << " arcsec) supplied: must be strictly positive.\n";
                    return (2);
                }
            } else {
                cerr << "Output maximum Great Circle Residual\nkeyword supplied with no corresponding "
                        "argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-minarc") {
            if (i + 1 < argc) {
                // There is still something to read;
                minarc = stod(argv[++i]);
                config.minarc = minarc;
                i++;
                if (!isnormal(minarc) && minarc != 0.0l) {
                    cerr << "Error: invalid minimum angular arc\n";
                    cerr << "(" << minarc << " arcsec) supplied.\n";
                    return (2);
                }
            } else {
                cerr << "Minimum angular arc\nkeyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-earth" || string(argv[i]) == "-e" || string(argv[i]) == "-Earth" ||
                   string(argv[i]) == "--earthfile" || string(argv[i]) == "--Earthfile" ||
                   string(argv[i]) == "--earth" || string(argv[i]) == "--Earth") {
            if (i + 1 < argc) {
                // There is still something to read;
                earthfile = argv[++i];
                i++;
            } else {
                cerr << "Earth file keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-mintrkpts" || string(argv[i]) == "-mpt" ||
                   string(argv[i]) == "-mintrackpts" || string(argv[i]) == "-minpts" ||
                   string(argv[i]) == "--minimumtrackletpoints" || string(argv[i]) == "--mintrackpoints" ||
                   string(argv[i]) == "--mintrackletpoints") {
            if (i + 1 < argc) {
                // There is still something to read;
                mintrkpts = stoi(argv[++i]);
                config.mintrkpts = mintrkpts;
                i++;
            } else {
                cerr << "Earth file keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-obscode" || string(argv[i]) == "-obs" || string(argv[i]) == "-oc" ||
                   string(argv[i]) == "-obscodes" || string(argv[i]) == "--obscode" ||
                   string(argv[i]) == "--obscodes" || string(argv[i]) == "--observatorycodes") {
            if (i + 1 < argc) {
                // There is still something to read;
                obscodefile = argv[++i];
                i++;
            } else {
                cerr << "Observatory code file keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else if (string(argv[i]) == "-colformat" || string(argv[i]) == "-format" ||
                   string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" ||
                   string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" ||
                   string(argv[i]) == "--cformat") {
            if (i + 1 < argc) {
                // There is still something to read;
                colformatfile = argv[++i];
                i++;
            } else {
                cerr << "Column format file keyword supplied with no corresponding argument\n";
                show_usage();
                return (1);
            }
        } else {
            cerr << "Warning: unrecognized keyword " << argv[i] << "\n";
            i++;
        }
    }

    if (indetfile.size() <= 0) {
        cerr << "Please supply an input detection file:\n\n";
        show_usage();
        return (1);
    }

    if (earthfile.size() <= 0) {
        cerr << "Please supply a heliocentric ephemeris file for the Earth:\n\n";
        show_usage();
        return (1);
    }

    if (obscodefile.size() <= 0) {
        cerr << "Please supply a observatory code file:\n\n";
        show_usage();
        return (1);
    }

    if (mintrkpts < 2) {
        mintrkpts = 2;
        config.mintrkpts = mintrkpts;
    }

    cout << "indet file " << indetfile << "\n";
    cout << "inimage file " << inimfile << "\n";
    cout << "column formatting file " << colformatfile << "\n";
    cout << "observatory code file " << obscodefile << "\n";
    cout << "output image file " << outimfile << "\n";
    cout << "pairfile file " << outpairfile << "\n";
    cout << "paired detection file " << pairdetfile << "\n";
    cout << "Heliocentric ephemeris file for the Earth: " << earthfile << "\n";
    cout << "image radius " << imrad << " degrees\n";
    cout << "max time interval " << maxtime * 24.0 << " hours\n";
    cout << "min time interval " << mintime * 24.0 << " hours\n";
    cout << "minvel " << minvel << " deg/day\n";
    cout << "maxvel " << maxvel << " deg/day\n";
    cout << "minimum number of points per tracklet " << mintrkpts << "\n";
    cout << "minarc " << minarc << " arcsec\n";
    cout << "maxGCR " << maxgcr << " arcsec\n";
    maxdist = maxtime * maxvel;
    config.maxdist = maxdist;

    // Read the column formatting file, if any
    if (colformatfile.size() > 0) {
        instream1.open(colformatfile);
        if (!instream1) {
            cerr << "ERROR: unable to open input file " << colformatfile << "\n";
            return (1);
        }
        colreadct = 0;
        while (!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct < cols_to_read) {
            instream1 >> stest;
            if (stest == "IDCOL") {
                instream1 >> idcol;
                if (!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
            } else if (stest == "MJDCOL") {
                instream1 >> mjdcol;
                if (!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
            } else if (stest == "RACOL") {
                instream1 >> racol;
                if (!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
            } else if (stest == "DECCOL") {
                instream1 >> deccol;
                if (!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
            } else if (stest == "MAGCOL") {
                instream1 >> magcol;
                if (!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
            } else if (stest == "BANDCOL") {
                instream1 >> bandcol;
                if (!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
            } else if (stest == "OBSCODECOL") {
                instream1 >> obscodecol;
                if (!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
            } else {
                cout << "WARNING: unrecognized string " << stest << " read from column formatting file\n";
            }
        }
        instream1.close();
        if (colreadct < cols_to_read) {
            cout << "WARNING: only " << colreadct << " column specifications, of " << cols_to_read
                 << " expected, were read from column format file " << colformatfile << ".\n";
        }
    }

    cout << "Column specifications:\n";
    cout << "IDCOL " << idcol << "\n";
    cout << "MJDCOL " << mjdcol << "\n";
    cout << "RACOL " << racol << "\n";
    cout << "DECCOL " << deccol << "\n";
    cout << "MAGCOL " << magcol << "\n";
    cout << "BANDCOL " << bandcol << "\n";
    cout << "OBSCODECOL " << obscodecol << "\n";

    // Read observatory code file
    observatory_list = readObscodeFile(obscodefile);
    if (observatory_list.size() == 0) {
        return (1);
    }

    // Read input detection file.
    detvec = readDetectionsFile(indetfile, idcol, mjdcol, racol, deccol,
                                magcol, bandcol, obscodecol);
    if (detvec.size() == 0) {
        return (1);
    }

    // Print out sorted detvec for debugging.
    // for(i=0;i<detvec.size();i++) {
    //    cout << "det# " << i << ": " << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << " "
    //    << detvec[i].idstring << " " << detvec[i].mag << " " << detvec[i].band << " " << detvec[i].obscode
    //    << " " << detvec[i].index << "\n";
    //  }

    // Get image information.
    if (inimfile.size() > 0) {
        img_log = readImageFile(config, detvec, inimfile);
        if (img_log.size() == 0) {
            return (1);
        }
    } else {
        // Make image file
        img_log = makeImageLogs(config, detvec);
    }

    // Read in Earth ephemerides
    std::tie(EarthMJD, Earthpos, Earthvel) = readEarthEphemerides(earthfile);
    cout << "Finished reading heliocentric ephemeris file " << earthfile << " for the Earth.\n";

    // Output for testing
    if (DEBUG >= 2) {
        // Test: print out time-sorted detection table.
        outstream1.open("testjunk01.txt");
        for (i = 0; i < detvec.size(); i++) {
            outstream1 << fixed << setprecision(6) << detvec[i].MJD << " " << detvec[i].RA << " "
                       << detvec[i].Dec << " " << detvec[i].mag << " " << detvec[i].band << " "
                       << detvec[i].obscode << "\n";
        }
        outstream1.close();
    }

    if (outimfile.size() > 0) {
        // Write and print image log table
        outstream1.open(outimfile);
        for (imct = 0; imct < img_log.size(); imct++) {
            //	  cout << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
            //	  cout << " " << img_log[imct].Dec << " " << img_log[imct].startind << " " <<
            //img_log[imct].endind << " "; 	  cout << img_log[imct].endind-img_log[imct].startind << "\n";
            outstream1 << fixed << setprecision(6) << img_log[imct].MJD << " " << img_log[imct].RA;
            outstream1 << fixed << setprecision(6) << " " << img_log[imct].Dec << " " << img_log[imct].obscode
                       << " ";
            outstream1 << img_log[imct].startind << " " << img_log[imct].endind << "\n";
        }
        outstream1.close();
    }

    // Compute the heliocentric positions for each detection at time of exposure.
    computeHelioPositions(detvec, img_log, observatory_list, EarthMJD, Earthpos);

    std::tie(pairdets, pairvec, indvecs) = buildTracklets(config, detvec, img_log);

    cout << "Writing paired detections file\n";
    outstream1.open(pairdetfile);
    outstream1 << "#MJD,RA,Dec,observerX,observerY,observerZ,stringID,mag,band,obscode,origindex\n";
    for (size_t i = 0; i < pairdets.size(); i++) {
        outstream1 << fixed << setprecision(7) << pairdets[i].MJD << "," << pairdets[i].RA << ","
                   << pairdets[i].Dec << ",";
        outstream1 << fixed << setprecision(3) << pairdets[i].x << "," << pairdets[i].y << ","
                   << pairdets[i].z << ",";
        outstream1 << fixed << setprecision(3) << pairdets[i].idstring << "," << pairdets[i].mag << ","
                   << pairdets[i].band << ",";
        outstream1 << pairdets[i].obscode << "," << pairdets[i].index << "\n";
    }
    outstream1.close();

    // This does some output for now; need to disentangle.
    refineTracklets(config, pairdets, indvecs, outpairfile);

    return 0;
}