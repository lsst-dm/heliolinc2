#include "makeTracklets.h"
#include "cmath"

int main(int argc, char *argv[]) {
    det_obsmag_indvec o1 = det_obsmag_indvec(0L, 0l, 0l, 0L, 0L, 0L, "null", 0l, "V", "I11", 0, {});
    vector<det_obsmag_indvec> detvec = {};
    vector<det_obsmag_indvec> pairdets = {};
    vector<det_obsmag_indvec> ppset = {};
    observatory obs1 = observatory("Ill", 0l, 0l, 0l);
    vector<observatory> observatory_list = {};
    img_log03 imlog = img_log03(0.0, 0.0, 0.0, "I11", 0, 0);
    vector<img_log03> img_log_tmp = {};
    vector<img_log03> img_log = {};
    longpair onepair = longpair(0, 0);
    vector<longpair> pairvec = {};
    point3d p3 = point3d(0, 0, 0);
    point3d p3avg = point3d(0, 0, 0);
    vector<point3LD> Earthpos;
    vector<point3LD> Earthvel;
    vector<point3LD> observer_heliopos;
    vector<long double> EarthMJD;
    point3LD outpos = point3LD(0, 0, 0);
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
    double maxvel = MAXVEL;  // Max angular velocity in deg/day
    double minvel = 0.0l;    // Min angular velocity in deg/day
    double minarc = 0.0l;    // Min total angular arc in arcseconds
    double angvel = 0.0l;
    double maxtime = MAXTIME;           // Max time interval a tracklet could span,
                                        // in days.
    double maxdist = MAXVEL * MAXTIME;  // Max angular distance a tracklet
                                        // could span, in degrees.
    double imrad = IMAGERAD;            // radius from image center to most distant corner (deg).
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
    long_index ppn = long_index(0, 0);
    vector<long_index> pair_partner_num = {};
    vector<long_index> tracklet_check = {};
    double dt, dtref, dx, dy;
    dt = dtref = dx = dy = 0.0l;
    xy_index xyind = xy_index(0.0, 0.0, 0);
    vector<xy_index> axyvec = {};
    double dist, pa;
    dist = pa = 0.0;
    int dettarg = 0;
    vector<double> timevec;
    vector<double> xvec;
    vector<double> yvec;
    vector<long> detindexvec;
    int biggest_tracklet = -1;
    int tracklet_size = 0;
    double slopex, slopey, interceptx, intercepty, worsterr;
    vector<double> fiterr = {};
    vector<double> fiterr2 = {};
    int worstpoint = -1;
    int istracklet = 0;
    int rp1, rp2, instep;
    rp1 = rp2 = instep = 0;
    double outra1, outra2, outdec1, outdec2;
    outra1 = outra2 = outdec1 = outdec2 = 0.0l;
    point3d_index p3di = point3d_index(0l, 0l, 0l, 0);
    vector<point3d_index> track_mrdi_vec;
    double mintime = IMAGETIMETOL / SOLARDAY;
    int trkptnum, istimedup = 1;
    double maxgcr = MAX_GCR;
    int idcol = IDCOL;
    int mjdcol = MJDCOL;
    int racol = RACOL;
    int deccol = DECCOL;
    int magcol = MAGCOL;
    int bandcol = BANDCOL;
    int obscodecol = OBSCODECOL;
    int colreadct = 0;
    ifstream instream1;
    ofstream outstream1;
    string stest;
    int mintrkpts = 2;

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
                i++;
                if (isnormal(maxtime) && maxtime > 0.0) {
                    maxtime /= 24.0;  // Convert from hours to days.
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

    if (mintrkpts < 2) mintrkpts = 2;

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

    // Read the column formatting file, if any
    if (colformatfile.size() > 0) {
        instream1.open(colformatfile);
        if (!instream1) {
            cerr << "ERROR: unable to open input file " << colformatfile << "\n";
            return (1);
        }
        colreadct = 0;
        while (!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct < COLS_TO_READ) {
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
        if (colreadct < COLS_TO_READ) {
            cout << "WARNING: only " << colreadct << " column specifications, of " << COLS_TO_READ
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
    instream1.open(obscodefile);
    if (!instream1) {
        cerr << "can't open input file " << obscodefile << "\n";
        return (1);
    }
    // Skip one-line header
    getline(instream1, lnfromfile);
    while (!instream1.eof() && !instream1.fail() && !instream1.bad()) {
        instream1 >> stest;
        stringncopy01(obscode, stest, MINSTRINGLEN);
        instream1 >> obslon;
        instream1 >> plxcos;
        instream1 >> plxsin;
        obs1 = observatory(obscode, obslon, plxcos, plxsin);
        observatory_list.push_back(obs1);
        // Skip the rest of the line
        getline(instream1, lnfromfile);
    }
    instream1.close();
    cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile
         << ":\n";

    if (DEBUG >= 2) {
        for (i = 0; i < observatory_list.size(); i++) {
            cout << observatory_list[i].obscode << " " << observatory_list[i].obslon << " "
                 << observatory_list[i].plxcos << " " << observatory_list[i].plxsin << "\n";
        }
    }

    // Read input detection file.
    instream1.open(indetfile);
    if (!instream1) {
        cerr << "can't open input file " << indetfile << "\n";
        return (1);
    }
    // Skip one-line header
    getline(instream1, lnfromfile);
    lct++;
    cout << lnfromfile << "\n";
    while (reachedeof == 0) {
        getline(instream1, lnfromfile);
        lct++;
        if (!instream1.eof() && !instream1.fail() && !instream1.bad())
            ;  // Read on.
        else if (instream1.eof())
            reachedeof = 1;  // End of file, fine.
        else if (instream1.fail())
            reachedeof = -1;  // Something wrong, warn
        else if (instream1.bad())
            reachedeof = -2;  // Worse problem, warn
        i = 0;
        j = 0;
        c = '0';
        while (i < lnfromfile.size() && lnfromfile.size() >= 30 && reachedeof == 0) {
            // Note check on line length: it is completely impossible for a
            // line containing all the required quantities at minimum plausible
            // precision to be less than 30 characters long.
            c = '0';
            stest = "";
            while (i < lnfromfile.size() && c != ',' && c != '\n' && c != EOF) {
                c = lnfromfile[i];
                if (c != ',' && c != '\n' && c != EOF) stest.push_back(c);
                i++;
            }
            // We just finished reading something
            j++;
            if (j == idcol) stringncopy01(idstring, stest, SHORTSTRINGLEN);
            if (j == mjdcol)
                MJD = stold(stest);
            else if (j == racol)
                RA = stold(stest);
            else if (j == deccol)
                Dec = stold(stest);
            else if (j == magcol)
                mag = stod(stest);
            else if (j == bandcol)
                stringncopy01(band, stest, MINSTRINGLEN);
            else if (j == obscodecol)
                stringncopy01(obscode, stest, MINSTRINGLEN);
            // cout<<"Column "<< j << " read as " << stest << ".\n";
        }
        if (reachedeof == 0 && lnfromfile.size() >= 30) {
            // cout<<"MJD, RA, Dec: " << MJD-floor(MJD) << " " << RA << " " << Dec << "\n";
            o1 = det_obsmag_indvec(MJD, RA, Dec, 0L, 0L, 0L, idstring, mag, band, obscode, -lct, {});
            detvec.push_back(o1);
        }
    }
    instream1.close();
    if (reachedeof == 1) {
        cout << "File read successfully to the end.\n";
    } else if (reachedeof == -1)
        cout << "Warning: file read failed\n";
    else if (reachedeof == -2)
        cout << "Warning: file possibly corrupted\n";
    else
        cout << "Warning: unknown file read problem\n";

    cout << "Last two obscodes: " << detvec[detvec.size() - 2].obscode << " and "
         << detvec[detvec.size() - 1].obscode << "\n";

    // time-sort the detection vector
    sort(detvec.begin(), detvec.end(), early_det_obsmag_indvec());

    cout << "Last two obscodes: " << detvec[detvec.size() - 2].obscode << " and "
         << detvec[detvec.size() - 1].obscode << "\n";

    // Print out sorted detvec for debugging.
    // for(i=0;i<detvec.size();i++) {
    //    cout << "det# " << i << ": " << detvec[i].MJD << " " << detvec[i].RA << " " << detvec[i].Dec << " "
    //    << detvec[i].idstring << " " << detvec[i].mag << " " << detvec[i].band << " " << detvec[i].obscode
    //    << " " << detvec[i].index << "\n";
    //  }

    // Get image information.
    if (inimfile.size() > 0) {
        // Read input image file: MJD, RA, Dec, obscode:
        instream1.open(inimfile);
        if (!instream1) {
            cerr << "can't open input file " << inimfile << "\n";
            return (1);
        }
        reachedeof = 0;
        while (reachedeof == 0) {
            getline(instream1, lnfromfile);
            if (!instream1.eof() && !instream1.fail() && !instream1.bad())
                ;  // Read on.
            else if (instream1.eof())
                reachedeof = 1;  // End of file, fine.
            else if (instream1.fail())
                reachedeof = -1;  // Something wrong, warn
            else if (instream1.bad())
                reachedeof = -2;  // Worse problem, warn
            i = 0;
            j = 0;
            c = '0';
            MJD = 0.0;
            while (i < lnfromfile.size() && reachedeof == 0) {
                stest = "";
                c = '0';
                while (i < lnfromfile.size() && c != ',' && c != ' ' && c != '\n' && c != EOF) {
                    // We allow the file to be delimited by comma or space.
                    c = lnfromfile[i];
                    if (c != ',' && c != ' ' && c != '\n' && c != EOF) stest.push_back(c);
                    i++;
                }
                // We just finished reading something
                j++;
                if (j == 1)
                    MJD = stod(stest);  // We assume we have MJD, RA, Dec, obscode
                else if (j == 2)
                    RA = stod(stest);
                else if (j == 3)
                    Dec = stod(stest);
                else if (j == 4)
                    stringncopy01(obscode, stest, MINSTRINGLEN);
            }
            if ((reachedeof == 0 || reachedeof == 1) && MJD > 0.0) {
                // Requirement of MJD>0.0 tests that we read a plausibly
                // valid line.
                imlog = img_log03(MJD, RA, Dec, obscode, 0, 0);
                img_log_tmp.push_back(imlog);
            }
        }
        if (reachedeof == 1) {
            cout << "File read successfully to the end.\n";
        } else if (reachedeof == -1)
            cout << "Warning: file read failed\n";
        else if (reachedeof == -2)
            cout << "Warning: file possibly corrupted\n";
        else
            cout << "Warning: unknown file read problem\n";
        // time-sort the image file
        sort(img_log_tmp.begin(), img_log_tmp.end(), early_imlg3());
        // find the indices in the time-sorted detection file
        // that correspond to the earliest and latest detections
        // on each image, and load these values into imglog02.
        detct = 0;
        for (imct = 0; imct < img_log_tmp.size(); imct++) {
            while (detct < detvec.size() &&
                   detvec[detct].MJD < img_log_tmp[imct].MJD - IMAGETIMETOL / SOLARDAY)
                detct++;  // Not on any image
            if (detct < detvec.size() &&
                fabs(detvec[detct].MJD - img_log_tmp[imct].MJD) <= IMAGETIMETOL / SOLARDAY) {
                // This should be the first detection on image imct.
                img_log_tmp[imct].startind = detct;
                while (detct < detvec.size() &&
                       fabs(detvec[detct].MJD - img_log_tmp[imct].MJD) <= IMAGETIMETOL / SOLARDAY)
                    detct++;  // Still on this same image
                // This should be the first detection on the next image
                img_log_tmp[imct].endind = detct;
            }
            if (img_log_tmp[imct].startind >= 0 && img_log_tmp[imct].endind > 0) {
                img_log.push_back(img_log_tmp[imct]);
            }
        }
        instream1.close();
    } else {
        // No input image file was supplied: we have to create one from
        // the sorted detection file.
        mjdnorm = 1.0;
        mjdmean = detvec[0].MJD;
        startind = 0;
        for (i = 1; i < detvec.size(); i++) {
            tdelt = detvec[i].MJD - detvec[i - 1].MJD;
            if (tdelt < IMAGETIMETOL / SOLARDAY &&
                stringnmatch01(detvec[i].obscode, detvec[i - 1].obscode, 3) == 0) {
                // This point corresponds to the same image as the previous one.
                mjdmean += detvec[i].MJD;
                mjdnorm += 1.0;
            } else {
                // Now we are considering a new image.
                // Calculate the meanmjd of the previous image, for which
                //  we have now seen all points.
                // Record the current detct i as the detection index just
                // after the end of the previous image
                endind = i;
                if (isnormal(mjdnorm))
                    mjdmean /= mjdnorm;
                else
                    mjdmean = 0.0;
                // Load it into the vector with mean MJD for all images,
                //  and increment image count.
                imlog = img_log03(mjdmean, 0.0, 0.0, detvec[endind - 1].obscode, startind, endind);
                img_log.push_back(imlog);
                // Set up for the next image, starting with detvec[i].MJD;
                mjdmean = detvec[i].MJD;
                mjdnorm = 1.0;
                startind = i;
            }
        }
        cout << "OK here, now for final image\n";
        fflush(stdout);
        // Account for the final image.
        if (isnormal(mjdnorm)) {
            endind = i;
            mjdmean /= mjdnorm;
            // Load it into the vector with mean MJD for all images,
            //  and increment image count.
            imlog = img_log03(mjdmean, 0.0, 0.0, detvec[endind - 1].obscode, startind, endind);
            img_log.push_back(imlog);
        }
        cout << "Done with final image\n";
        fflush(stdout);

        // We've now loaded the mean MJDs and the starting and ending
        // detection table indices for each image; it still remains to
        // get the mean RA and Dec.

        detnum = detvec.size();
        imnum = img_log.size();
        cout << img_log.size() << " unique images were identified.\n";
        cout << "Given our total of " << detvec.size() << " detections,\n";
        cout << "we have " << double(detvec.size()) / double(img_log.size())
             << " detections per image, on average\n";

        // Find the number of detections and the average RA, Dec on each image.
        // We perform the average after projection onto the unit circle, to
        // avoid wrapping issues.
        detct = imct = 0;
        while (imct < imnum && detct < detnum) {
            vector<det_obsmag_indvec> imobs = {};
            num_dets = 0;
            p3avg = point3d(0, 0, 0);
            while (detct < detnum && detvec[detct].MJD < img_log[imct].MJD + IMAGETIMETOL / SOLARDAY) {
                num_dets++;                      // Keep count of detections on this image
                imobs.push_back(detvec[detct]);  // Load vector of observations for this image
                p3 = celeproj01(detvec[detct].RA, detvec[detct].Dec);  // Project current detection
                p3avg.x += p3.x;
                p3avg.y += p3.y;
                p3avg.z += p3.z;  // Average projected coords
                detct++;
            }
            // If we got here, we must just have finished with an image.
            // Finish the average
            if (num_dets > 0) {
                p3avg.x /= double(num_dets);
                p3avg.y /= double(num_dets);
                p3avg.z /= double(num_dets);
                i = celedeproj01(p3avg, &img_log[imct].RA, &img_log[imct].Dec);
                if (i == 0)
                    ;  // All is well.
                else if (i == 1) {
                    cout << "Warning: vector of zeros fed to celedeproj01\n";
                    img_log[imct].RA = img_log[imct].Dec = 0.0;
                } else if (i == 2) {
                    cout << "Warning: impossible z value " << p3avg.z << " fed to celedeproj01\n";
                    img_log[imct].RA = img_log[imct].Dec = 0.0;
                } else {
                    cout << "Warning: unspecified failure from celedeproj01 with\n";
                    cout << "input " << p3avg.x << " " << p3avg.y << " " << p3avg.z << "\n";
                    img_log[imct].RA = img_log[imct].Dec = 0.0;
                }
            }
            if (DEBUG >= 1) {
                cout << "Image " << imct << " of " << img_log.size() << ": " << num_dets << " = "
                     << img_log[imct].endind - img_log[imct].startind;
                cout << " from " << img_log[imct].startind << " to " << img_log[imct].endind << " of "
                     << detvec.size() << ".\n";
                fflush(stdout);
            }
            imct++;
        }
    }

    EarthMJD = {};
    Earthpos = {};
    Earthvel = {};
    read_horizons_fileLD(earthfile, EarthMJD, Earthpos, Earthvel);
    cout << "Finished reading heliocentric ephemeris file " << earthfile << " for the Earth.\n";

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
    // Calculate observer's heliocentric position at the time of each image.
    observer_heliopos = {};
    for (imct = 0; imct < img_log.size(); imct++) {
        if (imct == 0 ||
            (imct > 0 && stringnmatch01(img_log[imct].obscode, img_log[imct - 1].obscode, 3) == 0)) {
            // Observatory has changed: get observatory coordinates for this image.
            status = obscode_lookup(observatory_list, img_log[imct].obscode, obslon, plxcos, plxsin);
            if (status > 0) {
                cerr << "ERROR: obscode_lookup failed for observatory code " << img_log[imct].obscode << "\n";
                return (3);
            }
        }
        observer_barycoords01LD(img_log[imct].MJD, 5, obslon, plxcos, plxsin, EarthMJD, Earthpos, outpos);
        observer_heliopos.push_back(outpos);
    }
    // assign observer heliocentric position to each detection.
    for (imct = 0; imct < img_log.size(); imct++) {
        for (i = img_log[imct].startind; i < img_log[imct].endind; i++) {
            detvec[i].x = observer_heliopos[imct].x;
            detvec[i].y = observer_heliopos[imct].y;
            detvec[i].z = observer_heliopos[imct].z;
            // Check that the obscodes match between image and detection
            if (stringnmatch01(detvec[i].obscode, img_log[imct].obscode, 3) != 0) {
                cout << "ERROR: obscode mismatch (" << detvec[i].obscode << " vs. " << img_log[imct].obscode
                     << " between image " << imct << " and detection " << i << "\n";
                return (4);
            }
            detvec[i].indvec = {};  // Just making sure the index vectors are empty at this point.
        }
    }

    makeTracklets(detvec, img_log, pairvec, pairdets, outpairfile, pairdetfile);
    return 0;
}