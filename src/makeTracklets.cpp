#include "makeTracklets.h"

void buildTracklets(
    MakeTrackletsConfig config,
    std::vector<det_obsmag_indvec> detvec,
    std::vector<img_log03> img_log,
    std::vector<longpair> &pairvec,
    std::vector<det_obsmag_indvec> &pairdets
) {
    longpair onepair = longpair(0, 0);
    int i = 0;
    int j = 0;
    int k = 0;
    int imct = 0;
    long detct = 0;
    int startind = 0;
    int endind = 0;
    long double MJD, RA, Dec;
    MJD = RA = Dec = 0.0L;
    double mag = 0l;
    double maxvel = config.maxvel;  // Max angular velocity in deg/day
    double minvel = 0.0l;    // Min angular velocity in deg/day
    double minarc = 0.0l;    // Min total angular arc in arcseconds
    double angvel = 0.0l;
    double maxtime = config.maxtime;           // Max time interval a tracklet could span,
                                        // in days.
    double maxdist = maxvel * maxtime;  // Max angular distance a tracklet
                                        // could span, in degrees.
    double imrad = config.imagerad;            // radius from image center to most distant corner (deg).
    double obslon = 289.26345L;
    double plxcos = 0.865020L;
    double plxsin = -0.500901L;
    xy_index xyind = xy_index(0.0, 0.0, 0);
    vector<xy_index> axyvec = {};
    double dist, pa;
    dist = pa = 0.0;
    int dettarg = 0;
    double mintime = config.imagetimetol / SOLARDAY;
    double maxgcr = config.maxgcr;
    int mintrkpts = config.mintrkpts;

    // PERFORM PAIRING
    long pdct = 0;    // count of detections that have been paired
    long pairct = 0;  // count of actual pairs
    // Loop over images for image A
    for (imct = 0; imct < img_log.size(); imct++) {
        int apct = 0;
        int adetct = 0;
        // See if there are any images that might match
        std::vector<int> imagematches = {};
        int imatchcount = 0;
        int imtarg = imct + 1;
        while (imtarg < img_log.size() && img_log[imtarg].MJD < img_log[imct].MJD + maxtime) {
            // See if the images are close enough on the sky.
            double timediff = img_log[imtarg].MJD - img_log[imct].MJD;
            if (!isnormal(timediff) || timediff < 0.0) {
                cerr << "WARNING: Negative time difference " << timediff << " encountered between images "
                     << imct << " and " << imtarg << "\n";
            }
            double imcendist = distradec01(img_log[imct].RA, img_log[imct].Dec, img_log[imtarg].RA, img_log[imtarg].Dec);
            if (imcendist < 2.0 * imrad + maxvel * timediff && timediff >= mintime) {
                if (DEBUG >= 1)
                    cout << "  pairs may exist between images " << imct << " and " << imtarg
                         << ": dist = " << imcendist << ", timediff = " << timediff << "\n";
                imagematches.push_back(imtarg);
            }
            imtarg++;
        }
        if (DEBUG >= 1)
            cout << "Looking for pairs for image " << imct << ": " << imagematches.size()
                 << " later images are worth searching\n";
        if (imagematches.size() > 0) {
            // Search is worth doing. Project all the detections
            // on image A.
            xyind = xy_index(0.0, 0.0, 0);
            axyvec = {};
            dist = pa = 0.0;
            dettarg = 0;
            for (detct = img_log[imct].startind; detct < img_log[imct].endind; detct++) {
                distradec02(img_log[imct].RA, img_log[imct].Dec, detvec[detct].RA, detvec[detct].Dec, &dist,
                            &pa);
                xyind = xy_index(dist * sin(pa / DEGPRAD), dist * cos(pa / DEGPRAD), detct);
                axyvec.push_back(xyind);
                if ((!isnormal(xyind.x) && xyind.x != 0) || (!isnormal(xyind.y) && xyind.y != 0)) {
                    cerr << "nan-producing input: ra1, dec1, ra2, dec2, dist, pa:\n";
                    cerr << img_log[imct].RA << " " << img_log[imct].Dec << " " << detvec[detct].RA << " "
                         << detvec[detct].Dec << " " << dist << " " << pa << " " << xyind.x << " " << xyind.x
                         << "\n";
                    // return(6);
                }
            }
            // Loop over images with potential matches (image B's)
            for (imatchcount = 0; imatchcount < imagematches.size(); imatchcount++) {
                imtarg = imagematches[imatchcount];
                double range = (img_log[imtarg].MJD - img_log[imct].MJD) * maxvel;
                vector<xy_index> bxyvec = {};
                // Project all detections on image B
                for (dettarg = img_log[imtarg].startind; dettarg < img_log[imtarg].endind; dettarg++) {
                    distradec02(img_log[imct].RA, img_log[imct].Dec, detvec[dettarg].RA, detvec[dettarg].Dec,
                                &dist, &pa);
                    xyind = xy_index(dist * sin(pa / DEGPRAD), dist * cos(pa / DEGPRAD), dettarg);
                    bxyvec.push_back(xyind);
                }
                // Create k-d tree of detections on image B (imtarg).
                int dim = 1;
                xy_index xyi = bxyvec[0];
                kdpoint root = kdpoint(xyi, -1, -1, dim);
                kdpoint lp1 = kdpoint(xyi, -1, -1, dim);
                kdpoint rp1 = kdpoint(xyi, -1, -1, dim);
                kdpoint kdtest = kdpoint(xyi, -1, -1, dim);
                vector<kdpoint> kdvec = {};
                long medpt;
                medpt = medindex(bxyvec, dim);
                root = kdpoint(bxyvec[medpt], -1, -1, 1);
                kdvec.push_back(root);
                kdtest = kdvec[0];
                kdtree01(bxyvec, dim, medpt, 0, kdvec);
                // Loop over detections on image A
                if (DEBUG >= 1)
                    cout << "Looking for pairs between " << axyvec.size() << " detections on image " << imct
                         << " and " << kdvec.size() << " on image " << imtarg << "\n";
                for (detct = 0; detct < axyvec.size(); detct++) {
                    vector<long> indexvec = {};
                    if ((isnormal(axyvec[detct].x) || axyvec[detct].x == 0) &&
                        (isnormal(axyvec[detct].y) || axyvec[detct].y == 0)) {
                        kdrange01(kdvec, axyvec[detct].x, axyvec[detct].y, range, indexvec);
                    } else {
                        cerr << "WARNING: detection " << detct << " on image " << imct
                             << " not normal: " << axyvec[detct].x << " " << axyvec[detct].y << "\n";
                    }
                    int matchnum = indexvec.size();
                    long matchpt = 0;
                    int matchct = 0;
                    if (matchnum > 0) {
                        // Record image A detection as paired, if not already recorded.
                        if (detvec[axyvec[detct].index].index < 0) {
                            // This detection has not yet been paired with any other.
                            detvec[axyvec[detct].index].index *=
                                    -1;  // Mark as paired by changing to positive sign.
                            pairdets.push_back(
                                    detvec[axyvec[detct].index]);  // Load into paired detection vector
                            pairdets[pdct].indvec = {};  // Make sure index vector is currently empty.
                            detvec[axyvec[detct].index].index =
                                    pdct;  // Re-assign index to apply to paired detection vector
                            pdct++;        // Increment count of paired detections
                            adetct++;
                            if (pdct != pairdets.size())
                                cerr << "\nWARNING: PAIRED DETECTION MISMATCH: " << pdct << " vs "
                                     << pairdets.size() << "\n\n";
                        }
                        // Record image B detections
                        for (matchct = 0; matchct < matchnum; matchct++) {
                            matchpt = indexvec[matchct];
                            if (detvec[kdvec[matchpt].point.index].index < 0) {
                                // This detection has not yet been paired with any other.
                                detvec[kdvec[matchpt].point.index].index *=
                                        -1;  // Mark as paired by changing to positive sign
                                pairdets.push_back(detvec[kdvec[matchpt].point.index]);  // Load into paired
                                                                                         // detection vector
                                pairdets[pdct].indvec = {};  // Make sure index vector is currently empty.
                                detvec[kdvec[matchpt].point.index].index =
                                        pdct;  // Re-assign index to apply to paired detection vector
                                pdct++;        // Increment count of paired detections
                                if (pdct != pairdets.size())
                                    cerr << "\nWARNING: PAIRED DETECTION MISMATCH: " << pdct << " vs "
                                         << pairdets.size() << "\n\n";
                            }
                            // Write index values for both components of the
                            // new pair to the pair vector, regardless of whether
                            // the index values are pre-existing or newly assigned.
                            onepair = longpair(detvec[axyvec[detct].index].index,
                                               detvec[kdvec[matchpt].point.index].index);
                            if (DEBUG >= 1)
                                cout << "Writing pair " << detvec[axyvec[detct].index].index << ", "
                                     << detvec[kdvec[matchpt].point.index].index << "\n";
                            pairvec.push_back(onepair);
                            pairct++;
                            apct++;
                            // Load index of each detection into the paired index vector of the other
                            pairdets[detvec[axyvec[detct].index].index].indvec.push_back(
                                    detvec[kdvec[matchpt].point.index].index);
                            pairdets[detvec[kdvec[matchpt].point.index].index].indvec.push_back(
                                    detvec[axyvec[detct].index].index);
                        }
                        // Close if-statement checking if image A detection was matched to anything.
                    }
                    // Close loop over detections on source image (image A)
                }
                // Close loop over image B candidates
            }
            // Close if-statement checking if any images could match image A
        }
        if (DEBUG >= 1)
            cout << "Image " << imct << ": found " << adetct << " newly paired detections and a total of "
                 << apct << " pairs.\n";
        // Close loop over images for image A
    }
    if (DEBUG >= 1) cout << "Test count of paired detections: " << pdct << " " << pairdets.size() << "\n";
    if (DEBUG >= 1) cout << "Test count of pairs: " << pairct << " " << pairvec.size() << "\n";

    return;
}

void refineTracklets(
    MakeTrackletsConfig config,
    std::vector<det_obsmag_indvec> &pairdets,
    string outpairfile
) {
    ofstream outstream1;
    vector<long_index> pair_partner_num = {};
    double pa, slopex, slopey, interceptx, intercepty, dist;
    double outra1, outra2, outdec1, outdec2;
    int worstpoint = -1;
    long pdct = 0;
    long double maxgcr = config.maxgcr;
    long double minarc = config.minarc;
    long double maxvel = config.maxvel;
    long double minvel = config.minvel;
    int mintrkpts = config.mintrkpts;

    // Load a vector storing the number of pair-partners found for each detection.
    for (size_t i = 0; i < pairdets.size(); i++) {
        long_index ppn = long_index(pairdets[i].indvec.size(), i);
        pair_partner_num.push_back(ppn);
    }
    // Sort the new vector by number of pair-partners
    sort(pair_partner_num.begin(), pair_partner_num.end(), lower_long_index());

    // Analyze paired detections in order of decreasing number of partners.
    // At the same time, write the output pair file, distinguishing
    // between real pairs comprising just two detections (indicated
    // with the letter P), and effective pairs that are really tracklets
    // created by fitting and averaging many points that lie along a
    // consistent trajectory.
    cout << "Constructing tracklets, and writing pairs to output file\n";
    outstream1.open(outpairfile);

    for (int i = pairdets.size() - 1; i >= 0; i--) {
        pdct = pair_partner_num[i].index;
        int istracklet = 0;  // Assume there is no tracklet unless one is confirmed to exist.
        if (pairdets[pdct].indvec.size() > mintrkpts - 1) {
            if (DEBUG >= 2) {
                cout << "Working on detection " << i << " = " << pdct << " with " << pair_partner_num[i].lelem
                     << " = " << pairdets[pdct].indvec.size() << " pair partners:\n";
                for (size_t j = 0; j < pairdets[pdct].indvec.size(); j++) {
                    cout << pairdets[pdct].indvec[j] << ", ";
                }
                cout << "\n";
            }
            // The corresponding detection is paired with more than one
            // other detection.
            // Project all of these pairs relative to detection pdct,
            // storing x,y projected coordinates in axyvec.
            std::vector<xy_index> axyvec = {};
            std::vector<det_obsmag_indvec> ppset = {};
            for (size_t j = 0; j < pairdets[pdct].indvec.size(); j++) {
                long detct = pairdets[pdct].indvec[j];
                if (pairdets[detct].indvec.size() > 0) {
                    // Detection detct hasn't already been allocated to a tracklet.
                    distradec02(pairdets[pdct].RA, pairdets[pdct].Dec, pairdets[detct].RA,
                                pairdets[detct].Dec, &dist, &pa);
                    dist *= 3600.0L;  // Convert distance from degrees to arcsec.
                    xy_index xyind = xy_index(dist * sin(pa / DEGPRAD), dist * cos(pa / DEGPRAD), detct);
                    axyvec.push_back(xyind);
                    ppset.push_back(pairdets[detct]);  // We need this vector (of type det_obsmag_indvec)
                                                       // mainly just to have some way to store the
                                                       // indices of mutually consistent pair partners
                                                       // on the next step
                }
            }
            if (DEBUG >= 2)
                cout << "Loaded axyvec and ppset vectors OK, with sizes " << axyvec.size() << " and "
                     << ppset.size() << "\n";
            if (axyvec.size() != ppset.size()) {
                cerr << "ERROR: vectors of projected and original\n";
                cerr << "pair partner candidates do not have the same length!\n";
                cerr << axyvec.size() << " != " << ppset.size() << "\n";
                // return (3);
            }
            // Perform n^2 search to find the largest cluster of consistent points.
            // Load all the pair-partners for detection pdct into a vector
            // (this is mainly just to be able to use the index vectors)
            for (size_t j = 0; j < axyvec.size(); j++) {
                double dtref = ppset[j].MJD - pairdets[pdct].MJD;
                if (dtref == 0) {
                    cerr << "ERROR: paired detections with no time separation!\n";
                    // return (4);
                }
                // Make sure corresponding index vector in ppset is empty
                ppset[j].indvec = {};
                // Count consistent pair partners
                if (DEBUG >= 2) cout << "Counting consistent pair partners\n";
                for (size_t k = 0; k < axyvec.size(); k++) {
                    if (j != k) {
                        double dt = ppset[k].MJD - pairdets[pdct].MJD;
                        double dx = axyvec[k].x - axyvec[j].x * (dt / dtref);
                        double dy = axyvec[k].y - axyvec[j].y * (dt / dtref);
                        dist = sqrt(dx * dx + dy * dy);
                        if (DEBUG >= 2)
                            cout << "Detection " << axyvec[j].index << ":" << axyvec[k].index
                                 << " dist = " << dist << "\n";
                        if (dist < 2.0 * maxgcr) {
                            ppset[j].indvec.push_back(k);
                        }
                    }
                }
            }
            // Find the largest set of pair-partners lying along a line.
            if (DEBUG >= 2) cout << "Find the largest set of pair-partners lying along a line.\n";
            int biggest_tracklet = -1;
            int tracklet_size = 0;
            if (DEBUG >= 2) cout << "size = " << ppset.size() << "\n";
            for (size_t j = 0; j < ppset.size(); j++) {
                if (DEBUG >= 2)
                    cout << j << ":" << ppset.size() - 1 << " size = " << ppset[j].indvec.size() << " ";
                if (ppset[j].indvec.size() + 2 > tracklet_size) {
                    tracklet_size = ppset[j].indvec.size() +
                                    2;  // We add one for pdct, one for j, to get actual tracklet size
                    biggest_tracklet = j;
                    if (DEBUG >= 2)
                        cout << "bt = " << biggest_tracklet << ", size = " << tracklet_size << "\n";
                } else if (DEBUG >= 2)
                    cout << "not the biggest\n";
            }
            if (DEBUG >= 2)
                cout << "Biggest tracklet is " << biggest_tracklet << ", which corresponds to "
                     << axyvec[biggest_tracklet].index << ", with size " << tracklet_size << "\n";
            istracklet = 0;  // Assume there is no tracklet until one is confirmed to exist.
            if (tracklet_size <= mintrkpts) {
                istracklet = 0;
            } else {
                // Perform linear fits to x and y vs time.
                // Load all the points from the biggest potential tracklet.
                std::vector<point3d_index> track_mrdi_vec = {};  // We need this vector purely so we can do a time-sort.
                                      // mrdi stands for MJD, RA, Dec, index
                // Load the reference point
                point3d_index p3di = point3d_index(0.0l, 0.0l, 0.0l, pdct);
                track_mrdi_vec.push_back(p3di);
                // Load anchor point corresponding to biggest_tracklet
                p3di = point3d_index(ppset[biggest_tracklet].MJD - pairdets[pdct].MJD,
                                     axyvec[biggest_tracklet].x, axyvec[biggest_tracklet].y,
                                     axyvec[biggest_tracklet].index);
                track_mrdi_vec.push_back(p3di);
                // Load the other points
                for (size_t j = 0; j < ppset[biggest_tracklet].indvec.size(); j++) {
                    p3di = point3d_index(ppset[ppset[biggest_tracklet].indvec[j]].MJD - pairdets[pdct].MJD,
                                         axyvec[ppset[biggest_tracklet].indvec[j]].x,
                                         axyvec[ppset[biggest_tracklet].indvec[j]].y,
                                         axyvec[ppset[biggest_tracklet].indvec[j]].index);
                    track_mrdi_vec.push_back(p3di);
                    // timevec.push_back(ppset[ppset[biggest_tracklet].indvec[j]].MJD - pairdets[pdct].MJD);
                    // xvec.push_back(axyvec[ppset[biggest_tracklet].indvec[j]].x);
                    // yvec.push_back(axyvec[ppset[biggest_tracklet].indvec[j]].y);
                    // detindexvec.push_back(axyvec[ppset[biggest_tracklet].indvec[j]].index);
                }
                // Sort track_mrdi_vec by time.
                sort(track_mrdi_vec.begin(), track_mrdi_vec.end(), lower_point3d_index_x());
                // Load time, x, y, and index vectors from sorted track_mrdi_vec.
                std::vector<double> timevec = {};
                std::vector<double> xvec = {};
                std::vector<double> yvec = {};
                std::vector<long> detindexvec = {};
                for (size_t j = 0; j < track_mrdi_vec.size(); j++) {
                    timevec.push_back(track_mrdi_vec[j].x);
                    xvec.push_back(track_mrdi_vec[j].y);
                    yvec.push_back(track_mrdi_vec[j].z);
                    detindexvec.push_back(track_mrdi_vec[j].index);
                }
                if (DEBUG >= 2) {
                    cout << "First iteration linear fit vectors:\n";
                    for (size_t j = 0; j < timevec.size(); j++) {
                        cout << detindexvec[j] << " " << timevec[j] << " " << xvec[j] << " " << yvec[j]
                             << "\n";
                    }
                }

                // Perform fit to projected x coordinate as a function of time
                linfituw01(timevec, xvec, slopex, interceptx);
                // Perform fit to projected y coordinate as a function of time
                linfituw01(timevec, yvec, slopey, intercepty);
                // Load vector of residuals
                std::vector<double> fiterr = {};
                for (size_t j = 0; j < timevec.size(); j++) {
                    fiterr.push_back(sqrt(DSQUARE(timevec[j] * slopex + interceptx - xvec[j]) +
                                          DSQUARE(timevec[j] * slopey + intercepty - yvec[j])));
                }
                // Ditch duplicate times, if there are any
                int istimedup = 1;  // Guilty until proven innocent
                while (istimedup == 1 && timevec.size() >= mintrkpts + 1) {
                    istimedup = 0;
                    int j = 1;
                    while (j < timevec.size() && istimedup == 0) {
                        if (fabs(timevec[j] - timevec[j - 1]) < config.imagetimetol / SOLARDAY) {
                            istimedup = 1;  // Point j and j-1 are time-duplicates.
                            // Mark for rejection whichever one has the largest fitting error
                            if (fiterr[j] >= fiterr[j - 1])
                                worstpoint = j;
                            else
                                worstpoint = j - 1;
                        }
                        j++;
                    }
                    if (istimedup == 1) {
                        // Reject the bad point
                        int trkptnum = timevec.size();
                        for (int j = worstpoint; j < trkptnum - 1; j++) {
                            timevec[j] = timevec[j + 1];
                            xvec[j] = xvec[j + 1];
                            yvec[j] = yvec[j + 1];
                            detindexvec[j] = detindexvec[j + 1];
                        }
                        trkptnum--;
                        timevec.resize(trkptnum);
                        xvec.resize(trkptnum);
                        yvec.resize(trkptnum);
                        detindexvec.resize(trkptnum);
                        // Re-do linear fit
                        // Perform fit to projected x coordinate as a function of time
                        linfituw01(timevec, xvec, slopex, interceptx);
                        // Perform fit to projected y coordinate as a function of time
                        linfituw01(timevec, yvec, slopey, intercepty);
                        // Load vector of residuals
                        fiterr = {};
                        for (j = 0; j < timevec.size(); j++) {
                            fiterr.push_back(sqrt(DSQUARE(timevec[j] * slopex + interceptx - xvec[j]) +
                                                  DSQUARE(timevec[j] * slopey + intercepty - yvec[j])));
                        }
                    }
                }
                // Find worst error.
                double worsterr = 0.0l;
                for (size_t j = 0; j < timevec.size(); j++) {
                    if (fiterr[j] > worsterr) {
                        worsterr = fiterr[j];
                        worstpoint = j;
                    }
                }
                // Reject successive points until either there are only three left
                // or the worst error drops below maxgcr.
                while (worsterr > maxgcr && timevec.size() > 3 && timevec.size() >= mintrkpts) {
                    // Reject the worst point
                    int trkptnum = timevec.size();
                    for (int j = worstpoint; j < trkptnum - 1; j++) {
                        timevec[j] = timevec[j + 1];
                        xvec[j] = xvec[j + 1];
                        yvec[j] = yvec[j + 1];
                        detindexvec[j] = detindexvec[j + 1];
                    }
                    trkptnum--;
                    timevec.resize(trkptnum);
                    xvec.resize(trkptnum);
                    yvec.resize(trkptnum);
                    detindexvec.resize(trkptnum);
                    // Perform fit to projected x coordinate as a function of time
                    linfituw01(timevec, xvec, slopex, interceptx);
                    // Perform fit to projected y coordinate as a function of time
                    linfituw01(timevec, yvec, slopey, intercepty);
                    // Load vector of residuals
                    fiterr = {};
                    for (size_t j = 0; j < timevec.size(); j++) {
                        fiterr.push_back(sqrt(DSQUARE(timevec[j] * slopex + interceptx - xvec[j]) +
                                              DSQUARE(timevec[j] * slopey + intercepty - yvec[j])));
                    }
                    // Find worst error.
                    worsterr = 0.0l;
                    for (size_t j = 0; j < timevec.size(); j++) {
                        if (fiterr[j] > worsterr) {
                            worsterr = fiterr[j];
                            worstpoint = j;
                        }
                    }
                }
                if (worsterr <= maxgcr && timevec.size() >= 3 && timevec.size() >= mintrkpts) {
                    // We succeeded in finding a tracklet with no time-duplicates, and
                    // no outliers beyond maxgcr. Prepare to write it to the pair file.
                    // Select points that will represent this tracklet.
                    int instep = (timevec.size() - 1) / 4;
                    int rp1 = instep;
                    int rp2 = timevec.size() - 1 - instep;
                    if (rp1 == rp2) {
                        cerr << "ERROR: both representative points for a tracklet are the same!\n";
                        cerr << "size, instep, rp1, rp2: " << timevec.size() << " " << instep << " " << rp1
                             << " " << rp2 << "\n";
                        // return (4);
                    }
                    // Calculate angular velocity in deg/day. The slope values
                    // correspond to velocities in arcsec/day.
                    double angvel = sqrt(slopex * slopex + slopey * slopey) / 3600.0l;

                    // Determine improved RA, Dec based on tracklet fit for the representative points
                    // Calculated projected x, y at rp1
                    double dx = timevec[rp1] * slopex + interceptx;
                    double dy = timevec[rp1] * slopey + intercepty;
                    // Calculate equivalent celestial position angle.
                    if (dx == 0l && dy >= 0l)
                        pa = 0.0l;
                    else if (dx == 0l && dy < 0l)
                        pa = M_PI;
                    else if (dx > 0l)
                        pa = M_PI / 2.0l - atan(dy / dx);
                    else if (dx < 0l)
                        pa = 3.0l * M_PI / 2.0l - atan(dy / dx);
                    else {
                        cerr << "ERROR: logical impossibility while trying to solve for PA\n";
                        cerr << "dx = " << dx << " dy = " << dy << "\n";
                    }
                    dist = sqrt(dx * dx + dy * dy) / 3600.0l;  // renders distance in degrees, not arcsec.
                    pa *= DEGPRAD;                             // position angle in degrees, not radians.
                    arc2cel01(pairdets[pdct].RA, pairdets[pdct].Dec, dist, pa, outra1, outdec1);
                    if (!isnormal(outra1)) {
                        cerr << "NAN WARNING: dx, dy, dist, pa: " << dx << " " << dy << " " << dist << " "
                             << pa << "\n";
                    }
                    // Calculated projected x, y at rp2
                    dx = timevec[rp2] * slopex + interceptx;
                    dy = timevec[rp2] * slopey + intercepty;
                    // Calculate equivalent celestial position angle.
                    if (dx == 0l && dy >= 0l)
                        pa = 0.0l;
                    else if (dx == 0l && dy < 0l)
                        pa = M_PI;
                    else if (dx > 0l)
                        pa = M_PI / 2.0l - atan(dy / dx);
                    else if (dx < 0l)
                        pa = 3.0l * M_PI / 2.0l - atan(dy / dx);
                    else {
                        cerr << "ERROR: logical impossibility while trying to solve for PA\n";
                        cerr << "dx = " << dx << " dy = " << dy << "\n";
                    }
                    dist = sqrt(dx * dx + dy * dy) / 3600.0l;  // renders distance in degrees, not arcsec.
                    pa *= DEGPRAD;                             // position angle in degrees, not radians.
                    arc2cel01(pairdets[pdct].RA, pairdets[pdct].Dec, dist, pa, outra2, outdec2);
                    // Calculate total angular arc
                    distradec02(outra1, outdec1, outra2, outdec2, &dist, &pa);
                    dist *= 3600.0l;
                    if (dist >= minarc && angvel >= minvel && angvel <= maxvel) {
                        // Write out representative pair, followed by RA, Dec and the total number of
                        // constituent points representative pair
                        outstream1 << "T " << detindexvec[rp1] << " " << detindexvec[rp2] << " ";
                        // RA1, Dec1, RA2, Dec2
                        outstream1 << fixed << setprecision(6) << outra1 << " " << outdec1 << " " << outra2
                                   << " " << outdec2 << " ";
                        // Number of points in final, refined tracklet.
                        outstream1 << detindexvec.size() << "\n";
                        // Now write out the detection indices for this full number of points,
                        // and wipe all the associated index vectors.
                        for (size_t j = 0; j < detindexvec.size(); j++) {
                            outstream1 << detindexvec[j] << "\n";
                            pairdets[detindexvec[j]].indvec = {};
                        }
                        istracklet = 1;
                        // Close if-statement confirming that a bona fide,
                        // aligned tracklet was found and written to the output file.
                    } else {
                        istracklet = 0;
                        cout << "A tracklet was rejected: arc = " << dist << " < " << minarc
                             << " or angvel = " << angvel << " < " << minvel << "\n";
                    }
                } else
                    istracklet = 0;
                // Close else-statement confirming there was a candidate for
                // being an aligned tracklet.
            }
            // Close if-statement checking that detection i has more than
            // one pair-partner, and hence COULD be part of a tracklet
        } else
            istracklet = 0;
        if (istracklet == 0 && mintrkpts == 2) {
            // Write out all the pairs as normal
            for (size_t j = 0; j < pairdets[pdct].indvec.size(); j++) {
                int k = pairdets[pdct].indvec[j];
                // Calculate angular arc and angular velocity
                distradec02(pairdets[pdct].RA, pairdets[pdct].Dec, pairdets[k].RA, pairdets[k].Dec, &dist,
                            &pa);
                double angvel = dist / fabs(pairdets[pdct].MJD - pairdets[k].MJD);  // Degrees per day
                dist *= 3600.0l;                                             // Arcseconds
                if (pairdets[k].indvec.size() > 0 && k > pdct && angvel >= minvel && dist >= minarc &&
                    angvel <= maxvel) {
                    outstream1 << "P " << pdct << " " << k << "\n";
                } else if (angvel < minvel || dist < minarc) {
                    cout << "A pair was rejected: arc = " << dist << " < " << minarc
                         << " or angvel = " << angvel << " < " << minvel << "\n";
                }
            }
        }
        // Close loop over all detections
    }
    outstream1.close();
}