#include "makeTracklets.h"

// like library function strncpy, but works the way I want it to.
// Copies first n characters of a string into a character array.
void stringncopy01(char *dest, string const& source, int n) {
    int i=0;
    int nchar = source.size();
    if(nchar<n) {
        // We have enough space to copy the whole thing.
        for(i=0;i<nchar;i++) {
            dest[i] = source[i];
        }
        dest[nchar]='\0';
    } else {
        // Not enough space: copy first n-1 characters.
        for(i=0;i<n;i++) {
            dest[i] = source[i];
        }
        dest[n-1]='\0';
    }
}

// like library function strncmp, but works the way I want it to.
// Compares first n characters of two character arrays
int stringnmatch01(char const* string1, char const* string2, int n) {
    int i=0;
    while(i<n && string1[i]!='\0' && string2[i]!='\0') {
        if(string1[i]<string2[i]) return(-1);
        else if(string1[i]>string2[i]) return(1);
        i++;
    }
    // If we get here without returning, the strings must have been equal
    return(0);
}

// Given an input state-vector ephemeris file downloaded directly
// from JPL Horizons, read it into position and velocity vectors.
// Note that the default unit convention is km for positions and
// km/sec for velocities. Note also that JPL state-vector
// ephemerides use dynamical TT, which is ahead of UT1 by about
// 70 seconds in 2022. This program does NOT correct TT to UT1,
// but programs making use of the ouput mjd, position, and velocity
// vectors might need to.
int read_horizons_fileLD(string infile, std::vector<long double> &mjdvec, std::vector<Point3LD> &pos, std::vector<Point3LD> &vel) {
    ifstream instream1 {infile};
    Point3LD pospoint = Point3LD(0.0,0.0,0.0);
    Point3LD velpoint = Point3LD(0.0,0.0,0.0);
    int reachedeof=0;
    int ondata=0;
    int i=0;
    char c = '0';
    int reachedend=0;
    string teststring, lnfromfile;
    long double x,y,z,vx,vy,vz,MJD;
    MJD=x=y=z=vz=vy=vz=0.0;

    if(!instream1) {
        cerr << "ERROR: can't open input file " << infile << "\n";
        return(1);
    }
    while(reachedeof==0 && !reachedend) {
        while(!ondata && !reachedend) {
            // See if this line contains the code for start-of-data.
            lnfromfile = "";
            teststring = "";
            getline(instream1,lnfromfile);
            if(instream1.eof()) reachedeof=1; //End of file, fine.
            else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
            else if(instream1.bad()) reachedeof=-2; //Worse problem, warn

            if(lnfromfile.size()>=5) {
                for(i=0;i<5;i++) {
                    teststring.push_back(lnfromfile[i]);
                }
                if(teststring == "$$SOE") ondata=1;
                else if(teststring == "$$EOE") reachedend=1;
            }
        }
        while(ondata && !reachedend && reachedeof==0) {
            lnfromfile = "";
            teststring = "";
            getline(instream1,lnfromfile);
            if(instream1.eof()) reachedeof=1; //End of file, fine.
            else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
            else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
            if(lnfromfile.size()>=5) {
                for(i=0;i<5;i++) {
                    teststring.push_back(lnfromfile[i]);
                }
                if(teststring == "$$EOE") reachedend=1;
            }
            if(!reachedend && reachedeof==0) {
                //Attempt to read entire four-line block.
                //First line has MJD
                teststring = "";
                c='0';
                i=0;
                while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!=' ' && c!='\n' && c!=EOF) {
                    c=lnfromfile[i];
                    if(c!=' ' && c!='=' && c!='\n' && c!=EOF) teststring.push_back(c);
                    if(c==EOF) reachedeof=1;
                    i++;
                }
                MJD=stold(teststring);
                //Next line has x,y,z positions
                lnfromfile = "";
                teststring = "";
                getline(instream1,lnfromfile);
                if(instream1.eof()) reachedeof=1; //End of file, fine.
                else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
                else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
                //Read to first equals sign
                c='0';
                i=0;
                while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
                    c=lnfromfile[i];
                    if(c==EOF) reachedeof=1;
                    i++;
                }
                //Read to next equals sign, loading into teststring to get X
                if(i<lnfromfile.size()) c=lnfromfile[i];
                while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
                    c=lnfromfile[i];
                    if(c==EOF) reachedeof=1;
                    if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
                    i++;
                }
                x = stold(teststring);
                teststring = "";
                //Read to next equals sign, loading into teststring to get Y
                if(i<lnfromfile.size()) c=lnfromfile[i];
                while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
                    c=lnfromfile[i];
                    if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
                    i++;
                }
                y = stold(teststring);
                teststring = "";
                //Read to next equals sign, loading into teststring to get Z
                if(i<lnfromfile.size()) c=lnfromfile[i];
                while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
                    c=lnfromfile[i];
                    if(c==EOF) reachedeof=1;
                    if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
                    i++;
                }
                z = stold(teststring);
                //Next line has x,y,z velocities
                lnfromfile = "";
                teststring = "";
                getline(instream1,lnfromfile);
                if(instream1.eof()) reachedeof=1; //End of file, fine.
                else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
                else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
                //Read to first equals sign
                c='0';
                i=0;
                while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
                    c=lnfromfile[i];
                    if(c==EOF) reachedeof=1;
                    i++;
                }
                //Read to next equals sign, loading into teststring to get XV
                if(i<lnfromfile.size()) c=lnfromfile[i];
                while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
                    c=lnfromfile[i];
                    if(c==EOF) reachedeof=1;
                    if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
                    i++;
                }
                vx = stold(teststring);
                teststring = "";
                //Read to next equals sign, loading into teststring to get VY
                if(i<lnfromfile.size()) c=lnfromfile[i];
                while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
                    c=lnfromfile[i];
                    if(c==EOF) reachedeof=1;
                    if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
                    i++;
                }
                vy = stold(teststring);
                teststring = "";
                //Read to next equals sign, loading into teststring to get VZ
                if(i<lnfromfile.size()) c=lnfromfile[i];
                while(i<lnfromfile.size() && reachedeof == 0 && c!='=' && c!='\n' && c!=EOF) {
                    c=lnfromfile[i];
                    if(c==EOF) reachedeof=1;
                    if(c!='=' && c!=' ' && c!='\n' && c!=EOF) teststring.push_back(c);
                    i++;
                }
                vz = stold(teststring);
                // Load output vectors
                pospoint = Point3LD(x,y,z);
                velpoint = Point3LD(vx,vy,vz);
                pos.push_back(pospoint);
                vel.push_back(velpoint);
                mjdvec.push_back(MJD-MJDOFF);
                // Next line is of no current interest: read and discard
                lnfromfile = "";
                teststring = "";
                getline(instream1,lnfromfile);
                if(instream1.eof()) reachedeof=1; //End of file, fine.
                else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
                else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
            }
        }
    }
    if(reachedeof==1 && ondata==1) {
        //Read file successfully to the end.
        return(0);
    } else if(reachedeof==1) {
        //Did not find any data
        return(1);
    } else return(reachedeof);
}

// Look up an observatory code from a list, and copy the
// coordinates to obslon, plxcos, and plxsin.
int obscode_lookup(std::vector<Observatory> const& observatory_list, char const* obscode, double &obslon, double &plxcos, double &plxsin) {
    int i=0;
    int nobs = observatory_list.size();
    for(i=0; i<nobs; i++) {
        // cout << "obscode_lookup comparing " << obscode << " with " << observatory_list[i].obscode << ":";
        if(stringnmatch01(observatory_list[i].obscode,obscode,3)==0) {
        // cout << " it\'s a match!\n";
        obslon = observatory_list[i].obslon;
        plxcos = observatory_list[i].plxcos;
        plxsin = observatory_list[i].plxsin;
        return(0);
        } // else cout << " not a match\n";
    }
    cerr << "ERROR: observatory " << obscode << " not found in list\n";
    return(1);
}

void make_LDvec(int nx, std::vector<long double> &ldvec) {
    int i = 0;
    ldvec = {};
    for (i = 0; i <= nx; i++) ldvec.push_back(0.0);
}

void make_LDmat(int nx, int ny, std::vector<std::vector<long double>> &ldmat) {
    int i = 0;
    int j = 0;
    std::vector<long double> tvec;
    ldmat = {};

    for (i = 0; i <= nx; i++) {
        tvec = {};
        for (j = 0; j <= ny; j++) tvec.push_back(0.0);
        ldmat.push_back(tvec);
    }
}

int solvematrix01LD(std::vector<std::vector<long double>> const& inmat, int eqnum, std::vector<long double> &outvec,
                    int verbose) {
    int eqhi, termhi, eqct, termct, i, j;
    long double max, pivot;
    std::vector<std::vector<long double>> newmat;
    std::vector<long double> coeffvec;
    std::vector<long double> outvec2;

    eqhi = termhi = eqct = termct = i = j = 0;
    max = pivot = 0.0;

    if (eqnum == 1) {
        if (inmat[0][1] != 0.0) {
            outvec[0] = -inmat[0][0] / inmat[0][1];
            return (0);
        } else {
            /*The coefficient for x1 was zero, so it is
              impossible to solve*/
            printf("ERROR: solvematrix01 fed a singular matrix!\n");
            outvec[0] = 0.0;
            return (1);
        }
    } else {
        make_LDmat(eqnum - 1, eqnum, newmat);
        make_LDvec(eqnum, coeffvec);
        make_LDvec(eqnum - 1, outvec2);
        /*REDUCE THE NUMBER OF EQUATIONS BY 1*/
        /*Find the coefficient with the largest absolute value*/
        eqhi = 0;
        termhi = 1;
        max = fabs(inmat[0][1]);
        for (eqct = 0; eqct < eqnum; eqct++) {
            for (termct = 1; termct < eqnum + 1; termct++) {
                if (max <= fabs(inmat[eqct][termct])) {
                    max = fabs(inmat[eqct][termct]);
                    eqhi = eqct;
                    termhi = termct;
                }
            }
        }
        pivot = inmat[eqhi][termhi];
        if (verbose >= 1)
            printf("At %Lf, coefficent %d of equation %d was the largest\n", pivot, termhi - 1, eqhi);
        if (max == 0.0) {
            printf("ERROR: solvematrix01 fed a singular matrix!\n");
            for (eqct = 0; eqct < eqnum; eqct++) outvec[eqct] = 0.0;
            return (1);
        }
        /*Solve equation eqhi for the x value corresponding to termhi*/
        j = 0;
        coeffvec[0] = inmat[eqhi][0] / pivot;
        for (termct = 1; termct < eqnum + 1; termct++) {
            if (termct != termhi) {
                j += 1;
                coeffvec[j] = inmat[eqhi][termct] / pivot;
            }
        }
        if (verbose >= 1) printf("Coefficient substitution vector:\n");
        if (verbose >= 1) printf("%Lf", coeffvec[0]);
        if (verbose >= 1)
            for (j = 1; j < eqnum; j++) printf(" %Lf", coeffvec[j]);
        if (verbose >= 1) printf("\n");
        /*Substitute this solution into the other equations,
          creating a new matrix with one fewer equations*/
        i = 0;
        for (eqct = 0; eqct < eqnum; eqct++) {
            if (eqct != eqhi) {
                j = 0;
                newmat[i][j] = inmat[eqct][0] - coeffvec[0] * inmat[eqct][termhi];
                for (termct = 1; termct < eqnum + 1; termct++) {
                    if (termct != termhi) {
                        j += 1;
                        newmat[i][j] = inmat[eqct][termct] - coeffvec[j] * inmat[eqct][termhi];
                    }
                }
                i += 1;
            }
        }
        if (verbose >= 1) printf("New reduced matrix:\n");
        for (i = 0; i < eqnum - 1; i++) {
            if (verbose >= 1) printf("%Lf", newmat[i][0]);
            if (verbose >= 1)
                for (j = 1; j < eqnum; j++) printf(" %Lf", newmat[i][j]);
            if (verbose >= 1) printf("\n");
        }
        /*Call solvematrix01 recursively on this new matrix*/
        if (solvematrix01LD(newmat, eqnum - 1, outvec2, verbose)) {
            printf("ERROR: recursive call of solvematrix01 failed\n");
            for (eqct = 0; eqct < eqnum; eqct++) outvec[eqct] = 0.0;
            return (1);
        }
        if (verbose >= 1) printf("Recursive result\n");
        if (verbose >= 1) printf("%Lf", outvec2[0]);
        if (verbose >= 1)
            for (i = 1; i < eqnum - 1; i++) printf(" %Lf", outvec2[i]);
        if (verbose >= 1) printf("\n");
        /*Load the solution for everything except the pivot*/
        i = 0;
        for (eqct = 0; eqct < eqnum; eqct++) {
            if (eqct != termhi - 1) {
                outvec[eqct] = outvec2[i];
                i += 1;
            }
        }
        /*Load the solution for the pivot*/
        outvec[termhi - 1] = -coeffvec[0];
        for (i = 0; i < eqnum - 1; i++) outvec[termhi - 1] -= coeffvec[i + 1] * outvec2[i];
    }
    return (0);
}

int perfectpoly01LD(std::vector<long double> const& x, std::vector<long double> const& y, std::vector<long double> &fitvec) {
    std::vector<std::vector<long double>> dmatrix;
    std::vector<long double> outvec;
    int i = 0;
    int j = 0;
    int k = 0;
    int status = 0;
    int npoints = x.size();
    if (y.size() != npoints) {
        cerr << "ERROR: x and y vectors in perfectpoly don't have the same number of points!\n";
        return (1);
    }
    if (npoints <= 1) {
        cerr << "ERROR: perfectpoly cannot fit just a single point!\n";
        return (2);
    }
    make_LDmat(npoints, npoints + 1, dmatrix);
    make_LDvec(npoints, outvec);
    // cout << "perfectpoly fitting matrix:\n";
    for (i = 0; i < npoints; i++) {
        dmatrix[i][0] = -y[i];
        // cout << dmatrix[i][0] << " ";
        for (j = 1; j <= npoints; j++) {
            dmatrix[i][j] = 1.0;
            for (k = 2; k <= j; k++) dmatrix[i][j] *= x[i];
            // cout << dmatrix[i][j] << " ";
        }
        // cout << "\n";
    }
    status = solvematrix01LD(dmatrix, npoints, fitvec, 0);
    return (status);
}

// Given a vector of MJD values and a vector of 3-D planet positions,
// use polynomial interpolation to obtain a precise estimate of the
// 3-D planet position at the time detmjd. It is assumed that the
// input time detmjd is in UT1 or some reasonable approximation
// thereof, while the planet ephemeris vectors are in dynamical TT.
// Hence, a correction is applied to the input time before the
// interpolation. If the calling function actually has time in TT
// already, planetpos01 should be called with the correction
// pre-subtracted from detmjd, so it will cancel out internally.
int planetpos01LD(long double detmjd, int polyorder, std::vector<long double> const& posmjd,
                  std::vector<Point3LD> const& planetpos, Point3LD &outpos) {
    int fitnum = polyorder + 1;
    int pointsbefore = fitnum - fitnum / 2;
    int pbf = 0;
    vector<long double> xvec;
    vector<long double> yvec;
    vector<long double> fitvec;
    long double tdelt = 0;
    long double sumvar = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    make_LDvec(fitnum, fitvec);

    // Convert input time from UT1 standard to the dynamical time (TT) used
    // by JPL Horizons for the state vectors.
    detmjd += TTDELTAT / SOLARDAY;
    // TT is ahead of UT1 because the Earth's rotation is slowing down.

    // Interpolate to find the planet's exact position at the time
    // of the detection.
    pbf = 0;
    i = posmjd.size();
    if (planetpos.size() != i) {
        cerr << "ERROR: planetpos01 finds time and position vectors\n";
        cerr << "to have different lengths\n";
        return (1);
    }
    while (i > 0 && pbf < pointsbefore) {
        i--;
        if (posmjd[i] < detmjd) pbf++;
    }
    pbf = i;
    xvec = {};
    yvec = {};
    tdelt = detmjd - posmjd[pbf];
    // Load vectors to fit x-coordinate of the planet's position.
    for (i = pbf; i < pbf + fitnum; i++) {
        xvec.push_back(posmjd[i] - posmjd[pbf]);
        yvec.push_back(planetpos[i].x);
    }
    // Solve for polynomial interpolation
    perfectpoly01LD(xvec, yvec, fitvec);
    // Calculate interpolated position.
    outpos.x = fitvec[0];
    for (j = 1; j < fitnum; j++) {
        sumvar = fitvec[j] * tdelt;
        for (k = 2; k <= j; k++) sumvar *= tdelt;
        outpos.x += sumvar;
    }
    // Load vector to fit y-coordinate
    yvec = {};
    for (i = pbf; i < pbf + fitnum; i++) yvec.push_back(planetpos[i].y);
    // Solve for polynomial interpolation
    perfectpoly01LD(xvec, yvec, fitvec);
    // Calculate interpolated position.
    outpos.y = fitvec[0];
    for (j = 1; j < fitnum; j++) {
        sumvar = fitvec[j] * tdelt;
        for (k = 2; k <= j; k++) sumvar *= tdelt;
        outpos.y += sumvar;
    }
    // Load vector to fit z-coordinate
    yvec = {};
    for (i = pbf; i < pbf + fitnum; i++) yvec.push_back(planetpos[i].z);
    // Solve for polynomial interpolation
    perfectpoly01LD(xvec, yvec, fitvec);
    // Calculate interpolated position.
    outpos.z = fitvec[0];
    for (j = 1; j < fitnum; j++) {
        sumvar = fitvec[j] * tdelt;
        for (k = 2; k <= j; k++) sumvar *= tdelt;
        outpos.z += sumvar;
    }
    return (0);
}

int celestial_to_statevecLD(long double RA, long double Dec, long double delta, Point3LD &baryvec) {
    long double x,y,z,theta,phi,thetapole,phipole;
    x = y = z = theta = phi = thetapole = phipole = 0.0;
    theta = Dec/DEGPRAD;
    phi = RA/DEGPRAD;
    thetapole = NEPDEC/DEGPRAD;

    z = sin(theta);
    x = -cos(theta)*sin(phi);
    y = cos(theta)*cos(phi);
    baryvec.z = delta*(z*sin(thetapole) + x*cos(thetapole));
    baryvec.y = delta*(z*cos(thetapole) - x*sin(thetapole));
    baryvec.x = delta*y;
    return(0);
}

int precess01aLD(long double ra1, long double dec1,long double mjd, long double *ra2, long double *dec2, int precesscon) {
    long double ndays,tds,zetaa,thetaa,zaa,ra4,dec4,cosra,sinra;

    /*time since standard epoch*/
    ndays = mjd-51544L; /*Number of days since Jan 1, 2000*/
    tds = ndays/36525.0L;

    /*cubic approximation to precession*/
    zetaa = ZET0 + ZET1*tds + ZET2*tds*tds + ZET3*tds*tds*tds + ZET4*tds*tds*tds*tds + ZET5*tds*tds*tds*tds*tds;
    zaa = Z0 + Z1*tds + Z2*tds*tds + Z3*tds*tds*tds + Z4*tds*tds*tds*tds + Z5*tds*tds*tds*tds*tds;
    thetaa = THET1*tds + THET2*tds*tds + THET3*tds*tds*tds + THET4*tds*tds*tds*tds + THET5*tds*tds*tds*tds*tds;

    /*transformation from arcseconds to radians*/
    zetaa*=(M_PI/648000.0L);
    zaa*=(M_PI/648000.0L);
    thetaa*=(M_PI/648000.0L);

    if(precesscon>=0) {
        /*Precess given J2000.0 coords to epoch of date*/

        /*get new declination*/
        if(dec1!=M_PI/2.0L) {
            /*printf("precess01 has normal declination case\n");*/
            dec4 = asin(cos(ra1+zetaa)*sin(thetaa)*cos(dec1) + cos(thetaa)*sin(dec1));
        } else {
            /*printf("precess01 has polar declination case\n");*/
            dec4 = asin(cos(thetaa));
        }
        /*if declination was obviously meant to be the pole, but
            it has gotten a little off by roundoff error, collapse
            it to the pole.*/
        if(fabs(dec4-M_PI/2.0L)<SMALLANG) {
            dec4 = M_PI/2.0L;
        }
        /*get new right ascension*/
        if(dec1!=M_PI/2.0L && dec4!=M_PI/2.0L) {
            /*printf("precess01 has normal right ascension case\n");*/
            cosra = (cos(ra1+zetaa)*cos(thetaa)*cos(dec1) - sin(thetaa)*sin(dec1))/cos(dec4);
            sinra = (sin(ra1+zetaa)*cos(dec1))/cos(dec4);
            if(sinra>=0.0) {
                ra4 = acos(cosra)+zaa;
            } else{
                ra4 = 2.0L*M_PI - acos(cosra)+zaa;
            }
        } else if(dec1==M_PI/2.0L && dec4!=M_PI/2.0L) {
            /*printf("precess01 has polar input right ascension case\n");*/
            ra4 = M_PI + zaa;
        } else if(dec4==M_PI/2.0L) {
            /*printf("precess01 has polar output right ascension case\n");*/
            ra4 = 0.0;
        } else {
            printf("IMPOSSIBLE CASE ERROR IN precess01\n");
        }
    } else {
        /*Deprecess given epoch of date coords to J2000.0*/

        /*get new declination*/
        if(dec1!=M_PI/2.0L) {
            /*printf("precess01 has normal declination case\n");*/
            dec4 = asin(-cos(ra1-zaa)*sin(thetaa)*cos(dec1) + cos(thetaa)*sin(dec1));
        } else {
            /*printf("precess01 has polar declination case\n");*/
            dec4 = asin(cos(thetaa));
        }
        /*if declination was obviously meant to be the pole, but
            it has gotten a little off by roundoff error, collapse
            it to the pole.*/
        if(fabs(dec4-M_PI/2.0L)<SMALLANG) {
            dec4 = M_PI/2.0L;
        }
        /*get new right ascension*/
        if(dec1!=M_PI/2.0L && dec4!=M_PI/2.0L) {
            /*printf("precess01 has normal right ascension case\n");*/
            cosra = (cos(ra1-zaa)*cos(thetaa)*cos(dec1) + sin(thetaa)*sin(dec1))/cos(dec4);
            sinra = (sin(ra1-zaa)*cos(dec1))/cos(dec4);
            if(sinra>=0.0) {
                ra4 = acos(cosra)-zetaa;
            } else{
                ra4 = 2.0L*M_PI - acos(cosra)-zetaa;
            }
        } else if(dec1==M_PI/2.0L && dec4!=M_PI/2.0L) {
            /*printf("precess01 has polar input right ascension case\n");*/
            ra4 = 2.0L*M_PI-zetaa; /*Note this could be wrong*/
                                /*There might be two solutions*/
        } else if(dec4==M_PI/2.0L) {
            /*printf("precess01 has polar output right ascension case\n");*/
            ra4 = 0.0;
        } else {
            printf("IMPOSSIBLE CASE ERROR IN precess01\n");
        }
    }
    *ra2 = ra4;
    *dec2 = dec4;
    return(1);
}

// Given the MJD of an observation, and a file giving barycentric state-vector
// coordinates for the Earth, the longitude and MPC latitude sin and cos terms
// for an observatory, calculate the observer's topocentric position
// in barycentric state vector coordinates.
// Note that the handling of Earth's rotation assumes that the
// input MJD is UT1, while the ephemeris vectors posmjd
// and planetpos are in dynamical TT. Hence, after calculating
// aspects related to Earth's rotation with detmjd as input,
// planetpos01 is called which internally converts the input
// UT1 into TT.
int observer_barycoords01LD(
    long double detmjd, int polyorder, long double lon, long double obscos,
    long double obssine, std::vector<long double> const& posmjd,
    std::vector<Point3LD> const& planetpos, Point3LD &outpos) {

    long double gmst=0;
    long double djdoff = detmjd-51544.5L;
    long double zenithRA=0.0;
    long double zenithDec=0.0;
    long double junkRA=0.0;
    long double junkDec=0.0;
    long double crad = sqrt(obscos*obscos + obssine*obssine)*EARTHEQUATRAD;
    Point3LD obs_from_geocen = Point3LD(0,0,0);
    Point3LD geocen_from_barycen = Point3LD(0,0,0);

    gmst = 18.697374558L + 24.06570982441908L*djdoff;
    // Add the longitude, converted to hours.
    // Note: at this point it stops being gmst.
    gmst += lon/15.0L;
    // Get a value between 0 and 24.0.
    while(gmst>=24.0L) gmst -= 24.0L;
    while(gmst<0.0L) gmst += 24.0L;
    // Convert to degrees
    zenithRA = gmst * 15.0L;
    // Get zenithDec
    if(obscos!=0.0L) {
        zenithDec = atan(obssine/obscos)*DEGPRAD;
    } else if(obssine>=0.0L) {
        zenithDec = 90.0L;
    } else {
        zenithDec=-90.0L;
    }
    // Now zenithRA and zenithDec are epoch-of-date coordinates.
    // If you want them in J2000.0, this is the place to convert them.
    int precesscon=-1; //Precess epoch-of-date to J2000.0
    junkRA = zenithRA/DEGPRAD;
    junkDec = zenithDec/DEGPRAD;
    precess01aLD(junkRA,junkDec,detmjd,&zenithRA,&zenithDec,precesscon);
    zenithRA*=DEGPRAD;
    zenithDec*=DEGPRAD;
    celestial_to_statevecLD(zenithRA,zenithDec,crad,obs_from_geocen);
    // crad is the distance from the geocenter to the observer, in AU.
    planetpos01LD(detmjd,polyorder,posmjd,planetpos,geocen_from_barycen);
    // cout << "obs_from_geocen: " << obs_from_geocen.x << " " << obs_from_geocen.y << " " << obs_from_geocen.z << " \n";
    // cout << "geocen_from_barycen: " << geocen_from_barycen.x << " " << geocen_from_barycen.y << " " << geocen_from_barycen.z << "\n";
    outpos.x = geocen_from_barycen.x + obs_from_geocen.x;
    outpos.y = geocen_from_barycen.y + obs_from_geocen.y;
    outpos.z = geocen_from_barycen.z + obs_from_geocen.z;
    return(0);
}

long medindex(std::vector<XYIndex> const& xyvec, int dim) {
    std::vector<XYIndex> xyv = xyvec; //Mutable copy of immutable input vector
    for(int i=0; i<xyv.size(); i++) xyv[i].index=i; //Redefine indices
    long medpt = xyv.size()/2;
    if(dim%2==1) sort(xyv.begin(), xyv.end(), compareXYIndexX);
    else sort(xyv.begin(), xyv.end(), compareXYIndexY);
    return(xyv[medpt].index);
}

// Given double precision RA, Dec in DEGREES, project
// onto the unit sphere and return an object of class point3d.
// Input coordinates are in degrees, input RA=0, Dec=0
// projects to x=1,y=0,z=0; then y increases for positive
// RA.
Point3D celeproj01(double RA, double Dec) {
    return(Point3D( cos(RA/DEGPRAD)*cos(Dec/DEGPRAD) , sin(RA/DEGPRAD)*cos(Dec/DEGPRAD), sin(Dec/DEGPRAD)));
}

// Given a 3-d point (class point3d), de-project it back to
// celestial coordinates IN DEGREES: i.e., reverse the process
// carried out by celeproj01.
int celedeproj01(Point3D p3, double *RA, double *Dec) {
    //Normalize the point
    double norm = sqrt(p3.x*p3.x + p3.y*p3.y + p3.z*p3.z);
    if(norm<=0.0) {
        *RA=0.0;
        *Dec=0.0;
        return(1);
    }
    double x = p3.x/norm;
    double y = p3.y/norm;
    double z = p3.z/norm;
    if(fabs(z)<=1.0) *Dec = asin(z)*DEGPRAD;
    else return(2);
    if(y==0 && x<0.0) {
        // y is zero and x is negative
        *RA = 180.0;
        return(0);
    }
    else if(y==0.0) {
        // y is zero and x is zero or positive
        *RA=0.0;
        return(0);
    }
    else if(y>0.0) {
        // y is strictly positive
        *RA = 90.0 - atan(x/y)*DEGPRAD;
        return(0);
    }
    else if(y<0.0) {
        // y is strictly negative
        *RA = 270.0 - atan(x/y)*DEGPRAD;
        return(0);
    }
    else {
        // Weird case, should be impossible
        return(3);
    }
};

// Given two pairs of RA, Dec coordinates, calculate
// their angular separation on the sky in degrees.
double distradec01(double RA1, double Dec1, double RA2, double Dec2) {
    double x1,y1,z1,x2,y2,z2,h;
    x1=cos(Dec1/DEGPRAD)*cos(RA1/DEGPRAD);
    y1=cos(Dec1/DEGPRAD)*sin(RA1/DEGPRAD);
    z1=sin(Dec1/DEGPRAD);
    x2=cos(Dec2/DEGPRAD)*cos(RA2/DEGPRAD);
    y2=cos(Dec2/DEGPRAD)*sin(RA2/DEGPRAD);
    z2=sin(Dec2/DEGPRAD);
    h=sqrt(DSQUARE(x1-x2)+DSQUARE(y1-y2)+DSQUARE(z1-z2));
    return(DEGPRAD*2.0*asin(h/2.0));
}

// calculates the celestial position angle from point 1 to point 2.
int distradec02(double ra1,double dec1,double ra2,double dec2,double *dist,double *pa) {
    double x1,y1,z1,x2,y2,z2,h,d,poleangle,celpa,sinepa,cosinepa,colat1,colat2;
    double arcsinepa=0.0l;

    // Handle trivial cases.
    if(fabs(ra1-ra2)/DEGPRAD < VSMALLANG) {
        // the two RA values are functionally identical
        if(fabs(dec1-dec2)/DEGPRAD < VSMALLANG) {
            // the two Dec values are functionally identical
            *dist=0.0l;
            *pa=0.0l; // Dummy value for zero distance.
            return(0);
        } else if(dec2<dec1) {
        *dist = dec1-dec2;
        *pa = 180.0l; // Due South
        return(0);
        } else {
            // dec2>dec1 by logical elimination
            *dist = dec2-dec1;
            *pa = 0.0l; // Due North
            return(0);
        }
    } else if(fabs(dec1-dec2)/DEGPRAD < VSMALLANG) {
        // the two Dec values are functionally identical,
        // although the two RA values are not.
        if(ra2<ra1) {
            *dist = ra1-ra2;
            *pa = 270.0l; // Due West.
            return(0);
        } else {
            // ra2>ra1 by logical elimination
            *dist = ra2-ra1;
            *pa = 90.0l; // Due East.
            return(0);
        }
    } else {
        // We have a non-trivial case
        // Calculate the distance d, in radians, between
        // the two points.
        x1=cos(dec1/DEGPRAD)*cos(ra1/DEGPRAD);
        y1=cos(dec1/DEGPRAD)*sin(ra1/DEGPRAD);
        z1=sin(dec1/DEGPRAD);
        x2=cos(dec2/DEGPRAD)*cos(ra2/DEGPRAD);
        y2=cos(dec2/DEGPRAD)*sin(ra2/DEGPRAD);
        z2=sin(dec2/DEGPRAD);
        h=sqrt(DSQUARE(x1-x2)+DSQUARE(y1-y2)+DSQUARE(z1-z2));
        if(h/2.0l <= 1.0l) {
            d=2.0*asin(h/2.0l);
        } else {
            cerr << "WARNING: distradec02 attempting to take arcsine of 1 + " << h/2.0l - 1.0l << "\n";
            d = M_PI/2.0l;
        }
        *dist = d*DEGPRAD;

        colat1 = M_PI/2.0 - dec1/DEGPRAD;
        colat2 = M_PI/2.0 - dec2/DEGPRAD;

        // Calculate the difference in RA, paying careful
        // attention to wrapping.
        cosinepa = (cos(colat2) - cos(d)*cos(colat1))/(sin(d)*sin(colat1));
        if(ra1<ra2 && (ra2-ra1)<=180.0) {
            // Simple case, point 2 is east of point 1,
            // so PA should be less than 180 degrees.
            poleangle = (ra2-ra1)/DEGPRAD;
            sinepa = sin(colat2)*sin(poleangle)/sin(d);
            if(sinepa<=1.0l) {
	            arcsinepa = asin(sinepa);
            } else {
	            cerr << "WARNING: distradec02 attempting to take the arcsine of 1 + " << sinepa-1.0l << "\n";
	            arcsinepa = M_PI/2.0l;
            }
            if(cosinepa>=0.0) celpa = arcsinepa;
            else celpa = M_PI - arcsinepa;
            *pa = celpa*DEGPRAD;
        } else if(ra1<ra2) {
            // Wrapped case with point 2 west of point 1
            // across zero RA: the PA should be greater
            // than 180 degrees.
            poleangle = (ra1+(double)360.0-ra2)/DEGPRAD;
            sinepa = sin(colat2)*sin(poleangle)/sin(d);
            if(sinepa<=1.0l) {
	            arcsinepa = asin(sinepa);
            } else {
	            cerr << "WARNING: distradec02 attempting to take the arcsine of 1 + " << sinepa-1.0l << "\n";
	            arcsinepa = M_PI/2.0l;
            }
            if(cosinepa>=0.0) celpa = arcsinepa;
            else celpa = M_PI - arcsinepa;
            *pa = (double)360.0 - celpa*DEGPRAD;
        } else if(ra1>ra2 && (ra1-ra2)<=180.0) {
            // Simple case, point 2 is west of point 1,
            // so PA should be greater than 180 degrees.
            poleangle = (ra1-ra2)/DEGPRAD;
            sinepa = sin(colat2)*sin(poleangle)/sin(d);
            if(sinepa<=1.0l) {
	            arcsinepa = asin(sinepa);
            } else {
	            cerr << "WARNING: distradec02 attempting to take the arcsine of 1 + " << sinepa-1.0l << "\n";
	            arcsinepa = M_PI/2.0l;
            }
            if(cosinepa>=0.0) celpa = arcsinepa;
            else celpa = M_PI - arcsinepa;
            *pa = (double)360.0 - celpa*DEGPRAD;
        } else if(ra1>ra2) {
            // Wrapped case with point 2 east of point 1
            // across zero RA: the PA should be less
            // than 180.0 degrees.
            poleangle = (ra2+(double)360.0-ra1)/DEGPRAD;
            sinepa = sin(colat2)*sin(poleangle)/sin(d);
            if(sinepa<=1.0l) {
	            arcsinepa = asin(sinepa);
            } else {
	            cerr << "WARNING: distradec02 attempting to take the arcsine of 1 + " << sinepa-1.0l << "\n";
	            arcsinepa = M_PI/2.0l;
            }
            if(cosinepa>=0.0) celpa = arcsinepa;
            else celpa = M_PI - arcsinepa;
            *pa = celpa*DEGPRAD;
        }
        return(0);
    }
    return(0);
}

int splitxy(std::vector<XYIndex> const& xyvec, int dim, long splitpoint,
            std::vector<XYIndex> & left, std::vector<XYIndex> & right) {
    long i=0;
    double xval = xyvec[splitpoint].x;
    double yval = xyvec[splitpoint].y;

    if(dim%2==1) {
        // Split on x
        for(i=0 ; i<xyvec.size() ; i++) {
            if(i!=splitpoint && xyvec[i].x<=xval) {
                left.push_back(xyvec[i]);
            } else if(i!=splitpoint) {
                right.push_back(xyvec[i]);
            }
        }
    } else {
        // Split on y
        for(i=0 ; i<xyvec.size() ; i++) {
            if(i!=splitpoint && xyvec[i].y<=yval) {
                left.push_back(xyvec[i]);
            } else if(i!=splitpoint) right.push_back(xyvec[i]);
        }
    }
    return(0);
}

int kdtree01(std::vector<XYIndex> const& xyvec, int dim, long rootptxy, long rootptkd, std::vector<KDPoint> &kdvec) {
    int status=0;
    int lmed=0;
    int rmed=0;
    int kdct = kdvec.size()-1;
    int i=0;
    long leftrootkd=-1;
    long rightrootkd=-1;
    XYIndex xyi = XYIndex(0.0,0.0,0);
    KDPoint root = kdvec[kdct];
    KDPoint lp = KDPoint(xyi,-1,-1,0);
    KDPoint rp = KDPoint(xyi,-1,-1,0);
    KDPoint kdtest = KDPoint(xyi,-1,-1,0);
    std::vector<XYIndex> leftvec = {};
    std::vector<XYIndex> rightvec = {};

    status = splitxy(xyvec,dim,rootptxy,leftvec,rightvec);

    if(dim==1) dim=2;
    else dim=1;

    if(leftvec.size()==1) {
        // Left branch is just a single leaf
        lp = KDPoint(leftvec[0],-1,-1,dim);
        kdvec.push_back(lp);
        kdct++;
        kdvec[rootptkd].left = kdct;
        kdtest = kdvec[kdct];
    } else if(leftvec.size()<=0) {
        // There is no left branch
        kdvec[rootptkd].left = -1;
    }
    if(rightvec.size()==1) {
        // Right branch is just a single leaf
        rp = KDPoint(rightvec[0],-1,-1,dim);
        kdvec.push_back(rp);
        kdct++;
        kdvec[rootptkd].right = kdct;
        kdtest = kdvec[kdct];
    } else if(rightvec.size()<=0) {
        // There is no right branch
        kdvec[rootptkd].right = -1;
    }

    if(leftvec.size()>1) {
        lmed = medindex(leftvec,dim);
        lp = KDPoint(leftvec[lmed],-1,-1,dim);
        kdvec.push_back(lp);
        kdct++;
        kdvec[rootptkd].left = kdct;
        leftrootkd = kdct;
        kdtest = kdvec[kdct];
    }

    if(rightvec.size()>1) {
        rmed = medindex(rightvec,dim);
        rp = KDPoint(rightvec[rmed],-1,-1,dim);
        kdvec.push_back(rp);
        kdct++;
        kdvec[rootptkd].right = kdct;
        rightrootkd = kdct;
        kdtest = kdvec[kdct];
    }
    // I moved these down out of the above loops, because I thought
    // that otherwise, a bunch of stuff might get pushed down by the
    // left loop that the right loop didn't know about.
    if(leftvec.size()>1 && leftrootkd>=0) kdtree01(leftvec,dim,lmed,leftrootkd,kdvec);
    else if(leftvec.size()>1 && leftrootkd<0) {
        cerr << "Error, kdtree01 finds leftroot less than zero with leftvec.size() = " << leftvec.size() << "\n";
    }
    if(rightvec.size()>1 && rightrootkd>=0) kdtree01(rightvec,dim,rmed,rightrootkd,kdvec);
    else if(rightvec.size()>1 && rightrootkd<0) {
        cerr << "Error, kdtree01 finds rightroot less than zero with rightvec.size() = " << rightvec.size() << "\n";
    }

    return(0);
}

// Given a k-d tree vector kdvec created by kdtree01,
// perform a range-query about the point x,y. Returns
// a vector indexing all of the points in the input k-d tree
// that lie within the specified range of the input coordinates.
// Assumes that kdvec[0] is the root of the k-d tree.
// NOTE THAT THIS IS FOR 2-D KD trees.
int kdrange01(std::vector<KDPoint> const& kdvec, double x, double y, double range, std::vector<long> &indexvec) {
    int branchct=0;
    double rng2 = range*range;
    int notdone=1;
    int kdveclen = kdvec.size();
    int dim=1;
    int currentpoint=0;
    int leftpoint=0;
    int rightpoint=0;
    int goleft=0;
    int goright=0;
    double xdiff=0.0;
    double ydiff=0.0;
    std::vector<long> checkit={};
    int i=0;
    int checknum=0;

    while(notdone>0) {
        // Climb to the top of the k-d tree, keeping track
        // of potentially interesting unexplored branches
        // in the vector checkit.
        while(leftpoint>=0 || rightpoint>=0) {
            // Previous step did not end on a leaf.
            leftpoint = kdvec[currentpoint].left;
            rightpoint = kdvec[currentpoint].right;
            dim = kdvec[currentpoint].dim;
            if(dim%2==1) {
                xdiff = kdvec[currentpoint].point.x - x;
                goright = (xdiff <= range); // possible hits lie to the left;
                goleft = (xdiff >= -range); // possible hits lie to the right;
                if(goleft && goright) {
                    // Current point might be within range.
                    ydiff = kdvec[currentpoint].point.y - y;
                    if(fabs(ydiff)<=range && (xdiff*xdiff + ydiff*ydiff)<=rng2) {
                        // Current point is within range. Add it to the output vector
                        indexvec.push_back(currentpoint);
                    }
                    if(leftpoint>=0) {
                        //Explore leftward first.
                        currentpoint = leftpoint;
                        if(rightpoint>=0) {
                            // Rightward branch will also be explored later
                            checknum++;
                            if(checknum>checkit.size()) {
                                checkit.push_back(rightpoint);
                            } else {
                                checkit[checknum-1] = rightpoint;
                            }
                        }
                    } else if(rightpoint>=0) {
                        // Leftward branch is a dead end: explore rightward branch
                        currentpoint = rightpoint;
                    }
                } else if(goleft) {
                    // Current point cannot be in range, but points that
                    // are in range may lie along the left branch.
                    if(leftpoint>=0) {
                        currentpoint = leftpoint;
                    } else rightpoint=-1; // Dead end, make sure while-loop exits.
                } else if(goright) {
                    // Current point cannot be in range, but points that
                    // are in range may lie along the right branch.
                    if(rightpoint>=0) {
                        currentpoint = rightpoint;
                    } else leftpoint=-1;  // Dead end, make sure while-loop exits.
                } else {
                    // Program concluded it should go neither left nor right.
                    // The likely cause is that it encountered a NAN. Give up on this point.
                    leftpoint=rightpoint=-1;
                    cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
                    cerr << "Input point " << x << ", " << y <<", target point " << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ".\n";
                }
                // Close x-dimension case
            } else if(dim%2==0) {
                ydiff = kdvec[currentpoint].point.y - y;
                goright = (ydiff <= range); // possible hits lie to the left;
                goleft = (ydiff >= -range); // possible hits lie to the right;
                if(goleft && goright) {
                    // Current point might be within range.
                    xdiff = kdvec[currentpoint].point.x - x;
                    if(fabs(ydiff)<=range && (xdiff*xdiff + ydiff*ydiff)<=rng2) {
                        // Current point is within range. Add it to the output vector
                        indexvec.push_back(currentpoint);
                    }
                    if(leftpoint>=0) {
                        //Explore leftward first.
                        currentpoint = leftpoint;
                        if(rightpoint>=0) {
                            // Rightward branch will also be explored later
                            checknum++;
                            if(checknum>checkit.size()) {
                                checkit.push_back(rightpoint);
                            } else {
                                checkit[checknum-1] = rightpoint;
                            }
                        }
                    } else if(rightpoint>=0) {
                        // Leftward branch is a dead end: explore rightward branch
                        currentpoint = rightpoint;
                    }
                } else if(goleft) {
                    // Current point cannot be in range, but points that
                    // are in range may lie along the left branch.
                    if(leftpoint>=0) {
                        currentpoint = leftpoint;
                    } else rightpoint = -1; // Dead end, make sure while loop exits.
                } else if(goright) {
                    // Current point cannot be in range, but points that
                    // are in range may lie along the right branch.
                    if(rightpoint>=0) {
                        currentpoint = rightpoint;
                    } else leftpoint=-1;  // Dead end, make sure while loop exits.
                } else {
                    // Program concluded it should go neither left nor right.
                    // The likely cause is that it encountered a NAN. Give up on this point.
                    leftpoint=rightpoint=-1;
                    cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
                    cerr << "Input point " << x << ", " << y <<", target point " << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ".\n";
                }
                // Note that we do not need to worry about the possiblity
                // that current point will get set to -1: i.e., we were
                // at a leaf or a one-sided branch. Such cases will
                // be caught by the while statement.
                // Close y-dimension case
            }
            // Close while-loop checking if we've hit a leaf.
        }
        // We have climbed up the tree to a leaf. Go backwards through
        // the checkit vector and see if there is anything to check.
        checknum=checkit.size();
        while(checknum>=1 && checkit[checknum-1]<0) checknum--;
        if(checknum<=0) {
            //There were no valid entries to check: we're done.
            notdone=0;
        } else {
            //Set currentpoint to the last valid entry in checkit
            currentpoint = checkit[checknum-1];
            //Mark this point as used.
            checkit[checknum-1]=-1;
            leftpoint=rightpoint=0;
        }
    }
    return(0);
}

// Simple and crude utility program, does an unweighted
// linear fit of the form y = mx * b, for m = slope, b = intercept
int linfituw01(std::vector<double> const& x, std::vector<double> const& y, double &slope, double &intercept) {
    int i;
    int pointnum = x.size();
    double delta,xal,yal,xty,xsq,nsum,rms,err,errmax;
    double siga,sigb;

    if(pointnum<=1) {
        cerr << "ERROR: linfituw01 CALLED WITH ONLY ONE POINT\n";
        return(1);
    }

    xal = yal = xty = xsq = nsum = 0.0;
    for(i=0;i<pointnum;i++) {
        xal += x[i];
        yal += y[i];
        xsq += x[i]*x[i];
        xty += x[i]*y[i];
        nsum += 1.0l;
    }
    delta = nsum*xsq - xal*xal;
    if(delta==0.0) {
        cerr << "ERROR: linfituw01 has non-finite slope\n";
        return(1);
    }
    intercept = (xsq*yal - xal*xty)/delta;
    slope = (nsum*xty - xal*yal)/delta;

    return(0);
}

// Given a central point and the arc distance and celestial
// position angle to second point, calculate the celestial
// coordinates of this second point. All input and output
// quantities are in degrees. Note that this is the reverse
// process of, e.g. distradec02, which finds
// the position angle and arc distance between two points
// on the celestial sphere. The current program finds the
// celestial coordinates of the second point, given the
// first point, and the arc distance and position angle.
// NOTE WELL: here the arc dist is in degrees.
int arc2cel01(double racenter, double deccenter, double dist, double pa, double &outra,double &outdec) {
    double colat1,tpa,rpa,arc,coscolat,colat2;
    double cosdra,deltaRA;

    tpa=pa;
    while(tpa>=360.0l) tpa-=360.0l;
    while(tpa<0.0l) tpa+=360.0l;

    // Handle trivial cases
    if(dist==0.0l) {
        outra = racenter;
        outdec = deccenter;
        return(0);
    } else if(dist==180.0l) {
        outra = racenter + 180.0l;
        if(outra >= 360.0l) outra -= 360.0l;
        outdec = -deccenter;
        return(0);
    } else if(deccenter==90.0l) {
        cerr << "WARNING: arc2cel01 called with starting point at north pole!\n";
        outra = tpa;
        outdec = 90.0l - dist;
        return(0);
    } else if(deccenter==-90.0l) {
        cerr << "WARNING: arc2cel01 called with starting point at south pole!\n";
        outra = tpa;
        outdec = -90.0l + dist;
        return(0);
    }

    colat1 = M_PI/(double)2.0l - deccenter/DEGPRAD;
    rpa = tpa/DEGPRAD;
    arc = dist/DEGPRAD;

    coscolat = cos(arc)*cos(colat1) + sin(arc)*sin(colat1)*cos(rpa);
    if(coscolat>1.0l) {
        cerr << "WARNING: arc2cel01 attempting to take arccos of 1 + " << coscolat-1.0l << "\n";
        colat2 = 0.0l;
    } else colat2 = acos(coscolat);

    outdec = 90.0l - colat2*DEGPRAD;
    if(sin(colat2)<=0.0l) {
        outra = 0.0l;
        return(0);
    }

    cosdra = (cos(arc) - cos(colat1)*cos(colat2)) / (sin(colat1)*sin(colat2));
    if(cosdra>1.0l) {
        cerr  << "WARNING: arc2cel01 attempting to take arccos of 1 + " << cosdra-1.0l << "\n";
        deltaRA = 0.0l;
    } else deltaRA = acos(cosdra)*DEGPRAD;

    // Direction of RA change
    if(tpa<=180.0l) {
        // Change is to the east
        outra = racenter + deltaRA;
        while(outra>=360.0l) outra-=360.0l;
        return(0);
    } else {
        // Change is to the west
        outra = racenter - deltaRA;
        while(outra<0.0l) outra+=360.0l;
        return(0);
    }
    return(0);
}

std::vector<Observatory> readObscodeFile(
    string obscodefile
) {
    std::vector<Observatory> observatory_list = {};
    char obscode[MINSTRINGLEN];
    long double obslon, plxcos, plxsin;

    ifstream instream;
    instream.open(obscodefile);
    if (!instream) {
        cerr << "can't open input file " << obscodefile << "\n";
        return observatory_list;
    }
    // Skip one-line header
    string lnfromfile;
    string stest;
    getline(instream, lnfromfile);
    while (!instream.eof() && !instream.fail() && !instream.bad()) {
        instream >> stest;
        stringncopy01(obscode, stest, MINSTRINGLEN);
        instream >> obslon;
        instream >> plxcos;
        instream >> plxsin;
        Observatory obs1 = Observatory(obscode, obslon, plxcos, plxsin);
        observatory_list.push_back(obs1);
        // Skip the rest of the line
        getline(instream, lnfromfile);
    }
    instream.close();
    cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile
         << ":\n";

    if (DEBUG >= 2) {
        for (size_t i = 0; i < observatory_list.size(); i++) {
            cout << observatory_list[i].obscode << " " << observatory_list[i].obslon << " "
                 << observatory_list[i].plxcos << " " << observatory_list[i].plxsin << "\n";
        }
    }
    return observatory_list;
}

std::vector<Detection> readDetectionsFile(
    string indetfile,
    int idcol,
    int mjdcol,
    int racol,
    int deccol,
    int magcol,
    int bandcol,
    int obscodecol
) {
    std::vector<Detection> detvec = {};
    long double MJD, RA, Dec;
    double mag;
    char idstring[SHORTSTRINGLEN];
    char band[MINSTRINGLEN];
    char obscode[MINSTRINGLEN];

    ifstream instream;
    instream.open(indetfile);
    if (!instream) {
        cerr << "can't open input file " << indetfile << "\n";
        return detvec;
    }
    // Skip one-line header
    string lnfromfile, stest;
    getline(instream, lnfromfile);
    long lct = 1;
    int reachedeof = 0;
    cout << lnfromfile << "\n";
    while (reachedeof == 0) {
        getline(instream, lnfromfile);
        lct++;
        if (!instream.eof() && !instream.fail() && !instream.bad())
            ;  // Read on.
        else if (instream.eof())
            reachedeof = 1;  // End of file, fine.
        else if (instream.fail())
            reachedeof = -1;  // Something wrong, warn
        else if (instream.bad())
            reachedeof = -2;  // Worse problem, warn
        int i = 0;
        int j = 0;
        char c = '0';
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
            Detection o1 = Detection(MJD, RA, Dec, 0L, 0L, 0L, idstring, mag, band, obscode, -lct, {});
            detvec.push_back(o1);
        }
    }
    instream.close();
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
    sort(detvec.begin(), detvec.end(), timeCompareDetections);

    cout << "Last two obscodes: " << detvec[detvec.size() - 2].obscode << " and "
         << detvec[detvec.size() - 1].obscode << "\n";

    return detvec;
}

std::vector<ImageLog> readImageFile(
    MakeTrackletsConfig const& config,
    std::vector<Detection> const& detvec,
    string inimfile
) {
    std::vector<ImageLog> img_log = {};
    std::vector<ImageLog> img_log_tmp = {};
    long double MJD, RA, Dec;
    char obscode[MINSTRINGLEN];

    ifstream instream;
    instream.open(inimfile);
    if (!instream) {
        cerr << "can't open input file " << inimfile << "\n";
        return img_log;
    }
    int reachedeof = 0;
    string lnfromfile;
    while (reachedeof == 0) {
        getline(instream, lnfromfile);
        if (!instream.eof() && !instream.fail() && !instream.bad())
            ;  // Read on.
        else if (instream.eof())
            reachedeof = 1;  // End of file, fine.
        else if (instream.fail())
            reachedeof = -1;  // Something wrong, warn
        else if (instream.bad())
            reachedeof = -2;  // Worse problem, warn
        int i = 0;
        int j = 0;
        char c = '0';
        MJD = 0.0;
        while (i < lnfromfile.size() && reachedeof == 0) {
            string stest = "";
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
            ImageLog imlog = ImageLog(MJD, RA, Dec, obscode, 0, 0);
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
    sort(img_log_tmp.begin(), img_log_tmp.end(), timeCompareImageLog);
    // find the indices in the time-sorted detection file
    // that correspond to the earliest and latest detections
    // on each image, and load these values into imglog02.
    long detct = 0;
    for (size_t imct = 0; imct < img_log_tmp.size(); imct++) {
        while (detct < detvec.size() &&
                detvec[detct].MJD < img_log_tmp[imct].MJD - config.imagetimetol / SOLARDAY)
            detct++;  // Not on any image
        if (detct < detvec.size() &&
            fabs(detvec[detct].MJD - img_log_tmp[imct].MJD) <= config.imagetimetol / SOLARDAY) {
            // This should be the first detection on image imct.
            img_log_tmp[imct].startind = detct;
            while (detct < detvec.size() &&
                    fabs(detvec[detct].MJD - img_log_tmp[imct].MJD) <= config.imagetimetol / SOLARDAY)
                detct++;  // Still on this same image
            // This should be the first detection on the next image
            img_log_tmp[imct].endind = detct;
        }
        if (img_log_tmp[imct].startind >= 0 && img_log_tmp[imct].endind > 0) {
            img_log.push_back(img_log_tmp[imct]);
        }
    }
    instream.close();
    return img_log;
};

std::vector<ImageLog> makeImageLogs(
    MakeTrackletsConfig const& config,
    std::vector<Detection> const& detvec
) {
    // No input image file was supplied: we have to create one from
    // the sorted detection file.
    std::vector<ImageLog> img_log = {};
    double mjdnorm = 1.0;
    double mjdmean = detvec[0].MJD;
    int i = 0;  // Used later (Can this be removed?)
    int startind = 0;
    int endind = 0;
    for (i = 1; i < detvec.size(); i++) {
        double tdelt = detvec[i].MJD - detvec[i - 1].MJD;
        if (tdelt < config.imagetimetol / SOLARDAY &&
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
            ImageLog imlog = ImageLog(mjdmean, 0.0, 0.0, detvec[endind - 1].obscode, startind, endind);
            img_log.push_back(imlog);
            // Set up for the next image, starting with detvec[i].MJD;
            mjdmean = detvec[i].MJD;
            mjdnorm = 1.0;
            startind = i;
        }
    }
    // Account for the final image.
    if (isnormal(mjdnorm)) {
        endind = i;
        mjdmean /= mjdnorm;
        // Load it into the vector with mean MJD for all images,
        //  and increment image count.
        ImageLog imlog = ImageLog(mjdmean, 0.0, 0.0, detvec[endind - 1].obscode, startind, endind);
        img_log.push_back(imlog);
    }

    // We've now loaded the mean MJDs and the starting and ending
    // detection table indices for each image; it still remains to
    // get the mean RA and Dec.

    long detnum = detvec.size();
    long imnum = img_log.size();
    cout << img_log.size() << " unique images were identified.\n";
    cout << "Given our total of " << detvec.size() << " detections,\n";
    cout << "we have " << double(detvec.size()) / double(img_log.size())
         << " detections per image, on average\n";

    // Find the number of detections and the average RA, Dec on each image.
    // We perform the average after projection onto the unit circle, to
    // avoid wrapping issues.
    long detct = 0;
    long imct = 0;
    while (imct < imnum && detct < detnum) {
        long num_dets = 0;
        vector<Detection> imobs = {};
        Point3D p3avg = Point3D(0.0, 0.0, 0.0);
        while (detct < detnum && detvec[detct].MJD < img_log[imct].MJD + config.imagetimetol / SOLARDAY) {
            num_dets++;                                  // Keep count of detections on this image
            Detection detc = detvec[detct];              // Current detection entry
            imobs.push_back(detc);                       // Load vector of observations for this image
            Point3D p3 = celeproj01(detc.RA, detc.Dec);  // Project current detection
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
            int i = celedeproj01(p3avg, &img_log[imct].RA, &img_log[imct].Dec);
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
        /*
        cout << "Image " << imct << " of " << img_log.size() << ": " << num_dets << " = "
             << img_log[imct].endind - img_log[imct].startind;
        cout << " from " << img_log[imct].startind << " to " << img_log[imct].endind << " of "
             << detvec.size() << ".\n";
        */
        fflush(stdout);
        imct++;
    }
    return img_log;
}

std::tuple<std::vector<long double>, std::vector<Point3LD>, std::vector<Point3LD>> readEarthEphemerides(
    string earthfile
) {
    std::vector<long double> earthMJD = {};
    std::vector<Point3LD> earthpos = {};
    std::vector<Point3LD> earthvel = {};
    read_horizons_fileLD(earthfile, earthMJD, earthpos, earthvel);
    return std::make_tuple(earthMJD, earthpos, earthvel);
}

void computeHelioPositions(
    std::vector<Detection> &detvec,
    std::vector<ImageLog> const& img_log,
    std::vector<Observatory> const& observatory_list,
    std::vector<long double> const& EarthMJD,
    std::vector<Point3LD> const& Earthpos
) {
    double obslon, plxcos, plxsin;

    // Calculate observer's heliocentric position at the time of each image.
    std::vector<Point3LD> observer_heliopos = {};
    for (size_t imct = 0; imct < img_log.size(); imct++) {
        if (imct == 0 ||
            (imct > 0 && stringnmatch01(img_log[imct].obscode, img_log[imct - 1].obscode, 3) == 0)) {
            // Observatory has changed: get observatory coordinates for this image.
            int status = obscode_lookup(observatory_list, img_log[imct].obscode, obslon, plxcos, plxsin);

            // Figure out what to do with this...
            /*if (status > 0) {
                cerr << "ERROR: obscode_lookup failed for observatory code " << img_log[imct].obscode << "\n";
                return (3);
            }*/
        }
        Point3LD outpos = Point3LD(0, 0, 0);
        observer_barycoords01LD(img_log[imct].MJD, 5, obslon, plxcos, plxsin, EarthMJD, Earthpos, outpos);
        observer_heliopos.push_back(outpos);
    }
    // assign observer heliocentric position to each detection.
    for (size_t imct = 0; imct < img_log.size(); imct++) {
        for (int i = img_log[imct].startind; i < img_log[imct].endind; i++) {
            detvec[i].x = observer_heliopos[imct].x;
            detvec[i].y = observer_heliopos[imct].y;
            detvec[i].z = observer_heliopos[imct].z;
            // Check that the obscodes match between image and detection
            if (stringnmatch01(detvec[i].obscode, img_log[imct].obscode, 3) != 0) {
                cout << "ERROR: obscode mismatch (" << detvec[i].obscode << " vs. " << img_log[imct].obscode
                     << " between image " << imct << " and detection " << i << "\n";
                // return (4); // fix this later
            }
            detvec[i].indvec = {};  // Just making sure the index vectors are empty at this point.
        }
    }
}

std::tuple<std::vector<Detection>, std::vector<LongPair>> buildTracklets(
    MakeTrackletsConfig const& config,
    std::vector<Detection> &detvec,
    std::vector<ImageLog> &img_log
) {
    std::vector<LongPair> pairvec = {};
    std::vector<Detection> pairdets = {};
    LongPair onepair = LongPair(0, 0);
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
    XYIndex xyind = XYIndex(0.0, 0.0, 0);
    vector<XYIndex> axyvec = {};
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
            xyind = XYIndex(0.0, 0.0, 0);
            axyvec = {};
            dist = pa = 0.0;
            dettarg = 0;
            for (detct = img_log[imct].startind; detct < img_log[imct].endind; detct++) {
                distradec02(img_log[imct].RA, img_log[imct].Dec, detvec[detct].RA, detvec[detct].Dec, &dist,
                            &pa);
                xyind = XYIndex(dist * sin(pa / DEGPRAD), dist * cos(pa / DEGPRAD), detct);
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
                vector<XYIndex> bxyvec = {};
                // Project all detections on image B
                for (dettarg = img_log[imtarg].startind; dettarg < img_log[imtarg].endind; dettarg++) {
                    distradec02(img_log[imct].RA, img_log[imct].Dec, detvec[dettarg].RA, detvec[dettarg].Dec,
                                &dist, &pa);
                    xyind = XYIndex(dist * sin(pa / DEGPRAD), dist * cos(pa / DEGPRAD), dettarg);
                    bxyvec.push_back(xyind);
                }
                // Create k-d tree of detections on image B (imtarg).
                int dim = 1;
                XYIndex xyi = bxyvec[0];
                KDPoint root = KDPoint(xyi, -1, -1, dim);
                KDPoint lp1 = KDPoint(xyi, -1, -1, dim);
                KDPoint rp1 = KDPoint(xyi, -1, -1, dim);
                KDPoint kdtest = KDPoint(xyi, -1, -1, dim);
                vector<KDPoint> kdvec = {};
                long medpt;
                medpt = medindex(bxyvec, dim);
                root = KDPoint(bxyvec[medpt], -1, -1, 1);
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
                            onepair = LongPair(detvec[axyvec[detct].index].index,
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

    return std::make_tuple(pairdets, pairvec);
}

void refineTracklets(
    MakeTrackletsConfig const& config,
    std::vector<Detection> &pairdets,
    string outpairfile
) {
    ofstream outstream1;
    vector<LongIndex> pair_partner_num = {};
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
        LongIndex ppn = LongIndex(pairdets[i].indvec.size(), i);
        pair_partner_num.push_back(ppn);
    }
    // Sort the new vector by number of pair-partners
    sort(pair_partner_num.begin(), pair_partner_num.end(), compareLongIndex);

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
            std::vector<XYIndex> axyvec = {};
            std::vector<Detection> ppset = {};
            for (size_t j = 0; j < pairdets[pdct].indvec.size(); j++) {
                long detct = pairdets[pdct].indvec[j];
                if (pairdets[detct].indvec.size() > 0) {
                    // Detection detct hasn't already been allocated to a tracklet.
                    distradec02(pairdets[pdct].RA, pairdets[pdct].Dec, pairdets[detct].RA,
                                pairdets[detct].Dec, &dist, &pa);
                    dist *= 3600.0L;  // Convert distance from degrees to arcsec.
                    XYIndex xyind = XYIndex(dist * sin(pa / DEGPRAD), dist * cos(pa / DEGPRAD), detct);
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
                std::vector<Point3DIndex> track_mrdi_vec = {};  // We need this vector purely so we can do a time-sort.
                                      // mrdi stands for MJD, RA, Dec, index
                // Load the reference point
                Point3DIndex p3di = Point3DIndex(0.0l, 0.0l, 0.0l, pdct);
                track_mrdi_vec.push_back(p3di);
                // Load anchor point corresponding to biggest_tracklet
                p3di = Point3DIndex(ppset[biggest_tracklet].MJD - pairdets[pdct].MJD,
                                     axyvec[biggest_tracklet].x, axyvec[biggest_tracklet].y,
                                     axyvec[biggest_tracklet].index);
                track_mrdi_vec.push_back(p3di);
                // Load the other points
                for (size_t j = 0; j < ppset[biggest_tracklet].indvec.size(); j++) {
                    p3di = Point3DIndex(ppset[ppset[biggest_tracklet].indvec[j]].MJD - pairdets[pdct].MJD,
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
                sort(track_mrdi_vec.begin(), track_mrdi_vec.end(), comparePoint3DIndexX);
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