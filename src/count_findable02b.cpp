// April 05, 2022: count_findable02b.cpp
// Hopefully a more polished and production-like version of
// count_findable01c.cpp. The important changes include calculating
// all the metrics even for objects ultimately not classified as
// findable.
//
// April 01, 2022: count_findable01c.cpp:
// Like count_findable01b.cpp, but also calculates the
// minimum, median, and maximum tracklet length for each object.
//
// February 23, 2022: count_findable01b.cpp:
// Like count_findable01a.cpp, but calculates more metrics,
// aimed at determining why some time-periods give heliolinc
// so much more trouble than others.
//
// Description of ancestor program count_findable01a.cpp
// Read a csv file of detections of astronomical objects, determine
// which ones were findable based on criteria involving the presence
// of two or more detections per night over at least NIGHT_MIN
// different nights. 

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define NUMPOS 3

#define IDCOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define MINOBSINTERVAL 1.0 // Minimum time-between-images in seconds
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME 1.5 // Default max inter-image time interval
                    // for tracklets, in hours (will be converted
                    // to days before use).
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)
#define MINPOINTS 6 // minimum number of points for a valid discovery.

 #define MINSPAN 1.0 // Temporal span must be at least this large (in days) for a bona fide cluster
#define MINDAYSTEPS 2 // A bona fide cluster must have at least this many intra-point
                      // time intervals greater than INTRANIGHTSTEP days.
#define INTRANIGHTSTEP 0.3 // Minimum interval in days between successive points
                           // in a tracklet, to enable them to be counted as being
                           // on separate nights.
#define MAXVELCALCSTEP 1.5 // Maximum interval in days to be used for estimating
                           // heliocentric velocity vector.
#define MINTRACKTIME 300.0 // Minimum tracklet timespan in seconds.
#define COLS_TO_READ 7

static void show_usage()
{
  cerr << "Usage: count_finddable02b -dets detfile -colformat column format file -maxtime max inter-image time interval (hr)/ \n";
  cerr << "-maxvel maximum angular velocity (deg/day) -mintracktime min tracklet time span (sec)\n";
  cerr << "-outfile output file\n";
  cerr << "\nor, at minimum\n\n";
  cerr << "count_findable02a -dets detfile -colformat column format file -outfile output file\n";
  cerr << "Note well that the minimum invocation will leave things\n";
  cerr << "set to defaults that may not be what you want.\n";
}
    
int main(int argc, char *argv[])
{
  det_OC_index o1 = det_OC_index(0L,0L,0L,0L,0L,0L,"null",0);
  vector <det_OC_index> detvec = {};
  vector <det_OC_index> trackvec ={};
  vector <det_OC_index> trackletvec ={};
  vector <long double> trackletmjd={};
  vector <long double> long_trackletmjd={};
  vector <long double> trackletarc={};
  vector <long double> tracklet_timespan={};
  long double minarc=0.0l;
  long double medarc=0.0l;
  long double maxarc=0.0l;
  long double min_timespan=0.0l;
  long double med_timespan=0.0l;
  long double max_timespan=0.0l;
  int validity = 0;
  double mjdmean = 0;
  string idstring;
  string lnfromfile;
  string stest;
  int i = 0;
  int j = 0;
  long detct=0;
  int reachedeof = 0;
  char c='0';
  long double MJD,RA,Dec;
  MJD = RA = Dec = 0.0L;
  double maxvel = MAXVEL; // Max angular velocity in deg/day
  double maxtime = MAXTIME; // Max time interval a tracklet could span,
                            // in days.
  double maxdist = MAXVEL*MAXTIME/24.0; // Max angular distance a tracklet
                                   // could span, in degrees.
  string indetfile;
  string outfile;
  string colformatfile;
  long double heliox = 0.0L;
  long double helioy = 0.0L;
  long double helioz = 0.0L;
  long lct=000000000000000;
  long double timespan = 0.0L;
  long double currentstep = 0.0L;
  double angdist=0l;
  double angspeed=0l;
  int tracknightcount=0;
  int long_tracknightcount=0;
  int trackletcount=0;
  int tracklet_firstpoint=0;
  long double timestep = 0L;
  long double maxtimestep = 0L;
  point3LD targvel = point3LD(0,0,0);
  point3LD targpos1 = point3LD(0,0,0);
  point3LD targpos2 = point3LD(0,0,0);
  int vp1=0;
  int vp2=0;
  long double r0,v0,E,lscalar,a,e,r1,r2,heliovel;
  r0=v0=E=lscalar=a=e=r1=r2=heliovel=0.0L;
  point3LD lvec = point3LD(0L,0L,0L);
  long double MGsun = GMSUN_KM3_SEC2;
  ofstream outstream1;
  double mintracktime = MINTRACKTIME; // Minimum tracklet time span in seconds.
  int idcol = IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int magcol, bandcol, obscodecol;
  magcol = bandcol = obscodecol = 0;
  ifstream instream1;
  int colreadct=0;
  
  if(argc!=13 && argc!=11 &&argc!=9 && argc!=5)
    {
      show_usage();
      return(1);
    }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-d" || string(argv[i]) == "-dets" || string(argv[i]) == "-det" || string(argv[i]) == "--dets" || string(argv[i]) == "--det" || string(argv[i]) == "--detection" || string(argv[i]) == "--detections") {
      if(i+1 < argc) {
	//There is still something to read;
	indetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format"  || string(argv[i]) == "-col" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" || string(argv[i]) == "--cformat") {
      if(i+1 < argc) {
	//There is still something to read;
	colformatfile=argv[++i];
	i++;
      }
      else {
	cerr << "Column format file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxtime") {
      if(i+1 < argc) {
	//There is still something to read;
        maxtime=stod(argv[++i]);
	i++;
	if(!isnormal(maxtime) || maxtime<=0.0) {
	  cerr << "Error: invalid maximum inter-image time interval\n";
	  cerr << "(" << maxtime << " hr) supplied: must be strictly positive.\n";
	  return(2);
	}      
      } else {
	cerr << "Output maximum inter-image time interval\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxvel") {
      if(i+1 < argc) {
	//There is still something to read;
        maxvel=stod(argv[++i]);
	i++;
	if(!isnormal(maxvel) || maxvel<=0.0) {
	  cerr << "Error: invalid maximum angular velocity\n";
	  cerr << "(" << maxvel << "deg/day) supplied: must be strictly positive.\n";
	  return(2);
	}
      }
      else {
	cerr << "Output maximum angular velocity\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-mintracktime") {
      if(i+1 < argc) {
	//There is still something to read;
        mintracktime=stod(argv[++i]);
	i++;
	if(!isnormal(mintracktime) && mintracktime!=0.0) {
	  cerr << "Error: invalid minimum tracklet timespan\n";
	  cerr << "(" << mintracktime << "sec) supplied.\n";
	  return(2);
	}
      }
      else {
	cerr << "Output minimum tracklet timespan\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile" || string(argv[i]) == "-of" || string(argv[i]) == "-out" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outputfile" || string(argv[i]) == "--compfile" || string(argv[i]) == "--outputcompfile") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else {
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }

  if(indetfile.size()<=0) {
    cerr << "Please supply an input detection file:\n\n";
    show_usage();
    return(1);
  }
  if(outfile.size()<=0) {
    cerr << "Please supply the name of an ouput file:\n\n";
    show_usage();
    return(1);
  }

  cout << "indet file " << indetfile << "\n";
  cout << "max time interval " << maxtime << "\n";
  cout << "maxvel " << maxvel << "\n";
  cout << "minimum tracklet timespan " << mintracktime << " seconds\n";
  maxtime/=24.0; /*Unit conversion from hours to days*/
  maxdist = maxtime*maxvel;
  cout << "output file: " << outfile << "\n";
  
  // Read the column formatting file, if any
  if(colformatfile.size()>0)
    {
      instream1.open(colformatfile);
      if(!instream1)  {
	cerr << "ERROR: unable to open input file " << colformatfile << "\n";
	return(1);
      }
      colreadct=0;
      while(!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct<COLS_TO_READ) {
	instream1 >> stest;
	if (stest == "IDCOL") {
	  instream1 >> idcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "MJDCOL") {
	  instream1 >> mjdcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "RACOL") {
	  instream1 >> racol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "DECCOL") {
	  instream1 >> deccol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "MAGCOL") {
	  instream1 >> magcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "BANDCOL") {
	  instream1 >> bandcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "OBSCODECOL") {
	  instream1 >> obscodecol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else {
	  cout << "WARNING: unrecognized string " << stest << " read from column formatting file\n";
	}
      }
      instream1.close();
      if(colreadct<COLS_TO_READ) {
	cout << "WARNING: only " << colreadct << " column specifications, of " << COLS_TO_READ << " expected, were read from column format file " << colformatfile << ".\n";
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

  instream1.open(indetfile);
  if(!instream1) {
    cerr << "can't open input file " << indetfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  lct++;
  cout << lnfromfile << "\n";
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    lct++;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    i=0;
    j = 0;
    c='0';
    while(i<long(lnfromfile.size()) && reachedeof == 0) {
      stest=""; // Note we used to declare stest here. Maybe should revert.
      c='0';
      while(i<long(lnfromfile.size()) && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==idcol) idstring=stest;
      if(j==mjdcol) MJD=stold(stest);
      else if(j==racol) RA=stold(stest);
      else if(j==deccol) Dec=stold(stest);
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0) {
      // cout<<"MJD, RA, Dec: " << MJD-floor(MJD) << " " << RA << " " << Dec << "\n";
      o1=det_OC_index(MJD,RA,Dec,heliox,helioy,helioz,idstring,-lct);
      detvec.push_back(o1);
    }
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  // Sort the detection vector by the id strings.
  sort(detvec.begin(), detvec.end(), stringsort_det_OC_index());

  //  for(detct=0;detct<long(detvec.size());detct++) {
  //  cout << detct << " " << detvec[detct].MJD << " " << detvec[detct].RA << " " << detvec[detct].Dec << " " << detvec[detct].idstring << "\n";
  // }
  // Pull out the sets of detections corresponding to
  // each simulated object.
  trackvec ={};
  trackvec.push_back(detvec[0]);
  idstring = detvec[0].idstring;
  detct=1;
  outstream1.open(outfile,ios_base::out);
  outstream1  << "idstring validity npts timespan tracknightcount r0 heliovel a e1 q minarc medarc maxarc min_timespan med_timespan max_timespan\n";
  while(detct<long(detvec.size())) {
    while(detct<long(detvec.size()) && detvec[detct].idstring == idstring) {
      // We are still on the same simulated object.
      trackvec.push_back(detvec[detct]);
      detct++;
    }
    if(detct<=long(detvec.size())) {
      // We must have just finished loading a set of points
      // corresponding to one simulated object.
      // Sort the track vector.
      sort(trackvec.begin(), trackvec.end(), early_det_OC_index());
      //for(i=0; i<trackvec.size(); i++) {
      //  cout << fixed << setprecision(6) << "trackpoint " << i << ": " << trackvec[i].MJD << " "<< trackvec[i].RA << " "<< trackvec[i].Dec << "\n";
      //}
      timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
      // Initialize some stuff to dummy values
      medarc = minarc = maxarc = med_timespan = min_timespan = max_timespan = 0.0L;
      r0 = -999.99L;
      heliovel = -999.99L;
      a = -999.999999L*1.495978700e8L;
      e = 0.0L;
      tracknightcount = -1;
      // Count the number of qualifying tracklets.
      trackletmjd={};
      trackletarc={};
      tracklet_timespan={};
      trackletvec={};
      trackletcount=0;
      tracklet_firstpoint=0;
      trackletvec.push_back(trackvec[tracklet_firstpoint]);
      // cout << "Preparing to probe track for object " << idstring << " with " << trackvec.size() << " points\n";
      for(i=1; i<long(trackvec.size()); i++) {
	currentstep = trackvec[i].MJD - trackvec[tracklet_firstpoint].MJD;
	angdist = distradec01(trackvec[tracklet_firstpoint].RA, trackvec[tracklet_firstpoint].Dec, trackvec[i].RA, trackvec[i].Dec); // Distance in degrees
	angspeed = angdist/currentstep; // Angular velocity in deg/day
	// cout << "firstpoint=" << tracklet_firstpoint << ", i=" << i << ", step dist speed = " << currentstep*86400.0 << " " << angdist*3600.0 << " " << angspeed << "\n";
	if(currentstep>MINOBSINTERVAL/SOLARDAY && currentstep<=maxtime && angdist<=maxdist && angspeed<=maxvel) {
	  // Detection i is part of an valid tracket starting with
	  // tracklet_firstpoint, safely contained within a single night,
	  // and meeting requirements on maximum separation and angular velocity.
	  trackletvec.push_back(trackvec[i]);
	  // cout << "Object " << idstring << ", tracklet " << trackletcount << ", adding point " << trackletvec.size()-1 << ": " << trackletvec[trackletvec.size()-1].MJD << " " <<  trackletvec[trackletvec.size()-1].RA << " " <<  trackletvec[trackletvec.size()-1].Dec << "\n"; 
	} else {
	  // Either we just finished a tracklet, or detection tracklet_firstpoint
	  // was in fact a singleton.
	  if(trackletvec.size()>1) {
	    // It was a valid tracklet.
	    //for(j=0;j<trackletvec.size();j++) {
	    //cout << "trackletpoint " << j << ": " << trackletvec[j].MJD << " "<< trackletvec[j].RA << " "<< trackletvec[j].Dec << "\n";
	    //}
	    // Calculate the mean MJD.
	    mjdmean=0.0;
	    for(j=0;j<long(trackletvec.size());j++) mjdmean += trackletvec[j].MJD;
	    mjdmean/=double(trackletvec.size());
	    trackletmjd.push_back(mjdmean);
	    // Calculate the angular arc
	    angdist = distradec01(trackletvec[0].RA, trackletvec[0].Dec, trackletvec[trackletvec.size()-1].RA, trackletvec[trackletvec.size()-1].Dec);
	    trackletarc.push_back(angdist*3600.0); // Angular arc spanned by tracklet, in arcsec.
	    tracklet_timespan.push_back((trackletvec[trackletvec.size()-1].MJD - trackletvec[0].MJD)*SOLARDAY); // Time interval spanned by the tracklet, in seconds
	    // Increment trackletcount
	    trackletcount+=1;
	    // Set up for a new tracklet
	    trackletvec={};
	    tracklet_firstpoint=i;
	    trackletvec.push_back(trackvec[tracklet_firstpoint]);
	  } else {
	    // It was a singleton. Set up for a new tracklet.
	    trackletvec={};
	    tracklet_firstpoint=i;
	    trackletvec.push_back(trackvec[tracklet_firstpoint]);
	  }
	}
      }
      // Done with looping over all points in the track.
      // Handle a possible final tracklet that may have been left hanging.
      if(trackletvec.size()>1) {
	// It was a valid tracklet.
	// Calculate the mean MJD.
	mjdmean=0.0;
	for(j=0;j<long(trackletvec.size());j++) mjdmean += trackletvec[j].MJD;
	mjdmean/=double(trackletvec.size());
	trackletmjd.push_back(mjdmean);
	// Calculate the angular arc
	angdist = distradec01(trackletvec[0].RA, trackletvec[0].Dec, trackletvec[trackletvec.size()-1].RA, trackletvec[trackletvec.size()-1].Dec);
	trackletarc.push_back(angdist*3600.0); // Angular arc spanned by tracklet, in arcsec.
	tracklet_timespan.push_back((trackletvec[trackletvec.size()-1].MJD - trackletvec[0].MJD)*SOLARDAY); // Time interval spanned by the tracklet, in seconds
	// Increment trackletcount
	trackletcount+=1;
      }
      if(trackletarc.size()>0) {
	ldmedian_minmax(trackletarc,medarc,minarc,maxarc);
      }
      if(tracklet_timespan.size()>0) {
	ldmedian_minmax(tracklet_timespan,med_timespan,min_timespan,max_timespan);
      }
      // Create a version of the trackletmjd vector that only
      // contains tracklets with timespan at least equal to mintracktime
      long_trackletmjd={};
      for(i=0;i<long(trackletmjd.size()); i++) {
	if(tracklet_timespan[i]>=mintracktime) {
	  // This tracklet is long enough: add it to the long tracklet vector
	  long_trackletmjd.push_back(trackletmjd[i]);
	}
      }
      tracknightcount=0;
      // See if the tracklets fall on different nights.
      for(i=1;i<long(trackletmjd.size()); i++) {
	currentstep = trackletmjd[i] - trackletmjd[i-1];
	if(currentstep >= INTRANIGHTSTEP) {
	  // The time interval between tracklets i-1 and i is long
	  // enough to suggest they occurred on different nights.
	  tracknightcount++;
	}
      }
      long_tracknightcount=0;
      // See if the long tracklets fall on different nights.
      for(i=1;i<long(long_trackletmjd.size()); i++) {
	currentstep = long_trackletmjd[i] - long_trackletmjd[i-1];
	if(currentstep >= INTRANIGHTSTEP) {
	  // The time interval between tracklets i-1 and i is long
	  // enough to suggest they occurred on different nights.
	  long_tracknightcount++;
	}
      }
      if(trackvec.size()==1) {
	// Cannot calculate a velocity, but at least get heliocentric distance.
	vp1=0;
	targpos1.x=trackvec[vp1].x;
	targpos1.y=trackvec[vp1].y;
	targpos1.z=trackvec[vp1].z;
	r0 = sqrt(dotprod3LD(targpos1,targpos1));
      } else if (trackvec.size()>1) {
	// Calculate the orbital semimajor axis and eccentricity
	// Find two points with suitable time-separation for
	// this purpose.
	vp1=vp2=0;
	maxtimestep=0L;
	for(i=0; i<long(trackvec.size())-1; i++) {
	  for(j=i+1;  j<long(trackvec.size()); j++) {
	    timestep = trackvec[j].MJD - trackvec[i].MJD;
	    if(timestep<=MAXVELCALCSTEP && timestep>maxtimestep) {
	      maxtimestep = timestep;
	      vp1=i;
	      vp2=j;
	    }
	  }
	}
	targpos1.x=trackvec[vp1].x;
	targpos1.y=trackvec[vp1].y;
	targpos1.z=trackvec[vp1].z;
	targpos2.x=trackvec[vp2].x;
	targpos2.y=trackvec[vp2].y;
	targpos2.z=trackvec[vp2].z;
	r1 = sqrt(dotprod3LD(targpos1,targpos1));
	r2 = sqrt(dotprod3LD(targpos2,targpos2));
	heliovel = (r2-r1)/maxtimestep/SOLARDAY;  // units should be km/sec.
	// Calculate heliocentric velocity vector.
	targvel.x = (targpos2.x-targpos1.x)/maxtimestep/SOLARDAY; // units should be km/sec.
	targvel.y = (targpos2.y-targpos1.y)/maxtimestep/SOLARDAY;
	targvel.z = (targpos2.z-targpos1.z)/maxtimestep/SOLARDAY;
	// average position vector
	targpos1.x = 0.5L*(targpos1.x + targpos2.x);
	targpos1.y = 0.5L*(targpos1.y + targpos2.y);
	targpos1.z = 0.5L*(targpos1.z + targpos2.z);

	r0 = sqrt(dotprod3LD(targpos1,targpos1));
	v0 = sqrt(dotprod3LD(targvel,targvel));
  
	// Calculate specific energy and angular momentum
	E = 0.5L*v0*v0 - MGsun/r0;
	lvec = crossprod3LD(targpos1,targvel);
	lscalar = sqrt(dotprod3LD(lvec,lvec));
		 
	// Calculate a and e: orbital semimajor axis and eccentricity.
	a = -MGsun*0.5L/E;
	e = sqrt(1.0L + 2.0L*E*lscalar*lscalar/MGsun/MGsun);
      }
      if(timespan >= MINSPAN && trackvec.size()>=MINPOINTS && long_tracknightcount >= MINDAYSTEPS) {
	// This simulated object had sufficient observations
	// to enable qualifying discovery.
	validity = 1;
      } else validity = 0;

      // Write stats to output file.
      // string ID, VALIDITY, npts, timespan, tracknightcount, heliodist (km), heliovel(km/sec), a (AU), e, peridist (AU)
      outstream1  << fixed << setprecision(6) << idstring << " " << validity << " " << trackvec.size() << " " << timespan << " " << tracknightcount << " ";
      outstream1  << fixed << setprecision(3) << r0 << " " << heliovel << " ";
      outstream1  << fixed << setprecision(6) << a/1.495978700e8L << " " << e << " " << (1.0L-e)*a/1.495978700e8L << " ";
      outstream1  << fixed << setprecision(3) << minarc << " " << medarc << " " << maxarc << " " << min_timespan << " " << med_timespan << " " << max_timespan << "\n";
      // Finished dealing with possible completed object.
      // Set up for the next object
      if(detct<long(detvec.size())) idstring = detvec[detct].idstring;
      trackvec ={};
      // Close if-statement checking if we just finished an object.
    }
  }
  return(0);
}
