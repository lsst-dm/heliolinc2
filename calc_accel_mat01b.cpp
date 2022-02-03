// January 17, 2022: calc_accel_mat01b.cpp:
// Like calc_accel_mat01a.cpp, but instead of printing a specific
// set of acceleration quantiles, prints out specified min and max
// quantiles (e.g. 5th to 95th percentiles) plus evenly sampled
// steps between them. Note that the steps are even in the actual
// acceleration, not even quantiles of the distribution (which is
// what calc_accel_mat01a gives). The units of the steps are GMsun/r^2.
// The even sampling in physical accleration rather than in quantiles
// of the distribution ensures we can configure the sampling of
// distance, velocity, acceleration quantiles so that we capture
// even rare objects with unusual orbits, rather than being sensitive
// only to objects in more common orbits.
//
// Description of ancestor program calc_accel_mat01a.cpp:
// Probe the search space of the heliolinc algorithm (that is,
// heliocentric distance and radial vecloity), and for each grid
// point in this space, calculate quantiles in acceleration for
// asteroids that appear at that point (that is, will at some
// point in their orbits have the corresponding heliocentric
// distance and radial velocity). Write a file providing the
// distance/velocity grid points and their associated heliocentric
// radial accelerations, for use as input in heliolinc runs using
// projectpairs04c.cpp.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: calc_accel_mat01b -orbfile orbit_file -outfile output_file -diststep dist0 diststep distnum -velstep vel0 velstep velnum -accelstep minquant maxquant accelstep -timestep timestep timenum\n";
}
    
int main(int argc, char *argv[])
{
  string orbfile,outfile,lnfromfile;
  ifstream instream1;
  ofstream outstream1;
  vector <heliogridpoint> gridvec;
  heliogridpoint hgp = heliogridpoint(0L,0L,{});
  long double diststart = 1.0L;
  long double velstart = -0.02L;
  long double diststep = 0.1L;
  long double velstep = -0.001L;
  int distnum = 91;
  int distct = 0;
  int velnum = 41;
  int velct = 0;
  int i=0;
  int j=0;
  long double timestep = 1.0; // Units are days.
  int tnum = 3650;
  int tct = 0;
  double minquant=0.0;
  double maxquant=1.0;
  double accelstep = 0.1;
  int accelstepnum=0;
  double minaccel=0.0l;
  double maxaccel=0.0l;
  vector <long double> semimaj;
  vector <long double> eccen;
  long double a=0L;
  long double e=0L;
  int orbitnum = 0;
  int orbitct = 0;
  long double MGsun = GMSUN_KM3_SEC2;
  long double E = 0.0L;
  long double lscalar = 0.0L;
  long double omega,Period,orbtime,tomega,psi,cospsi;
  omega = Period = orbtime = tomega = psi = cospsi = 0.0L;
  long double costheta,theta,r1,v1,vtan,vrad,ldc,lvc;
  costheta = theta = r1 = v1 = vtan = vrad = ldc = lvc = 0.0L;
  vector <long double> heliorad;
  vector <long double> heliovel;
  long double accel=0L;
  long double orboff=0L;
  vector <long double> accelvec;
  

  if(argc!=20) {
    show_usage();
    return(1);
  }
  
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-orbfile" || string(argv[i]) == "-orbitfile" || string(argv[i]) == "-ofile" || string(argv[i]) == "--orbfile" || string(argv[i]) == "--orbitfile" || string(argv[i]) == "--orbitf") {
      if(i+1 < argc) {
	//There is still something to read;
	orbfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input orbit file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile" || string(argv[i]) == "-of" || string(argv[i]) == "-outf" || string(argv[i]) == "--outputfile" || string(argv[i]) == "--ofile" || string(argv[i]) == "--fileout") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "output filename keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-diststep" || string(argv[i]) == "-ds" || string(argv[i]) == "-dstep" || string(argv[i]) == "-dist" || string(argv[i]) == "--diststep" || string(argv[i]) == "--distancestep" || string(argv[i]) == "--distance") {
      if(i+1 < argc) {
	//There is still something to read;
	diststart=stold(argv[++i]);
      }
      else {
	cerr << "Distance step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	diststep=stold(argv[++i]);
      }
      else {
	cerr << "Distance step keyword supplied with only one of the three required arguments\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	distnum=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Distance step keyword supplied with only two of the three required arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-velstep" || string(argv[i]) == "-vs" || string(argv[i]) == "-vstep" || string(argv[i]) == "-vel" || string(argv[i]) == "--velstep" || string(argv[i]) == "--velocitystep" || string(argv[i]) == "--velocity") {
      if(i+1 < argc) {
	//There is still something to read;
	velstart=stold(argv[++i]);
      }
      else {
	cerr << "Velocity step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	velstep=stold(argv[++i]);
      }
      else {
	cerr << "Velocity step keyword supplied with only one of the three required arguments\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	velnum=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Velocity step keyword supplied with only two of the three required arguments\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-accstep" || string(argv[i]) == "-acc"  || string(argv[i]) == "-accelstep" || string(argv[i]) == "-as" || string(argv[i]) == "--acceleration" || string(argv[i]) == "--accelstep" || string(argv[i]) == "--accelerationstep") {
      if(i+1 < argc) {
	//There is still something to read;
	minquant=stold(argv[++i]);
      }
      else {
	cerr << "output acceleration step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	maxquant=stold(argv[++i]);
      }
      else {
	cerr << "Acceleration step keyword supplied with only one of the three required arguments\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	accelstep=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Acceleration step keyword supplied with only two of the three required arguments\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-timestep" || string(argv[i]) == "-ts" || string(argv[i]) == "-tstep" || string(argv[i]) == "-time" || string(argv[i]) == "--timestep" || string(argv[i]) == "--time" || string(argv[i]) == "--tempstep") {
      if(i+1 < argc) {
	//There is still something to read;
	timestep=stold(argv[++i]);
      }
      else {
	cerr << "Time step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
	tnum=stoi(argv[++i]);
      }
      else {
	cerr << "Time step keyword supplied with only one of the two required arguments\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
    
  cout << "input orbint file " << orbfile << "\n";
  cout << "output file name " << outfile << "\n";
  cout << "distance sampling parameters " << diststart << " " << diststep << " " << distnum << "\n";
  cout << "velocity sampling parameters " << velstart << " " << velstep << " " << velnum << "\n";
  cout << "min acceleration quantile" << minquant << "\n";
  cout << "max acceleration quantile" << maxquant << "\n";
  cout << "acceleration step size" << accelstep << "\n";
  cout << "time step parameters " << timestep << " " << tnum << "\n";

  // Allocate grid vector with empty acceleration vectors
  gridvec={};
  for(distct=0; distct<distnum; distct++) {
    for(velct=0; velct<velnum; velct++) {
      hgp = heliogridpoint(diststart+diststep*(long double)distct, velstart+velstep*(long double)velct,{});
      gridvec.push_back(hgp);
      cout << "Loaded point " << distct*velnum+velct << " = " << gridvec.size()-1 << "\n";
    }
  }

  // Read input orbit file
  instream1.open(orbfile,ios_base::in);
  if(!instream1) {
    cerr << "can't open input file " << orbfile << "\n";
    return(1);
  }
  semimaj = {};
  eccen = {};
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    instream1 >> a >> e;
    if (!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      semimaj.push_back(a);
      eccen.push_back(e);
    }
    getline(instream1,lnfromfile);
  }
  instream1.close();
  orbitnum = semimaj.size();
  if(orbitnum != eccen.size()) {
    cerr << "ERROR: semimajor axis and eccentricity vectors do not have the same length\n";
    return(1);
  }
  cout << "Read " << orbitnum << " orbits from " << orbfile << "\n";
  
  outstream1.open(outfile,ios_base::out);

  for(orbitct=0; orbitct<orbitnum; orbitct++) {
    a = semimaj[orbitct]*AU_KM; // Convert from AU to km.
    e = eccen[orbitct]; // Just for convenience
    // Calculate specific energy E from semimajor axis a
    E = -MGsun*0.5L/a;
    // Calculate specific angular momentum from from semimajor axis a and eccentricity e
    lscalar = sqrt(MGsun*a*(1.0L-e*e));
    omega = sqrt(MGsun/(a*a*a)); // Mean angular frequency of orbit (rad/sec)
    Period = 2.0L*M_PI/omega;  // Orbital period in seconds.
    //cout << "a, e, E, ls, omega, Period: " << a << " " << e << " " << E << " " << lscalar << " " << omega << " " << Period << "\n";
    
    cout << "Orbit " << orbitct << ": a, e, Period: " << a/AU_KM << " " << e << " " << Period/SOLARDAY/365.2425 << "\n";
    
    // Calculate heliocentric distance and velocity for this orbit,
    // starting from perihelion at t=0 (using a random offset);
    heliorad = heliovel = {};
    orboff = (double(rand())/double(RAND_MAX))*Period;
    for(tct=0; tct<tnum; tct++) {
      orbtime = tct*timestep*SOLARDAY; // units of seconds
      // Add a random offset with a range of one period,
      // to avoid biasing statistics toward perihelion.
      orbtime += orboff;
      // Unwrap periods if necessary
      while(orbtime>=Period) orbtime -= Period;
      // Get input for Kepler Equation.
      tomega = orbtime*omega;
      // Solve Kepler Equation
      psi = kep_transcendental(tomega,e,KEPTRANSTOL);
      // Convert psi to orbital angle theta
      cospsi = cos(psi);
      if(1.0L - e*cospsi != 0.0L) {
	costheta = (cospsi - e)/(1.0L - e*cospsi);
	if(costheta >= -1.0L && costheta <= 1.0L) theta = acos(costheta);
	else if(costheta < -1.0L) {
	  //cout << "Warning: costheta = " << costheta << "\n";
	  theta = M_PI;
	} else {
	  //cout << "Warning: costheta = " << costheta << "\n";
	  theta = 0.0L;
	}
	if(psi>M_PI && theta<=M_PI) theta = 2.0L*M_PI - theta;
      } else {
	cerr << "Warning: e*cos(psi) = " << e*cospsi << " so 1 - e*cos(psi) = " << 1.0L - e*cospsi << "\n";
	theta = 0.0L;
      }
      while(theta<0.0L) theta += 2.0L*M_PI;
      while(theta>=2.0L*M_PI) theta -= 2.0L*M_PI;
      // Calculate r(t) from psi(t)
      r1 = a*(1.0L - e*cospsi);
      // Calculate v1 from r1 and the known energy
      v1 = sqrt((E +  MGsun/r1)*2.0L);
      vtan = lscalar/r1;
      if(vtan<=v1) vrad = sqrt(v1*v1 - vtan*vtan); // Should always apply, but check anyway...
      else vrad = 0.0L; // ...otherwise roundoff error could give us a NAN sometimes.
      if(psi>M_PI) vrad = -vrad; //Headed back toward perihelion.
      // Load calculated distance and velocity into vectors
      heliorad.push_back(r1);
      heliovel.push_back(vrad);
      //cout << "orbtime, tomega, psi, theta, r1, v1, vtan, vrad: " << orbtime << " " << tomega << " " << psi << " " << theta << " " << r1 << " " << v1 << " " << vtan << " " << vrad << "\n";
    }
    // Loop over vectors, calculating accelerations.
    for(tct=1; tct<tnum; tct++) {
      {
	// Bug-checking sanity output:
	//cout << "Vrad (km/sec) = " << heliovel[tct] << " = " << (heliorad[tct] - heliorad[tct-1])/(timestep*SOLARDAY) << "\n";
	// Get heliocentric distance in AU
	r1 = heliorad[tct]/AU_KM;
	// Get heliocentric velocity in AU/day
	vrad = heliovel[tct]*SOLARDAY/AU_KM;
	// Find which gridpoint we're looking at.
	ldc = (r1 - diststart+0.5*diststep)/diststep;
	lvc = (vrad - velstart+0.5*velstep)/velstep;
	distct = ldc;
	velct = lvc;
	if(distct>=0 && distct<distnum && velct>=0 && velct<velnum)
	  {
	    //cout << "dist = " << r1 << ", vel = " << vrad << ", matching gridpoint " << gridvec[distct*velnum+velct].dist << ", " << gridvec[distct*velnum+velct].vel << "\n";
	    // Calculate acceleration
	    accel = (heliovel[tct] - heliovel[tct-1])/(timestep*SOLARDAY); // Units are km/sec^2
	    accel /= -MGsun/(heliorad[tct]*heliorad[tct]); // Now it's relative to MGsun/r^2;
	    gridvec[distct*velnum+velct].acc.push_back(accel);
	    //cout << "Accel = " << accel << "\n";
	  }
      }	    
    }
  }

  // Loop over gridvec and output acceleration stats
  outstream1 << "#r(AU) rdot(AU/day) norm mean_accel mean_eccen mean_semimajor mean_perihelion\n";
  for(distct=0; distct<distnum; distct++) {
    for(velct=0; velct<velnum; velct++) {
      if(gridvec[distct*velnum+velct].acc.size() > 0) {
	// Load vector of accelerations
	accelvec={};
	for(i=0; i<gridvec[distct*velnum+velct].acc.size(); i++) {
	  accelvec.push_back(gridvec[distct*velnum+velct].acc[i]);
	}
	// Sort vector of accelerations.
	sort(accelvec.begin(), accelvec.end());
	// Extract requested quantiles;
	j=minquant*accelvec.size();
	minaccel = accelvec[j];
	j=maxquant*accelvec.size();
	maxaccel = accelvec[j];
	// Round to nearest step.
	minaccel = accelstep*round(minaccel/accelstep);
	maxaccel = accelstep*round(maxaccel/accelstep);
	accelstepnum = round((maxaccel-minaccel)/accelstep);
	for(j=0;j<=accelstepnum;j++) {
	  outstream1 << fixed << setprecision(6) << gridvec[distct*velnum+velct].dist << " " << gridvec[distct*velnum+velct].vel << " " << accelvec.size() <<  " " << minaccel+double(j)*accelstep << " -99.9 -99.9 -99.9\n";
	}
      }
    }
  }
  return(0);
}
