// September 06, 2022: sunelongcull01a.cpp:
// Given an input csv file with columns providing (at least) MJD, RA, and Dec,
// output a copy of the file, preserving a one-line header, that includes only
// detections within a specified range of solar elongation. For example,
// only solar elongations less than 90 degrees.
//
// Accepts an input column formatting file like those used by the maketrack
// programs, specifying in which columns of the input observation files the
// required quantities are to be found.
// Example column formatting file:
//
// MJDCOL 3
// RACOL 6
// DECCOL 8

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define MJDCOL 3
#define RACOL 6
#define DECCOL 8
      
static void show_usage()
{
  cerr << "Usage: sunelong_cull01a -in infile -sunelong minsunelong(deg) maxsunelong(deg) \n";
  cerr << "-colformat column format file -earth earthfile -out outfile\n";
}
    
int main(int argc, char *argv[])
{
  string infile,earthfile,colformatfile,outfile,stest,lnfromfile;
  double minsunelong,maxsunelong;
  int mjdcol,racol,deccol,i,j,c,reachedeof,lct;
  ifstream instream1;
  ofstream outstream1;
  double MJD,RA,Dec,oldMJD;
  oldMJD=0.0L;
  vector <point3d> Earthpos;
  vector <point3d> Earthvel;
  vector <double> EarthMJD;
  point3d outpos = point3d(0,0,0);
  point3d unitbary = point3d(0,0,0);
  double barydist,obsdot,opdistcos,opelong,sunelong;

  if(argc<7)
    {
      show_usage();
      return(1);
    }
  cerr << "Usage: sunelong_cull01a -in infile -sunelong minsunelong(deg) maxsunelong(deg) \n";
  cerr << "-colformat column format file -earth earthfile -out outfile\n";

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-in" || string(argv[i]) == "-infile" || string(argv[i]) == "-inputfile" || string(argv[i]) == "--in" || string(argv[i]) == "--infile" || string(argv[i]) == "--inputfile") {
      if(i+1 < argc) {
	//There is still something to read;
	infile=argv[++i];
	i++;
      }
      else {
	cerr << "Input file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-sunelong" || string(argv[i]) == "-se" || string(argv[i]) == "-selong" || string(argv[i]) == "-solelong" || string(argv[i]) == "--solarelongation" || string(argv[i]) == "--solar_elongation" || string(argv[i]) == "--sunelong" || string(argv[i]) == "--solelong") {
      if(i+1 < argc) {
	//There is still something to read;
        minsunelong=stod(argv[++i]);
	if(minsunelong<0.0) {
	  cerr << "Error: invalid minimum solar elongation (" << minsunelong << " deg) supplied.\n";
	  cerr << "Solar elongation must be non-negative!\n";
	  return(2);
	}
      }
      else {
	cerr << "Solar elongation keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
      if(i+1 < argc) {
	//There is still something to read;
        maxsunelong=stod(argv[++i]);
	i++;
	if(!isnormal(maxsunelong) || maxsunelong<=0.0) {
	  cerr << "Error: invalid minimum solar elongation (" << minsunelong << " deg) supplied.\n";
	  cerr << "Solar elongation must be strictly positive!\n";
	  return(2);
	}
      }
      else {
	cerr << "Solar elongation keyword supplied with only one of two required arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-earth" || string(argv[i]) == "-e" || string(argv[i]) == "-Earth" || string(argv[i]) == "--earthfile" || string(argv[i]) == "--Earthfile" || string(argv[i]) == "--earth" || string(argv[i]) == "--Earth") {
      if(i+1 < argc) {
	//There is still something to read;
	earthfile=argv[++i];
	i++;
      }
      else {
	cerr << "Earth file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" || string(argv[i]) == "--cformat") {
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
    }  else if(string(argv[i]) == "-outfile" || string(argv[i]) == "-of" || string(argv[i]) == "-out" || string(argv[i]) == "-output" || string(argv[i]) == "--outputfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--output") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Column format file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }

  if(infile.size()<=0) {
    cerr << "Please supply an input detection file:\n\n";
    show_usage();
    return(1);
  }

  if(earthfile.size()<=0) {
    cerr << "Please supply a heliocentric ephemeris file for the Earth:\n\n";
    show_usage();
    return(1);
  }
  
  cout << "input file " << infile << "\n";
  cout << "column formatting file " << colformatfile << "\n";
  cout << "output file " << outfile << "\n";
  cout << "Heliocentric ephemeris file for the Earth: " << earthfile << "\n";
  cout << "minimum solar elongation " << minsunelong << " deg\n";
  cout << "maximum solar elongation " << maxsunelong << " deg\n";
  
  // Read the column formatting file, if any
  if(colformatfile.size()>0)
    {
      instream1.open(colformatfile);
      if(!instream1)  {
	cerr << "ERROR: unable to open input file " << colformatfile << "\n";
	return(1);
      }
      while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
	instream1 >> stest;
	if(stest == "MJDCOL") {
	  instream1 >> mjdcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ;
	} else if(stest == "RACOL") {
	  instream1 >> racol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ;
	} else if(stest == "DECCOL") {
	  instream1 >> deccol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ;
	}
      }
      instream1.close();
    }

  cout << "Column specifications:\n";
  cout << "MJDCOL " << mjdcol << "\n";
  cout << "RACOL " << racol << "\n";
  cout << "DECCOL " << deccol << "\n";

  // Read input Earth ephemeris file
  EarthMJD={};
  Earthpos={};
  Earthvel={};
  read_horizons_csv(earthfile,EarthMJD,Earthpos,Earthvel);
  cout << "Finished reading heliocentric ephemeris file " << earthfile << " for the Earth.\n";

  
  // Read input detection file.
  instream1.open(infile);
  outstream1.open(outfile);
  if(!instream1) {
    cerr << "can't open input file " << infile << "\n";
    return(1);
  }
  if(!outstream1) {
    cerr << "can't open output file " << outfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  lct++;
  cout << lnfromfile << "\n";
  // Do copy it to the output file, though.
  outstream1 << lnfromfile << "\n";
  reachedeof=0;
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
    while(i<lnfromfile.size() && lnfromfile.size()>=30 && reachedeof == 0) {
      // Note check on line length: it is completely impossible for a
      // line containing all the required quantities at minimum plausible
      // precision to be less than 30 characters long.
      c='0';
      stest="";
      while(i<lnfromfile.size() && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==mjdcol) MJD=stold(stest);
      else if(j==racol) RA=stold(stest);
      else if(j==deccol) Dec=stold(stest);
    }
    if(reachedeof == 0 && lnfromfile.size()>=30) {
      if(MJD!=oldMJD) {
	// The time of this detection is not the same as the
	// one before: re-caculate the observer's position.
	planetpos01(MJD, 5, EarthMJD, Earthpos, outpos);
	oldMJD=MJD;
      }
      celestial_to_stateunit(RA,Dec,unitbary);      
      barydist = sqrt(dotprod3d(outpos,outpos));
      obsdot = dotprod3d(unitbary,outpos);
      opdistcos = obsdot/barydist;
      opelong = acos(opdistcos)*DEGPRAD;
      if(opelong < 0.0L) opelong = 90.0L - opelong;
      sunelong = 180.0L - opelong;
      if(sunelong>=minsunelong && sunelong<=maxsunelong) {
	// Write line to output file
	outstream1 << lnfromfile << "\n";
      }
    }
  }
  instream1.close();
  outstream1.close();
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  return(0);
}
