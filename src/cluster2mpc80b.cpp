// March 16, 2022: cluster2mpc80a.cpp
// Given input cluster and rms files produced by projectpairs06c,
// ancluster03c, or similar programs, and convert them into
// MPC 80-column format.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MAXCLUSTRMS 1.0e5
#define ONE_POINT_PER_IMAGE 1

static void show_usage()
{
  cerr << "Usage: cluster2mpc80b -clust clusterfile -prefix idprefix -outfile outfile\n";
}
    
int main(int argc, char *argv[])
{
  int i=1;
  int j=0;
  int linect=0;
  int linenum=0;
  string clusterfile,outfile,rating,stest;
  string lnfromfile;
  string w1,w2;
  ifstream instream1;
  ofstream outstream1;
  int ptct=0;
  double MJD,RA,Dec;
  MJD = RA = Dec = 0l;
  long i1,i2;
  int clusterct;
  int clustlinect;
  int badread=0;
  int startpoint=0;
  int endpoint=0;
  double mag = 0;
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  string idprefix;
  string str,signstring;
  char mpcid[8];
  int year,month;
  double day;
  int rahr,ramin;
  int decdeg,decmin;
  double rasec,decsec;
  int detfilelinect=0;
  int origind=0;
  int bandlen=0;
  
  if(argc<7) {
    show_usage();
    return(1);
  }
  
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-cf" || string(argv[i]) == "-c" || string(argv[i]) == "-clust" || string(argv[i]) == "--clusterfile" || string(argv[i]) == "--cfile" || string(argv[i]) == "-cfile" || string(argv[i]) == "--clustfile") {
      if(i+1 < argc) {
	//There is still something to read;
	clusterfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input cluster file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-id" || string(argv[i]) == "-idp"  || string(argv[i]) == "-prefix" || string(argv[i]) == "-pre" || string(argv[i]) == "-idpre" || string(argv[i]) == "--idprefix" || string(argv[i]) == "-stringid" || string(argv[i]) == "--stringidprefix" || string(argv[i]) == "--prefix") {
      if(i+1 < argc) {
	//There is still something to read;
	idprefix=argv[++i];
	i++;
      } else {
	cerr << "ID prefix keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-of" || string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--fileout" || string(argv[i]) == "--fileoutput") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }  else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
    
  cout << "input cluster file " << clusterfile << "\n";
  cout << "output id prefix " << idprefix << "\n";
  cout << "output cluster file " << outfile << "\n";

  if(idprefix.size()>6) {
    cerr << "ERROR: maximum length for idprefix is 6 characters\n";
    return(1);
  }


  // Open cluster file
  instream1.open(clusterfile,ios_base::in);
  if(!instream1) {
    cerr << "can't open cluster file " << clusterfile << "\n";
    return(2);
  }
  getline(instream1,lnfromfile);

  // Open output file
  outstream1.open(outfile,ios_base::out);

  clustlinect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the cluster file
    getline(instream1,lnfromfile);
    clustlinect++;
    badread=0;
    while(lnfromfile.size()<40 && !instream1.bad() && !instream1.fail() && !instream1.eof()) {
      cerr << "WARNING: line " << clustlinect << " of cluster file " << clusterfile << " is too short\n";
      // Read another line, maybe there's just a blank one.
      getline(instream1,lnfromfile);
      clustlinect++;
    }
    badread=0;
    if(lnfromfile.size()>40) {
      // Read all 11 quantities: ptct, MJD, RA, Dec, mag, band, obscode, i1, i2, and clusterct.
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) ptct = stoi(stest);
      else badread=1;
      cout << "ptct = " << ptct << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) MJD = stold(stest);
      else badread=1;
      cout << "MJD = " << MJD << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) RA = stold(stest);
      else badread=1;
      cout << "RA = " << RA << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Dec = stold(stest);
      else badread=1;
      cout << "Dec = " << Dec << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) ;
      else badread=1;
      cout << "idstring = " << stest << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) mag = stod(stest);
      else badread=1;
      cout << "mag = " << mag << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(band,stest,MINSTRINGLEN);
      else badread=1;
      cout << "band = " << band << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(obscode,stest,MINSTRINGLEN);
      else badread=1;
      cout << "obscode = " << obscode << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) i1 = stol(stest);
      else badread=1;
      cout << "i1 = " << i1 << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) i2 = stol(stest);
      else badread=1;
      cout << "i2 = " << i2 << "\n";
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) clusterct = stoi(stest);
      else badread=1;
      cout << "clusterct = " << clusterct << "\n";
      // All good. Write this detection to the output file.
      outstream1 << "     "; // 5 initial spaces
      stringstream ss;
      ss << idprefix;
      stest=intzero01i(clusterct,7-idprefix.size());
      ss << stest;
      str = ss.str();
      stringncopy01(mpcid,str,8);
      outstream1 << mpcid << "  C";
      mjd2mpcdate(MJD,year,month,day);
      outstream1 << year << " ";
      if(month<10) outstream1 << "0";
      outstream1 << month << " ";
      if(day<10.0l) outstream1  << fixed << setprecision(6) << "0";
      outstream1  << fixed << setprecision(6) << day;
      // Convert RA, Dec from decimal degrees to sexagesimal format.
      rahr = int(RA/15.0l);
      ramin = int(RA*4.0l - double(rahr)*60.0l);
      rasec = RA*240.0l - double(rahr)*3600.0l - double(ramin)*60.0l;
      if(Dec>=0) {
	signstring="+";
	Dec = Dec;
      } else {
	signstring="-";
	Dec = -Dec;
      }
      decdeg = int(Dec);
      decmin = int(Dec*60.0l - double(decdeg)*60.0l);
      decsec = Dec*3600.0l - double(decdeg)*3600.0l - double(decmin)*60.0l;
      // Write out RA and Dec with appropriate zero-padding.
      if(rahr<10) outstream1 << "0";
      outstream1 << rahr << " ";
      if(ramin<10) outstream1 << "0";
      outstream1 << ramin << " ";
      if(rasec<10.0l) outstream1 << "0";
      outstream1 << fixed << setprecision(3) << rasec << signstring;
      if(decdeg<10) outstream1 << "0";
      outstream1 << decdeg << " ";
      if(decmin<10) outstream1 << "0";
      outstream1 << decmin << " ";
      if(decsec<10.0l) outstream1 << "0";
      outstream1 << fixed << setprecision(2) << decsec << "         ";
      // Write out magnitude and band.
      outstream1 << fixed << setprecision(1) << mag << " " << band;
      // Add correct number of spaces after band.
      bandlen = j = 0;
      while(j<MINSTRINGLEN && band[j]!='\0') {
	bandlen++;
	j++;
      }		    
      for(j=0;j<7-bandlen;j++) outstream1 << " ";
      // Write out obscode
      outstream1 << obscode << "\n";
    }
    // If there was a file read error, abort.
    if(badread==1) {
      cerr << "ERROR reading line " << clustlinect << " of cluster file " << clusterfile << "\n";
      return(1);
    }
  }
  instream1.close();
  outstream1.close();

  return(0);
}
