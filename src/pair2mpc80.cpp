// March 17, 2023: pair2mpc80.cpp
// Given output from make_tracklets, turn it into an MPC80
// file, suitable for pasting into MPChecker or submitting
// to the MPC (if they don't require ADES yet).

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define TYPECODE_PAIR "P"
#define TYPECODE_TRACKLET "T"

static void show_usage()
{
  cerr << "Usage: pair2mpc80 -dets detfile -pairs pairfile -prefix idprefix -out outfile \n";
  cerr << " \n or, at minimum: \n";
  cerr << "pair2mpc80 -dets detfile -pairs pairfile\n";
  }
    
int main(int argc, char *argv[])
{
  det_obsmag_indvec o1 = det_obsmag_indvec(0L,0l,0l,0L,0L,0L,"null",0l,"V","I11",0,{});
  vector <det_obsmag_indvec> detvec = {};
  vector <long> tempind={};
  long double MJD,X,Y,Z;
  MJD = X = Y = Z = 0.0L;
  double RA,Dec,mag;
  RA = Dec = mag = 0.0l;
  char detid[SHORTSTRINGLEN];
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  string idprefix = "hl";
  string stest;
  long origind=0;
  long i1,i2,ipt;
  i1=i2=ipt=0;
  string typecode;
  string indetfile,inpairfile,lnfromfile;
  string outfile = "output.mpc";
  ifstream instream1;
  int default_prefix = 1;
  int default_outfile = 1;
  long detfilelinect;
  int badread=0;
  long startpoint=0;
  long endpoint=0;
  long tracklet_ct=0;
  string str,signstring;
  char mpcid[8];
  int year,month;
  year = month = 0;
  double day=0.0l;
  int rahr,ramin;
  rahr = ramin = 0;
  int decdeg,decmin;
  decdeg = decmin = 0;
  double rasec,decsec;
  rasec = decsec = 0.0l;
  long i=0;
  long j=0;
  long bandlen=0;
  int trackpointnum=0;
  int trackpointct=0;
  
  
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
    } else if(string(argv[i]) == "-p" || string(argv[i]) == "-pair" || string(argv[i]) == "-pairs" || string(argv[i]) == "--pairs" || string(argv[i]) == "--pair" || string(argv[i]) == "--pairfile" || string(argv[i]) == "--pairsfile") {
      if(i+1 < argc) {
	//There is still something to read;
	inpairfile=argv[++i];
	i++;
      }
      else {
	cerr << "Pair file keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-id" || string(argv[i]) == "-idp"  || string(argv[i]) == "-prefix" || string(argv[i]) == "-pre" || string(argv[i]) == "-idpre" || string(argv[i]) == "--idprefix" || string(argv[i]) == "-stringid" || string(argv[i]) == "--stringidprefix" || string(argv[i]) == "--prefix") {
      if(i+1 < argc) {
	//There is still something to read;
	idprefix=argv[++i];
	default_prefix = 0;
	i++;
      } else {
	cerr << "ID prefix keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outpair" || string(argv[i]) == "--outpairs") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	default_outfile = 0;
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
  
  if(argc<5)
    {
      cerr << "Too few arguments even for minimalist invocation:\n";
      show_usage();
      return(1);
    }
  
  cout.precision(17);  
  cout << "input detection file " << indetfile << "\n";
  cout << "input pair file " << inpairfile << "\n";
  if(default_prefix==0) cout << "ID prefix " << idprefix << "\n";
  else cout << "Using default ID prefix " << idprefix << "\n";
  if(default_outfile==0) cout << "output file " << outfile << "\n";
  else cout << "Using default output file name " << outfile << "\n";

  // Read input detection file.
  instream1.open(indetfile,ios_base::in);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << indetfile << "\n";
    return(1);
  }
  detvec={};
  // Skip header line
  getline(instream1,lnfromfile);
  cout << "Header line from input paired detection file " << indetfile << ":\n";
  cout << lnfromfile << "\n";
  // Read body of the file
  detfilelinect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the paired detections file, and load an object of class det_obsmag_indvec
    getline(instream1,lnfromfile);
    detfilelinect++;
    badread=0;
    if(lnfromfile.size()>60) {
      // Read MJD, RA, Dec, observer x, y, z
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) MJD = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) RA = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Dec = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) X = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Y = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) Z = stold(stest);
      else badread=1;
      // Read the IDstring
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(detid,stest,SHORTSTRINGLEN);
      else badread=1;
      // Read the magnitude
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) mag = stod(stest);
      else badread=1;
      // Read the band and observatory code
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(band,stest,MINSTRINGLEN);
      else badread=1;
       startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) stringncopy01(obscode,stest,MINSTRINGLEN);
      else badread=1;
      // Read the original detection index
       startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) origind = stol(stest);
      else badread=1;

      // If there was a file read error, abort.
      if(badread==1) {
	cerr << "ERROR reading line " << detfilelinect << " of paired detection file" << indetfile << "\n";
	return(1);
      }
      // If we reach this point, the line was read OK. Write it to detvec.
      o1=det_obsmag_indvec(MJD,RA,Dec,X,Y,Z,detid,mag,band,obscode,origind,{});
      detvec.push_back(o1);
    } else if(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      cerr << "WARNING: line " << detfilelinect << " of paired detection file " << indetfile << " was too short\n";
    }
  }
  instream1.close();
  cout << detvec.size() << " detection records read from " << indetfile << ".\n";
  
  // Read input image pair file
  instream1.open(inpairfile);
  ofstream outstream1 {outfile};
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << inpairfile << "\n";
    return(1);
  }
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    typecode = ""; // Wipe previously read typecode
    // Read the current type code
    instream1 >> typecode;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad() && typecode == TYPECODE_PAIR) {
      // Read a single, isolated pair
      instream1 >> i1 >> i2;
      tempind.push_back(i1);
      tempind.push_back(i2);
    } else if(!instream1.eof() && !instream1.fail() && !instream1.bad() && typecode == TYPECODE_TRACKLET) {
      // Read the representative pair
      instream1 >> i1 >> i2;
      // Read revised RA, Dec for representative detection 1.
      instream1 >> RA >> Dec;
      instream1 >> RA >> Dec;
      // Read number of points in the tracklet
      instream1 >> trackpointnum;
      cout << "Reading points of a " << trackpointnum << "-point tracklet\n";
      for(trackpointct=0; trackpointct<trackpointnum; trackpointct++) {
	instream1 >> ipt;
	tempind.push_back(ipt);
      }
    } else if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      cerr << "ERROR: unrecognized pair type code " << typecode << "\n";
      return(1);
    }
    // Loop over stored indices and write output
    for(i=0;i<long(tempind.size());i++) {
      outstream1 << "     "; // 5 initial spaces
      stringstream ss;
      ss << idprefix;
      stest=intzero01i(tracklet_ct,7-idprefix.size());
      ss << stest;
      str = ss.str();
      stringncopy01(mpcid,str,8);
      outstream1 << mpcid << "  C";
      mjd2mpcdate(detvec[tempind[i]].MJD,year,month,day);
      outstream1 << year << " ";
      if(month<10) outstream1 << "0";
      outstream1 << month << " ";
      if(day<10.0l) outstream1  << fixed << setprecision(6) << "0";
      outstream1  << fixed << setprecision(6) << day;
      // Convert RA, Dec from decimal degrees to sexagesimal format.
      rahr = int(detvec[tempind[i]].RA/15.0l);
      ramin = int(detvec[tempind[i]].RA*4.0l - double(rahr)*60.0l);
      rasec = detvec[tempind[i]].RA*240.0l - double(rahr)*3600.0l - double(ramin)*60.0l;
      if(detvec[tempind[i]].Dec>=0) {
	signstring="+";
	Dec = detvec[tempind[i]].Dec;
      } else {
	signstring="-";
	Dec = -detvec[tempind[i]].Dec;
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
      outstream1 << fixed << setprecision(1) << detvec[tempind[i]].mag << " " << detvec[tempind[i]].band;
      // Add correct number of spaces after band.
      bandlen = j = 0;
      while(j<MINSTRINGLEN && detvec[tempind[i]].band[j]!='\0') {
	bandlen++;
	j++;
      }		    
      for(j=0;j<7-bandlen;j++) outstream1 << " ";
      // Write out obscode
      outstream1 << detvec[tempind[i]].obscode << "\n";
    }
    tempind={};
    tracklet_ct++;
    if(tracklet_ct>=100000) {
      cerr << "Warning: maximum number of distinct tracklets exceeded\n";
      cerr << "Output file will be incomplete\n";
       instream1.close();
       outstream1.close();
       return(1);
    }
  }
  instream1.close();
  outstream1.close();
  return(0);
}
