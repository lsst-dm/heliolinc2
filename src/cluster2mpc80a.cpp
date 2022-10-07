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
  cerr << "Usage: cluster2mpc80a -pairdet detection_file -clust clusterfile -rms rmsfile -prefix idprefix -outfile outfile\n";
}
    
int main(int argc, char *argv[])
{
  int i=1;
  int j=0;
  int linect=0;
  int linenum=0;
  string pairdetfile,clusterfile,rmsfile,outfile,rating,stest;
  string lnfromfile;
  string w1,w2;
  vector <det_obsmag_indvec> detvec;
  det_obsmag_indvec dsv = det_obsmag_indvec(0l,0l,0l,0L,0L,0L,"null",0.0,"V","I11",1,{});
  double posrms,velrms,totrms;
  posrms = velrms = totrms = 0l;
  double heliodist, heliovel, helioacc;
  heliodist = heliovel = helioacc = 0l;
  double x,y,z,vx,vy,vz;
  x = y = z = vz = vy = vz = 0l;
  vector <string> inputlines;
  ifstream instream1;
  ifstream instream2;
  ofstream outstream1;
  vector <double> mjdvec;
  float clustmetric=0;
  int ptnum=0;
  int pairnum=0;
  int ptct=0;
  int nightobs=0;
  double timespan=0l;
  int obsnights=0;
  int globalct = 0;
  int ijunk1 = 0;
  int ijunk2 = 0;
  long i1 = 0;
  long i2 = 0;
  double MJD,RA,Dec;
  MJD = RA = Dec = 0l;
  int rmslinect,clustlinect;
  int badread=0;
  int startpoint=0;
  int endpoint=0;
  char detid[SHORTSTRINGLEN];
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
  
  if(argc<9) {
    show_usage();
    return(1);
  }
  
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pd" || string(argv[i]) == "-pdet" || string(argv[i]) == "--pairdet" || string(argv[i]) == "--paireddetections" || string(argv[i]) == "--pairdetfile" || string(argv[i]) == "--pairdetections") {
      if(i+1 < argc) {
	//There is still something to read;
	pairdetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input pairdet file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-cf" || string(argv[i]) == "-c" || string(argv[i]) == "-clust" || string(argv[i]) == "--clusterfile" || string(argv[i]) == "--cfile" || string(argv[i]) == "-cfile" || string(argv[i]) == "--clustfile") {
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
    }  else if(string(argv[i]) == "-rf" || string(argv[i]) == "-r" || string(argv[i]) == "-rms" || string(argv[i]) == "-rmsfile" || string(argv[i]) == "--rmsfile" || string(argv[i]) == "-rfile" || string(argv[i]) == "--rfile") {
      if(i+1 < argc) {
	//There is still something to read;
	rmsfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input RMS file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }   else if(string(argv[i]) == "-id" || string(argv[i]) == "-idp"  || string(argv[i]) == "-prefix" || string(argv[i]) == "-pre" || string(argv[i]) == "-idpre" || string(argv[i]) == "--idprefix" || string(argv[i]) == "-stringid" || string(argv[i]) == "--stringidprefix" || string(argv[i]) == "--prefix") {
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
    
  cout << "input paired detection file " << pairdetfile << "\n";
  cout << "input cluster file " << clusterfile << "\n";
  cout << "input rms file " << rmsfile << "\n";
  cout << "output id prefix " << idprefix << "\n";
  cout << "output cluster file " << outfile << "\n";

  if(idprefix.size()>6) {
    cerr << "ERROR: maximum length for idprefix is 6 characters\n";
    return(1);
  }

  // Read paired detection file
  instream1.open(pairdetfile,ios_base::in);
  if(!instream1) {
    cerr << "can't open input file " << pairdetfile << "\n";
    return(1);
  }
  detvec={};
  // Skip header line
  getline(instream1,lnfromfile);
  cout << "Header line from input paired detection file " << pairdetfile << ":\n";
  cout << lnfromfile << "\n";
  // Read body of the file
  detfilelinect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the paired detections file, and load an object of class det_obsmag_indvec
    getline(instream1,lnfromfile);
    detfilelinect++;
    badread=0;
    //while(instream1 >> MJD >> RA >> Dec >> X >> Y >> Z >> detid >> mag >> band >> obscode >> origind) {
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
      if(endpoint>0) x = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) y = stold(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) z = stold(stest);
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
	cerr << "ERROR reading line " << detfilelinect << " of paired detection file" << pairdetfile << "\n";
	return(1);
      }
      // If we reach this point, the line was read OK. Write it to detvec.
      dsv=det_obsmag_indvec(MJD,RA,Dec,x,y,z,detid,mag,band,obscode,origind,{});
      dsv.indvec = {};
      detvec.push_back(dsv);
    } else if(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      cerr << "WARNING: line " << detfilelinect << " of paired detection file " << pairdetfile << " was too short\n";
    }
  }
  instream1.close();
  cout << "Successfully read " << detvec.size() << " detections from " << pairdetfile << "\n";

  // Open cluster file
  instream1.open(clusterfile,ios_base::in);
  instream2.open(rmsfile,ios_base::in);
  if(!instream1) {
    cerr << "can't open cluster file " << clusterfile << "\n";
    return(2);
  } else cout << "Reading cluster file " << clusterfile << "\n";
  if(!instream2) {
    cerr << "can't open rms file " << rmsfile << "\n";
    return(2);
  } else cout << "Reading rms file " << rmsfile << "\n";
  // Skip header lines
  getline(instream1,lnfromfile);
  getline(instream2,lnfromfile);
  rmslinect=clustlinect=0;

  // Open output file
  outstream1.open(outfile,ios_base::out);

  rmslinect=clustlinect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof() && !instream2.bad() && !instream2.fail() && !instream2.eof()) {
    // Read a line from the rms file
    getline(instream2,lnfromfile);
    rmslinect++;
    badread=0;
    if(lnfromfile.size()>40) {
      // Read cluster index number;
      startpoint=0;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) globalct = stoi(stest);
      else badread=1;
      // Read three RMS values
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) posrms = stof(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) velrms = stof(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) totrms = stof(stest);
      else badread=1;

      // Read the integer pairnum
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) pairnum = stoi(stest);
      else badread=1;
	
      // Read the double timespan
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) timespan = stod(stest);
      else badread=1;
	
      // Read the integers ptnum and obsnights
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) ptnum = stoi(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) obsnights = stoi(stest);
      else badread=1;
	
      // Read the float clustmetric
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) clustmetric = stof(stest);
      else badread=1;
	
      // read the string rating
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) rating=stest;

      // Read three heliocentric hypothesis parameters.
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) heliodist = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) heliovel = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) helioacc = stod(stest);
      else badread=1;
	
      // Read the six elements of the mean state vector
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) x = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) y = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) z = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) vx = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) vy = stod(stest);
      else badread=1;
      startpoint = endpoint+1;
      if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      if(endpoint>0) vz = stod(stest);
      else badread=1;
     
      // If there was a file read error, abort.
      if(badread==1) {
	cerr << "ERROR reading line " << rmslinect << " of rms file " << rmsfile << "\n";
	return(1);
      }
      
      // Read the associated lines from the cluster file, but retain only the
      // indices to the original input file
      for(i=0;i<ptnum;i++) {
	if (!instream1.bad() && !instream1.fail() && !instream1.eof()) {
	  // Read a line from the cluster file.
	  getline(instream1,lnfromfile);
	  clustlinect++;
	  while(lnfromfile.size()<40 && !instream1.bad() && !instream1.fail() && !instream1.eof()) {
	    cerr << "WARNING: line " << clustlinect << " of cluster file " << clusterfile << " is too short\n";
	    // Read another line, maybe there's just a blank one.
	    getline(instream1,lnfromfile);
	    clustlinect++;
	  }
	  badread=0;
	  if(lnfromfile.size()>40) {
	    // Read and discard the first nine quantities: ptct, MJD, RA, Dec, mag, band, obscode, and i1.
	    startpoint=0;
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    // Read the essential quantity: the index to the original file
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) {
	      startpoint = endpoint+1;
	      i1 = stoi(stest);
	    } else badread=1;
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) startpoint = endpoint+1;
	    else badread=1;
	    // Read clusterct, and check it against the index read from the rms file
	    if(badread==0) endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	    if(endpoint>0) {
	      startpoint = endpoint+1;
	      i2 = stoi(stest);
	      if(i2 != globalct) {
		cerr << "ERROR: cluster count mismatch at rms line " << rmslinect << ", cluster file line " << clustlinect << "\n";
		return(1);
	      } else {
		// All good. Write this detection to the output file.
		outstream1 << "     "; // 5 initial spaces
		stringstream ss;
		ss << idprefix;
		stest=intzero01i(globalct,7-idprefix.size());
		ss << stest;
		str = ss.str();
		stringncopy01(mpcid,str,8);
		outstream1 << mpcid << "  C";
		mjd2mpcdate(detvec[i1].MJD,year,month,day);
		outstream1 << year << " ";
		if(month<10) outstream1 << "0";
		outstream1 << month << " ";
		if(day<10.0l) outstream1  << fixed << setprecision(6) << "0";
		outstream1  << fixed << setprecision(6) << day;
		// Convert RA, Dec from decimal degrees to sexagesimal format.
		rahr = int(detvec[i1].RA/15.0l);
		ramin = int(detvec[i1].RA*4.0l - double(rahr)*60.0l);
		rasec = detvec[i1].RA*240.0l - double(rahr)*3600.0l - double(ramin)*60.0l;
		if(detvec[i1].Dec>=0) {
		  signstring="+";
		  Dec = detvec[i1].Dec;
		} else {
		  signstring="-";
		  Dec = -detvec[i1].Dec;
		}
		decdeg = int(Dec);
		decmin = int(Dec*60.0l - double(decdeg)*60.0l);
		decsec = Dec*3600.0l - double(decdeg)*3600.0l - double(decmin)*60.0l;
	      }
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
	      outstream1 << fixed << setprecision(1) << detvec[i1].mag << " " << detvec[i1].band;
	      // Add correct number of spaces after band.
	      bandlen = j = 0;
	      while(j<MINSTRINGLEN && detvec[i1].band[j]!='\0') {
		bandlen++;
		j++;
	      }		    
	      for(j=0;j<7-bandlen;j++) outstream1 << " ";
	      // Write out obscode
	      outstream1 << detvec[i1].obscode << "\n";
	    } else badread=1;
	    // If there was a file read error, abort.
	    if(badread==1) {
	      cerr << "ERROR reading line " << clustlinect << " of cluster file " << clusterfile << "\n";
	      return(1);
	    } 
	  }
	  // Close if-statement checking if we're still reading valid lines from cluster file.
	} else {
	  cerr << "WARNING: unsuccessful read at line " << clustlinect << " of cluster file " << clusterfile << "\n";
	}
	// Close loop over points in the cluster
      }
    }
  }
  instream1.close();
  instream2.close();
  outstream1.close();

  return(0);
}
