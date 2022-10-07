// September 13, 2022: helio_input_select.cpp
// Given a file of astronomical detections suitable for input
// into the make_tracklets programs and a file of
// linked detections produced by heliolinc or heliovane,
// create a new astronomical detection file, identical in
// format to the original input file, consisting only of
// detections that WERE (-keeplink 1) or WERE NOT (-keeplink 0)
// recorded as linked in the linked detection file.
//
// helio_input_select relies on the fact that in the comprehensive
// output files (linked detection files) produced both by heliolinc
// itself and by the post-processing codes link_refine and
// link_refine_multisite, the line number **in the original input
// file** is recorded in column 10 of the csv files. Note well
// that an awk-style counting convention is used for these line
// numbers: the header line is line 1, so the first data line is
// line number 2.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define LINKFILE_INDEX_COL 10

static void show_usage()
{
  cerr << "Usage: helio_input_select -infile input_detections_file -linked linked_detection_file -keeplink 0 -outfile output_file\n";
}
    
int main(int argc, char *argv[])
{
  int keeplink=0;
  int i=0;
  int startpoint=0;
  int endpoint=0;
  int linect=0;
  int badread=0;
  string infile,linkfile,outfile,stest,lnfromfile;
  int ldnum=0;
  int ldct=0;
  ifstream instream1;
  ofstream outstream1;
  long origind=0;
  vector <long> indvec,tempvec;

  if(argc!=9) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-if" || string(argv[i]) == "-in" || string(argv[i]) == "--infile" || string(argv[i]) == "--inputfile" || string(argv[i]) == "--input_detections_file" || string(argv[i]) == "--indetfile") {
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
    } else if(string(argv[i]) == "-linked" || string(argv[i]) == "-linkfile" || string(argv[i]) == "-lf" || string(argv[i]) == "--linked" || string(argv[i]) == "--linkfile" || string(argv[i]) == "--linked_detection_file" || string(argv[i]) == "-ldf") {
      if(i+1 < argc) {
	//There is still something to read;
	linkfile=argv[++i];
	i++;
      }
      else {
	cerr << "Linked detection file supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-keeplink" || string(argv[i]) == "-kl" || string(argv[i]) == "-keep" || string(argv[i]) == "-klink" || string(argv[i]) == "--keeplink" || string(argv[i]) == "--kl") {
      if(i+1 < argc) {
	//There is still something to read;
	keeplink=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Keep/reject keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-of" || string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--fileout" || string(argv[i]) == "--fileoutput") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  cout << "input raw detection file " << infile << "\n";
  cout << "input linked detection file " << linkfile << "\n";
  if(keeplink<=0) {
    cout << "keeplink = " << keeplink << ": previously linked detections will be rejected\n";
    keeplink=0;
  } else if(keeplink>=1) {
    cout << "keeplink = " << keeplink << ": detections not previously linked will be rejected\n";
    keeplink=1;
  }
  cout << "output file to be written " << outfile << "\n";

  // Read linked detection file
  instream1.open(linkfile,ios_base::in);
  if(!instream1) {
    cerr << "can't open input file " << linkfile << "\n";
    return(1);
  }
  tempvec={};
  // Skip header line
  getline(instream1,lnfromfile);
  cout << "Header line from input linked detection file " << linkfile << ":\n";
  cout << lnfromfile << "\n";
  // Read body of the file
  linect=0;
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    // Read a line from the linked detections file, and save the original file index
    getline(instream1,lnfromfile);
    linect++;
    badread=0;
    if(lnfromfile.size()>60 && !instream1.bad() && !instream1.fail() && !instream1.eof()) {
      i=startpoint=0;
      while(i<LINKFILE_INDEX_COL) {
	endpoint = get_csv_string01(lnfromfile,stest,startpoint);
	startpoint = endpoint+1;
	i++;
      }
      if(endpoint>0) {
	origind = stol(stest);
	tempvec.push_back(origind);
      } else badread=1;
      if(!instream1.bad() && !instream1.fail() && !instream1.eof() && badread==1) {
	cerr << "ERROR reading line " << linect << " of linked detection file " << linkfile << "\n";
	return(1);
      }
    } else if (!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      cerr << "WARNING: line " << linect << " of linked detection file " << linkfile << " was too short\n";
    }
  }
  instream1.close();
  ldnum = tempvec.size();
  cout << "Successfully read " << ldnum << " detections from " << linkfile << "\n";
  // Sort tempvec
  sort(tempvec.begin(),tempvec.end());
  // Copy unique-ified version of tempvec into indvec

  indvec={};
  indvec.push_back(tempvec[0]);
  for(ldct=1;ldct<ldnum;ldct++) {
    if(tempvec[ldct]>indvec[indvec.size()-1]) indvec.push_back(tempvec[ldct]);
  }
  cout << "Of these " << ldnum << " detections read, ";
  ldnum = indvec.size();
  cout << ldnum << " are unique\n";
  
  // Read the input raw detection file
  instream1.open(infile,ios_base::in);
  outstream1.open(outfile,ios_base::out);
  linect=0;
  if(!instream1) {
    cerr << "can't open input file " << infile << "\n";
    return(1);
  }
  // Read header line
  getline(instream1,lnfromfile);
  // Copy it to output file
  outstream1 << lnfromfile << "\n";
  if(keeplink==1) {
    // Read body of input file, and write out lines whose numbers
    // correspond to the indices recorded in indvec
    getline(instream1,lnfromfile);
    ldct=0;
    linect=2; // Original-file line indexing counts lines starting from 1
              // and includes the header, so the first data line is 2.
              // This makes it match awk NR.
    while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
      if(linect==indvec[ldct]) {
	// Output this line, and advance both linect and ldct
	outstream1 << lnfromfile << "\n";
	ldct++;
	getline(instream1,lnfromfile);
	linect++;
      } else if(linect<indvec[ldct]) {
	// Skip to the next line
	getline(instream1,lnfromfile);
	linect++;
      } else if(linect>indvec[ldct]) {
	// Move to the next entry in indvec.
	ldct++;
      } else {
	cerr << "ERROR: logically excluded case linect = " << linect << ", indvec[ldct] = " << indvec[ldct] << ", ldct = " << ldct << "\n";
	return(1);
      }
    }
  } else {
    // Read body of input file, and write out lines whose numbers
    // DO NOT correspond to the indices recorded in indvec
    getline(instream1,lnfromfile);
    ldct=0;
    linect=2; // Original-file line indexing counts lines starting from 1
              // and includes the header, so the first data line is 2.
              // This makes it match awk NR.
    while(!instream1.bad() && !instream1.fail() && !instream1.eof() && ldct<ldnum) {
      if(linect==indvec[ldct]) {
	// Advance both linect and ldct
	ldct++;
	getline(instream1,lnfromfile);
	linect++;
      } else if(linect<indvec[ldct]) {
	// Output this line and then go on to the next one
	outstream1 << lnfromfile << "\n";
	getline(instream1,lnfromfile);
	linect++;
      } else if(linect>indvec[ldct]) {
	// Move to the next entry in indvec.
	ldct++;
      } else {
	cerr << "ERROR: logically excluded case linect = " << linect << ", indvec[ldct] = " << indvec[ldct] << ", ldct = " << ldct << "\n";
	return(1);
      }
    }
  }
  instream1.close();
  outstream1.close();
  return(0);
}

	
