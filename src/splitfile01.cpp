// May 04, 2023: splitfile01.cpp: Given an input file, split
// it into many different subfiles with a fixed maximum number
// of lines. For example, give a file with 473 lines and a maximum
// output size of 20 lines, will produce 23 output files with
// 20 lines each, and then one final file with 13 lines, since
// 473 = 23*20 + 13. Can be set to preserve a header line.

#include "solarsyst_dyn_geo01.h"
#include "cmath"


static void show_usage()
{
  cerr << "Usage: splitfile01 -infile input_file -linesout max_number_of_lines -header 1=header,0=no_header -outroot rootname_for_output_files -suffix suffix_for_output\n";
}

int main(int argc, char *argv[])
{
  string infile;
  string outroot="splittemp";
  string suffix=".txt";
  string outname,lnfromfile,headerline;
  int linesout=100;
  int useheader=0;
  ifstream instream1;
  ofstream outstream1;
  int readcount=0;
  int filecount=0;
  int i=0;
  string numstring;
  vector <string> newlines;

  if(argc<3) {
    show_usage();
    return(1);
  }
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-infile" || string(argv[i]) == "-in" || string(argv[i]) == "-input") {
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
    } else if(string(argv[i]) == "-linesout" || string(argv[i]) == "-linenum" || string(argv[i]) == "-lnum" || string(argv[i]) == "-lines" || string(argv[i]) == "-maxlines" || string(argv[i]) == "-outlines" || string(argv[i]) == "--linesout" || string(argv[i]) == "--output_lines") {
      if(i+1 < argc) {
	//There is still something to read;
	linesout=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Output file line number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-header" || string(argv[i]) == "-head" || string(argv[i]) == "-keepheader" || string(argv[i]) == "-isheader" || string(argv[i]) == "--header" ) {
      if(i+1 < argc) {
	//There is still something to read;
	useheader=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Header use keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outroot" || string(argv[i]) == "-outname") {
      if(i+1 < argc) {
	//There is still something to read;
	outroot=argv[++i];
	i++;
      }
      else {
	cerr << "Output file root name keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-suffix" || string(argv[i]) == "-postfix" || string(argv[i]) == "-sf" || string(argv[i]) == "-suf" || string(argv[i]) == "-pf") {
      if(i+1 < argc) {
	suffix=argv[++i];
	i++;
      }
      else {
	cerr << "Output file suffix keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  // Catch required parameters if missing
  if(infile.size()<=0) {
    cerr << "\nERROR: input file name is required\n";
    show_usage();
    return(1);
  }
  if(linesout<=0) {
    cerr << "ERROR: max line number for output files must be positive: " << linesout << " was supplied\n";
    show_usage();
    return(1);
  }
    
  cout << "input file " << infile << "\n";
  cout << "root name for output files: " << outroot << "\n";
  cout << "suffix for output files: " << suffix << "\n";
  cout << "output files will have " << linesout << " or fewer lines\n";
  if(useheader==1) {
    cout << "A one-line header from the input files\nwill be retained and used for all the output files.\n";
  } else {
    useheader=0;
    cout << "Input files are assumed to have no headers\n";
  }
  
  // Read input file
  instream1.open(infile);
  readcount=filecount=0;
  newlines = {};
  if(useheader==1) getline(instream1,headerline);
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    while(!instream1.eof() && !instream1.fail() && !instream1.bad() && readcount<linesout) {
      getline(instream1,lnfromfile);
      newlines.push_back(lnfromfile);
      readcount++;
    }
    if(newlines.size()>0) {
      // Construct output file name
      if(filecount<=999) {
	numstring = to_string(filecount);
	while(numstring.size()<3) {
	  numstring = "0" + numstring;
	}
	outname = outroot + numstring + suffix;       
	// Write output file.
	cout << "Writing output file called " << outname << " with " << long(newlines.size()) << " data lines\n";
	outstream1.open(outname);
	if(useheader) outstream1 << headerline << "\n";
	for(i=0;i<long(newlines.size());i++) {
	  outstream1 << newlines[i] << "\n";
	}
	outstream1.close();
	newlines={};
	readcount=0;
	filecount++;
      } else {
	cerr << "ERROR: splitfile01 cannot write more than 1000 output files\n";
	return(1);
      }
    }
  }
  
  instream1.close();
  return(0);
}
