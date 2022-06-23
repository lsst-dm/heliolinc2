// March 18, 2022: file_lineselect01.cpp
// Given a data file, and another file containing simply a list
// of line numbers, create a new output file containing only
// the lines of the data file corresponding to the line numbers,
// (-choose select) or only those that do not (-choose cull).

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MAXCLUSTRMS 1.0e5
#define ONE_POINT_PER_IMAGE 1


static void show_usage()
{
  cerr << "Usage: file_lineselect01 -data datafile -lines linefile -choose select/cull -out outfile\n";
}
    
int main(int argc, char *argv[])
{
  int i=1;
  string datafile,linefile,choice,outfile,lnfromfile;
  int datanum,datact,lnum,lct,linect,thisline;
  datanum=datact=lnum=lct=linect=thisline=0;
  vector <int> lnumvec;
  vector <string> inputlines;
  ifstream instream1;
  ofstream outstream1;
  
  if(argc<9) {
    show_usage();
    return(1);
  }
  
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-d" || string(argv[i]) == "-dat" || string(argv[i]) == "-data" || string(argv[i]) == "--datafile") {
      if(i+1 < argc) {
	//There is still something to read;
        datafile=argv[++i];
	i++;
      }
      else {
	cerr << "Data file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-l" || string(argv[i]) == "-lf" || string(argv[i]) == "-line"  || string(argv[i]) == "-lines"  || string(argv[i]) == "--lines" || string(argv[i]) == "--linefile" || string(argv[i]) == "-linefile") {
      if(i+1 < argc) {
	//There is still something to read;
	linefile=argv[++i];
	i++;
      }
      else {
	cerr << "Line file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-c" || string(argv[i]) == "-ch" || string(argv[i]) == "-choose" || string(argv[i]) == "-choice" || string(argv[i]) == "--choose" || string(argv[i]) == "--choice") {
      if(i+1 < argc) {
	//There is still something to read;
	choice=argv[++i];
	i++;
      }
      else {
	cerr << "Selection choice keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-of" || string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--fileout") {
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
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  if(choice!="select" && choice!="cull") {
    cout << "Invalid election choice keyword" << choice << "\n";
    cout << "must be either \"select\" or \"cull\"\n";
    show_usage();
    return(1);
  }
 
  cout << "input data file " << datafile << "\n";
  cout << "input line file " << linefile << "\n";
  cout << "Choice to select or cull the lines listed in the line file " << choice << "\n";
  cout << "output file " << outfile << "\n";

  // Read original input file into a string vector
  instream1.open(datafile,ios_base::in);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << datafile << "\n";
    return(1);
  }
  datact=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      inputlines.push_back(lnfromfile);
      datact++;
    }
  }
  datanum = datact;
  cout << "File " << datafile << " appears to have " << datanum << " lines\n";
  instream1.close();

  // Open line file
  instream1.open(linefile,ios_base::in);
  if(!instream1) {
    cerr << "can't open line file " << linefile << "\n";
    return(2);
  } else cout << "Reading line file " << linefile << "\n";
  lct=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    instream1 >> linect;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      lnumvec.push_back(linect);
      lct++;
    }
  }
  lnum = lct;
  cout << "File " << linefile << " appears to have " << lnum << " lines\n";
  instream1.close();
  
  // Open output file
  outstream1.open(outfile,ios_base::out);
  lct=0;
  datact=0;
  if(choice=="select") {
    // Print all the lines whose numbers are listed in the linefile
    for(lct=0;lct<lnum;lct++) {
      thisline = lnumvec[lct]-1;
      if(thisline>=0 && thisline<datanum) {
	outstream1 << inputlines[thisline] << "\n";
      } else {
	cerr << "WARNING: out-of-range line " << thisline << " cannot be printed\n";
      }
    }
  } else if(choice=="cull") {
    // Print all the lines whose numbers are NOT listed in the linefile
    lct=0;
    datact=0;
    while(datact<datanum) {
      if(lct<lnum) thisline = lnumvec[lct]-1; // We assume input line numbers are indexed from 1, not zero.
      if(datact<thisline) {
	while(datact<thisline && datact<datanum) {
	  outstream1 << inputlines[datact] << "\n";
	  datact++;
	}
	if(datact==thisline) datact++; 
      } else if(datact==thisline) {
	datact++; // Skip past the line
      } else if(datact>thisline && lct>=lnum) {
	// We must have gone past the last line in the line file
	// print all the rest of the lines in the data file.
	while(datact<datanum) {
	  outstream1 << inputlines[datact] << "\n";
	  datact++;
	}
      }      
      lct++;
    }
  } else {
    cout << "Invalid election choice keyword" << choice << "\n";
    cout << "must be either \"select\" or \"cull\"\n";
    show_usage();
    return(1);
  }
 
  
  outstream1.close();
  return(0);
}
