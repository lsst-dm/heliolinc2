// February 04, 2022: stringcm01.cpp:
// Given two files, each with a string column,
// write a new file consisting either of the
// lines from the two files where the strings matched,
// or the lines from the first file that had no counterpart
// in the second.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: stringcm01 file1 col1 file2 col2 match outfile\nOR:\nstringcm01 file1 col1 file2 col2 unmatch outfile\n";
}
    
int main(int argc, char *argv[])
{
  string file1,file2,mstring,s1,lnfromfile,idstring,outfile;
  int col1,col2;
  long linect1,linenum1,linect2,linenum2;
  linect1 = linect2 = linenum1 = linenum2 = 0;
  string_index si1 = string_index("",0);
  vector <string_index> lines01;
  vector <string_index> lines02;
  vector <string_index> stringid01;
  vector <string_index> stringid02;
  ifstream instream1;
  ofstream outstream1;
  int i=0;
  
  if(argc!=7)
    {
      show_usage();
      return(1);
    }

  file1=argv[1];
  col1 = stoi(argv[2]);
  file2=argv[3];
  col2 = stoi(argv[4]);
  mstring = argv[5];
  outfile = argv[6];

  cout << "input files: " << file1 << ", column " << col1 << "; " << file2 << ", column " << col2 <<"\n";
  cout << "output file: " << outfile << "\n";
  if(mstring=="match") {
    cout << "Matching lines will be output\n";
  } else if(mstring=="unmatch") {
    cout << "Lines from file " << file1 << " with no match in file " << file2 << " will be output.\n";
  } else {
    cout << "ERROR: unrecognized match specifier " << mstring << "\n";
    return(2);
  }
  
  // Read file 1 for lines.
  instream1.open(file1);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << file1 << "\n";
    return(1);
  }
  linect1=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      si1 = string_index(lnfromfile,linect1);
      lines01.push_back(si1);
      linect1++;
    }
  }
  linenum1 = linect1;
  cout << "File " << file1 << " appears to have " << linenum1 << " lines\n";
  instream1.close();
  // Read file 1 for string ID column
  instream1.open(file1);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << file1 << "\n";
    return(1);
  }
  linect1=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    // Skip lines up to the reference column
    for(i=1;i<col1;i++) instream1 >> s1;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      // Read reference column
      instream1 >> idstring;
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
	si1 = string_index(idstring,linect1);
	stringid01.push_back(si1);
	linect1++;
	// Read and discard remainder of the line
	getline(instream1,lnfromfile);
      }
    }
  }
  instream1.close();
  cout << "Vector sizes from file1: " << lines01.size() << " and " << stringid01.size() << "\n"; 
  
  // Read file 2 for lines.
  instream1.open(file2);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << file2 << "\n";
    return(1);
  }
  linect2=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      si1 = string_index(lnfromfile,linect2);
      lines02.push_back(si1);
      linect2++;
    }
  }
  linenum2 = linect2;
  cout << "File " << file2 << " appears to have " << linenum2 << " lines\n";
  instream1.close();
  
  // Read file 2 for string ID column
  instream1.open(file2);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << file2 << "\n";
    return(1);
  }
  linect2=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    // Skip lines up to the reference column
    for(i=1;i<col2;i++) instream1 >> s1;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      // Read reference column
      instream1 >> idstring;
      if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
	si1 = string_index(idstring,linect2);
	stringid02.push_back(si1);
	linect2++;
	// Read and discard remainder of the line
	getline(instream1,lnfromfile);
      }
    }
  }
  cout << "Vector sizes from file2: " << lines02.size() << " and " << stringid02.size() << "\n"; 
  instream1.close();

  // Sort both string ID vectors
  sort(stringid01.begin(), stringid01.end(), lower_string_index());
  sort(stringid02.begin(), stringid02.end(), lower_string_index());

  outstream1.open(outfile);
  if(mstring=="match") {
    // Identify and output matching lines
    linect1=linect2=0;
    while(linect1<linenum1 && linect2<linenum2) {
      if(stringid01[linect1].selem==stringid02[linect2].selem) {
	outstream1 << lines01[stringid01[linect1].index].selem << " | " << lines02[stringid02[linect2].index].selem << "\n";
	linect1++;
	linect2++;
      } else if(stringid01[linect1].selem < stringid02[linect2].selem) {
	// Increment linect1 to move forward and find a file1 entry that
	// will match the current one in file2.
	linect1++;
      } else if(stringid01[linect1].selem > stringid02[linect2].selem) {
	// Increment linect2 to move forward and find a file2 entry that
	// will match the current one in file1.
	linect2++;
      } else {
	cout << "WARNING: logically impossible case comparing " << stringid01[linect1].selem << " and " << stringid02[linect2].selem << "\n";
      }
    }
  } else if(mstring=="unmatch") {
    // Identify and output lines in the first file that
    // have no counterpart in the second.
    linect1=linect2=0;
    while(linect1<linenum1 && linect2<linenum2) {
      if(stringid01[linect1].selem==stringid02[linect2].selem) {
	// These lines match, so we are not interested in
	// either one. Move forward in both files
	linect1++;
	linect2++;
      } else if(stringid01[linect1].selem < stringid02[linect2].selem) {
	// File2 is already past this element of file1: there cannot
	// be a match in file2. Write out the line from file 1, and
	// advance to the next line
	outstream1 << lines01[stringid01[linect1].index].selem << "\n";
	linect1++;
      } else if(stringid01[linect1].selem > stringid02[linect2].selem) {
	// Increment linect2 to move forward and see if there is a file2 entry that
	// will match the current one in file1.
	linect2++;
      } else {
	cout << "WARNING: logically impossible case comparing " << stringid01[linect1].selem << " and " << stringid02[linect2].selem << "\n";
      }
    }
  }
  outstream1.close();
  return(0);
}
