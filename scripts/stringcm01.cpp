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
  int badlines01=0;
  int badlines02=0;
  int stringnum=0;
  vector <string> outstrings;
  int dupline1=0;
  int dupline2=0;
  
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
  
  // Read file 1.
  instream1.open(file1);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << file1 << "\n";
    return(1);
  }
  linect1=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      outstrings={};
      stringnum = stringline01(lnfromfile, outstrings);
      if(stringnum>=col1) {
	// Load whole line into vector lines01
	si1 = string_index(lnfromfile,linect1);
	lines01.push_back(si1);
	// Load string ID column into vector stringid01
	idstring = outstrings[col1-1];
	si1 = string_index(idstring,linect1);
	stringid01.push_back(si1);
	// Increment line count
	linect1++;
      } else badlines01++;
    }
  }
  instream1.close();

  linenum1 = linect1;
  cout << "File " << file1 << " appears to have " << linenum1 << " good lines\n";
  if(badlines01>0) {
    cerr << "WARNING: File " << file1 << " also has " << badlines01 << " bad lines, which were not read\n";
  }

  cout << "Vector sizes from file1: " << lines01.size() << " and " << stringid01.size() << "\n"; 
  
  // Read file 2.
  instream1.open(file2);
  if(!instream1)  {
    cerr << "ERROR: unable to open input file " << file2 << "\n";
    return(1);
  }
  linect2=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      outstrings={};
      stringnum = stringline01(lnfromfile, outstrings);
      if(stringnum>=col2) {
	si1 = string_index(lnfromfile,linect2);
	lines02.push_back(si1);
	// Load string ID column into vector stringid02
	idstring = outstrings[col2-1];
	si1 = string_index(idstring,linect2);
	stringid02.push_back(si1);
	// Increment line count
	linect2++;
      } else badlines02++;
    }
  }
  instream1.close();
  linenum2 = linect2;
  cout << "File " << file2 << " appears to have " << linenum2 << " good lines\n";
  if(badlines02>0) {
    cerr << "WARNING: File " << file2 << " also has " << badlines02 << " bad lines, which were not read\n";
  }
  cout << "Vector sizes from file2: " << lines02.size() << " and " << stringid02.size() << "\n"; 

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
	// See if we have repeated lines in either file.
	dupline1=dupline2=0;
	if(linect1+1<linenum1 && stringid01[linect1+1].selem==stringid01[linect1].selem) dupline1=1;
	if(linect2+1<linenum2 && stringid02[linect2+1].selem==stringid02[linect2].selem) dupline2=1;
	if(dupline1==1 && dupline2==0) {
	  // Increment linect1 only, so that all of the duplicated
	  // lines in file1 can be matched with the same, non-duplicated line
	  // in file2.
	  linect1++;
	} else if(dupline1==0 && dupline2==1) {
	  // Incement linect2 only, so that all of the duplicated lines in file2
	  // can be matched with the same, non-duplicated line in file1.
	  linect2++;
	} else {
	  // Either neither file has duplicated lines or both do: increment
	  // linect in both files
	  linect1++;
	  linect2++;
	}
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
	// either one.
	// See if we have repeated lines in either file.
	dupline1=dupline2=0;
	if(linect1+1<linenum1 && stringid01[linect1+1].selem==stringid01[linect1].selem) dupline1=1;
	if(linect2+1<linenum2 && stringid02[linect2+1].selem==stringid02[linect2].selem) dupline2=1;
	if(dupline1==1 && dupline2==0) {
	  // Increment linect1 only, so that all of the duplicated
	  // lines in file1 can be matched with the same, non-duplicated line
	  // in file2.
	  linect1++;
	} else if(dupline1==0 && dupline2==1) {
	  // Incement linect2 only, so that all of the duplicated lines in file2
	  // can be matched with the same, non-duplicated line in file1.
	  linect2++;
	} else {
	  // Either neither file has duplicated lines or both do: increment
	  // linect in both files
	  linect1++;
	  linect2++;
	}
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
    if(linect1<linenum1 && linect2>=linenum2) {
      // See if there are additional unmatched lines in file 1
      while(linect1<linenum1) {
	// Since there are no lines left in file2, all remaining
	// lines in file 1 must be unmatched.
	// Write them all out.
	outstream1 << lines01[stringid01[linect1].index].selem << "\n";
	linect1++;
      }
    }
  }
  outstream1.close();
  return(0);
}
