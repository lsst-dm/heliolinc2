// November 10, 2023: split_hlfile.cpp:
// Splits heliolinc output files, for use in parallel runs of
// link_refine_Herget_univar.

#include "solarsyst_dyn_geo01.h"
#include "cmath"


static void show_usage()
{
  cerr << "Usage: split_hlfile -sum summary_file -clust2det clust2detfile -splitnum number of parts -root root name for output cluster file\n";
}

int main(int argc, char *argv[])
{
  string clusterlist,clusterlist2;
  string sumfile,clust2detfile;
  string sumfile_part,clust2detfile_part;
  vector <hlclust> sumvec;
  vector <hlclust> sumvecmain;
  vector  <longpair> clust2det;
  vector  <longpair> clust2detmain;
  string imfile, pairdetfile,stest;
  string outsumfile = "LRHsumfile_test.csv";
  string outclust2detfile = "LRHclust2detfile_test.csv";
  ifstream instream1;
  ofstream outstream1;
  ofstream outstream2;
  ofstream outstream3;
  long i=0;
  int charnum;
  long clustnum=0;
  int status=0;
  long clustct=0;
  long splitnum=10;
  long splitct=0;
  long subclustnum=1000;
  long subclustct=0;
  long indexoffset=0;
  long clust2detnum=10000;
  long clust2detct=0;
  string rootname = "hltest";
  string sumsuffix = ".csv";
  string clust2detsuffix = ".csv";
  string sumroot,clust2detroot;
  string numsuffix;
  char nsfx[64];
  
  if(argc<9) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-sum" || string(argv[i]) == "-summary" || string(argv[i]) == "-rms" || string(argv[i]) == "--summary" || string(argv[i]) == "-sumfile" || string(argv[i]) == "--summaryfile" || string(argv[i]) == "--summary" || string(argv[i]) == "--sum") {
      if(i+1 < argc) {
	//There is still something to read;
	sumfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input summary file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clust2det" || string(argv[i]) == "-c2d" || string(argv[i]) == "-clust2detfile" || string(argv[i]) == "--clust2detfile" || string(argv[i]) == "--clust2det" || string(argv[i]) == "--cluster_to_detection") {
      if(i+1 < argc) {
	//There is still something to read;
	clust2detfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input cluster-to-detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-root" || string(argv[i]) == "-rootname" || string(argv[i]) == "-fileroot" || string(argv[i]) == "-rt" || string(argv[i]) == "--rootname" || string(argv[i]) == "--root") {
      if(i+1 < argc) {
	//There is still something to read;
	rootname=argv[++i];
	i++;
      }
      else {
	cerr << "Input root name eyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-splitnum" || string(argv[i]) == "-num" || string(argv[i]) == "-split" || string(argv[i]) == "-parts" || string(argv[i]) == "--splitnum" || string(argv[i]) == "--partnum" || string(argv[i]) == "--split" || string(argv[i]) == "--parts") {
      if(i+1 < argc) {
	//There is still something to read;
	splitnum=stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Split number keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  // Catch required parameters if missing
  if(sumfile.size()<=0) {
    cout << "\nERROR: input summary file is required\n";
    show_usage();
    return(1);
  } else if(clust2detfile.size()<=0) {
    cout << "\nERROR: input cluster-to-detection file is required\n";
    show_usage();
    return(1);
  } 

  cout << "input summary file " << sumfile << "\n";
  cout << "input cluster-to-detection file " << clust2detfile << "\n";
  cout << "Files will be split into " << splitnum << "parts\n";
  cout << "Root name for output lists: " << rootname << "\n";

  charnum = sumfile.size();
  sumroot = sumfile.substr(0,charnum-4);
  sumsuffix = sumfile.substr(charnum-4,4);
  charnum = clust2detfile.size();
  clust2detroot = clust2detfile.substr(0,charnum-4);
  clust2detsuffix = clust2detfile.substr(charnum-4,4);

  cout << "Root names for summary and clust2det files: " << sumroot << " and " << clust2detroot << "\n";
  cout << "Suffixes for summary and clust2det files: " << sumsuffix << " and " << clust2detsuffix << "\n";
  
  // Read cluster summary file
  sumvecmain={};
  status=read_clustersum_file(sumfile, sumvecmain, 0);
  if(status!=0) {
    cerr << "ERROR: could not successfully read input cluster summary file " << sumfile << "\n";
    cerr << "read_clustersum_file returned status = " << status << ".\n";
    return(1);
  }
  clustnum = sumvecmain.size();
  cout << "Read " << clustnum << " data lines from cluster summary file " << sumfile << "\n";

  clust2detmain={};
  // Read cluster-to-detection file
  status=read_longpair_file(clust2detfile, clust2detmain, 0);
  if(status!=0) {
    cerr << "ERROR: could not successfully read cluster-to-detection file " << clust2detfile << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
    return(1);
  }
  clust2detnum = clust2detmain.size();
  cout << "Read " << clust2detnum << " data lines from cluster-to-detection file " << clust2detfile << "\n";

  subclustnum = clustnum/splitnum + 1;
  clustct=0;
  clust2detct=0;
  clusterlist2 = "clusterlist_" + rootname + "_main";
  outstream3.open(clusterlist2);
  
  for(splitct=0; splitct<splitnum; splitct++) {
    i=sprintf(nsfx,"%ld",splitct);
    numsuffix = nsfx;
    if(i<=1) numsuffix = "00" + numsuffix;
    else if(i==2) numsuffix = "0" + numsuffix;
    sumfile_part = sumroot + numsuffix + sumsuffix;
    clust2detfile_part = clust2detroot + numsuffix + clust2detsuffix;

    clusterlist = "clusterlist_" + rootname + numsuffix;
    outstream1.open(clusterlist);
    outstream1 << sumfile_part << " " << clust2detfile_part << "\n";
    outstream1.close();

    outstream3 << "LRH" << sumfile_part << " LRH" << clust2detfile_part << "\n";

    outstream1.open(sumfile_part);
    outstream2.open(clust2detfile_part);
    indexoffset=clustct;
    cout << "Writing output cluster-summary file " << sumfile_part << "\n";
    cout << "And output clust2det file " << clust2detfile_part << "\n";
    outstream1 << "#clusternum,posRMS,velRMS,totRMS,astromRMS,pairnum,timespan,uniquepoints,obsnights,metric,rating,heliohyp0,heliohyp1,heliohyp2,posX,posY,posZ,velX,velY,velZ,orbit_a,orbit_e,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count\n";
    outstream2 << "#clusternum,detnum\n";
    
    for(subclustct=0;subclustct<subclustnum;subclustct++) {
      if(clustct<clustnum) {
	outstream1 << fixed << setprecision(3) << sumvecmain[clustct].clusternum-indexoffset << "," << sumvecmain[clustct].posRMS << "," << sumvecmain[clustct].velRMS << "," << sumvecmain[clustct].totRMS << ",";
	outstream1 << fixed << setprecision(4) << sumvecmain[clustct].astromRMS << ",";
	outstream1 << fixed << setprecision(6) << sumvecmain[clustct].pairnum << "," << sumvecmain[clustct].timespan << "," << sumvecmain[clustct].uniquepoints << "," << sumvecmain[clustct].obsnights << "," << sumvecmain[clustct].metric << "," << sumvecmain[clustct].rating << ",";
	outstream1 << fixed << setprecision(6) << sumvecmain[clustct].heliohyp0 << "," << sumvecmain[clustct].heliohyp1 << "," << sumvecmain[clustct].heliohyp2 << ",";
	outstream1 << fixed << setprecision(1) << sumvecmain[clustct].posX << "," << sumvecmain[clustct].posY << "," << sumvecmain[clustct].posZ << ",";
	outstream1 << fixed << setprecision(4) << sumvecmain[clustct].velX << "," << sumvecmain[clustct].velY << "," << sumvecmain[clustct].velZ << ",";
	outstream1 << fixed << setprecision(6) << sumvecmain[clustct].orbit_a << "," << sumvecmain[clustct].orbit_e << "," << sumvecmain[clustct].orbit_MJD << ",";
	outstream1 << fixed << setprecision(1) << sumvecmain[clustct].orbitX << "," << sumvecmain[clustct].orbitY << "," << sumvecmain[clustct].orbitZ << ",";
	outstream1 << fixed << setprecision(4) << sumvecmain[clustct].orbitVX << "," << sumvecmain[clustct].orbitVY << "," << sumvecmain[clustct].orbitVZ << "," << sumvecmain[clustct].orbit_eval_count << "\n";
	while(clust2detmain[clust2detct].i1==sumvecmain[clustct].clusternum && clust2detct<clust2detnum) {
	  outstream2 << clust2detmain[clust2detct].i1-indexoffset << "," << clust2detmain[clust2detct].i2 << "\n";
	  clust2detct++;
	}
	clustct++;
      }
    }
    outstream1.close();
    outstream2.close();
  }
  outstream3.close();
  
  return(0);
}
