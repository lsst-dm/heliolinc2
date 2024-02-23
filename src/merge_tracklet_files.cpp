// February 22, 2024: merge_tracklet_files.cpp:
// Given output from two different runs of make_tracklets_new
// and/or make_trailed_tracklets, merge the files.
// Note well that this code only operates on three of the four
// output files of the make_tracklets programs: the fourth,
// the output image catalog, must have the same order and
// number of entries for both the two make_tracklets outputs
// being combined.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: merge_tracklet_files -inputs1 paired_detection_file#1 tracklet_file#1 tracklet-to-detection_file#1 -inputs2 paired_detection_file#2 tracklet_file#2 tracklet-to-detection_file#2 -outputs output_paired_detection_file output_tracklet_file output_tracklet-to-detection_file\n";
}
    
int main(int argc, char *argv[])
{
  string pairdetfile_in1,trackletfile_in1,trk2detfile_in1;
  string pairdetfile_in2,trackletfile_in2,trk2detfile_in2;
  string pairdetfile_out,trackletfile_out,trk2detfile_out;
  vector <hldet> pairdets_in1;
  vector <tracklet> tracklets_in1;
  vector <longpair> trk2det_in1;
  vector <hldet> pairdets_in2;
  vector <tracklet> tracklets_in2;
  vector <longpair> trk2det_in2;
  vector <hldet> pairdets;
  vector <tracklet> tracklets;
  vector <longpair> trk2det;
  tracklet onetrack;
  longpair onepair;
  int verbose=0;
  ofstream outstream1;
  int status=0;

  if(argc!=13) {
    cout << "Need argc=13, got " << argc << "\n";
    show_usage();
    return(1);
  }
  
  int i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-inputs1" || string(argv[i]) == "-input1" || string(argv[i]) == "-in1" || string(argv[i]) == "-inputs01") {
      if(i+3 < argc) {
	//There is still something to read;
        pairdetfile_in1=argv[++i];
        trackletfile_in1=argv[++i];
        trk2detfile_in1=argv[++i];
	i++;
      }
      else {
	cerr << "Input1 keyword supplied with fewer than the required three arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-inputs2" || string(argv[i]) == "-input2" || string(argv[i]) == "-in2" || string(argv[i]) == "-inputs02") {
      if(i+3 < argc) {
	//There is still something to read;
        pairdetfile_in2=argv[++i];
        trackletfile_in2=argv[++i];
        trk2detfile_in2=argv[++i];
	i++;
      }
      else {
	cerr << "Input2 keyword supplied with fewer than the required three arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outputs" || string(argv[i]) == "-output" || string(argv[i]) == "-outfile" || string(argv[i]) == "-outfiles" || string(argv[i]) == "-out") {
      if(i+3 < argc) {
	//There is still something to read;
        pairdetfile_out=argv[++i];
        trackletfile_out=argv[++i];
        trk2detfile_out=argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with fewer than the required three arguments\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  // Read first set of input files
  status=read_pairdet_file(pairdetfile_in1, pairdets_in1 , verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << pairdetfile_in1 << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << pairdets_in1.size() << " data lines from paired detection file " << pairdetfile_in1 << "\n";
    
  status=read_tracklet_file(trackletfile_in1, tracklets_in1, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read tracklet file " << trackletfile_in1 << "\n";
    cerr << "read_tracklet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << tracklets_in1.size() << " data lines from tracklet file " << trackletfile_in1 << "\n";
  
  status=read_longpair_file(trk2detfile_in1, trk2det_in1, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read trk2det file " << trk2detfile_in1 << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << trk2det_in1.size() << " data lines from trk2det file " << trk2detfile_in1 << "\n";
  
  // Read second set of input files
  status=read_pairdet_file(pairdetfile_in2, pairdets_in2 , verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << pairdetfile_in2 << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << pairdets_in2.size() << " data lines from paired detection file " << pairdetfile_in2 << "\n";
    
  status=read_tracklet_file(trackletfile_in2, tracklets_in2, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read tracklet file " << trackletfile_in2 << "\n";
    cerr << "read_tracklet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << tracklets_in2.size() << " data lines from tracklet file " << trackletfile_in2 << "\n";
  
  status=read_longpair_file(trk2detfile_in2, trk2det_in2, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read trk2det file " << trk2detfile_in2 << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << trk2det_in2.size() << " data lines from trk2det file " << trk2detfile_in2 << "\n";

  // Copy first set of input files directly to ouput files:
  pairdets = pairdets_in1;
  tracklets = tracklets_in1;
  trk2det = trk2det_in1;

  long numdets = pairdets.size();
  long numtrack = tracklets.size();

  // Concatenate the new files onto the old
  for(long j=0; j<long(pairdets_in2.size()); j++) {
    pairdets.push_back(pairdets_in2[j]);
  }
  for(long j=0; j<long(tracklets_in2.size()); j++) {
    onetrack = tracklets_in2[j];
    onetrack.trk_ID += numtrack;
    tracklets.push_back(onetrack);
  }
  for(long j=0; j<long(trk2det_in2.size()); j++) {
    onepair = trk2det_in2[j];
    onepair.i1 += numtrack;
    onepair.i2 += numdets;
    trk2det.push_back(onepair);
  }
  
  // Write paired detection file
  cout << "Writing paired detection file with " << pairdets.size() << " lines\n";
  outstream1.open(pairdetfile_out);
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(i=0;i<long(pairdets.size());i++) {
    outstream1 << fixed << setprecision(7) << pairdets[i].MJD << "," << pairdets[i].RA << "," << pairdets[i].Dec << ",";
    outstream1 << fixed << setprecision(4) << pairdets[i].mag << ",";
    outstream1 << fixed << setprecision(2) << pairdets[i].trail_len << "," << pairdets[i].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << pairdets[i].sigmag << ",";
    outstream1 << fixed << setprecision(3) << pairdets[i].sig_across << "," << pairdets[i].sig_along << ",";
    outstream1 << pairdets[i].image << "," << pairdets[i].idstring << "," << pairdets[i].band << ",";
    outstream1 << pairdets[i].obscode << "," << pairdets[i].known_obj << ","; 
    outstream1 << pairdets[i].det_qual << "," << pairdets[i].index << "\n"; 
  }
  outstream1.close();

  // Write tracklet file
  cout << "Writing tracklet file with " << tracklets.size() << " lines\n";
  outstream1.open(trackletfile_out);
  outstream1 << "#Image1,RA1,Dec1,Image2,RA2,Dec2,npts,trk_ID\n";
  for(i=0;i<long(tracklets.size());i++) {
    outstream1 << fixed << setprecision(7) << tracklets[i].Img1 << "," << tracklets[i].RA1 << "," << tracklets[i].Dec1 << ",";
    outstream1 << fixed << setprecision(7) << tracklets[i].Img2 << "," << tracklets[i].RA2 << "," << tracklets[i].Dec2 << ",";
    outstream1 << tracklets[i].npts << "," << tracklets[i].trk_ID << "\n"; 
  }
  outstream1.close();

  // Write trk2det file
  cout << "Writing trk2det file with " << trk2det.size() << " lines\n";
  outstream1.open(trk2detfile_out);
  outstream1 << "#trk_ID,detnum\n";
  for(i=0;i<long(trk2det.size());i++) {
    outstream1 << trk2det[i].i1 << "," << trk2det[i].i2 << "\n"; 
  }
  outstream1.close();

  return(0);
}
