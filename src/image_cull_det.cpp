// September 20, 2022: image_cull_deg.cpp:
// Uses much of the machinery from make_tracklets to cull an input
// data set in a manner that establishes a maximum number
// of detections per image. These are supposed to be the 'best'
// detections based on a metric supplied in the input file.
#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define NUMPOS 3

#define IDCOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define MAGCOL 5
#define BANDCOL 6
#define OBSCODECOL 7
#define COLS_TO_READ 7
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME (1.5/24.0) // Default max inter-image time interval
                           // for tracklets, in days.
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)
#define MAX_GCR 0.5 // Default maximum Great Circle Residual allowed for a valid tracklet
#define DEBUG 0
#define DEBUGA 0

      
static void show_usage()
{
  cerr << "Usage: image_cull_det -dets detfile -out outfile -colformat column format file ";
  cerr << "-meritcol column number holding the merit fuction -percut percentile in image detections to use for a maximum\n\nOR\n\n";
  cerr << "image_cull_det -dets detfile -out outfile -colformat column format file ";
  cerr << "-meritcol column number holding the merit fuction -maxdet maximum number of detections allowed per image\n";
}
    
int main(int argc, char *argv[])
{
  det_obsmag_indvec o1 = det_obsmag_indvec(0L,0l,0l,0L,0L,0L,"null",0l,"V","I11",0,{});
  vector <det_obsmag_indvec> detvec = {};
  vector <det_obsmag_indvec> ppset = {};
  char idstring[SHORTSTRINGLEN];
  char band[MINSTRINGLEN];
  char obscode[MINSTRINGLEN];
  string lnfromfile;
  int status = 0;
  int i = 0;
  int j = 0;
  int imct=0;
  int imctp=0;
  int imnum=0;
  long detnum=0;
  long num_dets=0;
  long detct=0;
  int startind=0;
  int endind=0;
  int reachedeof = 0;
  char c='0';
  long double MJD,RA,Dec,merit;
  MJD = RA = Dec = merit = 0.0L;
  double mag = 0l;
  string indetfile;
  string outfile;
  string colformatfile;
  long lct=0;
  img_log03 imlog = img_log03(0.0,0.0,0.0,"I11",0,0);
  vector <img_log03> img_log_tmp = {};
  vector <img_log03> img_log = {};
  long_index ppn = long_index(0,0);
  vector <long_index> pair_partner_num={};
  vector <long_index> tracklet_check={};
  int idcol=IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int magcol = MAGCOL;
  int bandcol = BANDCOL;
  int obscodecol = OBSCODECOL;
  int meritcol = 0;
  int colreadct=0;
  double maxdet=0l;
  double percut=0l;
  ifstream instream1;
  ofstream outstream1;
  string stest;
  vector <long> dets_per_image;
  vector <long> keep_index;
  ldouble_index ldi = ldouble_index(0,0.0L);
  vector <ldouble_index> metricvec;
  int small_image_ct=0;
  int big_image_ct=0;
  int good_detnum,good_detct;
  long double mjdmean,mjdnorm,tdelt;
  
  if(argc!=11)
    {
      show_usage();
      return(1);
    }
  
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
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-of") {
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
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format"  || string(argv[i]) == "-col" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" || string(argv[i]) == "--cformat") {
      if(i+1 < argc) {
	//There is still something to read;
	colformatfile=argv[++i];
	i++;
      }
      else {
	cerr << "Column format file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-meritcol") {
      if(i+1 < argc) {
	//There is still something to read;
	meritcol=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Merit column keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-percut") {
      if(i+1 < argc) {
	//There is still something to read;
	percut=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Percentile cutoff keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxdet") {
      if(i+1 < argc) {
	//There is still something to read;
	maxdet=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Maximum detections keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }

  if(indetfile.size()<=0) {
    cerr << "Please supply an input detection file:\n\n";
    show_usage();
    return(1);
  }
  if(meritcol<=0) {
    cerr << "Please specify the column number for the merit function\n";
    cerr << "in your input file (must be an integer greater than 0)\n";
    show_usage();
    return(1);
  }
  if(maxdet>0.0) cout << "A maximum of " << maxdet << " detections will be retained for each image\n";
  else if(percut>0.0) {
    if(percut>1.0l) percut/=100.0l;
    cout << "The number of detections per image will be thresholded at the\n";
    cout << percut*100.0 << "percentile for the raw distribution of detections per image\n";
  } else if(maxdet<=0.0l && percut<=0.0l) {
    cerr << "Error: you must specify either a maximum detection number\n";
    cerr << "or percentile threshold for detections per image\n";
    show_usage();
    return(1);
  }
  cout << "\nInput detection file is called " << indetfile << "\n\n";

  if(colformatfile.size()<=0) {
    cerr << "ERROR: no column formatting file specified\n";
    show_usage();
    return(1);
  }
  cout << "Output file will be called " << outfile << "\n";
  
  // Read the column formatting file, if any
  if(colformatfile.size()>0)
    {
      instream1.open(colformatfile);
      if(!instream1)  {
	cerr << "ERROR: unable to open input file " << colformatfile << "\n";
	return(1);
      }
      colreadct=0;
      while(!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct<COLS_TO_READ) {
	instream1 >> stest;
	if (stest == "IDCOL") {
	  instream1 >> idcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "MJDCOL") {
	  instream1 >> mjdcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "RACOL") {
	  instream1 >> racol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "DECCOL") {
	  instream1 >> deccol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "MAGCOL") {
	  instream1 >> magcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "BANDCOL") {
	  instream1 >> bandcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "OBSCODECOL") {
	  instream1 >> obscodecol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else {
	  cout << "WARNING: unrecognized string " << stest << " read from column formatting file\n";
	}
      }
      instream1.close();
      if(colreadct<COLS_TO_READ) {
	cout << "WARNING: only " << colreadct << " column specifications, of " << COLS_TO_READ << " expected, were read from column format file " << colformatfile << ".\n";
      }
    }

  cout << "Column specifications:\n";
  cout << "IDCOL " << idcol << "\n";
  cout << "MJDCOL " << mjdcol << "\n";
  cout << "RACOL " << racol << "\n";
  cout << "DECCOL " << deccol << "\n";
  cout << "MAGCOL " << magcol << "\n";
  cout << "BANDCOL " << bandcol << "\n";
  cout << "OBSCODECOL " << obscodecol << "\n";
  cout << "meritcol " << meritcol << "\n";
  
  // Read input detection file.
  instream1.open(indetfile);
  if(!instream1) {
    cerr << "can't open input file " << indetfile << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  lct++;
  //cout << lnfromfile << "\n";
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    lct++;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    i=0;
    j = 0;
    c='0';
    while(i<lnfromfile.size() && lnfromfile.size()>=30 && reachedeof == 0) {
      // Note check on line length: it is completely impossible for a
      // line containing all the required quantities at minimum plausible
      // precision to be less than 30 characters long.
      c='0';
      stest="";
      while(i<lnfromfile.size() && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==idcol) stringncopy01(idstring,stest,SHORTSTRINGLEN);
      if(j==mjdcol) MJD=stold(stest);
      else if(j==racol) RA=stold(stest);
      else if(j==deccol) Dec=stold(stest);
      else if(j==magcol) mag=stod(stest);
      else if(j==bandcol) stringncopy01(band,stest,MINSTRINGLEN);
      else if(j==obscodecol) stringncopy01(obscode,stest,MINSTRINGLEN);
      else if(j==meritcol) merit=stold(stest);
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0 && lnfromfile.size()>=30) {
      // cout<<"MJD, RA, Dec: " << MJD-floor(MJD) << " " << RA << " " << Dec << "\n";
      o1=det_obsmag_indvec(MJD,RA,Dec,merit,0L,0L,idstring,mag,band,obscode,lct,{});
      detvec.push_back(o1);
    }
  }
  instream1.close();
  if(reachedeof==1) { 
    cout << "Input file " << indetfile << " read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";

  // time-sort the detection vector
  sort(detvec.begin(), detvec.end(), early_det_obsmag_indvec());
  
  // Determine the number of detections per image
  dets_per_image = {};
  mjdnorm = 1.0;
  mjdmean = detvec[0].MJD;
  startind=0;
  for(i=1;i<detvec.size();i++) {
    tdelt = detvec[i].MJD - detvec[i-1].MJD;
    if(tdelt < IMAGETIMETOL/SOLARDAY && stringnmatch01(detvec[i].obscode,detvec[i-1].obscode,3)==0) {
      //This point corresponds to the same image as the previous one.
      mjdmean += detvec[i].MJD;
      mjdnorm += 1.0;
    }
    else {
      //Now we are considering a new image.
      //Calculate the meanmjd of the previous image, for which
      // we have now seen all points.
      //Record the current detct i as the detection index just
      //after the end of the previous image
      endind=i;
      if(isnormal(mjdnorm)) mjdmean /= mjdnorm;
      else mjdmean = 0.0;
      //Load it into the vector with mean MJD for all images,
      // and increment image count.
      imlog = img_log03(mjdmean,0.0,0.0,detvec[endind-1].obscode,startind,endind);
      dets_per_image.push_back(endind-startind);
      img_log.push_back(imlog);
      // Set up for the next image, starting with detvec[i].MJD;
      mjdmean = detvec[i].MJD;
      mjdnorm = 1.0;
      startind=i;
    }
  }
  // Account for the final image.
  if(isnormal(mjdnorm)) {
    endind=i;
    mjdmean /= mjdnorm;
    //Load it into the vector with mean MJD for all images,
    // and increment image count.
    imlog = img_log03(mjdmean,0.0,0.0,detvec[endind-1].obscode,startind,endind);
    dets_per_image.push_back(endind-startind);
    img_log.push_back(imlog);
  }

  detnum = detvec.size();
  imnum = img_log.size();
  cout << img_log.size() << " unique images were identified.\n";
  cout << "Given our total of " << detvec.size() << " detections,\n";
  cout << "we have " << double(detvec.size())/double(img_log.size()) << " detections per image, on average\n";

  sort(dets_per_image.begin(),dets_per_image.end());
  j = dets_per_image.size()/2;
  cout << "Of " << dets_per_image.size() << " images, the median number of detections per image is " << dets_per_image[j] << "\n";
  cout << "The minimum number of detections on an image is " << dets_per_image[0] << "\n";
  cout << "and the maximum number of detections on an image is " << dets_per_image[dets_per_image.size()-1] << "\n";
  if(maxdet<=0.0) {
    j=percut*dets_per_image.size();
    maxdet = dets_per_image[j];
    cout << "The " << percut*100.0l << " percentile in detections per image is " << maxdet << "\n";
  }

  // Pass through the detection catalog again, writing out
  // indices of detections to be kept.
  keep_index={};
  small_image_ct = big_image_ct = 0;
  for(imct=0;imct<imnum;imct++) {
    num_dets = img_log[imct].endind - img_log[imct].startind;
    if(num_dets<=maxdet) {
      // There aren't too many detections on this image: keep them all.
      small_image_ct++;
      for(i=img_log[imct].startind;i<img_log[imct].endind;i++) {
	keep_index.push_back(detvec[i].index);
      }
    } else {
      // There are too many detections on this image. We will retain only
      // the maxdet best ones.
      big_image_ct++;
      metricvec={};
      for(i=img_log[imct].startind;i<img_log[imct].endind;i++) {
	ldi = ldouble_index(detvec[i].x,detvec[i].index);
	metricvec.push_back(ldi);
      }
      sort(metricvec.begin(), metricvec.end(), lower_ldouble_index());
      for(i=metricvec.size()-1;i>=metricvec.size()-maxdet;i--) {
	keep_index.push_back(metricvec[i].index);
      }
    }
  }
  // Now the vector keep_index should have the indices of all the
  // detections good enough to keep.
  cout << "Of a total " << imnum << " images, " << small_image_ct << " had acceptable numbers of detections,\n";
  cout << "while " << big_image_ct << " had excessive numbers of detections\n";
  cout << "Of " << detvec.size() << " detections, " << keep_index.size() << " will be retained and the rest rejected\n";

  // Sort keep_index
  sort(keep_index.begin(),keep_index.end());
  
  // Read input detection file again, outputting all the good lines.
  instream1.open(indetfile);
  outstream1.open(outfile);
  if(!instream1) {
    cerr << "can't open input file " << indetfile << "\n";
    return(1);
  }
  // Read one-line header
  good_detnum = keep_index.size();
  reachedeof=0;
  good_detct = 0;
  lct=0;
  getline(instream1,lnfromfile);
  // Copy header to output file
  outstream1 << lnfromfile << "\n";
  lct++;
  while(reachedeof==0 && good_detct < good_detnum) {
    if(lct<keep_index[good_detct]) {
      getline(instream1,lnfromfile);
      lct++;
    } else if(lct==keep_index[good_detct]) {
      outstream1 << lnfromfile << "\n";
      getline(instream1,lnfromfile);
      lct++;
      good_detct++;
    }
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
  }
  instream1.close();
  outstream1.close();
  if(reachedeof==1 || good_detct==good_detnum) { 
    cout << "Input file " << indetfile << " read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  return(0);
}
