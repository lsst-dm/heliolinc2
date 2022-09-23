# Heliolinc3D C++ Implementation #

## Obtain the following files from github: ##

#### Source code: ####

solarsyst_dyn_geo01.cpp <br>
solarsyst_dyn_geo01.h <br>
make_tracklets.cpp <br>
heliolinc.cpp <br>
link_refine.cpp <br>
link_refine_multisite.cpp <br>

#### General-use input files: ####


Earth1day2020s_02a.txt <br>
ObsCodes.txt <br>
colformat_LSST_01.txt <br>
accelmat_mb08a_sp04.txt <br>

#### Test data files: ####

LSST_raw_input_data01a.csv <br>
LSST_pairdets_01_check.csv <br>
LSST_pairs_01_check <br>
LSST_linkref_all01_check.csv <br>
LSST_linkref_summary01_check.csv <br>


## Compile the heliolinc suite ##

#### Use these suggested compile commands, or their equivalents for your system: ####

`c++ -O3 -o heliolinc heliolinc.cpp solarsyst_dyn_geo01.cpp -std=c++11` <br>
`c++ -O3 -o make_tracklets make_tracklets.cpp solarsyst_dyn_geo01.cpp -std=c++11` <br>
`c++ -O3 -o link_refine_multisite link_refine_multisite.cpp solarsyst_dyn_geo01.cpp -std=c++11` <br>
`c++ -O3 -o link_refine link_refine.cpp solarsyst_dyn_geo01.cpp -std=c++11` <br>


## What the programs do: ##

**make_tracklets:** Perform image-based pairing to create input pair/tracklet files for heliolinc

**heliolinc:** link together pairs/tracklets produced by make_trackelts into candidate asteroid discoveries

**link_refine and link_refine_multisite:** post-process linkages produced by heliolinc to produce a de-duplicated, non-overlapping set of candidate discoveries with the highest likelihood of being real.

## Testing your installation: ##

### Testing make_tracklets ###

If you've downloaded the recommended test data files, you can test your installation immediately. Here is the minimalist invocation for `make_tracklets`, the first program in the suite:

`make_tracklets -dets LSST_raw_input_data01a.csv -earth Earth1day2020s_02a.txt -obscode ObsCodes.txt -colformat colformat_LSST_01.txt`

**This run should take about 10 seconds.**

If you got a final message stating, "Constructing tracklets, and writing pairs to output file," the code probably executed correctly. It will have produced output files called `pairdetfile01.csv` and `outpairfile01`. The names of these output files can be set with command-line options, but in the minimalist invocation they default to the names given here.

`pairdetfile01.csv` is the **paired detection file**, which is just a reformatted version of the input detection catalog limited only to detections that formed pairs or longer tracklets.

`outpairfile01` is the **pair file**, which records the pairs and longer tracklets that were found using integer indices that specify their position in the paired detection file.

If the minimalist invocation failed, check that you have the input files used in the example:

LSST_raw_input_data01a.csv <br>
Earth1day2020s_02a.txt <br>
ObsCodes.txt <br>
colformat_LSST_01.txt <br>

If the minimalist invocation succeeded, you can try a new run in which you specify the output names:

`make_tracklets -dets LSST_raw_input_data01a.csv -pairdets LSST_pairdets_01.csv -pairs LSST_pairs_01 -earth Earth1day2020s_02a.txt -obscode ObsCodes.txt -colformat colformat_LSST_01.txt`

If this run succeeds, you can compare the output files to the 'check' files you downloaded along with the source code:

`diff LSST_pairdets_01.csv LSST_pairdets_01_check.csv` <br>
`diff LSST_pairs_01 LSST_pairs_01_check`

Both diffs should be clean if everything has gone well. If they are not clean, you can check the word count and sizes of the files to see if the differences are innocuous (e.g., roundoff error in the last digit) or significant.

### Testing heliolinc ###

As a first test of heliolinc, try this invocation:

`heliolinc -dets LSST_pairdets_01.csv`

This is well short of the minimal executable invocation, but it serves a useful purpose: since you have supplied a detection file (LSST_pairdets_01.csv) and no reference time, heliolinc automatically calculates the time at the center of the spread of time values given in the detection file (that is, the average of the earliest time and the latest time), and suggests this as a reference time in its output error message. The time is required to be in Modified Julian Days (MJD), hence the relevant part of the output error message should be:

> ERROR: input positive-valued reference MJD is required <br>
> Suggested value is 60608.63 <br>
> based on your input detection catalog LSST_pairdets_01.csv <br>

If you don't see this output (as part of a longer error message including suggested usage) check to see if the input file `LSST_pairdets_01.csv`, which is supposed to have been produced when you tested `make_tracklets`, actually exists and matches the comparison file `LSST_pairdets_01_check.csv`.

If the error message came out as expected, try the actual minimum invocation of heliolinc using the suggested reference time (i.e. referece MJD):

`heliolinc -dets LSST_pairdets_01.csv -pairs LSST_pairs_01 -mjd 60608.63 -obspos Earth1day2020s_02a.txt -heliodist accelmat_mb08a_sp04.txt -out LSST_hl_all01.csv -outsum LSST_hl_summary01.csv`

After printing some descriptive lines about its input and configuration (designed to help users catch errors with the input files or other arguments), heliolinc should start reporting its progress like this:

> 122776 detection records read from LSST_pairdets_01.csv. <br>
> Writing to LSST_hl_all01.csv <br>
> Read 122776 detections and 5948 pairs. <br>
> Working on grid point 0: r = 2.8400 AU, v = 13.8517 km/sec, dv/dt = -0.4000 GMsun/r^2 <br>
> 5948 input pairs/tracklets led to 5832 physically reasonable state vectors <br>
> Identified 1960 candidate linkages <br>
> Working on grid point 1: r = 2.8400 AU, v = 13.8517 km/sec, dv/dt = -0.8000 GMsun/r^2 <br>
> 5948 input pairs/tracklets led to 5831 physically reasonable state vectors <br>
> Identified 1956 candidate linkages <br>
> Working on grid point 2: r = 2.8800 AU, v = -13.8517 km/sec, dv/dt = -0.4000 GMsun/r^2 <br>
> 5948 input pairs/tracklets led to 5426 physically reasonable state vectors <br>
> Identified 1815 candidate linkages <br>

It will go through a total of 200 grid points, which are specified by the input file `accelmat_mb08a_sp04.txt`. Each grid point constitutes a hypothesis about the target asteroids' heliocentric distance as a function of time, and the result is supposed to be the correct linkage of any real asteroids whose motion matches the hypothesis, provided the input data contains a sufficient number of detections of each such asteroid. To turn off heliolinc's progress reporting, add the flag `-verbose -1` to the command-line arguments.

This should finish in about five minutes.

We have not provided comparison files for checking the heliolinc output, because the output files are a bit large to post on github. However, you can check the validity of the input by seeing if the wordcounts match:

`wc  LSST_hl_all01.csv`
> 48482740   48482740 3875752557 LSST_hl_all01.csv


`wc LSST_hl_summary01.csv`
> 388242   388242 62674788 LSST_hl_summary01.csv

### Testing link_refine ###

The file `LSST_hl_summary01.csv` output by heliolinc in the the previous test summarizes each of the candidate linkages with a single line. The line count given above therefore indicates that heliolinc identified 388241 distinct candidate linkages (subtracting 1 for the header line). However, `LSST_raw_input_data01a.csv`, the test data file we provided as input for `make_tracklets`, contains detections only for 1000 simulated asteroids. Hence, there are almost 400 times as many candidate linkages as there are actual distince asteroids in this example.

There are several reasons why there are more candidates than real objects:

* All detections of the same asteroid might be correctly linked under more than one heliocentric hypothesis.
* Multiple non-identical subsets of the detections for a given asteroid might lead to multiple distinct linkages.
* Spurious linkages consisting of detections from more than one asteroid can exist.

In the heliolinc suite, the program `link_refine` exists to deal with these cases and cull the vast number of linkages produced by heliolinc down to an optimized, non-redundant list. It does this by identifying sets of mutually overlapping candidate linkages, choosing the best one, and rejecting all the inferior linkages that overlap (share detections) with it. The criteria used to define the *best* linkages will be descried in more detail below.

For the present, in order to test your implementaion of `link_refine`, it is necessary first to create a text file listing the output from `heliolinc` that will be analyzed. The reason `link_refine` requires a list of files is that it can be used to analyze output from multiple executions of `heliolinc` at one time. Each line of this text file gives the names of both the output files from a single execution of heliolinc: first the comprehensive output file (e.g., `LSST_hl_all01.csv`) and then the summary output file (e.g., `LSST_hl_summary01.csv`). For this test, your link file list will have only one line since you are analyzing results from just one execution of `heliolinc`. You can construct the list file by hand in a text editor such as emacs, or automatically with a command such as:

`printf "LSST_hl_all01.csv LSST_hl_summary01.csv\n" > LSST_lflist01`

Its contents should be simply the single line:

> `LSST_hl_all01.csv LSST_hl_summary01.csv`

You can then test `link_refine` with the invocation:

`link_refine -pairdet LSST_pairdets_01.csv -lflist LSST_lflist01 -outfile LSST_linkref_all01.csv -outsum LSST_linkref_summary01.csv`

This should finish in about two minutes.

The output files `LSST_linkref_all01.csv` and `LSST_linkref_summary01.csv` are formatted identically to the comprehensive and summary files, respectively, that were originally outputted by `heliolinc`. This allows `link_refine` to be run recursively -- that is, its output files can be cycled back as inputs to a new execution. Reasons why you might want to do this will be discussed below.

For now, you can test the successful execution using the comparison files `LSST_linkref_all01_check.csv` and `LSST_linkref_summary01_check.csv`

`diff LSST_linkref_all01.csv LSST_linkref_all01_check.csv` <br>
`diff LSST_linkref_summary01.csv LSST_linkref_summary01_check.csv`

Both diffs should be clean. If they are not, you can compare wordcounts and mumerical values to see if the differences are benign -- e.g., roundoff errors affecting the least significant digits -- or substantive.

Assuming the ouputs match the supplied check files, your wordcounts will indicate that `link_refine` has culled down the 388241 candidate linkages that were input into a pure and complete list of exactly 1000 simulated real asteroids:

`wc LSST_linkref_summary01.csv` <br>
> 1001   1001 152585 LSST_linkref_summary01.csv

Where the one extra line is the header.

If you had input data from multiple observing sites, instead of just a single site as in this example, you would want to use the program `link_refine_multisite` instead of `link_refine`. This not because `link_refine_multisite` is more sophisticated or powerful -- actually the opposite is true -- but `link_refine` has an powerful feature that significantly improves the results on single-site data but that (for well-understood reasons) does not correctly analyze linkages with detections from multiple observing sites. Hence, the less sophisticated program `link_refine_multisite` has to be used if you have observations from multiple observing sites.

This completes the test run of the regular heliolinc suite of programs. The format of the output files is csv (comma separated values). Most of the files have single-line headers that should give you some idea of what the columns mean. If you want to manipulate these files using `awk`, you just inform `awk` that the values are separated by commas rather than whitespace by writing `awk -v FS=","` followed by the sequence of instructions you'd use if the values were separated by whitespace.

For example, if you wanted to find which cluster (that is, distinct linkage) spans the largest amount of time, you could look at the header of the summary file `LSST_linkref_summary01.csv` and see that 'clusternum' is the first column and 'timespan' is the sixth. Then you could write:

`awk -v FS="," '{printf "%s %s\n", $1, $6}' LSST_linkref_summary01.csv | sort -k2 -g`

This generates a list of just the cluster (linkage) number and the time-span in days, sorted by time-span so the one with the longest time span will be at the bottom. You will see that this is **cluster number 611**, with observations spanning 14.944 days. To pull out all the observations corresponding to this particular cluster, you could use `awk` on the comprehensive output file `LSST_linkref_all01.csv` first looking at the header line to see that the cluster number is the 11th and last column. Then, you could write:

`awk -v FS="," 'NR==1 || $11==611 {printf "%s\n", $0}' LSST_linkref_all01.csv`

You get the header line for context (that's what the `NR==1` part does), plus what turns out to be 22 detections heliolinc has decided belong to this particular asteroid, with their observation times (MJD), celestial coordinates (RA and Dec), magnitude, photometric band, and a few other thing that will be explained below.

For now, if all the tests have been successful, you probably want to run the heliolinc suite on your own data. In order to accomplish that, all you need to do is understand the input files and translate your data into the required format. 

## Understanding the File Formats ##

### Original Input Data ###

To input new astronomical data into the heliolinc suite, you just need to creat a single file: the input data file for `make_tracklets`. The format required is a csv file with a one-line header followed by one line for each detection, with columns containing at least the following quantities:

* A string identifier for each detection.
* The time the detection was made (e.g. mid-exposure) in Modified Julian Days (MJD)
* The right ascension (RA) of the detection, in decimal degrees.
* The declination (Dec) of the detection, in decimal degrees.
* The magnitude of the detection.
* The photometric band in which the observation was made (a single character such as u, g, r, i, etc).
* The 3-character observatory code of the observing site (e.g. 695 for Kitt Peak or I11 for Gemini South)

Any number of additional columns is allowed, and the required columns can be in any order. How does `make_tracklets` find the required data when the file-format is so flexible? It requires a column format file, introduce on the command line by the `-colformat` keyword, which simply tells it which column has each of the seven required quantities. For example, here is the full contents of the column format file `colformat_LSST_01.txt` used in our example invocation above:

`IDCOL 1` <br>
`MJDCOL 3` <br>
`RACOL 6` <br>
`DECCOL 8` <br>
`MAGCOL 32` <br>
`BANDCOL 26` <br>
`OBSCODECOL 38` <br>

That's all there is to it. To run `make_tracklets` on your own data, create a csv file with any columns you want (as long as it includes the required seven) and then create a column format file that tells `make_tracklets` where to find the seven things it needs.

### Input Image File ###

If you have easy access to information about the time (MJD) of mid-exposure and the boresight pointing for each image in your input data, you can (and should) supply this information to `make_tracklets` in the form of a 4-column file that has **no header** and is **not a csv** but has columns separated by spaces. The name of this file is supplied to `make_tracklets` using the command line keyword `-imgs`. The four columns in the file are MJD, RA, Dec, and observatory code. The RA and Dec must be in decimal degrees. The observatory code must be a three-character string, and it and the MJD **must match** the corresponding entries in the input detection data file. Here is an example of the first few lines in an acceptable input image file. 

`59790.705516 146.323066 -34.685595 M22` <br>
`59790.705967 150.855833 -29.401107 M22` <br>
`59790.706458 155.796972 -18.751285 M22` <br>
`59790.706908 154.952426 -24.263799 M22` <br>

If you do not have easy access to this information, don't worry: `make_tracklets` reconstructs it internally from the input detection data. If you want `make_tracklets` to write out its internally-reconstructed image log in the required format, supply the name of an output image file to be written using the command line keyword `-outimgs`. The output file has two additional columns recording the first and last entries (actually, one after the last entry) in the output paired detection catalog (see below) that come from that image. For example:

`59790.705516 146.323066 -34.685595 M22 0 529` <br>
`59790.705967 150.855833 -29.401107 M22 529 955` <br>
`59790.706458 155.796972 -18.751285 M22 955 1434` <br>
`59790.706908 154.952426 -24.263799 M22 1434 1768` <br>

Because `make_tracklets` simply ignores input image columns beyond the first four, the output image file written by one invocation of `make_tracklets` can be used as an input image file for future invocations on detections from the same images. A case where you might want to do this is if the output image file was produced by `make_tracklets` running on all available data, and you were re-running `make_tracklets` after aggressively culling the input data down to a much smaller number of detections. Because `make_tracklets` estimates the boresight RA and Dec of a given image by simply averaging over all the detections taken at the same time and observatory code, the old values based on a larger number of detections will be more accurate that what can be reconstructed from the smaller, culled set of detections.

### Observer Location Files ###

##### Executive summary: just use the ones we gave you. #####

The remaining files required by `make_tracklets` are the ones used to determine the observer's position in the solar system at the time of each observation. In our example invocation, these are `Earth1day2020s_02a.txt` and `ObsCodes.txt`, which are introduced with the command line keywords `-earth` and `-obscode`, respectively. The files used in the example are broadly applicable, and it's unlikely you'll need different files. The `-earth` file specifies the Earth's position relative to the sun as a function of time, while the `-obscode` (observatory code) file allows `make_tracklets` to look up the location of the observatory on the Earth's surface based on the three-character observatory codes specified in the input detection file. 

Cases where you might need different files include observations before December 2019 or after December 2030; or made from a new observatory that doesn't have an existing observing code. You can get a new `-earth` file using the JPL Horizons ephemeris generator (`https://ssd.jpl.nasa.gov/horizons/app.html`). Set Ephmeris Type to `Vector Table`; Target Body to `Earth`; and Coordinate Center to `Sun (body center)` (try entering 'sol' or '@10' in the search box). Set Time Specification to whatever range of dates you want, at 1-day sampling, and leave Table Settings at *defaults*. Click `Generate Ephemeris` and then `Download Results`; save the resulting file; and feed it directly into `make_tracklets` using the command line keyword `-earth`. For a new observatory code, your best option is probably to go to `https://minorplanetcenter.net/iau/lists/ObsCodesF.html`, find the lines for the new observatories you need, and just add them to the existing file `ObsCodes.txt`. Alternatively, you can copy the whole file, but be careful to do the following before you feed it to `make_tracklets`:

* Leave a one-line header
* Remove all the lines that have columns 2-4 empty (these are space-based observatories, which `make_tracklets` cannot handle)

### Output files from make_tracklets ###

The optional output image file introduced with the `-outimgs` command line keyword has already been described above. The primary output files produced by `make_tracklets` are the catalog of paired detections ('paired detection file') and the file specifying the pairs and tracklets ('pair file').

The default name of the output **paired detection file** is `pairdetfile01.csv`. A different name can be specified using the command line keyword `-pairdets`. This file echoes the data in the input detection catalog, but in a standardized form expected by `heliolinc`, and with two important changes. First, the observer's position in Cartesian ecliptic coordinates relative to the center of the sun at the time of each detection is calculated and recorded. Second, only detections that were included in some viable pair or tracklet are included -- hence, the paired detection file does not necessarily have as many entries as the input detection catalog. Here is a description of the columns:

| Column | Name        | Description                                                                |
| -------| ----------- | -------------------------------------------------------------------------- |
1        | `MJD`       |  Modified Julian Day			                                    |
2 	 | `RA`        |  right ascension in decimal degrees					    |
3 	 | `Dec`       |  declination in decimal degrees					    |
4 	 | `observerX` |  observer's Cartesian x coordinate relative to the heliocenter, in km	    |
5 	 | `observerY` |  observer's Cartesian y coordinate relative to the heliocenter, in km	    |
6 	 | `observerZ` | observer's Cartesian z coordinate relative to the heliocenter, in km	    |
7 	 | `stringID`  |  string identifier copied from input detection catalog		   	    |
8 	 | `mag`       |  magnitude								    |
9 	 | `band`      | photometric band, single character: u, g, r, i, etc			    |
10	 | `obscode`   |  Observatory code, three-character string    	 			    |
11 	 | `origindex` | Line number of corresponding entry in the original input detection catalog |

The `origindex` column is not explicitly used by `heliolinc` or `link_refine`, but it is retained into the final outputs of both programs, enabling easy mapping of the final linkages back to lines in the original detection catalog. This is likely to be useful in the (quite probable) case that you had (and wrote into your input detection catalog) a bunch of interesting data besides the seven columns required by `heliolinc`. You can then use the preserved line numbers to recover this original, detailed information for every linked detection. Note that `origindex` literally records the line numbers of your original input data file, starting with the header as line number 1. Hence, the first data line is `origindex=2`.

The default name of the output **pair file** is `outpairfile01`. A different name can be specified using the command line keyword `-pairs`. The format of the pair file is designed to store potentially many millions of pairs and tracklets in a text file with relatively small on-disk size. Hence, the precise floating-point values of the MJD, RA, and Dec are not provided whenever it can be avoided -- instead, integer indices to the paired detection file are used.

Making things more complicated is the need to represent both pairs (always consisting of just two detections) and tracklets (consisting of more than two detections) in the same file. A pair is represented simply by a 'P' followed by the indices of the two detections in the paired detection file. Here are some example lines indicating simple pairs:

`P 119468 119469` <br>
`P 119466 119467` <br>
`P 119464 119465` <br>
`P 119462 119463` <br>

Here, for example, the first line says that detections number 119468 and 119469 in the paired detection file could plausibly be two detections of the same object, and hence they form a pair.

A tracklet with more than two points could be represented as a pair with just its first and last detection, but this would throw away the opportunity for greater astrometric precision from the larger number of detections. Hence, when `make_tracklets` identifies a viable tracklet, it performs a least-squares fit modeling the object's motion as a Great Circle at constant angular velocity (the tracklet is discarded if the residuals from this fit are too large). Otherwise, a representative point is chosen near the beginning of the tracklet and another near the end, and the Great Circle fit is evaluated at the times corresponding to these representative points. This produces new RA, Dec positions that represent the on-sky motion of the tracklet as a whole more accurately than the original RA, Dec at the representative points alone. Since these new RA, Dec positions exist nowhere in the paired detection file, it is necessary to write them to the pair file -- although the MJDs can still be indicated by indexing the paired detection file.

Here is an example of a set of lines from the pair file describing a tracklet with 6 points:

`T 1969 2161 358.432630 -7.439334 358.424803 -7.439827 6` <br>
`1968` <br>
`1969` <br>
`2150` <br>
`2156` <br>
`2161` <br>
`2168` <br>

The initial 'T' informs `heliolinc` that it is reading a tracklet with more than two points. The next two integer values are the representative points: detections number 1968 and 2161 in the paired detection file. There follow the RA and Dec from the Great Circle fit for the first representative point, then the RA and Dec from the fit to the other representative point; and finally the number of points in the tracklet: in this case, six, all of which are listed on the six lines that follow. For its internal analysis, `heliolinc` uses only the Great Circle fit RA and Dec for the representative points, and the corresponding MJDs looked up from the paired detection file. However, it reads the indices of the other detections in the tracklet in order to include this information in its output files.

### Input Files for heliolinc ###

The most important `heliolinc` input files are the ones generated by `make_tracklets` and described above. It also requires an Earth position file, but we have already discussed this since it is also a required input for `make_tracklets`. The only `heliolinc` input file that has yet to be described is the one containing the hypotheses about heliocentric radial motion, which is introduced by the command line keyword `-heliodist`, and was called `accelmat_mb08a_sp04.txt` in our example invocation. This file has four columns (or more; additional columns are ignored). It has a one-line header and is **not a csv** but rather a space-delimited file like the optional image files that can be read and written by `make_tracklets`. The four required columns are as follows:

| Column | Name        | Description                                                                |
| -------| ----------- | -------------------------------------------------------------------------- |
1        | `r(AU)`       |  Distance from the sun at the reference time		                                    |
2 	 | `rdot(AU/day)`|  Radial velocity at the reference time			    |
3 	 | `norm`        |  Positive integer = use this point; 0 = skip this point		    |
4 	 | `mean_accel`  |  Time-derivative of radial velocity, divided by (-GMsun/r^2)	    |

The original meaning of the `norm` column comes from when we generated files of this type by integrating actual asteroid orbits: then `norm` meant the number of times (in a very long integration) that any real asteroid was instantaneously found with the corresponding heliocentric distance and radial velocity (within a tolerance defined by the grid spacing).

The solar gravity always imposes, on any object in the solar system, a vector acceleration of GMsun/r^2 toward the center of the sun. The `mean_accel` column in the hypothesis file is not the vector acceleration but rather the time-derivative of the heliocentric radial velocity (equivalently, the second time-derivative of the heliocentric distance). We are inclined to call this the **radial acceleration**, and it is generally not equal to the vector acceleration. For example, in a circular orbit the radial velocity and radial acceleration are both always zero, even though the vector acceleration is still GMsun/r^2 toward the heliocenter. This is why we use units of -GMsun/r^2 for the radial acceleration term `mean_accel`: the value is 0 for a circular orbit and 1.0 for an object with zero angular momentum (Keplerian eccentricity exactly 1.0) that is falling directly toward the center of the sun.

Realistic values for actual bound orbits in the solar system are typically between -1.0 and +1.0, with a distribution that is highly dependent on the heliocentric distance and radial velocity. This is why we used integration of actual asteroid orbits to generate heliocentric hypothesis files.




REQUIRED:
indetfile=argv[++i];
	earthfile=argv[++i];
	obscodefile=argv[++i];

NOT REQUIRED, DEFAULT NONE:
inimfile=argv[++i];
outimfile=argv[++i];

NOT REQUIRED, DEFAULT EXISTS:
outpairfile=argv[++i];
pairdetfile=argv[++i];
        imrad=stod(argv[++i]);
        maxtime=stod(argv[++i]);
        mintime=stod(argv[++i]);
        minvel=stod(argv[++i]);
        maxvel=stod(argv[++i]);
        maxgcr=stod(argv[++i]);
        minarc=stod(argv[++i]);
	mintrkpts=stoi(argv[++i]);
	colformatfile=argv[++i];
