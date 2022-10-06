# Heliolinc3D C++ Implementation #
## Why download heliolinc? ##

The motivation for this implementation of the heliolinc algorithm is to enable asteroid discovery using the LSST survey strategy of taking just two images of each field per night, rather than the usual practice of taking four images per field per night. This simple change in survey strategy requires a paradigm shift in asteroid detection. Surveys taking the usual four images per night can identify candidate asteroid discoveries based on just a single night's data. With only two images per night, the attempt to do this would result in overwhelming numbers of false positives. Identifying a reasonable discovery candidate requires **linking** multiple detection pairs across multiple nights of data. The new C++ implementation of heliolinc has proven capability to link simulated **and real** asteroid detections across multiple nights, fast enough that ingesting detection catalogs spanning two weeks of output for major surveys -- including detections from more than one observing site -- is computationally tractable. 

This does not mean that heliolinc has rendered obsolete the typical survey strategy of taking four images per field per night. This strategy continues to have significant advantages over a two image per night strategy, even with heliolinc. Relative to LSST, the currently-operating four-image surveys have much better sensitivity to small asteroids that have brief, close encounters with Earth. These include small asteroids passing within the 0.01 AU of Earth, and very small asteroids on their 'final plunge' toward impact with the Earth -- several or which have been discovered by currently operating surveys. The high-performance implementation of the heliolinc algorithm that we present here should **not** be interpreted as a reason for ongoing asteroids surveys to abandon their highly successful four-image strategy.

Instead, heliolinc has something to offer to surveys using the conventional four-image strategy -- even though it was invented to enable a completely different strategy. What heliolinc can offer to a four-image survey is the likelihood of successful linking and discovery of objects detected fewer than four times on a give night. This can happen for many reasons: the survey might not acquire four images of a particular field because of weather or other contingencies, or the asteroid might not be successfully detected on all four images because of varying image quality, superposition on a star, rotational brightness variations, or a host of other reasons. Such 'missing' detections are most likely for the faintest objects, hence heliolinc has the potential to extend a survey's sensitivity toward fainter magnitudes.

Additionally, heliolinc has the potential to enable remarkable asteroid discoveries from data sets not acquired with the objective of finding asteroids at all.


The `heliolinc` suite described herein has proven ability to ingest up to a month of detection catalogs from major surveys. 


## Obtain the following files from GitHub: ##

#### Source code: ####

```
solarsyst_dyn_geo01.cpp
solarsyst_dyn_geo01.h
make_tracklets.cpp
heliolinc.cpp
link_refine.cpp
link_refine_multisite.cpp
```

#### General-use input files: ####

```
Earth1day2020s_02a.txt
ObsCodes.txt
colformat_LSST_01.txt
accelmat_mb08a_sp04.txt
```

#### Test data files: ####

```
LSST_raw_input_data01a.csv
LSST_pairdets_01_check.csv
LSST_pairs_01_check
LSST_linkref_all01_check.csv
LSST_linkref_summary01_check.csv
```

## Compile the heliolinc suite ##

#### Use these suggested compile commands, or their equivalents for your system: ####

```
c++ -O3 -o heliolinc heliolinc.cpp solarsyst_dyn_geo01.cpp -std=c++11
c++ -O3 -o make_tracklets make_tracklets.cpp solarsyst_dyn_geo01.cpp -std=c++11
c++ -O3 -o link_refine_multisite link_refine_multisite.cpp solarsyst_dyn_geo01.cpp -std=c++11
c++ -O3 -o link_refine link_refine.cpp solarsyst_dyn_geo01.cpp -std=c++11
```

## What the programs do: ##

**make_tracklets:** Perform image-based pairing to create input pair/tracklet files for heliolinc

**heliolinc:** link together pairs/tracklets produced by make_tracklets into candidate asteroid discoveries

**link_refine and link_refine_multisite:** post-process linkages produced by heliolinc to produce a de-duplicated, non-overlapping set of candidate discoveries with the highest likelihood of being real.

## Testing your installation: ##

### Testing make_tracklets ###

If you've downloaded the recommended test data files, you can test your installation immediately. Here is the minimalist invocation for `make_tracklets`, the first program in the suite:

`./make_tracklets -dets LSST_raw_input_data01a.csv -earth Earth1day2020s_02a.txt -obscode ObsCodes.txt -colformat colformat_LSST_01.txt`

**This run should take about 10 seconds.**

If you got a final message stating, "Constructing tracklets, and writing pairs to output file," the code probably executed correctly. It will have produced output files called `pairdetfile01.csv` and `outpairfile01`. The names of these output files can be set with command-line options, but in the minimalist invocation they default to the names given here.

`pairdetfile01.csv` is the **paired detection file**, which is just a reformatted version of the input detection catalog limited only to detections that formed pairs or longer tracklets.

`outpairfile01` is the **pair file**, which records the pairs and longer tracklets that were found using integer indices that specify their position in the paired detection file.

If the minimalist invocation failed, check that you have the input files used in the example:

```
LSST_raw_input_data01a.csv
Earth1day2020s_02a.txt
ObsCodes.txt
colformat_LSST_01.txt
```
If the minimalist invocation succeeded, you can try a new run in which you specify the output names:

`./make_tracklets -dets LSST_raw_input_data01a.csv -pairdets LSST_pairdets_01.csv -pairs LSST_pairs_01 -earth Earth1day2020s_02a.txt -obscode ObsCodes.txt -colformat colformat_LSST_01.txt`

If this run succeeds, you can compare the output files to the 'check' files you downloaded along with the source code:

```
diff LSST_pairdets_01.csv LSST_pairdets_01_check.csv
diff LSST_pairs_01 LSST_pairs_01_check
```

Both diffs should be clean if everything has gone well. If they are not clean, you can check the word count and sizes of the files to see if the differences are innocuous (e.g., roundoff error in the last digit) or significant.

### Testing heliolinc ###

As a first test of heliolinc, try this invocation:

`./heliolinc -dets LSST_pairdets_01.csv`

This is well short of the minimal executable invocation, but it serves a useful purpose: since you have supplied a detection file (LSST_pairdets_01.csv) and no reference time, heliolinc automatically calculates the time at the center of the spread of time values given in the detection file (that is, the average of the earliest time and the latest time), and suggests this as a reference time in its output error message. The time is required to be in Modified Julian Days (MJD), hence the relevant part of the output error message should be:

```
ERROR: input positive-valued reference MJD is required
Suggested value is 60608.63
based on your input detection catalog LSST_pairdets_01.csv
```
If you don't see this output (as part of a longer error message including suggested usage) check to see if the input file `LSST_pairdets_01.csv`, which is supposed to have been produced when you tested `make_tracklets`, actually exists and matches the comparison file `LSST_pairdets_01_check.csv`.

If the error message came out as expected, try the actual minimum invocation of heliolinc using the suggested reference time (i.e. reference MJD):

`./heliolinc -dets LSST_pairdets_01.csv -pairs LSST_pairs_01 -mjd 60608.63 -obspos Earth1day2020s_02a.txt -heliodist accelmat_mb08a_sp04.txt -out LSST_hl_all01.csv -outsum LSST_hl_summary01.csv`

After printing some descriptive lines about its input and configuration (designed to help users catch errors with the input files or other arguments), heliolinc should start reporting its progress like this:
```
122776 detection records read from LSST_pairdets_01.csv.
Writing to LSST_hl_all01.csv
Read 122776 detections and 5948 pairs.
Working on grid point 0: r = 2.8400 AU, v = 13.8517 km/sec, dv/dt = -0.4000 GMsun/r^2
5948 input pairs/tracklets led to 5832 physically reasonable state vectors
Identified 1960 candidate linkages
Working on grid point 1: r = 2.8400 AU, v = 13.8517 km/sec, dv/dt = -0.8000 GMsun/r^2
5948 input pairs/tracklets led to 5831 physically reasonable state vectors
Identified 1956 candidate linkages
Working on grid point 2: r = 2.8800 AU, v = -13.8517 km/sec, dv/dt = -0.4000 GMsun/r^2
5948 input pairs/tracklets led to 5426 physically reasonable state vectors
Identified 1815 candidate linkages
```
It will go through a total of 200 grid points, which are specified by the input file `accelmat_mb08a_sp04.txt`. Each grid point constitutes a hypothesis about the target asteroids' heliocentric distance as a function of time, and the result is supposed to be the correct linkage of any real asteroids whose motion matches the hypothesis, provided the input data contains a sufficient number of detections of each such asteroid. To turn off heliolinc's progress reporting, add the flag `-verbose -1` to the command-line arguments.

This should finish in about five minutes.

We have not provided comparison files for checking the heliolinc output, because the output files are a bit large to post on GitHub. However, you can check the validity of the output by seeing if the word-counts match:

`wc  LSST_hl_all01.csv`
> 48482740   48482740 3875752557 LSST_hl_all01.csv


`wc LSST_hl_summary01.csv`
> 388242   388242 62674788 LSST_hl_summary01.csv

### Testing link_refine ###

The file `LSST_hl_summary01.csv` output by heliolinc in the the previous test summarizes each of the candidate linkages with a single line. The line count given above therefore indicates that heliolinc identified 388241 distinct candidate linkages (subtracting 1 for the header line). However, `LSST_raw_input_data01a.csv`, the test data file we provided as input for `make_tracklets`, contains detections only for 1000 simulated asteroids. Hence, there are almost 400 times as many candidate linkages as there are actual distinct asteroids in this example.

There are several reasons why there are more candidates than real objects:

* All detections of the same asteroid might be correctly linked under more than one heliocentric hypothesis.
* Multiple non-identical subsets of the detections for a given asteroid might lead to multiple distinct linkages.
* Spurious linkages consisting of detections from more than one asteroid can exist.

In the heliolinc suite, the program `link_refine` exists to deal with these cases and cull the vast number of linkages produced by heliolinc down to an optimized, non-redundant list. It does this by identifying sets of mutually overlapping candidate linkages, choosing the best one, and rejecting all the inferior linkages that overlap (share detections) with it. The criteria used to define the *best* linkages will be descried in more detail below.

For the present, in order to test your implementation of `link_refine`, it is necessary first to create a text file listing the output from `heliolinc` that will be analyzed. The reason `link_refine` requires a list of files is that it can be used to analyze output from multiple executions of `heliolinc` at one time. Each line of this text file gives the names of both the output files from a single execution of heliolinc: first the comprehensive output file (e.g., `LSST_hl_all01.csv`) and then the summary output file (e.g., `LSST_hl_summary01.csv`). For this test, your link file list will have only one line since you are analyzing results from just one execution of `heliolinc`. You can construct the list file by hand in a text editor such as emacs, or automatically with a command such as:

`printf "LSST_hl_all01.csv LSST_hl_summary01.csv\n" > LSST_lflist01`

Its contents should be simply the single line:

> `LSST_hl_all01.csv LSST_hl_summary01.csv`

You can then test `link_refine` with the invocation:

`link_refine -pairdet LSST_pairdets_01.csv -lflist LSST_lflist01 -outfile LSST_linkref_all01.csv -outsum LSST_linkref_summary01.csv`

This should finish in about two minutes.

The output files `LSST_linkref_all01.csv` and `LSST_linkref_summary01.csv` are formatted identically to the comprehensive and summary files, respectively, that were originally outputted by `heliolinc`. This allows `link_refine` to be run recursively -- that is, its output files can be cycled back as inputs to a new execution. Reasons why you might want to do this will be discussed below.

For now, you can test the successful execution using the comparison files `LSST_linkref_all01_check.csv` and `LSST_linkref_summary01_check.csv`
```
diff LSST_linkref_all01.csv LSST_linkref_all01_check.csv
diff LSST_linkref_summary01.csv LSST_linkref_summary01_check.csv
```
Both diffs should be clean. If they are not, you can compare word-counts and numerical values to see if the differences are benign -- e.g., roundoff errors affecting the least significant digits -- or substantive.

Assuming the outputs match the supplied check files, your word-counts will indicate that `link_refine` has culled down the 388241 candidate linkages that were input into a pure and complete list of exactly 1000 simulated real asteroids:

`wc LSST_linkref_summary01.csv`
> 1001   1001 152585 LSST_linkref_summary01.csv

Where the one extra line is the header.

If you had input data from multiple observing sites, instead of just a single site as in this example, you would want to use the program `link_refine_multisite` instead of `link_refine`. This not because `link_refine_multisite` is more sophisticated or powerful -- actually the opposite is true -- but `link_refine` has an powerful feature that significantly improves the results on single-site data but that (for well-understood reasons) does not correctly analyze linkages with detections from multiple observing sites. Hence, the less sophisticated program `link_refine_multisite` has to be used if you have observations from multiple observing sites.

This completes the test run of the regular heliolinc suite of programs. The format of the output files is csv (comma separated values). Most of the files have single-line headers that should give you some idea of what the columns mean. If you want to manipulate these files using `awk`, you just inform `awk` that the values are separated by commas rather than white space by writing `awk -v FS=","` followed by the sequence of instructions you'd use if the values were separated by white space.

For example, if you wanted to find which cluster (that is, distinct linkage) spans the largest amount of time, you could look at the header of the summary file `LSST_linkref_summary01.csv` and see that 'clusternum' is the first column and 'timespan' is the sixth. Then you could write:

`awk -v FS="," '{printf "%s %s\n", $1, $6}' LSST_linkref_summary01.csv | sort -k2 -g`

This generates a list of just the cluster (linkage) number and the time-span in days, sorted by time-span so the one with the longest time span will be at the bottom. You will see that this is **cluster number 611**, with observations spanning 14.944 days. To pull out all the observations corresponding to this particular cluster, you could use `awk` on the comprehensive output file `LSST_linkref_all01.csv` first looking at the header line to see that the cluster number is the 11th and last column. Then, you could write:

`awk -v FS="," 'NR==1 || $11==611 {printf "%s\n", $0}' LSST_linkref_all01.csv`

You get the header line for context (that's what the `NR==1` part does), plus what turns out to be 22 detections heliolinc has decided belong to this particular asteroid, with their observation times (MJD), celestial coordinates (RA and Dec), magnitude, photometric band, and a few other thing that will be explained below.

For now, if all the tests have been successful, you probably want to run the heliolinc suite on your own data. In order to accomplish that, all you need to do is understand the input files and translate your data into the required format. 

## Understanding the File Formats ##

### Original Input Data ###

To input new astronomical data into the heliolinc suite, you just need to create a single file: the input data file for `make_tracklets`. The format required is a csv file with a one-line header followed by one line for each detection, with columns containing at least the following quantities:

* A string identifier for each detection.
* The time the detection was made (e.g. mid-exposure) in Modified Julian Days (MJD)
* The right ascension (RA) of the detection, in decimal degrees.
* The declination (Dec) of the detection, in decimal degrees.
* The magnitude of the detection.
* The photometric band in which the observation was made (a single character such as u, g, r, i, etc).
* The 3-character observatory code of the observing site (e.g. 695 for Kitt Peak or I11 for Gemini South)

Any number of additional columns is allowed, and the required columns can be in any order. The columns that hold the required data are communicated to `make_tracklets` through a column format file, introduced on the command line by the `-colformat` keyword, which simply tells it which column has each of the seven required quantities. For example, the full contents of the column format file `colformat_LSST_01.txt` used in our example invocation above are:
```
IDCOL 1
MJDCOL 3
RACOL 6
DECCOL 8
MAGCOL 32
BANDCOL 26
OBSCODECOL 38
```
That's all there is to it. To run `make_tracklets` on your own data, create a csv file with any columns you want (as long as it includes the required seven) and then create a column format file that tells `make_tracklets` where to find the seven things it needs. Note that it starts counting columns from 1, not from 0.

### Input Image File ###

If you have easy access to information about the time (MJD) of mid-exposure and the bore-sight pointing for each image in your input data, you can (and should) supply this information to `make_tracklets` in the form of a 4-column file that has **no header** and is **not a csv** but has columns separated by spaces. The name of this file is supplied to `make_tracklets` using the command line keyword `-imgs`. The four columns in the file are MJD, RA, Dec, and observatory code. The RA and Dec must be in decimal degrees. The observatory code must be a three-character string, and it and the MJD **must match** the corresponding entries in the input detection data file. Here is an example of the first few lines in an acceptable input image file. 
```
59790.705516 146.323066 -34.685595 M22
59790.705967 150.855833 -29.401107 M22
59790.706458 155.796972 -18.751285 M22
59790.706908 154.952426 -24.263799 M22
```
If you do not have easy access to this information, don't worry: `make_tracklets` reconstructs it internally from the input detection data. If you want `make_tracklets` to write out its internally-reconstructed image log in the required format, supply the name of an output image file to be written using the command line keyword `-outimgs`. The output file has two additional columns recording the first and last entries (actually, one after the last entry) in the output paired detection catalog (see below) that come from that image. For example:
```
59790.705516 146.323066 -34.685595 M22 0 529
59790.705967 150.855833 -29.401107 M22 529 955
59790.706458 155.796972 -18.751285 M22 955 1434
59790.706908 154.952426 -24.263799 M22 1434 1768
```
Because `make_tracklets` simply ignores input image columns beyond the first four, the output image file written by one invocation of `make_tracklets` can be used as an input image file for future invocations on detections from the same images. A case where you might want to do this is if the output image file was produced by `make_tracklets` running on all available data, and you were re-running `make_tracklets` after aggressively culling the input data down to a much smaller number of detections. Because `make_tracklets` estimates the bore-sight RA and Dec of a given image by simply averaging over all the detections taken at the same time and observatory code, the old values based on a larger number of detections will be more accurate that what can be reconstructed from the smaller, culled set of detections.

### Observer Location Files ###

##### Executive summary: just use the ones we gave you. #####

The remaining files required by `make_tracklets` are the ones used to determine the observer's position in the solar system at the time of each observation. In our example invocation, these are `Earth1day2020s_02a.txt` and `ObsCodes.txt`, which are introduced with the command line keywords `-earth` and `-obscode`, respectively. The files used in the example are broadly applicable, and it's unlikely you'll need different files. The `-earth` file specifies the Earth's position relative to the sun as a function of time, while the `-obscode` (observatory code) file allows `make_tracklets` to look up the location of the observatory on the Earth's surface based on the three-character observatory codes specified in the input detection file. 

Cases where you might need different files include observations before December 2019 or after December 2030; or made from a new observatory that doesn't have an existing observing code. You can get a new `-earth` file using the JPL Horizons ephemeris generator (`https://ssd.jpl.nasa.gov/horizons/app.html`). Set Ephemeris Type to `Vector Table`; Target Body to `Earth`; and Coordinate Center to `Sun (body center)` (try entering 'sol' or '@10' in the search box). Set Time Specification to whatever range of dates you want, at 1-day sampling, and leave Table Settings at *defaults*. Click `Generate Ephemeris` and then `Download Results`; save the resulting file; and feed it directly into `make_tracklets` using the command line keyword `-earth`.

For a new observatory code, your best option is probably to go to `https://minorplanetcenter.net/iau/lists/ObsCodesF.html`, find the lines for the new observatories you need, and just add them to the existing file `ObsCodes.txt`. Alternatively, you can copy the whole file, but be careful to do the following before you feed it to `make_tracklets`:

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
```
P 119468 119469
P 119466 119467
P 119464 119465
P 119462 119463
```
Here, for example, the first line says that detections number 119468 and 119469 in the paired detection file could plausibly be two detections of the same object, and hence they form a pair.

A tracklet with more than two points could be represented as a pair with just its first and last detection, but this would throw away the opportunity for greater astrometric precision from the larger number of detections. Hence, when `make_tracklets` identifies a viable tracklet, it performs a least-squares fit modeling the object's motion as a Great Circle at constant angular velocity (the tracklet is discarded if the residuals from this fit are too large). Otherwise, a representative point is chosen near the beginning of the tracklet and another near the end, and the Great Circle fit is evaluated at the times corresponding to these representative points. This produces new RA, Dec positions that represent the on-sky motion of the tracklet as a whole more accurately than the original RA, Dec at the representative points alone. Since these new RA, Dec positions exist nowhere in the paired detection file, it is necessary to write them to the pair file -- although the MJDs can still be indicated by indexing the paired detection file.

Here is an example of a set of lines from the pair file describing a tracklet with 6 points:
```
T 1969 2161 358.432630 -7.439334 358.424803 -7.439827 6
1968
1969
2150
2156
2161
2168
```
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

Besides the very short example file `accelmat_mb08a_sp04.txt` provided for test runs, we have included the following additional heliocentric hypothesis files:
```
heliohypo_mb01
heliohypo_mb02
heliohypo_mb03
heliohypo_mb04
heliohypo_mb05
heliohypo_neo01
heliohypo_neo02
```

The 'mb' suffix indicates they probe mostly main-belt asteroids, while 'neo' indicates they probe mostly NEOs. They are listed in order of increasing number of hypotheses probed: for example, `heliohypo_mb01` has only 1091 data lines (plus the header) while `heliohypo_mb05` has 53590.

For processing large input data files with many hypothesis, we recommend breaking the files into smaller pieces for embarrassingly parallel runs. For example, the file heliohypo_mb05 could be broken into 27 pieces with 2000 hypotheses each (except the last, which would have only 1590. **Don't forget that each individual hypothesis file needs its own one-line header**. Then heliolinc could be run with a series of commands like this:
```
./heliolinc -dets big_paired_detection_file.csv -pairs big_pair_file -mjd 60608.63 -obspos Earth1day2020s_02a.txt -heliodist heliohypo_mb05_part01 -out big_all_part01.csv -outsum big_summary_part01.csv &
./heliolinc -dets big_paired_detection_file.csv -pairs big_pair_file -mjd 60608.63 -obspos Earth1day2020s_02a.txt -heliodist heliohypo_mb05_part02 -out big_all_part02.csv -outsum big_summary_part02.csv &
./heliolinc -dets big_paired_detection_file.csv -pairs big_pair_file -mjd 60608.63 -obspos Earth1day2020s_02a.txt -heliodist heliohypo_mb05_part03 -out big_all_part03.csv -outsum big_summary_part03.csv &
...
...
./heliolinc -dets big_paired_detection_file.csv -pairs big_pair_file -mjd 60608.63 -obspos Earth1day2020s_02a.txt -heliodist heliohypo_mb05_part27 -out big_all_part27.csv -outsum big_summary_part27.csv &
```

This creates a list of output files from `big_all_part01.csv` to `big_all_part27.csv` and `big_summary_part01.csv` to `big_summary_part27.csv`. All of these files can be fed into the program `link_refine` at one time: `link_refine` was designed to take file lists rather than individual file names for this very reason. If the files are very large, feeding them all into `link_refine` together might exceed the available memory. In this case, `link_refine` can be run individually on each pair of files (e.g., `big_all_part01.csv` and `big_summary_part01.csv`, then `big_all_part02.csv` and `big_summary_part02.csv`, etc.), and then `link_refine` can be run again on the whole list of smaller output files from each of the individual runs. This is made possible by the fact that `link_refine` uses input and output files of identical format: therefore it can re-ingest its own output whenever needed.

You should feel free to make your own heliocentric hypothesis files rather than using the ones we have provided. In some cases very coarse sampling using a very small number of hypotheses can produce good results. For example, you could set all the accelerations (4th column) to zero -- or even do the same to the velocities (2nd column), and probe only fixed heliocentric distances. At the other end of the scale, you might want to probe heliocentric distances of less than 1AU (that is, NEOs interior to Earth's orbit), which are not explored by any of the included files. Good results in this regime require fine sampling: we have probed over one million hypotheses in tests. To save time, make sure you figure out the approximate solar escape velocity for the range of heliocentric distances you are exploring, and avoid probing velocities outside this range. Remember that the velocity units are AU/day.

This discussion of the heliocentric hypothesis files, which are input to heliolinc using the command-line keyword -heliodist, concludes the description of file formats used by programs in the heliolinc C++ suite.

## Taking Control: understanding the optional arguments ##

All the programs in the heliolinc C++ suite allow many optional arguments beyond the minimal invocations used above. In many cases these optional arguments are the key to getting the specific behavior you want from the programs. 

### Arguments for make_tracklets ###

The table below describes all of the arguments 

| command line keyword | type      | Required? | units   | Default | Description                                            |
| -------------------- | --------- | --------- | ------- | ------- | --------------------------------------------- |
 -dets                 | file name | yes       | NA      |  none   | csv-formatted input file containing the detection catalog |
 -earth                | file name | yes       | km, km/sec  |  none    | heliocentric ephemeris for Earth, from JPL Horizons        |
 -obscode              | file name | yes       | NA      |  none   | maps MPC observatory codes to Earth latitude, longitude, and elevation |
 <nobr>-colformat</nobr>            | file name | no        | NA      | see below | tells what columns in the input detection catalog hold each of the required values |
 -imgs                 | file name | no        | NA      | none    | input space delimited file giving MJD, RA, Dec for each image |
 -outimgs              | file name | no	   | NA      | none    | output image catalog file |
 -pairs		       | file name | no        | NA      | outpairfile01 | output file specifying the pairs and tracklets that were found, using indices to the output paired detection catalog |
 -pairdets	       | file name | no        | NA      | pairdetfile01.csv | output paired detection catalog |
 -imrad 	       | float     | no        | degrees | 2.5     | center-to-corner size of images |
 -maxtime              | float     | no        | hours   | 1.5     | maximum time between detections in a pair |
 -mintime              | float     | no        | hours   | 0.000278     | minimum time between detections in a pair |
 -maxGCR	       | float     | no	   | arcsec  | 0.5     | maximum residual from best Great Circle fit for tracklets with more than two points |
 <nobr>-mintrkpts</nobr>           | integer   | no        | NA      | 2       | minimum points for a valid tracklet (using 3 or greater will reject pairs) |
 -minvel               | float     | no        | deg/day | 0.0     | minimum angular velocity for valid pairs or tracklets |
 -maxvel               | float     | no        | deg/day | 1.5     | maximum angular velocity for valid pairs or tracklets |
 -minarc               | float     | no        | arcsec  | 0.0     | minimum angular arc for valid pairs or tracklets      |
 
#### A few more details about previously discussed arguments ####

The column formatting file is technically not required, because `make_tracklets` defaults to the following column specifications:
```
IDCOL 1
MJDCOL 2
RACOL 3
DECCOL 4
MAGCOL 5
BANDCOL 6
OBSCODECOL 7
```

However, in the very likely event that you have constructed your input detection catalog with additional useful data beyond the minimum required columns and/or with a different column order than the above, you will need to create a new column formatting file (using the above as a template) and supply it using the `-colformat` command line keyword.

As previously discussed, the input image file and the output image file are both optional. When no input image file is supplied `make_tracklets` will construct one internally, which will be written out if an output image filename is provided, and discarded otherwise. We recommend always supplying a name for the output image file, since this file is easy to interpret and can be a helpful sanity check. 

The names of the output pair file and the output paired detection file, though technically optional because of the default names describe above, should always be specified.

#### Purpose and use of the other arguments ####

Here, we will attempt to describe how you can use the arguments that haven't yet been discussed to fine-tune the behavior of `make_tracklets` and obtain the desired results for your particular science objective.

**Image radius** `-imrad`: This is nominally the center-to-corner distance, in degrees, for each of your images. For example, with a square field 1 degree on a side, it would be 0.707 degrees. The image radius is used when `make_tracklets` loops over all images, considering each in turn as 'Image A' and constructing a list of possible 'Image B' images whose detections might form pairs or tracklet with detections on Image A. A viable Image B must have a central RA and Dec within twice the image radius (plus an additional margin accounting for the movement of the objects) of the central RA and Dec of image A. Hence, setting the image radius too small can prevent the formation of valid pairs, while setting it too large wastes compute time probing 'candidate' pairs on images that are really too far apart on the sky. Note that if no input image file is provided, `make_tracklets` estimates the central RA and Dec of each image through an average of detections made at the same time from the same observing site. If the average number of detections per image is small, the estimated image center may be quite inaccurate. The possible bad effects of this can be mitigated by using a nominal image radius larger than the true value.

**Maximum time interval** `-maxtime`: This parameter again relates to identifying 'Image B' candidates whose detections might pair with those of a given Image A. If the images are separated by more than the specified maximum time, in hours, no pairs will be sought. This helps avoid excessive numbers of pairs, which can greatly increase the computational load both for `make_tracklets` itself, and for `heliolinc` when it is run on the paired detection files produced by `make_tracklets`.

**Minimum time interval** `-mintime`: This parameter defaults to 1 second (0.000278 hours). It can be raised to a larger value to avoid making tracklets spanning so little time that, for many objects, negligible motion will be detected. Typically, `heliolinc` cannot link objects with negligible motion because the velocity is essentially undetermined -- so eliminating such objects to begin with could reduce computation time. We have typically not used this option, but it might be valuable when targeting very slow-moving objects. In many cases, the minimum arc length (see below) is likely to be a better way to achieve the same effect.

**Maximum Great Circle Residual** `-maxGCR`: This parameter takes effect only in cases where `make_tracklets` detects overlapping pairs that might indicate more than two detections of the same object, and hence attempts to construct a 'tracklet' with more than two detections. It attempts to fit such a candidate tracklet using a Great Circle trajectory on the sky, at constant angular velocity. If the RMS residual from this Great Circle fit is less than the value specified (in units of arcseconds), the tracklet is considered valid and the pairs are merged. If the residual is too high, the tracklet is considered invalid and its constituent pairs are written to the output file individually. Note that unlimited sharing of detections between different pairs is allowed, but tracklets are entirely exclusive. In other words, if a set of detections are found to make up a valid tracklet, none of them will be included in any other pair (or tracklet). The maximum Great Circle residual defaults to 0.5 arcsec. For surveys with 1-2 arcsecond pixels, larger values up to 2 arcseconds can be reasonable. If for some reason it is desired that **only pairs** should be generated (no tracklets with more than two points), one way of ensuring this is to set the maximum Great Circle residual to an unreasonably small value such as 0.0001 arcseconds.

**Minimum number of points in a tracklet** `-mintrkpts`: This defaults to 2, and values smaller than 2 are of course meaningless. Any value larger than 2 disallows pairs, and will in general greatly reduce the size of the output 'pair file' and hence the runtime of `heliolinc` when this file is used as input. In the context of LSST, the `heliolinc` suite of programs is explicitly designed to enable multi-night linking of pairs, and setting the `-mintrkpts` to 3 or 4 is self-defeating. For surveys or data sets that have more than two images per field per night, however, such settings can greatly reduce the number of candidate objects, the resulting `heliolinc` runtimes, and the purity of the output.

**Minimum angular velocity** `-minvel`: This is the minimum angular velocity, in degrees per day, for a valid pair or tracklet. It defaults to zero, but can be used to exclude extremely slow-moving candidate detections. One reason you might want to do this is if you find that vast numbers of spurious slow-moving candidates are being produced due to false-positive detections caused by stationary stars.  As with the minimum time interval, however, the minimum arc length (see below) may provide a better way to achieve the same effect.

**Maximum angular velocity** `-maxvel`: This defaults to 1.5 degrees per day, and should be set with care because it has an extremely strong effect on the number of pairs and tracklets produced. For a given detection, the area of sky containing other detections that could be pair-partners after a given interval of time scales with the square of the maximum angular velocity. Hence, setting large values of the `-maxvel` can result in a huge increase in spurious pairs. `heliolinc` can handle this, up to a point, but runtimes and file sizes may be unnecessarily increased. An important consideration is at what angular velocity the detections would be noticeably trailed. If trailed-source information is available for the input data, it might make sense to run `make_tracklets` two or more times with different target regimes. For example, run it once with `-maxvel` set to the lowest angular velocity where detections are expected to be clearly trailed -- and then create a new input catalog culled down to **only** the trailed detections, and run `make_tracklets` again with a much faster maximum angular velocity. Note that the margin for object motion, mentioned in the discussion of image radius above, is equal to the maximum angular velocity times the time separation of two images being considered.

**Minimum angular arc** `-minarc`: This defaults to zero, and has units of arcseconds. It refers to the minimum angular separation for the endpoints of a valid pair or tracklet. It is very useful for rejecting spurious pairs/tracklets caused by stationary stars. A value similar to that used for `-maxGCR` -- e.g., two or three times the astrometric precision on the faintest detectable objects -- will often be a sensible setting.
