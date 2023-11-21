SOURCE CODE:
solarsyst_dyn_geo01.cpp
solarsyst_dyn_geo01.h
make_tracklets_new.cpp
heliolinc_danby.cpp
link_refine_Herget_univar.cpp


FILE LIST:

Earth1day2020s_02a.csv
ObsCodesNew.txt
radhyp_test01.txt
testin_mt01.csv
outim_test01_compare.txt
clust2det_test01_compare.csv
LRHclust2det_test01_compare.csv
LRHsum_test01_compare.csv
pairdet_test01_compare.csv
sum_test01_compare.csv
tracklets_test01_compare.csv
trk2det_test01_compare.csv
colformat_LSST_01.txt

COMPILE:

c++ -O3 -Wall -o make_tracklets_new make_tracklets_new.cpp solarsyst_dyn_geo01.cpp -fopenmp -std=c++11
c++ -O3 -Wall -o heliolinc_danby heliolinc_danby.cpp solarsyst_dyn_geo01.cpp -fopenmp -std=c++11
c++ -O3 -Wall -o link_refine_Herget_univar link_refine_Herget_univar.cpp solarsyst_dyn_geo01.cpp -fopenmp -std=c++11


RUN MAKE_TRACKLETS:

time make_tracklets_new -dets testin_mt01.csv -outim outim_test01.txt -pairdets pairdet_test01.csv -tracklets tracklets_test01.csv -trk2det trk2det_test01.csv -colformat colformat_LSST_01.txt -imrad 2.0 -maxtime 2.0 -maxGCR 1.5 -maxvel 1.5 -minarc 1.0 -earth Earth1day2020s_02a.csv -obscode ObsCodesNew.txt

Writing paired detection file with 273 lines
Writing tracklet file with 8 lines
Writing trk2det file with 271 lines

real	0m0.057s
user	0m0.045s
sys	0m0.001s

TEST OUTPUT:

diff outim_test01.txt outim_test01_compare.txt 
diff pairdet_test01.csv pairdet_test01_compare.csv
diff tracklets_test01.csv tracklets_test01_compare.csv
diff trk2det_test01.csv trk2det_test01_compare.csv


RUN HELIOLINC:

time heliolinc_danby -imgs outim_test01.txt -pairdets pairdet_test01.csv -tracklets tracklets_test01.csv -trk2det trk2det_test01.csv -mjd 60607.74 -obspos Earth1day2020s_02a.csv -heliodist radhyp_test01.txt -outsum sum_test01.csv -clust2det clust2det_test01.csv

Writing 2 lines to output cluster-summary file sum_test01.csv
Writing 542 lines to output clust2det file clust2det_test01.csv

real	0m0.010s
user	0m0.010s
sys	0m0.000s

TEST OUTPUT:

diff sum_test01.csv sum_test01_compare.csv
diff clust2det_test01.csv clust2det_test01_compare.csv


RUN LINK_REFINE:

printf "sum_test01.csv clust2det_test01.csv\n" > clusterlist_test01

time link_refine_Herget_univar -imgs outim_test01.txt -pairdets pairdet_test01.csv -lflist clusterlist_test01 -mjd 60607.74 -outsum LRHsum_test01.csv -clust2det LRHclust2det_test01.csv

Writing 1 lines to output cluster-summary file LRHsum_test01.csv
Writing 271 lines to output clust2det file LRHclust2det_test01.csv

real	0m0.069s
user	0m0.057s
sys	0m0.001s

TEST OUTPUT:

diff LRHsum_test01.csv LRHsum_test01_compare.csv
diff LRHclust2det_test01.csv LRHclust2det_test01_compare.csv

