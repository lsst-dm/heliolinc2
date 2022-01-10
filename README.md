Creation and testing of modular pieces for C++ implementation of heliolinc3D.

maketrack02b.cpp: Perform image-based pairing to create tracklets.

readpairs01a.cpp: Read pairing output files to realize tracklets.

maketrack03a.cpp: Updated tracklet creation, stores observer barycentric
                  positions and string IDs enabling the mapping of created
		  tracklets to specific simulated objects.

projectpairs04b.cpp: First complete implementation of heliolinc in C++.
		     Requires input files produced by maketrack03a.

projectpairs04c.cpp: Experimental improvement on projectpairs04b, using
		     integerized versions of state vectors to speed searching.

Suggested compile commands:
c++ -O3 -o maketrack03a maketrack03a.cpp solarsyst_dyn_geo01.cpp -std=c++11
c++ -O3 -o projectpairs04b projectpairs04b.cpp solarsyst_dyn_geo01.cpp -std=c++11
c++ -O3 -o projectpairs04c projectpairs04c.cpp solarsyst_dyn_geo01.cpp -std=c++11

