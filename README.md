# HelioLinC2

### Algorithm:
Based on HelioLinC (Holman et al. 2018) we transform topocentric observations to heliocentric states assuming a distance and radial velocity. The resulting 3D positions are collected into tracklets. Tracklets contain at least two observations and can, thus, be used to create velocity vectors. A tracklet + velocity vector is called an "arrow". Arrows are propagated to a common epoch using spiceypy's 2body propagator, and then clustered using dbscan.




