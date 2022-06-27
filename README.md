# HelioLinC3D

### Algorithm:
Based on HelioLinC (Holman et al. 2018) we transform topocentric observations to heliocentric states assuming a distance and radial velocity. The resulting 3D positions are collected into tracklets. Tracklets contain at least two observations and can, thus, be used to create velocity vectors. A tracklet + velocity vector is called an "arrow". Arrows are propagated to a common epoch using spiceypy's 2body propagator, and then clustered using dbscan.

### Getting Started:
To get started, please create a conda environment using the environment.yml file

`conda env create -f environment.yml`

This will create a conda environment named "heliolinc3d" that contains all necessary dependencies.
The environment can be activated via

`conda activate heliolinc3d`


### Install HelioLinC3D:
The package is pip installable.
Activate the corresponding conda environment. 
`conda activate heliolinc3d`

Use pip to install.
`pip install .`

### Run the demo notebook:
In order to run the demo notebook please activate the conda environment named heliolinc3d

`conda activate heliolinc3d`

and add it as a Jupyter notebook kernel to be used with the demo .ipynb files

`python -m ipykernel install --user --name heliolinc3d`

The last step only needs to be done once.

Open jupyterhub. Open the jupyter notebook demo. Select the 'heliolinc3d' kernel from the list of available kernels in the 'Change Kernel' menu.








