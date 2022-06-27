from setuptools import setup

setup(
   name="HelioLinC3D",
   license="BSD 3-Clause License",
   author="Siegfried Eggl",
   author_email="eggl@illinois.edu",
   long_description=open("README.md").read(),
   long_description_content_type="text/markdown",
   url="https://github.com/lsst_dm/heliolinc2",
   packages=["heliolinc3d"],
   package_dir={"heliolinc2": "heliolinc3d"},
   package_data={"heliolinc2": ["tests/data*.tar.gz"]},
   use_scm_version=True,
   setup_requires=["pytest-runner", "setuptools_scm"],
   tests_require=["pytest"],
   install_requires=[
       "astropy",
       "astroquery",
       "joblib",
       "numba",
       "numpy",
       "scipy",
       "scikit-learn",
       "spiceypy",
       "pandas",
       "pytest",
   ],
)
