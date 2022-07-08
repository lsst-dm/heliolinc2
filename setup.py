from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "hela",
        ["hela/hela.cc", "src/makeTracklets.cpp"],
        cxx_std=17,
        include_dirs=["include/", "pybind11/include"]
    ),
]

setup(
    name="hela",
    version="0.0.1",
    description="HelioLINC Advanced",
    ext_modules=ext_modules,
)
