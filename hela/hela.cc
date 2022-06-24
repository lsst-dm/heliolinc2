#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "makeTracklets.h"

namespace py = pybind11;

PYBIND11_MODULE(hela, m) {
    m.doc() = "HelioLINC Advanced (hela)";

    py::class_<MakeTrackletsConfig>(m, "MakeTrackletsConfig")
        .def(py::init<>())
        .def_readwrite("mintrkpts", &MakeTrackletsConfig::mintrkpts)
        .def_readwrite("imagetimetol", &MakeTrackletsConfig::imagetimetol)
        .def_readwrite("maxvel", &MakeTrackletsConfig::maxvel)
        .def_readwrite("minvel", &MakeTrackletsConfig::minvel)
        .def_readwrite("minarc", &MakeTrackletsConfig::minarc)
        .def_readwrite("maxtime", &MakeTrackletsConfig::maxtime)
        .def_readwrite("mintime", &MakeTrackletsConfig::mintime)
        .def_readwrite("angvel", &MakeTrackletsConfig::angvel)
        .def_readwrite("maxdist", &MakeTrackletsConfig::maxdist)
        .def_readwrite("imagerad", &MakeTrackletsConfig::imagerad)
        .def_readwrite("maxgcr", &MakeTrackletsConfig::maxgcr);
}