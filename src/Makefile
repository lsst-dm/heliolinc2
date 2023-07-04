#
# Makefile for "HelioLinc Advanced"
#
# To build all programs (in parallel), run:      make -j4
# To install to ../bin, run:                     sudo make install
# To install to (say) /usr/local, run:           PREFIX=/usr/local make install
#

PREFIX ?= ..
CXXFLAGS ?= -O3

###

PROGRAMS_PATH = $(PREFIX)/bin
PROGRAMS = heliolinc make_tracklets link_refine link_refine_multisite link_refine_Herget cluster2mpc80a pair2mpc80 make_tracklets_new heliolinc_new link_refine_Herget_new parse_clust2det_new parse_clust2det_MPC80 calc_accel_mat02 heliolinc_omp link_refine_omp

LIB = libhela.a
LIB_SOURCES = solarsyst_dyn_geo01.cpp

SOURCES = $(LIB_SOURCES) $(PROGRAMS:%=%.cpp)

.PHONY: all
all: $(PROGRAMS)

.PHONY: install
install: $(PROGRAMS)
	mkdir -p $(PREFIX)/bin
	cp -p $(PROGRAMS) $(PREFIX)/bin

$(PROGRAMS): %: %.o $(LIB)
	$(CXX) $< -L. -lhela -fopenmp $(OUTPUT_OPTION)

$(LIB): $(LIB_SOURCES:%.cpp=%.o)
	$(AR) -rc $@ $^

clean:
	rm -f $(PROGRAMS) $(LIB) $(SOURCES:%.cpp=%.o)
	rm -rf .deps

###

DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -fopenmp -std=c++11 -I../include -c $(OUTPUT_OPTION) $<
	@ $(POSTCOMPILE)

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(SOURCES:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))
