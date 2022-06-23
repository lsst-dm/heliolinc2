#
# Makefile for "HelioLinc Advanced"
#
# Builds heliolinc library (libhela.so) and binaries.
#

PREFIX = ..
PROGRAMS_PATH = $(PREFIX)/scripts
PROGRAMS = test_makeTracklets

LIB = libhela.so
LIB_SOURCES = solarsyst_dyn_geo01.cpp makeTracklets.cpp

SOURCES = $(LIB_SOURCES) $(PROGRAMS:%=%.cpp)

.PHONY: all
all: $(PROGRAMS)

.PHONY: install
install: $(LIB) $(PROGRAMS)
	mkdir -p $(PREFIX)/{lib,bin}
	cp -p $(LIB) $(PREFIX)/lib
	cp -p $(PROGRAMS) $(PREFIX)/bin

$(PROGRAMS): %: %.o $(LIB)
	$(CXX) -L. -lhela -Wl,-rpath=$$ORIGIN/../lib -o $@ $<

$(LIB): $(LIB_SOURCES:%.cpp=%.o)
	$(CXX) -shared -fPIC $(LDFLAGS) $(OUTPUT_OPTION) $^

clean:
	rm -f $(PROGRAMS) $(SOURCES:%.cpp=%.o)
	rm -rf .deps

DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -fPIC -std=gnu++17 -I../include -c $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(SOURCES:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))
