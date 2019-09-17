# main34.py is a part of the PYTHIA event generator.
# Copyright (C) 2019 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.

# Author: Philip Ilten, March 2016.

# An example where the hard process (p p -> mu+ mu-) is automatically
# produced externally with MadGraph 5, read in, and the remainder of
# the event is then produced by Pythia (MPI, showers, hadronization,
# and decays). A comparison is made between events produced with
# Pythia at LO, MadGraph 5 at LO, and aMC@NLO at NLO.

# For this example to run, MadGraph 5 must be installed and the
# command "exe" (set by default as "mg5_aMC") must be available via
# the command line. Additionally, GZIP support must be enabled via
# the "--with-gzip" configuration option(s). Note that this example has
# only been tested with MadGraph 5 version 2.3.3; due to rapid
# MadGraph development, this example may not work with other
# versions. For more details on the LHAMadgraph class see the
# comments of Pythia8Plugins/LHAMadgraph.h.

# To set the path to the Pythia 8 Python interface do either (in a
# shell prompt):
#      export PYTHONPATH=$(PREFIX_LIB):$PYTHONPATH
# or the following which sets the path from within Python.
import sys
cfg = open('Makefile.inc')
lib = '../lib'
for line in cfg:
    if line.startswith('PREFIX_LIB='): lib = line[11:-1]; break
sys.path.insert(0, lib)
import pythia8

#==========================================================================

# A simple method to run Pythia, analyze the events, and fill a histogram.

def run(pythia, hist, nEvent):
  pythia.readString("Random:setSeed = on")
  pythia.readString("Random:seed = 1")
  pythia.init()
  for iEvent in range(0, nEvent):
      if not pythia.next(): continue
      iMu1 = 0; iMu2 = 0
      for prt in pythia.event:
          if not iMu1 and prt.id() == 13:  iMu1 = prt.index()
          if not iMu2 and prt.id() == -13: iMu2 = prt.index()
          if iMu1 and iMu2:
              iMu1 = pythia.event[iMu1].iBotCopyId()
              iMu2 = pythia.event[iMu2].iBotCopyId()
              hist.fill((pythia.event[iMu1].p() + pythia.event[iMu2].p()).pT())
              break
  pythia.stat()

#==========================================================================

# The name of the MadGraph5_aMC@NLO executable.
# You must prepend this string with the path to the executable
# on your local installation, or otherwise make it available.
exe = "mg5_aMC"

# Create the histograms.
pyPtZ = pythia8.Hist("Pythia dN/dpTZ", 100, 0., 100.)
mgPtZ = pythia8.Hist("MadGraph dN/dpTZ", 100, 0., 100.)
amPtZ = pythia8.Hist("aMC@NLO dN/dpTZ", 100, 0., 100.)

# Produce leading-order events with Pythia.
pythia = pythia8.Pythia()
pythia.readString("Beams:eCM = 13000.")
pythia.readString("WeakSingleBoson:ffbar2gmZ = on")
pythia.readString("23:onMode = off")
pythia.readString("23:onIfMatch = -13 13")
pythia.readString("PhaseSpace:mHatMin = 80.")
run(pythia, pyPtZ, 100)

# Produce leading-order events with MadGraph 5.
pythia = pythia8.Pythia()
madgraph = pythia8.LHAupMadgraph(pythia, True, "madgraphrun", exe)
madgraph.readString("generate p p > mu+ mu-");
# Note the need for a blank character before "set".
madgraph.readString(" set ebeam1 6500")
madgraph.readString(" set ebeam2 6500")
madgraph.readString(" set mmll 80")
pythia.setLHAupPtr(madgraph)
run(pythia, mgPtZ, 1000)

# Produce next-to-leading-order events with aMC@NLO.
pythia = pythia8.Pythia()
amcatnlo = pythia8.LHAupMadgraph(pythia, True, "amcatnlorun", exe)
amcatnlo.readString("generate p p > mu+ mu- [QCD]")
# Note the need for a blank character before "set".
amcatnlo.readString(" set ebeam1 6500")
amcatnlo.readString(" set ebeam2 6500")
amcatnlo.readString(" set mll 80")
pythia.setLHAupPtr(amcatnlo);
run(pythia, amPtZ, 1000)

# Print the histograms.
print(pyPtZ)
print(mgPtZ)
print(amPtZ)
