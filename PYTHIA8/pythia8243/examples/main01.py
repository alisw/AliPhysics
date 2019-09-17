# main01.py is a part of the PYTHIA event generator.
# Copyright (C) 2019 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
#
# This is a simple test program. It fits on one slide in a talk.  It
# studies the charged multiplicity distribution at the LHC. To set the
# path to the Pythia 8 Python interface do either (in a shell prompt):
#      export PYTHONPATH=$(PREFIX_LIB):$PYTHONPATH
# or the following which sets the path from within Python.
#
# Use "python-config --include" to find the include directory and
# then configure Pythia ""--with-python-include=*".

import sys
cfg = open("Makefile.inc")
lib = "../lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Import the Pythia module.
import pythia8
pythia = pythia8.Pythia()
pythia.readString("Beams:eCM = 8000.")
pythia.readString("HardQCD:all = on")
pythia.readString("PhaseSpace:pTHatMin = 20.")
pythia.init()
mult = pythia8.Hist("charged multiplicity", 100, -0.5, 799.5)
# Begin event loop. Generate event. Skip if error. List first one.
for iEvent in range(0, 100):
    if not pythia.next(): continue
    # Find number of all final charged particles and fill histogram.
    nCharged = 0
    for prt in pythia.event:
        if prt.isFinal() and prt.isCharged(): nCharged += 1
    mult.fill(nCharged)
# End of event loop. Statistics. Histogram. Done.
pythia.stat();
print(mult)
