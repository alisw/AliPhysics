# main10.py is a part of the PYTHIA event generator.
# Copyright (C) 2019 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.

# Example how you can use UserHooks to trace pT spectrum through
# program, and veto undesirable jet multiplicities. To set the path to
# the Pythia 8 Python interface do either (in a shell prompt):
#      export PYTHONPATH=$(PREFIX_LIB):$PYTHONPATH
# or the following which sets the path from within Python.
import sys
cfg = open("Makefile.inc")
lib = "../lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)
import pythia8

#==========================================================================

# Put histograms here to make them global, so they can be used both
# in MyUserHooks and in the main program.

pTtrial   = pythia8.Hist("trial pT spectrum", 100, 0., 400.)
pTselect  = pythia8.Hist("selected pT spectrum (before veto)", 100, 0., 400.)
pTaccept  = pythia8.Hist("accepted pT spectrum (after veto)", 100, 0., 400.)
nPartonsB = pythia8.Hist("number of partons before veto", 20, -0.5, 19.5)
nJets     = pythia8.Hist("number of jets before veto", 20, -0.5, 19.5)
nPartonsA = pythia8.Hist("number of partons after veto", 20, -0.5, 19.5)
nFSRatISR = pythia8.Hist("number of FSR emissions at first ISR emission",
                         20, -0.5, 19.5)

#==========================================================================

# Write own derived UserHooks class.

class MyUserHooks(pythia8.UserHooks):

    # Constructor creates anti-kT jet finder with (-1, R, pTmin, etaMax).
    def __init__(self):
        pythia8.UserHooks.__init__(self)
        self.slowJet = pythia8.SlowJet(-1, 0.7, 10., 5.)
        self.pTHat   = 0.

    # Allow process cross section to be modified...
    def canModifySigma(self): return True

    # ...which gives access to the event at the trial level, before selection.
    def multiplySigmaBy(self, sigmaProcessPtr, phaseSpacePtr, inEvent):

        # All events should be 2 -> 2, but kill them if not.
        if sigmaProcessPtr.nFinal() != 2: return 0.

        # Extract the pT for 2 -> 2 processes in the event generation chain
        # (inEvent = false for initialization).
        if inEvent:
            self.pTHat = phaseSpacePtr.pTHat()
            # Fill histogram of pT spectrum.
            pTtrial.fill(self.pTHat)

        # Here we do not modify 2 -> 2 cross sections.
        return 1.

    # Allow a veto for the interleaved evolution in pT.
    def canVetoPT(self): return True

    # Do the veto test at a pT scale of 5 GeV.
    def scaleVetoPT(self): return 5.

    # Access the event in the interleaved evolution.
    def doVetoPT(self, iPos, event):

        # iPos <= 3 for interleaved evolution; skip others.
        if iPos > 3: return False

        # Fill histogram of pT spectrum at this stage.
        pTselect.fill(self.pTHat)

        # Extract a copy of the partons in the hardest system.
        self.subEvent(event)
        nPartonsB.fill(self.workEvent.size())

        # Find number of jets with given conditions.
        self.slowJet.analyze(event);
        nJet = self.slowJet.sizeJet()
        nJets.fill(nJet)

        # Veto events which do not have exactly three jets.
        if nJet != 3: return True

        # Statistics of survivors.
        nPartonsA.fill(self.workEvent.size())
        pTaccept.fill(self.pTHat)

        # Do not veto events that got this far.
        return False

    # Allow a veto after (by default) first step.
    def canVetoStep(self): return True

    # Access the event in the interleaved evolution after first step.
    def doVetoStep(self, iPos, nISR, nFSR, event):

        # Only want to study what happens at first ISR emission
        if iPos == 2 and nISR == 1: nFSRatISR.fill(nFSR)

        # Not intending to veto any events here.
        return False

#==========================================================================

# Generator.
pythia = pythia8.Pythia()

#  Process selection. No need to study hadron level.
pythia.readString("HardQCD:all = on")
pythia.readString("PhaseSpace:pTHatMin = 50.")
pythia.readString("HadronLevel:all = off")

# Set up to do a user veto and send it in.
myUserHooks = MyUserHooks()
pythia.setUserHooksPtr(myUserHooks)

# Tevatron initialization.
pythia.readString("Beams:idB = -2212")
pythia.readString("Beams:eCM = 1960.")
pythia.init()

# Generate events.
for iEvent in range(0, 1000): pythia.next();

# Statistics. Histograms.
pythia.stat()
print(pTtrial)
print(pTselect)
print(pTaccept)
print(nPartonsB)
print(nJets)
print(nPartonsA)
print(nFSRatISR)
