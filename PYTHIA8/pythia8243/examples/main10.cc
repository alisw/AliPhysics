// main10.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Example how you can use UserHooks to trace pT spectrum through program,
// and veto undesirable jet multiplicities.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Write own derived UserHooks class.

class MyUserHooks : public UserHooks {

public:

  // Constructor creates anti-kT jet finder with (-1, R, pTmin, etaMax).
  MyUserHooks() { slowJet = new SlowJet(-1, 0.7, 10., 5.);
    // Histograms for analysis inside User Hooks.
    pTtrial = new Hist("trial pT spectrum", 100, 0., 400.);
    pTselect = new Hist("selected pT spectrum (before veto)", 100, 0., 400.);
    pTaccept = new Hist("accepted pT spectrum (after veto)", 100, 0., 400.);
    nPartonsB = new Hist("number of partons before veto", 20, -0.5, 19.5);
    nJets = new Hist("number of jets before veto", 20, -0.5, 19.5);
    nPartonsA = new Hist("number of partons after veto", 20, -0.5, 19.5);
    nFSRatISR = new Hist("number of FSR emissions at first ISR emission",
        20, -0.5, 19.5);
  }

  // Destructor deletes anti-kT jet finder and prints histograms.
  ~MyUserHooks() {delete slowJet; cout << *pTtrial << *pTselect << *pTaccept
       << *nPartonsB << *nJets << *nPartonsA << *nFSRatISR; }

  // Allow process cross section to be modified...
  virtual bool canModifySigma() {return true;}

  // ...which gives access to the event at the trial level, before selection.
  virtual double multiplySigmaBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool inEvent) {

    // All events should be 2 -> 2, but kill them if not.
    if (sigmaProcessPtr->nFinal() != 2) return 0.;

    // Extract the pT for 2 -> 2 processes in the event generation chain
    // (inEvent = false for initialization).
    if (inEvent) {
      pTHat = phaseSpacePtr->pTHat();
      // Fill histogram of pT spectrum.
      pTtrial->fill( pTHat );
    }

    // Here we do not modify 2 -> 2 cross sections.
    return 1.;
  }

  // Allow a veto for the interleaved evolution in pT.
  virtual bool canVetoPT() {return true;}

  // Do the veto test at a pT scale of 5 GeV.
  virtual double scaleVetoPT() {return 5.;}

  // Access the event in the interleaved evolution.
  virtual bool doVetoPT(int iPos, const Event& event) {

    // iPos <= 3 for interleaved evolution; skip others.
    if (iPos > 3) return false;

    // Fill histogram of pT spectrum at this stage.
    pTselect->fill(pTHat);

    // Extract a copy of the partons in the hardest system.
    subEvent(event);
    nPartonsB->fill( workEvent.size() );

    // Find number of jets with given conditions.
    slowJet->analyze(event);
    int nJet = slowJet->sizeJet();
    nJets->fill( nJet );

    // Veto events which do not have exactly three jets.
    if (nJet != 3) return true;

    // Statistics of survivors.
    nPartonsA->fill( workEvent.size() );
    pTaccept->fill(pTHat);

    // Do not veto events that got this far.
    return false;

  }

  // Allow a veto after (by default) first step.
  virtual bool canVetoStep() {return true;}

  // Access the event in the interleaved evolution after first step.
  virtual bool doVetoStep( int iPos, int nISR, int nFSR, const Event& ) {

    // Only want to study what happens at first ISR emission
    if (iPos == 2 && nISR == 1) nFSRatISR->fill( nFSR );

    // Not intending to veto any events here.
    return false;
  }

private:

  // The anti-kT (or kT, C/A) jet finder.
  SlowJet* slowJet;

  // Save the pThat scale.
  double pTHat;

  // The list of histograms.
  Hist *pTtrial, *pTselect, *pTaccept, *nPartonsB, *nJets, *nPartonsA,
       *nFSRatISR;

};

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  //  Process selection. No need to study hadron level.
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 50.");
  pythia.readString("HadronLevel:all = off");

  // Set up to do a user veto and send it in.
  MyUserHooks* myUserHooks = new MyUserHooks();
  pythia.setUserHooksPtr( myUserHooks);

  // Tevatron initialization.
  pythia.readString("Beams:idB = -2212");
  pythia.readString("Beams:eCM = 1960.");
  pythia.init();

  // Begin event loop.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {

    // Generate events.
    pythia.next();

  // End of event loop.
  }

  // Statistics.
  pythia.stat();

  // Done.
  delete myUserHooks;
  return 0;
}
