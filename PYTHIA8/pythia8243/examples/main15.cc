// main15.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how either
// (a) B decays (sections marked "Repeated decays:), or
// (b) all hadronization (sections marked "Repeated hadronization:")
// could be repeated a number of times for each event,
// to improve statistics when this could be a problem.
// Option (a) is faster than (b), but less generic.

// Note 1: the compartmentalization of hadronization in forceHadronLevel
// from the rest of the event processing somewhat limits the ways the
// program can retry in case of problems, and so an occasional abort
// may occur more easily than normally.

// Note 2: for simple cases, where it is only one particle that is to be
// decayed repeatedly, the event[i].undoDecay() method is handy.
// When used for several particles, remember that the position of
// some particles may be moved by the undoDecay step.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Main switches: redo B decays only or redo all hadronization, but not both.
  bool redoBDecays = false;
  bool redoHadrons = true;
  if (redoHadrons) redoBDecays = false;

  // Number of events. Number to list redone events.
  int nEvent = 100;
  int nListRedo = 1;

  // Number of times decays/hadronization should be redone for each event.
  int nRepeat = 10;
  if (!redoBDecays && !redoHadrons) nRepeat = 1;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Simulate b production above given pTmin scale.
  // Warning: these processes do not catch all possible production modes.
  // You would need to use HardQCD:all or even SoftQCD:nonDiffractive for that.
  pythia.readString("HardQCD:gg2bbbar = on");
  pythia.readString("HardQCD:qqbar2bbbar = on");
  pythia.readString("PhaseSpace:pTHatMin = 50.");

  // Repeated decays: list of weakly decaying B hadrons.
  // Note: this list is overkill; some will never be produced.
  int bCodes[28] = {511, 521, 531, 541, 5122, 5132, 5142, 5232, 5242,
    5332, 5342, 5412, 5414, 5422, 5424, 5432, 5434, 5442, 5444, 5512,
    5514, 5522, 5524, 5532, 5534, 5542, 5544, 5544 };
  int nCodes = 28;

  // Repeated decays: location of B handrons.
  vector<int> iBHad;
  int nBHad = 0;

  // Repeated hadronization: spare copy of event.
  Event savedEvent;

  // Repeated hadronization: switch off normal HadronLevel call.
  if (redoHadrons) pythia.readString("HadronLevel:all = off");

  // Initialize for LHC energies; default 14 TeV
  pythia.init();

  // Histogram invariant mass of muon pairs.
  Hist nBperEvent("number of b quarks in an event", 10, -0.5, 9.5);
  Hist nSameEvent("number of times same event is used", 10, -0.5, 9.5);
  Hist oppSignMass("mass of opposite-sign muon pair", 100, 0.0, 100.0);
  Hist sameSignMass("mass of same-sign muon pair", 100, 0.0, 100.0);

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Repeated decays: switch off decays of weakly decaying B hadrons.
    // (More compact solution than repeated readString(..).)
    if (redoBDecays) for (int iC = 0; iC < nCodes; ++iC)
      pythia.particleData.mayDecay( bCodes[iC], false);

    // Generate event. Skip it if error.
    if (!pythia.next()) continue;

    // Find and histogram number of b quarks.
    int nBquark = 0;
    int stat;
    for (int i = 0; i < event.size(); ++i) {
      stat = event[i].statusAbs();
      if (event[i].idAbs() == 5 && (stat == 62 || stat == 63)) ++nBquark;
    }
    nBperEvent.fill( nBquark );

    // Repeated decays: find all locations where B hadrons are stored.
    if (redoBDecays) {
      iBHad.resize(0);
      for (int i = 0; i < event.size(); ++i) {
        int idAbs = event[i].idAbs();
        for (int iC = 0; iC < 28; ++iC)
        if (idAbs == bCodes[iC]) {
          iBHad.push_back(i);
          break;
        }
      }

      // Repeated decays: check that #b = #B.
      nBHad = iBHad.size();
      if (nBquark != nBHad) cout << " Warning: " << nBquark
        << " b quarks but " << nBHad << " B hadrons" << endl;

      // Repeated decays: store size of current event.
      event.saveSize();

      // Repeated decays: switch back on weakly decaying B hadrons.
      for (int iC = 0; iC < nCodes; ++iC)
        pythia.particleData.mayDecay( bCodes[iC], true);

    //  Repeated hadronization: copy event into spare position.
    } else if (redoHadrons) {
      savedEvent = event;
    }

    // Begin loop over rounds of decays / hadronization for same event.
    int nWithPair = 0;
    for (int iRepeat = 0; iRepeat < nRepeat; ++iRepeat) {

      // Repeated decays: remove B decay products from previous round.
      if (redoBDecays) {
        if (iRepeat > 0) {
          event.restoreSize();

          // Repeated decays: mark decayed B hadrons as undecayed.
          for (int iB = 0; iB < nBHad; ++iB) event[ iBHad[iB] ].statusPos();
        }

        // Repeated decays: do decays of B hadrons, sequentially for products.
        // Note: modeDecays does not work for bottomonium (or heavier) states,
        // since there decays like Upsilon -> g g g also need hadronization.
        // Also, there is no provision for Bose-Einstein effects.
        if (!pythia.moreDecays()) continue;


      // Repeated hadronization: restore saved event record.
      } else if (redoHadrons) {
        if (iRepeat > 0) event = savedEvent;

        // Repeated hadronization: do HadronLevel (repeatedly).
        // Note: argument false needed owing to bug in junction search??
        if (!pythia.forceHadronLevel(false)) continue;
      }

      // List last repetition of first few events.
      if ( (redoBDecays || redoHadrons) && iEvent < nListRedo
        && iRepeat == nRepeat - 1) event.list();

      // Look for muons among decay products (also from charm/tau/...).
      vector<int> iMuNeg, iMuPos;
      for (int i = 0; i < event.size(); ++i) {
        int id = event[i].id();
        if (id ==  13) iMuNeg.push_back(i);
        if (id == -13) iMuPos.push_back(i);
      }

      // Check whether pair(s) present.
      int nMuNeg = iMuNeg.size();
      int nMuPos = iMuPos.size();
      if (nMuNeg + nMuPos > 1) {
        ++nWithPair;

        // Fill masses of opposite-sign pairs.
        for (int iN = 0; iN < nMuNeg; ++iN)
        for (int iP = 0; iP < nMuPos; ++iP)
          oppSignMass.fill(
            (event[iMuNeg[iN]].p() + event[iMuPos[iP]].p()).mCalc() );

        // Fill masses of same-sign pairs.
        for (int i1 = 0; i1 < nMuNeg - 1; ++i1)
        for (int i2 = i1 + 1; i2 < nMuNeg; ++i2)
          sameSignMass.fill(
            (event[iMuNeg[i1]].p() + event[iMuNeg[i2]].p()).mCalc() );
        for (int i1 = 0; i1 < nMuPos - 1; ++i1)
        for (int i2 = i1 + 1; i2 < nMuPos; ++i2)
          sameSignMass.fill(
            (event[iMuPos[i1]].p() + event[iMuPos[i2]].p()).mCalc() );

      // Finished analysis of current round.
      }

    // End of loop over many rounds. fill number of rounds with pairs.
    }
    nSameEvent.fill( nWithPair );

  // End of event loop.
  }

  // Statistics. Histograms.
  pythia.stat();
  cout << nBperEvent << nSameEvent << oppSignMass << sameSignMass << endl;

  // Done.
  return 0;
}
