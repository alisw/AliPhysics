// main102.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program provides a demonstration of the flavour rope model supplied
// in the Rope Hadronization framework. It produces four histograms
// showing the ratio of respectively K^0_s, Lambda_0, Cascade and Omega^-
// to pions as function of event multiplicity.
// No kind of Levy-Tsallis fitting, triggering or similar
// is done is this simple demonstration analysis. It should therefore not be
// taken as anything but a proof of concept.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {
  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 7000.");
  pythia.readString("SoftQCD:nonDiffractive = on");
  // Enabling flavour ropes, setting model parameters.
  // The model is still untuned. These parameter values
  // are choosen for illustrative purposes.
  pythia.readString("Ropewalk:RopeHadronization = on");
  pythia.readString("Ropewalk:doShoving = off");
  pythia.readString("Ropewalk:doFlavour = on");
  pythia.readString("Ropewalk:r0 = 0.5");
  pythia.readString("Ropewalk:m0 = 0.2");
  pythia.readString("Ropewalk:beta = 0.1");
  // Enabling setting of vertex information.
  pythia.readString("PartonVertex:setVertex = on");
  // Prevent unstable particles from decaying.
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.init();

  // Histograms.
  Hist pion("pions (mult)", 50, 10., 135.);
  Hist kaon("kaons (mult)", 50, 10., 135.);
  Hist lambda("lambdas (mult)", 50, 10., 135.);
  Hist xi("xi (mult)", 50, 10., 135.);
  Hist omega("omega (mult)", 50, 10., 135.);

  // Note: High statistics is needed to fill the high multiplicity
  // end of the histograms, especially for Omega.
  const int nEvent = 100000;
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    // Counters for particle species.
    int nCharged = 0, nPions = 0, nKaons = 0;
    int nLambdas = 0, nXis = 0, nOmegas = 0;
    // Event short notation.
    Event& event = pythia.event;
    for (int i = 0; i < event.size(); ++i){
        Particle& p = event[i];
        // Apply simple, particle level, cuts.
        if(p.isFinal() && abs(p.eta()) < 2.5 && p.pT() > 0.1 ){
          if(p.isCharged()) ++nCharged;
          int absid = abs(p.id());
          if(absid == 211) ++nPions;
          else if(absid == 310) ++nKaons;
          else if(absid == 3122) ++nLambdas;
          else if(absid == 3312) ++nXis;
          else if(absid == 3334) ++nOmegas;
        }
     }
     // Discard events with event multiplicity less than 10.
     if(nCharged < 10) continue;
     // Fill histograms.
     pion.fill( double(nCharged), double(nPions) );
     kaon.fill( double(nCharged), double(nKaons) );
     lambda.fill( double(nCharged), double(nLambdas) );
     xi.fill( double(nCharged), double(nXis) );
     omega.fill( double(nCharged), double(nOmegas) );
  // End of event loop.
  }
  // Construct ratio histograms.
  Hist kp = kaon / pion;
  kp.title("kaon / pion (multiplicity)");
  Hist lp = lambda / pion;
  lp.title("lambda / pion (multiplicity)");
  Hist xp = xi / pion;
  xp.title("xi / pion (multiplicity)");
  Hist op = omega / pion;
  op.title("omega / pion (multiplicity)");
  // Statistics. Histograms. Done.
  pythia.stat();
  cout << kp << lp << xp << op;
  return 0;
}
