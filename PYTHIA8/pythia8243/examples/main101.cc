// main101.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program provides a demonstration of the string shoving model supplied
// in the Rope Hadronization framework. It produces four histograms
// shoving the emergence of a "ridge" in two-particle correlations at high
// event multiplicity.
// No kind of background subtraction, event mixing, triggering or similar
// is done is this simple demonstration analysis. It should therefore not be
// taken as anything but a proof of concept.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {
  // Generator. Process selection.
  Pythia pythia;
  pythia.readString("Beams:eCM = 7000.");
  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("Next:numberShowEvent = 0");
  // Enabling string shoving, setting model parameters.
  // The model is still untuned. These parameter values
  // are choosen for illustrative purposes.
  pythia.readString("Ropewalk:RopeHadronization = on");
  pythia.readString("Ropewalk:doShoving = on");
  pythia.readString("Ropewalk:doFlavour = off");
  pythia.readString("Ropewalk:rCutOff = 10.0");
  pythia.readString("Ropewalk:limitMom = on");
  pythia.readString("Ropewalk:pTcut = 2.0");
  pythia.readString("Ropewalk:r0 = 0.41");
  pythia.readString("Ropewalk:m0 = 0.2");
  pythia.readString("Ropewalk:gAmplitude = 10.0");
  pythia.readString("Ropewalk:gExponent = 1.0");
  pythia.readString("Ropewalk:deltat = 0.1");
  pythia.readString("Ropewalk:tShove = 1.");
  pythia.readString("Ropewalk:deltay = 0.1");
  pythia.readString("Ropewalk:tInit = 1.5");
  // Enabling setting of vertex information.
  pythia.readString("PartonVertex:setVertex = on");
  pythia.readString("PartonVertex:protonRadius = 0.7");
  pythia.readString("PartonVertex:emissionWidth = 0.1");
  pythia.init();
  // Histograms.
  Hist deltaPhi1("dPhi, 0 < Nch < 20", 16, -M_PI/2., 3.);
  Hist deltaPhi2("dPhi, 20 < Nch < 40", 16, -M_PI/2., 3.);
  Hist deltaPhi3("dPhi, 40 < Nch < 60", 16, -M_PI/2., 3.);
  Hist deltaPhi4("dPhi, 60 < Nch < 120", 16, -M_PI/2., 3.);
  // Note: High statistics is needed to fill the high multiplicity
  // histogram.
  const int nEvent = 100000;
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    // Event short notation.
    Event& event = pythia.event;
    // First we find the particles we need for the analysis,
    // as well as event multiplicity.
    vector<Particle*> parts;
    int mult = 0;
    for (int i = 0; i < event.size(); ++i){
        // Particle short notation
        Particle& p = event[i];
        // Apply simple, particle level, cuts.
        if(p.isFinal() && p.isCharged() && abs(p.eta()) < 2.5 &&
          p.pT() > 0.5){
                ++mult;
                if(p.pT() > 1.0 && p.pT() < 3.0)
                  parts.push_back(&p);

        }
    }
    // We discard events outside multiplicity bounds.
    int np = parts.size();
    if(mult < 2) continue;
    // We loop over all particle pairs.
    for (int i = 0; i < np; ++i)
      for(int j = 0; j < np; ++j){
        // Skip if same particle.
        if( i == j) continue;
        // The distance in eta between the pair.
        double dEta = abs(parts[i]->eta() - parts[j]->eta());
        if(dEta < 4 && dEta > 2){
          // Calculate the phase difference.
          double dPhi = parts[i]->phi() - parts[j]->phi();
          while (dPhi < -M_PI/2) dPhi += 2*M_PI;
          while (dPhi > 3*M_PI/2) dPhi -= 2*M_PI;
          if(mult <= 20)
            deltaPhi1.fill(dPhi);
          else if(mult <= 40)
            deltaPhi2.fill(dPhi);
          else if(mult <= 60)
            deltaPhi3.fill(dPhi);
          else
            deltaPhi4.fill(dPhi);
        }
      }
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << deltaPhi1 << deltaPhi2 << deltaPhi3 << deltaPhi4;
  return 0;
}
