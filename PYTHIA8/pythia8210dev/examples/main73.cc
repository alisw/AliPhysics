// main73.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Example how to compare "parton-level" and "hadron-level" properties.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//--------------------------------------------------------------------------

// Generic routine to extract the particles that existed right before
// the hadronization machinery was invoked.

void getPartonLevelEvent( Event& event, Event& partonLevelEvent) {

  // Copy over all particles that existed right before hadronization.
  partonLevelEvent.reset();
  for (int i = 0; i < event.size(); ++i)
  if (event[i].isFinalPartonLevel()) {
    int iNew = partonLevelEvent.append( event[i] );

    // Set copied properties more appropriately: positive status,
    // original location as "mother", and with no daughters.
    partonLevelEvent[iNew].statusPos();
    partonLevelEvent[iNew].mothers( i, i);
    partonLevelEvent[iNew].daughters( 0, 0);
  }

}

//--------------------------------------------------------------------------

// Generic routine to extract the particles that exist after the
// hadronization machinery. Normally not needed, since SlowJet
// contains the standard possibilities preprogrammed, but this
// method illustrates further discrimination.

void getHadronLevelEvent( Event& event, Event& hadronLevelEvent) {

  // Iterate over all final particles.
  hadronLevelEvent.reset();
  for (int i = 0; i < event.size(); ++i) {
    bool accept = false;
    if (event[i].isFinal()) accept = true;

    // Example 1: reject neutrinos (standard option).
    int idAbs = event[i].idAbs();
    if (idAbs == 12 || idAbs == 14 || idAbs == 16) accept = false;

    // Example 2: reject particles with pT < 0.1 GeV (new possibility).
    if (event[i].pT() < 0.1) accept = false;

    // Copy over accepted particles, with original location as "mother".
    if (accept) {
      int iNew = hadronLevelEvent.append( event[i] );
      hadronLevelEvent[iNew].mothers( i, i);
    }
  }

}

//--------------------------------------------------------------------------

int main() {

  // Number of events, generated and listed ones.
  int nEvent    = 1000;
  int nListEvts = 1;
  int nListJets = 5;

  // Generator. LHC process and output selection. Initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 200.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();

  // Parton and Hadron Level event records. Remeber to initalize.
  Event partonLevelEvent;
  partonLevelEvent.init("Parton Level event record", &pythia.particleData);
  Event hadronLevelEvent;
  hadronLevelEvent.init("Hadron Level event record", &pythia.particleData);

  //  Parameters for the jet finders. Need select = 1 to catch partons.
  double radius   = 0.5;
  double pTjetMin = 10.;
  double etaMax   = 4.;
  int select      = 1;

  // Set up anti-kT clustering, comparing parton and hadron levels.
  SlowJet antiKTpartons( -1, radius, pTjetMin, etaMax, select);
  SlowJet antiKThadrons( -1, radius, pTjetMin, etaMax, select);

  // Histograms.
  Hist nJetsP("number of jets, parton level", 50, -0.5, 49.5);
  Hist nJetsH("number of jets, hadron level", 50, -0.5, 49.5);
  Hist pTallP("pT for jets, parton level", 100, 0., 500.);
  Hist pTallH("pT for jets, hadron level", 100, 0., 500.);
  Hist pThardP("pT for hardest jet, parton level", 100, 0., 500.);
  Hist pThardH("pT for hardest jet, hadron level", 100, 0., 500.);
  Hist pTdiff("pT for hardest jet, hadron - parton", 100, -100., 100.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Construct parton and hadron level event.
    getPartonLevelEvent( pythia.event, partonLevelEvent);
    getHadronLevelEvent( pythia.event, hadronLevelEvent);

    // List first few events.
    if (iEvent < nListEvts) {
      pythia.event.list();
      partonLevelEvent.list();
      hadronLevelEvent.list();
    }

    // Analyze jet properties and list first few analyses.
    antiKTpartons.analyze( partonLevelEvent );
    antiKThadrons.analyze( hadronLevelEvent );
    if (iEvent < nListJets) {
      antiKTpartons.list();
      antiKThadrons.list();
    }

    // Fill jet properties distributions.
    nJetsP.fill( antiKTpartons.sizeJet() );
    nJetsH.fill( antiKThadrons.sizeJet() );
    for (int i = 0; i < antiKTpartons.sizeJet(); ++i)
      pTallP.fill( antiKTpartons.pT(i) );
    for (int i = 0; i < antiKThadrons.sizeJet(); ++i)
      pTallH.fill( antiKThadrons.pT(i) );
    if ( antiKTpartons.sizeJet() > 0 && antiKThadrons.sizeJet() > 0) {
      pThardP.fill( antiKTpartons.pT(0) );
      pThardH.fill( antiKThadrons.pT(0) );
      pTdiff.fill( antiKThadrons.pT(0) - antiKTpartons.pT(0) );
    }

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();
  cout << nJetsP << nJetsH << pTallP << pTallH
       <<  pThardP << pThardH << pTdiff;

  // Done.
  return 0;
}
