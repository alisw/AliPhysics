// main74.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Example how to use the modified Mass Drop Tagger on Pythia jets.
// Note: to run this you must install and link the FastJet Contrib add-ons.

// Pythia include and namespace.
#include "Pythia8/Pythia.h"
using namespace Pythia8;

// FastJet include and namespace.
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"
using namespace fastjet;

int main() {

  // Number of events.
  int nEvent = 1000;

  // Set up Pythia generation of Z + jet; Z -> hadrons; m_Z restricted.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("WeakBosonAndParton:qqbar2gmZg = on");
  pythia.readString("WeakBosonAndParton:qg2gmZq = on");
  pythia.readString("PhaseSpace:pTHatMin = 400.");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");
  pythia.readString("23:mMin = 70.");
  pythia.readString("23:mMax = 120.");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();

  // Detector size, anti-kT radius, and modified mass-drop tagger z.
  double etaMax = 5.;
  double radius = 1.;
  double z_cut  = 0.04;

  // Set up FastJet jet finders and modified mass-drop tagger.
  JetDefinition jetDefAKT( antikt_algorithm, radius);
  JetDefinition jetDefCA( cambridge_algorithm, JetDefinition::max_allowable_R);
  contrib::ModifiedMassDropTagger mMDT(z_cut);

  // Histograms for Z mass: truth, before and after mass drop.
  Hist rZjet( "R separation true vs. reconstructed Z", 100, 0., 1.);
  Hist mTrue(   "Z0 mass as generated", 100, 0., 200.);
  Hist mBefDrop( "Z0 mass before mMDT", 100, 0., 200.);
  Hist mAftDrop( "Z0 mass after mMDT",  100, 0., 200.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Store final visible central particle four-momenta as start
    // configuration. Also find last copy 0f Z0, i.e. right before decay.
    vector<PseudoJet> particles;
    int iZ = 0;
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].isFinal() && event[i].isVisible()
        && abs(event[i].eta()) < etaMax) particles.push_back( PseudoJet(
        event[i].px(), event[i].py(), event[i].pz(), event[i].e() ) );
      if (event[i].id() == 23) iZ = i;
    }

    // Run Fastjet anti-kT algorithm and sort jets in pT order.
    ClusterSequence clustSeq1( particles, jetDefAKT );
    vector<PseudoJet> sortedJets = sorted_by_pt( clustSeq1.inclusive_jets() );

    // Z should be close to either of two hardest jets (in R space).
    if (sortedJets.size() < 2) continue;
    double y0Z   = sortedJets[0].rap() - event[iZ].y();
    double phi0Z = abs(sortedJets[0].phi_std() - event[iZ].phi());
    if (phi0Z > M_PI) phi0Z = 2. * M_PI - phi0Z;
    double r0Z   = sqrt( pow2(y0Z) + pow2(phi0Z) );
    double y1Z   = sortedJets[1].rap() - event[iZ].y();
    double phi1Z = abs(sortedJets[1].phi_std() - event[iZ].phi());
    if (phi1Z > M_PI) phi1Z = 2. * M_PI - phi1Z;
    double r1Z   = sqrt( pow2(y1Z) + pow2(phi1Z) );
    if (min( r0Z, r1Z) > 1.) continue;
    int iJet     = (r1Z > r0Z) ? 0 : 1;

    // Extract Z0-associated jet and run C/A on it. Should give one jet.
    vector<PseudoJet> constituents = sortedJets[iJet].constituents();
    ClusterSequence clustSeq2( constituents, jetDefCA );
    vector<PseudoJet> subJets = sorted_by_pt( clustSeq2.inclusive_jets() );
    if (subJets.size() > 1) continue;

    // Use modified mass-drop tagger to clean up jet.
    PseudoJet reclusteredJet = subJets[0];
    PseudoJet taggedJet = mMDT(reclusteredJet);

    // Fill histograms.
    rZjet.fill( min( r0Z, r1Z) );
    mTrue.fill( event[iZ].m() );
    mBefDrop.fill( reclusteredJet.m() );
    mAftDrop.fill( taggedJet.m() );
  }

  // End of event loop. Statistics. Histograms. Done.
  pythia.stat();
  cout << rZjet << mTrue << mBefDrop << mAftDrop;
  return 0;
}
