// main75.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program to study jets in Dark Matter production.

#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;

int main() {

  // Generator. Process selection. Initialization. Event shorthand.
  Pythia pythia;
  pythia.readFile("main75.cmnd");
  pythia.init();
  Event& process = pythia.process;

  // Fastjet analysis - select algorithm and parameters.
  double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition( fastjet::kt_algorithm, Rparam,
           recombScheme, strategy);

  // Fastjet input.
  std::vector <fastjet::PseudoJet> fjInputs;

  // Histograms. Error counter.
  Hist pTj("dN/dpTj", 100, 0., 100.);
  Hist mRec("mRec", 100, 0., 1000.);
  int iErr = 0;

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) {
      if (++iErr < 100) continue;
      else {
        cout << "Too many errors" << endl;
        break;
      }
    }

    // Invariant mass of DM system.
    Vec4 mRes = process[5].p() + process[6].p();
    mRec.fill(mRes.mCalc());

    // Keep track of missing ET.
    Vec4 missingETvec;
    fjInputs.clear();

    // Loop over event record to decide what to pass to FastJet.
    for (int i = 0; i < pythia.event.size(); ++i) {
      // Final state only.
      if (!pythia.event[i].isFinal()) continue;

      // No neutrinos or DM.
      if ( pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14
        || pythia.event[i].idAbs() == 16 || pythia.event[i].idAbs() == 52)
        continue;

      // Only |eta| < 3.6.
      if (abs(pythia.event[i].eta()) > 3.6) continue;

      // Missing ET.
      missingETvec += pythia.event[i].p();

      // Store as input to Fastjet.
      fjInputs.push_back( fastjet::PseudoJet( pythia.event[i].px(),
        pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e() ) );
    }

    // Check that event contains analyzable particles.
    if (fjInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }

    // Run Fastjet algorithm.
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV).
    inclusiveJets = clustSeq.inclusive_jets(20.0);
    sortedJets    = sorted_by_pt(inclusiveJets);

    if(sortedJets.size() < 1) {
      // cout << " No jets found in event " << iEvent << endl;
      continue;
    }
    pTj.fill( sortedJets[0].pt() );

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  cout << pTj;

  // Done.
  return 0;
}
