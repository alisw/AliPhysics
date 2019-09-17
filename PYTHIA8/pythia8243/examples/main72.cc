// main72.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It compares SlowJet, FJcore and FastJet, showing that they
// find the same jets.

#include "Pythia8/Pythia.h"

// The FastJet3.h header enables automatic initialisation of
// fastjet::PseudoJet objects from Pythia8 Particle and Vec4 objects,
// as well as advanced features such as access to (a copy of)
// the original Pythia 8 Particle directly from the PseudoJet,
// and fastjet selectors that make use of the Particle properties.
// See the extensive comments in the header file for further details
// and examples.
#include "Pythia8Plugins/FastJet3.h"

using namespace Pythia8;

int main() {

  // Number of events, generated and listed ones (for jets).
  int nEvent    = 1000;
  int nListJets = 3;

  // Select common parameters for SlowJet and FastJet analyses.
  int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
  double R       = 0.7;    // Jet size.
  double pTMin   = 5.0;    // Min jet pT.
  double etaMax  = 5.0;    // Pseudorapidity range of detector.
  int    select  = 2;      // Which particles are included?
  int    massSet = 2;      // Which mass are they assumed to have?

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Process selection.
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 200.");

  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // LHC initialization.
  pythia.readString("Beams:eCM = 14000.");
  pythia.init();

  // Set up SlowJet jet finder in native mode.
  SlowJet slowJet( power, R, pTMin, etaMax, select, massSet, 0, false);

  // Set up SlowJet jet finder as a wrapper to the fjcore package.
  // Note that this is now the default SlowJet constructor choice.
  SlowJet fjCore( power, R, pTMin, etaMax, select, massSet, 0, true);

  // Set up FastJet jet finder.
  //   one can use either explicitly use antikt, cambridge, etc., or
  //   just use genkt_algorithm with specification of power
  //fastjet::JetAlgorithm algorithm;
  //if (power == -1)      algorithm = fastjet::antikt_algorithm;
  //if (power ==  0)      algorithm = fastjet::cambridge_algorithm;
  //if (power ==  1)      algorithm = fastjet::kt_algorithm;
  //fastjet::JetDefinition jetDef(algorithm, R);
  // there's no need for a pointer to the jetDef (it's a fairly small object)
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);
  std::vector <fastjet::PseudoJet> fjInputs;

  // Histograms.
  Hist nJetsS("number of jets (SlowJet)", 100, -0.5, 99.5);
  Hist nJetsC("number of jets, FJcore - SlowJet ", 99, -49.5, 49.5);
  Hist nJetsF("number of jets, FastJet - SlowJet", 99, -49.5, 49.5);
  Hist pTjetsS("pT for jets (SlowJet)", 100, 0., 500.);
  Hist pTjetsC("pT for jets, FJcore - SlowJet ", 100, -10., 10.);
  Hist pTjetsF("pT for jets, FastJet - SlowJet", 100, -10., 10.);
  Hist etaJetsS("eta for jets (SlowJet)", 100, -5., 5.);
  Hist phiJetsS("phi for jets (SlowJet)", 100, -M_PI, M_PI);
  Hist RdistC("R distance FJcore to SlowJet ", 100, 0., 1.);
  Hist RdistF("R distance FastJet to SlowJet", 100, 0., 1.);
  Hist distJets("R distance between jets", 100, 0., 10.);
  Hist pTdiff("pT difference between consecutive jets", 100, -100., 400.);
  Hist nAna("multiplicity of analyzed event", 100, -0.5, 999.5);
  Hist tGen("generation time as fn of multiplicity", 100, -0.5, 999.5);
  Hist tSlow("SlowJet time as fn of multiplicity", 100, -0.5, 999.5);
  Hist tCore("FJcore time as fn of multiplicity ", 100, -0.5, 999.5);
  Hist tFast("FastJet time as fn of multiplicity", 100, -0.5, 999.5);
  Hist tSlowGen("SlowJet/generation time as fn of multiplicity",
    100, -0.5, 999.5);
  Hist tCoreGen("FJcore/generation time as fn of multiplicity",
    100, -0.5, 999.5);
  Hist tFastGen("FastJet/generation time as fn of multiplicity",
    100, -0.5, 999.5);
  Hist tFastCore("FastJet/FJcore time as fn of multiplicity",
    100, -0.5, 999.5);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    clock_t befGen = clock();
    if (!pythia.next()) continue;
    clock_t aftGen = clock();

    // Begin SlowJet analysis of jet properties. List first few.
    clock_t befSlow = clock();
    slowJet.analyze( pythia.event );
    clock_t aftSlow = clock();
    if (iEvent < nListJets) slowJet.list();

    // Fill inclusive SlowJet jet distributions.
    int nSlow = slowJet.sizeJet();
    nJetsS.fill( nSlow );
    for (int i = 0; i < nSlow; ++i) {
      pTjetsS.fill( slowJet.pT(i) );
      etaJetsS.fill( slowJet.y(i) );
      phiJetsS.fill( slowJet.phi(i) );
    }

    // Fill distance between SlowJet jets.
    for (int i = 0; i < nSlow - 1; ++i)
    for (int j = i + 1; j < nSlow; ++j) {
      double dY = slowJet.y(i)  - slowJet.y(j);
      double dPhi = abs( slowJet.phi(i) - slowJet.phi(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dY) + pow2(dPhi) );
      distJets.fill( dR );
    }

    // Fill pT-difference between SlowJet jets (to check ordering of list).
    for (int i = 1; i < nSlow; ++i)
      pTdiff.fill( slowJet.pT(i-1)- slowJet.pT(i) );

    // Begin FJcore analysis of jet properties. List first few.
    clock_t befCore = clock();
    fjCore.analyze( pythia.event );
    clock_t aftCore = clock();
    if (iEvent < nListJets) fjCore.list();

    // Fill distribution of fjCore jets relative to SlowJet ones.
    int nCore = fjCore.sizeJet();
    nJetsC.fill( nCore - nSlow);
    if (nCore == nSlow) {
      for (int i = 0; i < nCore; ++i) {
        pTjetsC.fill( fjCore.pT(i) - slowJet.pT(i) );
        double dist2 = pow2( fjCore.y(i) - slowJet.y(i))
          + pow2( fjCore.phi(i) - slowJet.phi(i) );
        RdistC.fill( sqrt(dist2) );
      }
    }

    // Begin FastJet analysis: extract particles from event record.
    clock_t befFast = clock();
    fjInputs.resize(0);
    Vec4   pTemp;
    double mTemp;
    int nAnalyze = 0;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

      // Require visible/charged particles inside detector.
      if      (select > 2 &&  event[i].isNeutral() ) continue;
      else if (select == 2 && !event[i].isVisible() ) continue;
      if (etaMax < 20. && abs(event[i].eta()) > etaMax) continue;

      // Create a PseudoJet from the complete Pythia particle.
      fastjet::PseudoJet particleTemp = event[i];

      // Optionally modify mass and energy.
      pTemp = event[i].p();
      mTemp = event[i].m();
      if (massSet < 2) {
        mTemp = (massSet == 0 || event[i].id() == 22) ? 0. : 0.13957;
        pTemp.e( sqrt(pTemp.pAbs2() + mTemp*mTemp) );
        particleTemp.reset_momentum( pTemp.px(), pTemp.py(),
           pTemp.pz(), pTemp.e() );
      }

      // Store acceptable particles as input to Fastjet.
      // Conversion to PseudoJet is performed automatically
      // with the help of the code in FastJet3.h.
      fjInputs.push_back( particleTemp);
      ++nAnalyze;
    }

    // Run Fastjet algorithm and sort jets in pT order.
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets = clustSeq.inclusive_jets(pTMin);
    sortedJets    = sorted_by_pt(inclusiveJets);
    clock_t aftFast = clock();

    // List first few FastJet jets and some info about them.
    // Note: the final few columns are illustrative of what information
    // can be extracted, but does not exhaust the possibilities.
    if (iEvent < nListJets) {
      cout << "\n --------  FastJet jets, p = " << setw(2) << power
           << "  --------------------------------------------------\n\n "
           << "  i         pT        y      phi  mult chgmult photons"
           << "      hardest  pT in neutral " << endl
           << "                                                       "
           << "  constituent        hadrons " << endl;
      for (int i = 0; i < int(sortedJets.size()); ++i) {
        vector<fastjet::PseudoJet> constituents
          = sortedJets[i].constituents();
        fastjet::PseudoJet hardest
          = fastjet::SelectorNHardest(1)(constituents)[0];
        vector<fastjet::PseudoJet> neutral_hadrons
          = ( fastjet::SelectorIsHadron()
           && fastjet::SelectorIsNeutral())(constituents);
        double neutral_hadrons_pt = join(neutral_hadrons).perp();
        cout << setw(4) << i << fixed << setprecision(3) << setw(11)
             << sortedJets[i].perp() << setw(9)  << sortedJets[i].rap()
             << setw(9) << sortedJets[i].phi_std()
             << setw(6) << constituents.size()
             << setw(8) << fastjet::SelectorIsCharged().count(constituents)
             << setw(8) << fastjet::SelectorId(22).count(constituents)
             << setw(13) << hardest.user_info<Particle>().name()
             << "     " << setw(10) << neutral_hadrons_pt << endl;
      }
      cout << "\n --------  End FastJet Listing  ------------------"
           << "---------------------------------" << endl;
    }

    // Fill distribution of FastJet jets relative to SlowJet ones.
    int nFast = sortedJets.size();
    nJetsF.fill( nFast - nSlow);
    if (nFast == nSlow) {
      for (int i = 0; i < nFast; ++i) {
        pTjetsF.fill( sortedJets[i].perp() - slowJet.pT(i) );
        double dist2 = pow2( sortedJets[i].rap() - slowJet.y(i))
          + pow2( sortedJets[i].phi_std() - slowJet.phi(i));
        RdistF.fill( sqrt(dist2) );
      }
    }

    // Comparison of time consumption by analyzed multiplicity.
    nAna.fill( nAnalyze);
    tGen.fill( nAnalyze, aftGen - befGen);
    tSlow.fill( nAnalyze, aftSlow - befSlow);
    tCore.fill( nAnalyze, aftCore - befCore);
    tFast.fill( nAnalyze, aftFast - befFast);

  // End of event loop.
  }

  // Statistics. Histograms.
  pythia.stat();
  tSlowGen  = tSlow / tGen;
  tCoreGen  = tCore / tGen;
  tFastGen  = tFast / tGen;
  tFastCore = tFast / tCore;
  tGen     /= nAna;
  tSlow    /= nAna;
  tCore    /= nAna;
  tFast    /= nAna;
  cout << nJetsS << nJetsC << nJetsF << pTjetsS << pTjetsC << pTjetsF
       << etaJetsS << phiJetsS << RdistC << RdistF << distJets << pTdiff
       << nAna << tGen << tSlow << tCore << tFast << tSlowGen << tCoreGen
       << tFastGen << tFastCore;

  // Done.
  return 0;
}
