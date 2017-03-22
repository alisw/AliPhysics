// main29.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Omnibus version for colour reconnection (CR) effect studies.

// Links to two UserHooks that, along with the internal models,
// implement all the models used for the top mass study in
// S. Argyropoulos and T. Sjostrand,
// arXiv:1407.6653 [hep-ph] (LU TP 14-23, DESY 14-134, MCnet-14-15)

// Warning: some small modifications have been made when collecting
// the models, but nothing intended to change the behaviour.
// Note: the move model is also available with ColourReconnection:mode = 2,
// while the ColourReconnection:mode = 1 model has not been used here.
// Note: the new models tend to be slower than the default CR scenario,
// since they have to probe many more reconnection possibilities.

// Important: the top mass shift analysis encoded here is very primitive,
// does not perform well at all, and should not be taken seriously.
// The important part is that you see how the different scenarios
// should be set up to operate as intended.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/ColourReconnectionHooks.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events to generate.
  // Warning: much statistics is needed for significant results,
  // so this is just an appetizer. Anyway, the reconstruction is
  // pretty lousy, so not useful for real studies.
  int nEvent = 1000;

  // Target t and W masses.
  double mT = 173.3;
  double mW = 80.385;

  // Set up anti-kT jet finder.
  double Rjet = 0.5;
  double pTjetMin = 20.;
  SlowJet sJet( -1, Rjet, pTjetMin);

  // Loop over different reconnection scenarios.
  for (int mLoop = 0; mLoop < 14; ++mLoop) {
  cout << "\n\n ================================ Now begin mLoop = "
       << mLoop << " ================================\n" << endl;

    // Generator at 8 TeV LHC.
    Pythia pythia;
    Event& event = pythia.event;
    pythia.readString("Beams:eCM = 8000.");

    // q qbar, g g -> t tbar.
    pythia.readString("Top:qqbar2ttbar = on");
    pythia.readString("Top:gg2ttbar = on");

    // Colour reconnection setups.
    UserHooks* myUserHooks;

    // Tuning parameters for CR scenarios have tune 4C as a starting point.
    pythia.readString("Tune:pp = 5");

    // No reconnection at all.
    if (mLoop == 0) {
      pythia.readString("ColourReconnection:reconnect = off");
      pythia.readString("PartonLevel:earlyResDec = off");
      pythia.readString("MultipartonInteractions:pT0Ref = 2.30");

    // Standard reconnection, but top decay products unaffected.
    } else if (mLoop == 1) {
      pythia.readString("ColourReconnection:reconnect = on");
      pythia.readString("PartonLevel:earlyResDec = off");

    // Standard reconnection, including top decay products.
    } else if (mLoop == 2) {
      pythia.readString("ColourReconnection:reconnect = on");
      pythia.readString("PartonLevel:earlyResDec = on");

    // New gluon swap and move scenarios, including top decay products.
    // (Note: the move scenario is also implemented internally.)
    } else if (mLoop >= 3 && mLoop <= 8) {
      pythia.readString("ColourReconnection:reconnect = off");
      pythia.readString("PartonLevel:earlyResDec = off");
      // Swap (mode = 1) or move (2), and flip (1, 2) or not (0).
      int mode = (mLoop <= 5) ? 1 : 2;
      int flip = ( mLoop - 3 * mode) % 3;
      // Possibilities to vary effects by further parameters.
      double dLamCut   = 0.;
      double fracGluon = 1.;
      if ( mode == 1 ) {
        if ( flip > 0 )
             pythia.readString("MultipartonInteractions:pT0Ref = 2.20");
        else pythia.readString("MultipartonInteractions:pT0Ref = 2.30");
      }
      else {
        if ( flip > 0 )
             pythia.readString("MultipartonInteractions:pT0Ref = 2.15");
        else pythia.readString("MultipartonInteractions:pT0Ref = 2.25");
      }
      myUserHooks = new MBReconUserHooks(mode, flip, dLamCut, fracGluon);
      pythia.setUserHooksPtr( myUserHooks);

    // New scenaros that do top reconnections separately from normal one.
    // =  9: reconnect with random background gluon;
    // = 10: reconnect with nearest (smallest-mass) background gluon;
    // = 11: reconnect with furthest (largest-mass) background gluon;
    // = 12: reconnect with smallest (with sign) lambda measure shift;
    // = 13: reconnect only if reduced lamda, and then to most reduction.
    } else if (mLoop >= 9 && mLoop <= 13) {
      pythia.readString("ColourReconnection:reconnect = on");
      pythia.readString("PartonLevel:earlyResDec = off");
      // Possibility with reduced reconnection strength.
      double strength = ( mLoop == 13 ) ? 1. : 0.075;
      myUserHooks = new TopReconUserHooks(mLoop - 8, strength);
      pythia.setUserHooksPtr( myUserHooks);
    }

    // Simplify generation. For tryout only.
    //pythia.readString("ProcessLevel:resonanceDecays = off");
    //pythia.readString("PartonLevel:ISR = off");
    //pythia.readString("PartonLevel:FSR = off");
    //pythia.readString("PartonLevel:MPI = off");
    //pythia.readString("BeamRemnants:primordialKT = off");
    //pythia.readString("HadronLevel:all = off");

    // Top and W masses. Semileptonic top decay chosen by W decay.
    // (One of two charge states, so properly ought to symmetrize.)
    // Trick: only allow decay to stable tau, standing in for e and mu
    // as well, but the tau is easy to remove before jet finding.
    pythia.readString("6:m0 = 173.3");
    pythia.readString("24:m0 = 80.385");
    pythia.readString("24:onPosIfAny = 1 2 3 4 5");
    pythia.readString("24:onNegIfAny = 15");
    pythia.readString("24:offIfAny = 11 13");
    pythia.readString("15:mayDecay = off");

    // Reduce printout.
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberCount = 100000");

    // Initialize.
    pythia.init();

    // Histograms for current scenario.
    Hist nRecH(  "number of top reconnections",  100, -0.5, 99.5);
    Hist nchH(   "charged multiplicity",         100,  -1., 799.);
    Hist nJetH(  "jet multiplicity",              20, -0.5, 19.5);
    Hist mWH(    "reconstructed W mass",         100,  40., 140.);
    Hist mTH(    "reconstructed t mass",         100, 120., 220.);
    Hist mWerrH( "reconstructed W mass error",   100, -10.,  10.);
    Hist mTerrH( "reconstructed t mass error",   100, -20.,  20.);
    Hist pTTH(   "reconstructed pT_t",            11,   0., 275.);
    Hist mTpTH(  "reconstructed delta-m_t(pT_t)", 11,   0., 275.);

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Charged multiplicity.
      int nch = 0;
      for (int i = 0; i < event.size(); ++i)
        if (event[i].isFinal() && event[i].isCharged()) ++nch;
      nchH.fill(nch);

      // Remove tau leptons. (Recall: they were put stable, so simple.)
      for (int i = 0; i < event.size(); ++i)
        if (event[i].idAbs() == 15) event[i].statusNeg();

      // Find number of jets. At least four to keep going.
      sJet.analyze(event);
      int nJet = sJet.sizeJet();
      nJetH.fill( nJet);
      if (nJet < 4) continue;

      // Find two jets that form mass closest to mW.
      int i1min     = 0;
      int i2min     = 0;
      double m12min = 0.;
      double diff   = 1e10;
      for (int i1 = 0; i1 < nJet - 1; ++i1)
      for (int i2 = i1 + 1; i2 < nJet; ++i2) {
        double m12 = (sJet.p(i1) + sJet.p(i2)).mCalc();
        if (abs(m12 - mW) < diff) {
          i1min  = i1;
          i2min  = i2;
          m12min = m12;
          diff   = abs(m12 - mW);
        }
      }
      mWH.fill( m12min);
      mWerrH.fill( m12min - mW);

      // Only keep going if within +-5 GeV.
      if (abs(m12min - mW) > 5.) continue;

      // Find third jet that forms mass closest to mT.
      int i3min      = 0;
      double m123min = 0.;
      diff           = 1e10;
      for (int i3 = 0; i3 < nJet; ++i3)
      if (i3 != i1min && i3 != i2min) {
        double m123 = (sJet.p(i1min) + sJet.p(i2min) + sJet.p(i3)).mCalc();
        if (abs(m123 - mT) < diff) {
          i3min   = i3;
          m123min = m123;
          diff    = abs(m123 - mT);
        }
      }
      mTH.fill( m123min);
      mTerrH.fill( m123min - mT);

      // Only keep going if within +-20 GeV.
      if (abs(m123min - mT) > 20.) continue;

      // Study top pT and dependence of top mass error.
      double pTT = (sJet.p(i1min) + sJet.p(i2min) + sJet.p(i3min)).pT();
      if (pTT > 250.) pTT = 260.;
      pTTH.fill( pTT);
      mTpTH.fill( pTT, m123min - mT);

    // End of event loop. Statistics. Histograms.
    }
    pythia.stat();
    mTpTH /= pTTH;
    cout <<  nchH << nJetH << mWH << mTH << mWerrH
         << mTerrH << pTTH << mTpTH;

    // End loop over top colour reconnection scenarios.
    if (mLoop >=3) delete myUserHooks;
  }

  // Done.
  return 0;
}
