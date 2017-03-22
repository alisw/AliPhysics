// main06.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It studies event properties of LEP1 events.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;

  // Allow no substructure in e+- beams: normal for corrected LEP data.
  pythia.readString("PDF:lepton = off");
  // Process selection.
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  // Switch off all Z0 decays and then switch back on those to quarks.
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");

  // LEP1 initialization at Z0 mass.
  pythia.readString("Beams:idA =  11");
  pythia.readString("Beams:idB = -11");
  double mZ = pythia.particleData.m0(23);
  pythia.settings.parm("Beams:eCM", mZ);
  pythia.init();

  // Histograms.
  Hist nCharge("charged multiplicity", 100, -0.5, 99.5);
  Hist spheri("Sphericity", 100, 0., 1.);
  Hist linea("Linearity", 100, 0., 1.);
  Hist thrust("thrust", 100, 0.5, 1.);
  Hist oblateness("oblateness", 100, 0., 1.);
  Hist sAxis("cos(theta_Sphericity)", 100, -1., 1.);
  Hist lAxis("cos(theta_Linearity)", 100, -1., 1.);
  Hist tAxis("cos(theta_Thrust)", 100, -1., 1.);
  Hist nLund("Lund jet multiplicity", 40, -0.5, 39.5);
  Hist nJade("Jade jet multiplicity", 40, -0.5, 39.5);
  Hist nDurham("Durham jet multiplicity", 40, -0.5, 39.5);
  Hist eDifLund("Lund e_i - e_{i+1}", 100, -5.,45.);
  Hist eDifJade("Jade e_i - e_{i+1}", 100, -5.,45.);
  Hist eDifDurham("Durham e_i - e_{i+1}", 100, -5.,45.);

  // Set up Sphericity, "Linearity", Thrust and cluster jet analyses.
  Sphericity sph;
  Sphericity lin(1.);
  Thrust thr;
  ClusterJet lund("Lund");
  ClusterJet jade("Jade");
  ClusterJet durham("Durham");

  // Begin event loop. Generate event. Skip if error. List first few.
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    if (!pythia.next()) continue;

    // Find and histogram charged multiplicity.
    int nCh = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) ++nCh;
    nCharge.fill( nCh );

    // Find and histogram sphericity.
    if (sph.analyze( pythia.event )) {
      if (iEvent < 3) sph.list();
      spheri.fill( sph.sphericity() );
      sAxis.fill( sph.eventAxis(1).pz() );
      double e1 = sph.eigenValue(1);
      double e2 = sph.eigenValue(2);
      double e3 = sph.eigenValue(3);
      if (e2 > e1 || e3 > e2) cout << "eigenvalues out of order: "
      << e1 << "  " << e2 << "  " << e3 << endl;
    }

    // Find and histogram linearized sphericity.
    if (lin.analyze( pythia.event )) {
      if (iEvent < 3) lin.list();
      linea.fill( lin.sphericity() );
      lAxis.fill( lin.eventAxis(1).pz() );
      double e1 = lin.eigenValue(1);
      double e2 = lin.eigenValue(2);
      double e3 = lin.eigenValue(3);
      if (e2 > e1 || e3 > e2) cout << "eigenvalues out of order: "
      << e1 << "  " << e2 << "  " << e3 << endl;
    }

    // Find and histogram thrust.
    if (thr.analyze( pythia.event )) {
      if (iEvent < 3) thr.list();
      thrust.fill( thr.thrust() );
      oblateness.fill( thr.oblateness() );
      tAxis.fill( thr.eventAxis(1).pz() );
      if ( abs(thr.eventAxis(1).pAbs() - 1.) > 1e-8
        || abs(thr.eventAxis(2).pAbs() - 1.) > 1e-8
        || abs(thr.eventAxis(3).pAbs() - 1.) > 1e-8
        || abs(thr.eventAxis(1) * thr.eventAxis(2)) > 1e-8
        || abs(thr.eventAxis(1) * thr.eventAxis(3)) > 1e-8
        || abs(thr.eventAxis(2) * thr.eventAxis(3)) > 1e-8 ) {
        cout << " suspicious Thrust eigenvectors " << endl;
        thr.list();
      }
    }

    // Find and histogram cluster jets: Lund, Jade and Durham distance.
    if (lund.analyze( pythia.event, 0.01, 0.)) {
      if (iEvent < 3) lund.list();
      nLund.fill( lund.size() );
      for (int j = 0; j < lund.size() - 1; ++j)
        eDifLund.fill( lund.p(j).e() - lund.p(j+1).e() );
    }
    if (jade.analyze( pythia.event, 0.01, 0.)) {
      if (iEvent < 3) jade.list();
      nJade.fill( jade.size() );
      for (int j = 0; j < jade.size() - 1; ++j)
        eDifJade.fill( jade.p(j).e() - jade.p(j+1).e() );
    }
    if (durham.analyze( pythia.event, 0.01, 0.)) {
      if (iEvent < 3) durham.list();
      nDurham.fill( durham.size() );
      for (int j = 0; j < durham.size() - 1; ++j)
        eDifDurham.fill( durham.p(j).e() - durham.p(j+1).e() );
    }

  // End of event loop. Statistics. Output histograms.
  }
  pythia.stat();
  cout << nCharge << spheri << linea << thrust << oblateness
       << sAxis << lAxis << tAxis
       << nLund << nJade << nDurham
       << eDifLund << eDifJade << eDifDurham;

  // Done.
  return 0;
}
