// main70.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Main program to demonstrate how to define a photon flux and use that
// to generate charged-particle pT spectra in photo-production processes.
// Author: Ilkka Helenius, September 2017.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

// Photon flux from leptons, corresponds to internal Lepton2gamma.

class Lepton2gamma2 : public PDF {

public:

  // Constructor.
  Lepton2gamma2(int idBeamIn) : PDF(idBeamIn) {}

  // Update the photon flux.
  void xfUpdate(int , double x, double Q2) {
    xgamma = 0.5 * 0.007297353080 / M_PI * (1. + pow2(1. - x)) / Q2;
  }
};

// Photon flux from lead-ions. Integrated over impact parameters > 2*r_Pb.
// Suitable for photo-nuclear processes but not for photon-photon.
// This should be considered as an experimental setup and used with caution.

class Nucleus2gamma : public PDF {

public:

  // Constructor.
  Nucleus2gamma(int idBeamIn) : PDF(idBeamIn) {}

  // Update the photon flux.
  void xfUpdate(int , double x, double ) {

    // Minimum impact parameter (~2*radius) [fm].
    double bmin = 2 * 6.636;

    // Charge of the nucleus.
    double z = 82.;

    // Per-nucleon mass for lead.
    double m2 = pow2(0.9314);
    double alphaEM = 0.007297353080;
    double hbarc = 0.197;
    double xi = x * sqrt(m2) * bmin / hbarc;
    double bK0 = besselK0(xi);
    double bK1 = besselK1(xi);
    double intB = xi * bK1 * bK0 - 0.5 * pow2(xi) * ( pow2(bK1) - pow2(bK0) );
    xgamma = 2. * alphaEM * pow2(z) / M_PI * intB;
  }

};

// The main program.

int main() {

  // Generator.
  Pythia pythia;

  // Decrease the output.
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount = 1000");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 1");
  pythia.readString("Next:numberShowEvent = 1");

  // Two possible process to consider here 1=ep at HERA, 2=UPC at LHC.
  int process = 1;

  // Pointer to externally defined photon flux.
  PDF* photonFlux = 0;

  // Beam parameters.
  pythia.readString("Beams:idA = -11");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("PDF:lepton2gamma = on");
  pythia.readString("PDF:lepton2gammaSet = 2");

  // Photoproduction at HERA.
  // NOTE: Same results more effectively could be obtained with
  // pythia.readString("PDF:lepton2gammaSet = 1"), this is for demonstration
  // purposed only.
  if (process == 1) {
    pythia.readString("Beams:frameType = 2");
    pythia.readString("Beams:eA = 27.5");
    pythia.readString("Beams:eB = 820.");
    photonFlux = new Lepton2gamma2(-11);

  // Experimental UPC generation in PbPb at LHC.
  // Photon flux only from leptons but here need just the photon flux.
  // Since the sampling is optimized for leptons this is not very efficient.
  } else if (process == 2) {
    pythia.readString("Beams:eCM = 5020.");
    // Use nuclear PDF for the hard process generation in the proton side.
    pythia.readString("PDF:useHardNPDFB = on");
    // Modify the minimum impact parameter to match the flux defined above.
    pythia.readString("PDF:gammaFluxApprox2bMin = 13.272");
    // Optimized sampling for photon flux from nuclei.
    pythia.readString("PDF:lepton2gammaApprox = 2");
    // Do not sample virtuality since use b-integrated flux here.
    pythia.readString("Photon:sampleQ2 = off");
    photonFlux = new Nucleus2gamma(-11);
  }

  // Set the external photon flux for beam A.
  pythia.setPhotonFluxPtr(photonFlux, 0);

  // Switch relevant processes on.
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhotonParton:all = on");

  // Limit partonic pThat.
  pythia.readString("PhaseSpace:pTHatMin = 10.0");

  // Use optimized pT0ref for photon-hadron.
  pythia.readString("MultipartonInteractions:pT0Ref = 3.0");

  // Parameters for histograms.
  double pTmin = 0.0;
  double pTmax = 40.0;
  int nBinsPT  = 40;

  // Initialize the histograms.
  Hist pTch("Charged hadron pT distribution", nBinsPT, pTmin, pTmax);

  // Initialize the generator.
  pythia.init();

  // Number of events.
  int nEvent = 10000;

  // Begin event loop. Skip if fails.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate next event.
    if (!pythia.next()) continue;

    // Event weight.
    double weight = pythia.info.weight();

    // Loop over event record.
    for (int i = 0; i < pythia.event.size(); ++i){
      if ( pythia.event[i].isFinal() && pythia.event[i].isCharged() ) {

        // Store the pT value.
        pTch.fill(pythia.event[i].pT(), weight );
      }
    }

  } // End of event loop.

  // Delete photon flux pointer.
  delete photonFlux;

  // Show statistics.
  pythia.stat();

  // Normalize to cross section [mb].
  double sigmaNorm = pythia.info.sigmaGen() / pythia.info.weightSum();
  double pTbin     = (pTmax - pTmin) / (1. * nBinsPT);

  pTch *= sigmaNorm / pTbin;

  cout << pTch;

  // Done.
  return 0;
}
