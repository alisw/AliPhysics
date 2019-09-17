// main69.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Main program to generate charged hadron spectra from photon-initiated
// hard processes, by combining sub-runs with direct or resolved photons
// or by generating all with contributions in a single run.

// In case of photon-photon interactions four different contributions are
// present:              ProcessType:
//  - resolved-resolved  1
//  - resolved-direct    2
//  - direct-resolved    3
//  - direct-direct      4
// Events can be generated either with photon beams or with photons emitted
// from lepton beams.

// In case of photon-proton interaction two contributions are present
//  - resolved
//  - direct
// When the photon is from beam A the relevant contributions are
// set with "Photon:ProcessType" values 1 (resolved) and 3 (direct) as the
// convention follows the photon-photon case above.
// Also lepton->photon + proton can be generated.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;

  // Decrease the output.
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Shorthand for some public members of pythia (also static ones).
  Settings& settings = pythia.settings;
  Info& info         = pythia.info;

  // Photon-proton collisions.
  bool photonProton         = false;

  // Generate photon-photon events in leptonic or photon beams.
  bool photonsFromElectrons = false;

  // Each contributions separately or in a one combined run.
  bool automaticMix         = true;

  // Optionally use different PDFs from LHAPDF for hard process.
  // Requires linkin with LHAPDF5.
  // pythia.readString("PDF:useHard = on");
  // pythia.readString("PDF:GammaHardSet = LHAPDF5:SASG.LHgrid/5");

  // Beam parameters.
  pythia.readString("Beams:eCM = 200.");

  // Set up beam particles for electron -> photon + proton.
  if ( photonProton) {
    if ( photonsFromElectrons) {
      pythia.readString("Beams:idA = 11");
      pythia.readString("Beams:idB = 2212");
      pythia.readString("PDF:lepton2gamma = on");

    // Set up beam particles for photon + proton.
    } else {
      pythia.readString("Beams:idA = 22");
      pythia.readString("Beams:idB = 2212");
    }

  // Set up beam particles for photon-photon in e+e-.
  } else if ( photonsFromElectrons) {
    pythia.readString("Beams:idA = -11");
    pythia.readString("Beams:idB =  11");
    pythia.readString("PDF:lepton2gamma = on");

  // Set up beam particles for photon-photon.
  } else {
    pythia.readString("Beams:idA = 22");
    pythia.readString("Beams:idB = 22");
  }

  // Cuts on photon virtuality and invariant mass of gamma-gamma/hadron pair.
  if ( photonsFromElectrons) {
    pythia.readString("Photon:Q2max = 1.0");
    pythia.readString("Photon:Wmin  = 10.0");
  }

  // For photon-proton increase pT0Ref (for better agreement with HERA data).
  // Photon-photon has a new default pT0 parametrization tuned to LEP data.
  if ( photonProton)
    pythia.readString("MultipartonInteractions:pT0Ref = 3.00");

  // Limit partonic pThat.
  settings.parm("PhaseSpace:pTHatMin", 5.0);

  // Reset statistics after each subrun.
  pythia.readString("Stat:reset = on");

  // Parameters for histograms.
  double pTmin = 0.0;
  double pTmax = 40.0;
  int nBinsPT  = 40;

  // Initialize the histograms.
  Hist pTtot("Total charged hadron pT distribution", nBinsPT, pTmin, pTmax);
  Hist pTresres("Resolved-resolved contribution", nBinsPT, pTmin, pTmax);
  Hist pTresdir("Resolved-direct contribution", nBinsPT, pTmin, pTmax);
  Hist pTdirres("Direct-resolved contribution", nBinsPT, pTmin, pTmax);
  Hist pTdirdir("Direct-direct contribution", nBinsPT, pTmin, pTmax);
  Hist pTiRun("Contribution from Run i", nBinsPT, pTmin, pTmax);

  // Initialize hard QCD processes with 0, 1, or 2 initial photons.
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhotonParton:all = on");
  if ( !photonProton) {
    pythia.readString("PhotonCollision:gmgm2qqbar = on");
    pythia.readString("PhotonCollision:gmgm2ccbar = on");
    pythia.readString("PhotonCollision:gmgm2bbbar = on");
  }

  // Number of runs.
  int nRuns = photonProton ? 2 : 4;
  if (automaticMix) nRuns = 1;

  // Number of events per run.
  int nEvent = 10000;

  // Loop over relevant processes.
  for ( int iRun = 1; iRun < nRuns + 1; ++iRun) {

    // Turn of MPIs for processes with unresolved photons.
    if (iRun == 2) pythia.readString("PartonLevel:MPI = off");

    // For photon+proton direct contribution with processType = 3.
    if (photonProton && iRun == 2) iRun = 3;

    // Set the type of gamma-gamma process:
    // 0 = mix of all below,
    // 1 = resolved-resolved,
    // 2 = resolved-direct,
    // 3 = direct-resolved,
    // 4 = direct-direct.
    if (automaticMix) settings.mode("Photon:ProcessType", 0);
    else              settings.mode("Photon:ProcessType", iRun);

    // Initialize the generator.
    pythia.init();

    // Clear the histogram.
    pTiRun.null();

    // Begin event loop. Skip if fails.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate next event.
      if (!pythia.next()) continue;

      // List the first process and event for each run.
      if (iEvent == 0) {
        pythia.process.list();
        pythia.event.list();
      }

      // Possible event weights.
      double weight = info.weight();

      // Loop over event record and find charged final state particles.
      for (int i = 0; i < pythia.event.size(); ++i){
        if ( pythia.event[i].isFinal() && pythia.event[i].isCharged() ) {

          // Store the pT value.
          double pTch = pythia.event[i].pT();

          // Fill the correct histogram depending on the process type.
          if (automaticMix) {
            pTtot.fill(pTch, weight);
            if (info.photonMode() == 1) pTresres.fill(pTch, weight);
            if (info.photonMode() == 2) pTresdir.fill(pTch, weight);
            if (info.photonMode() == 3) pTdirres.fill(pTch, weight);
            if (info.photonMode() == 4) pTdirdir.fill(pTch, weight);
          } else {
            pTiRun.fill(pTch, weight);
          }
        }
      }
    } // End of event loop.

    // Show statistics after each run.
    pythia.stat();

    // Normalize to cross section [mb].
    double sigmaNorm = info.sigmaGen() / info.weightSum();
    double pTBin     = (pTmax - pTmin) / (1. * nBinsPT);

    // For mix of all contributions normalize with total cross section.
    if (automaticMix) {
      pTtot    *= sigmaNorm / pTBin;
      pTresres *= sigmaNorm / pTBin;
      pTresdir *= sigmaNorm / pTBin;
      pTdirres *= sigmaNorm / pTBin;
      pTdirdir *= sigmaNorm / pTBin;

    // For each contribution normalize with cross section for the given run.
    } else {
      pTiRun *= sigmaNorm / pTBin;
      if (iRun == 1) pTresres = pTiRun;
      if (iRun == 2) pTresdir = pTiRun;
      if (iRun == 3) pTdirres = pTiRun;
      if (iRun == 4) pTdirdir = pTiRun;
      pTtot += pTiRun;
    }

  // End of loop over runs.
  }

  // Print histograms.
  cout << pTresres << pTresdir << pTdirres << pTdirdir << pTtot;

  // Done.
  return 0;
}
