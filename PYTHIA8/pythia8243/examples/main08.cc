// main08.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates methods to emphasize generation at high pT.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Different modes are illustrated for setting the pT ranges.
  // 1 : Hardcoded in the main program.
  // 2 : Using the Main:subrun keyword in a separate command file.
  // A third method instead biases selection continuously.
  // 3 : Bias high-pT selection by a pT^4 factor.
  // Matching also to low-pT processes is more complicated.
  // 4 : Matching between low- and high-pT. (No diffraction.)
  // 5: As 4, but bias high-pT selection by a pT^4 factor.
  int mode = 5;

  // Number of events to generate per bin.
  int nEvent = 10000;

  // One does not need complete events to study pThard spectrum only.
  bool completeEvents = false;

  // Optionally minimize output (almost) to final results.
  bool smallOutput = true;

  // Book histograms.
  int nRange = 100;
  double pTrange = (mode < 4) ? 1000. : 100.;
  Hist pTraw("pTHat distribution, unweighted", nRange, 0., pTrange);
  Hist pTnorm("pTHat distribution, weighted", nRange, 0., pTrange);
  Hist pTpow3("pTHat distribution, pT3*weighted", nRange, 0., pTrange);
  Hist pTpow5("pTHat distribution, pT5*weighted", nRange, 0., pTrange);
  Hist pTnormPart("pTHat distribution, weighted", nRange, 0., pTrange);
  Hist pTpow3Part("pTHat distribution, pT3*weighted", nRange, 0., pTrange);
  Hist pTpow5Part("pTHat distribution, pT5*weighted", nRange, 0., pTrange);

  // Generator.
  Pythia pythia;

  // Shorthand for some public members of pythia (also static ones).
  Settings& settings = pythia.settings;
  Info& info = pythia.info;

  // Optionally limit output to minimal one.
  if (smallOutput) {
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Next:numberCount = 1000000000");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
  }

  // Number of bins to use. In mode 2 read from main08.cmnd file.
  int nBin = 5;
  if (mode == 2) {
    pythia.readFile("main08.cmnd");
    nBin = pythia.mode("Main:numberOfSubruns");
  }
  else if (mode == 3) nBin = 1;
  else if (mode == 4) nBin = 4;
  else if (mode == 5) nBin = 2;

  // Mode 1: set up five pT bins - last one open-ended.
  double pTlimit[6] = {100., 150., 250., 400., 600., 0.};

  // Modes 4 & 5: set up pT bins for range [0, 100]. The lowest bin
  // is generated with soft processes, to regularize pT -> 0 blowup.
  // Warning: if pTlimitLow[1] is picked too low there will be a
  // visible discontinuity, since soft processes are generated with
  // dampening and "Sudakov" for pT -> 0, while hard processes are not.
  double pTlimitLow[6] = {0., 20., 40., 70., 100.};
  double pTlimitTwo[3] = {0., 20., 100.};

  // Loop over number of bins, i.e. number of subruns.
  for (int iBin = 0; iBin < nBin; ++iBin) {

    // Normally HardQCD, but in two cases nonDiffractive.
    // Need MPI on in nonDiffractive to get first interaction, but not else.
    if (mode > 3 && iBin == 0) {
      pythia.readString("HardQCD:all = off");
      pythia.readString("SoftQCD:nonDiffractive = on");
      if (!completeEvents) {
      pythia.readString("PartonLevel:all = on");
        pythia.readString("PartonLevel:ISR = off");
        pythia.readString("PartonLevel:FSR = off");
        pythia.readString("HadronLevel:all = off");
      }
    } else {
      pythia.readString("HardQCD:all = on");
      pythia.readString("SoftQCD:nonDiffractive = off");
      if (!completeEvents) pythia.readString("PartonLevel:all = off");
    }

    // Mode 1: hardcoded here. Use settings.parm for non-string input.
    if (mode == 1) {
      settings.parm("PhaseSpace:pTHatMin", pTlimit[iBin]);
      settings.parm("PhaseSpace:pTHatMax", pTlimit[iBin + 1]);
    }

    // Mode 2: subruns stored in the main08.cmnd file.
    else if (mode == 2) pythia.readFile("main08.cmnd", iBin);

    // Mode 3: The whole range in one step, but pT-weighted.
    else if (mode == 3) {
      settings.parm("PhaseSpace:pTHatMin", pTlimit[0]);
      settings.parm("PhaseSpace:pTHatMax", 0.);
      pythia.readString("PhaseSpace:bias2Selection = on");
      pythia.readString("PhaseSpace:bias2SelectionPow = 4.");
      pythia.readString("PhaseSpace:bias2SelectionRef = 100.");
    }

    // Mode 4: hardcoded here. Use settings.parm for non-string input.
    else if (mode == 4) {
      settings.parm("PhaseSpace:pTHatMin", pTlimitLow[iBin]);
      settings.parm("PhaseSpace:pTHatMax", pTlimitLow[iBin + 1]);
    }

    // Mode 5: hardcoded here. Use settings.parm for non-string input.
    // Hard processes in one step, but pT-weighted.
    else if (mode == 5) {
      settings.parm("PhaseSpace:pTHatMin", pTlimitTwo[iBin]);
      settings.parm("PhaseSpace:pTHatMax", pTlimitTwo[iBin + 1]);
      if (iBin == 1) {
        pythia.readString("PhaseSpace:bias2Selection = on");
        pythia.readString("PhaseSpace:bias2SelectionPow = 4.");
        pythia.readString("PhaseSpace:bias2SelectionRef = 20.");
      }
    }

    // Initialize for LHC at 14 TeV.
    pythia.readString("Beams:eCM = 14000.");
    pythia.init();

    // Reset local histograms (that need to be rescaled before added).
    pTnormPart.null();
    pTpow3Part.null();
    pTpow5Part.null();

    // Begin event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events. Skip if failure.
      if (!pythia.next()) continue;

      // Soft events have no upper pT limit. They therefore overlap
      // with hard events, and the overlap must be removed by hand.
      // No overlap for elastic/diffraction, which is only part of soft.
      double pTHat  = info.pTHat();
      if (mode > 3 && iBin == 0 && info.isNonDiffractive()
        && pTHat > pTlimitLow[1]) continue;

      // Fill hard scale of event.
      double weight = info.weight();
      pTraw.fill( pTHat );
      pTnormPart.fill( pTHat, weight);
      pTpow3Part.fill( pTHat, weight * pow3(pTHat) );
      pTpow5Part.fill( pTHat, weight * pow5(pTHat) );

    // End of event loop. Statistics.
    }
    if (!smallOutput) pythia.stat();

    // Normalize to cross section for each case, and add to sum.
    double sigmaNorm = (info.sigmaGen() / info.weightSum())
                     * (nRange / pTrange);
    pTnorm += sigmaNorm * pTnormPart;
    pTpow3 += sigmaNorm * pTpow3Part;
    pTpow5 += sigmaNorm * pTpow5Part;

  // End of pT-bin loop.
  }

  // Output histograms.
  cout << pTraw << pTnorm << pTpow3 << pTpow5;

  // Done.
  return 0;
}
