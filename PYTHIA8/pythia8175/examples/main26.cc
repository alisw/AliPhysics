// main26.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a test program for the extra dimensions processes.
// Author: Stefan Ask (Stefan DOT Ask AT cern DOT ch)
// Documentation: S. Ask et al., arXiv:0809.4750 and arXiv:0912.4233

#include "Pythia.h"

using namespace Pythia8; 

// The main program.
int main() {

  // Test cases
  // 1  = Jet + G       (real G emission) 
  // 2  = Jet + U       (real U emission) 
  // 3  = Z + G         (real G emission) 
  // 4  = Z + U         (real U emission)
  // 5  = gamma gamma   (LED G* exchange)
  // 6  = l lbar        (LED U* exchange). 
  //      Note: charged leptons only!
  // 7  = Z_KK/gamma_KK (TEV ED resonance) 
  // 8  = G*            (RS resonance, SM on the TeV brane)
  // 9  = kk-gluon*     (RS resonance)
  int nTest = 1;  

  // Number of events to generate. Max number of errors.
  int nEvent     = 1000;      
  int nAbort     = 50;         

  // Pythia generator.
  Pythia pythia;

  // PYTHIA paramters:
  pythia.readString("PhaseSpace:showViolation = off");

  // Test case parameters
  if (nTest == 1) { 
    pythia.readString("ExtraDimensionsLED:monojet = on");
    pythia.readString("ExtraDimensionsLED:n = 4");
    pythia.readString("ExtraDimensionsLED:MD = 4000.");
    pythia.readString("ExtraDimensionsLED:CutOffmode = 3");
    pythia.readString("ExtraDimensionsLED:t = 2");
    pythia.readString("5000039:m0 = 2500.");
    pythia.readString("5000039:mWidth = 1500.");
    pythia.readString("5000039:mMin = 1.");
    pythia.readString("5000039:mMax = 13990.");
    pythia.readString("PhaseSpace:pTHatMin = 700.");
  } else if (nTest == 2){ 
    pythia.readString("ExtraDimensionsUnpart:gg2Ug = off");
    pythia.readString("ExtraDimensionsUnpart:qg2Uq = on");
    pythia.readString("ExtraDimensionsUnpart:qqbar2Ug = on");
    pythia.readString("ExtraDimensionsUnpart:spinU = 1");
    pythia.readString("ExtraDimensionsUnpart:dU = 1.2");
    pythia.readString("ExtraDimensionsUnpart:LambdaU = 1000");
    pythia.readString("ExtraDimensionsUnpart:lambda = 1.0");
    pythia.readString("ExtraDimensionsUnpart:CutOffmode = 0");
    pythia.readString("5000039:m0 = 300.");
    pythia.readString("5000039:mWidth = 500.");
    pythia.readString("5000039:mMin = 1.");
    pythia.readString("5000039:mMax = 13990.");
    pythia.readString("PhaseSpace:pTHatMin = 700.");
  } else if (nTest == 3){
    pythia.readString("ExtraDimensionsLED:ffbar2GZ = on");
    pythia.readString("ExtraDimensionsLED:n = 6");
    pythia.readString("ExtraDimensionsLED:MD = 2000.");
    pythia.readString("ExtraDimensionsLED:CutOffmode = 1");
    pythia.readString("5000039:m0 = 3000.");
    pythia.readString("5000039:mWidth = 1500.");
    pythia.readString("5000039:mMin = 1.");
    pythia.readString("5000039:mMax = 13990.");
    pythia.readString("PhaseSpace:pTHatMin = 50.");
  } else if (nTest == 4){ 
    pythia.readString("ExtraDimensionsUnpart:ffbar2UZ = on");
    pythia.readString("ExtraDimensionsUnpart:spinU = 1");
    pythia.readString("ExtraDimensionsUnpart:dU = 2.0");
    pythia.readString("ExtraDimensionsUnpart:LambdaU = 1000");
    pythia.readString("ExtraDimensionsUnpart:lambda = 1.000");
    pythia.readString("ExtraDimensionsUnpart:CutOffmode = 0");
    pythia.readString("5000039:m0 = 500.");
    pythia.readString("5000039:mWidth = 1000.");
    pythia.readString("5000039:mMin = 1.");
    pythia.readString("5000039:mMax = 13990.");
    pythia.readString("PhaseSpace:pTHatMin = 50.");
  } else if (nTest == 5){ 
    pythia.readString("ExtraDimensionsLED:ffbar2gammagamma = on");
    pythia.readString("ExtraDimensionsLED:gg2gammagamma = on");
    pythia.readString("ExtraDimensionsLED:LambdaT = 3300.");
    pythia.readString("PhaseSpace:mHatMin = 800.");
  } else if (nTest == 6){ 
    pythia.readString("ExtraDimensionsUnpart:ffbar2llbar = on");
    pythia.readString("ExtraDimensionsUnpart:gg2llbar = off");
    pythia.readString("ExtraDimensionsUnpart:spinU = 1");
    pythia.readString("ExtraDimensionsUnpart:dU = 1.3");
    pythia.readString("ExtraDimensionsUnpart:LambdaU = 1000");
    pythia.readString("ExtraDimensionsUnpart:lambda = 1.0");
    pythia.readString("ExtraDimensionsUnpart:gXX = 0");
    pythia.readString("ExtraDimensionsUnpart:gXY = 0");
    pythia.readString("PhaseSpace:mHatMin = 300.");
  } else if (nTest == 7){
    pythia.readString("ExtraDimensionsTEV:ffbar2mu+mu- = on");
    pythia.readString("ExtraDimensionsTEV:gmZmode = 3"); 
    pythia.readString("ExtraDimensionsTEV:nMax = 100"); 
    pythia.readString("ExtraDimensionsTEV:mStar = 4000");
    pythia.readString("PhaseSpace:mHatMin = 1000");
    pythia.readString("PhaseSpace:mHatMax = 6000");
    pythia.readString("5000023:isResonance = false");
  } else if (nTest == 8){
    pythia.readString("ExtraDimensionsG*:all = on");
  } else if (nTest == 9){
    pythia.readString("ExtraDimensionsG*:qqbar2KKgluon* = on");
    pythia.readString("ExtraDimensionsG*:KKintMode = 2");
    pythia.readString("ExtraDimensionsG*:KKgqR = -0.2");
    pythia.readString("ExtraDimensionsG*:KKgqL = -0.2");
    pythia.readString("ExtraDimensionsG*:KKgbR = -0.2");
    pythia.readString("ExtraDimensionsG*:KKgbL = 1.0");
    pythia.readString("ExtraDimensionsG*:KKgtR = 5.0");
    pythia.readString("ExtraDimensionsG*:KKgtL = 1.0");
    pythia.readString("5100021:m0 = 2000");
  }

  // Switch off sophisticated tau treatment: not yet matched to SUSY.
  pythia.readString("ParticleDecays:sophisticatedTau = 0");

  // Initialization for LHC.
  pythia.readString("Beams:eCM = 14000.");      
  pythia.init();

  // Validation histograms
  Hist hEtjet("dN/dETjet: monojet check", 100, 0., 7000.);
  Hist hMass("dN/m: graviton mass spectrum", 100, 0., 7000.);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    
    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      std::cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    // Checked particle index
    int tmp_monojet = -1;

    // Particle loop
    for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {

      // From hard process (inital = 21, interm.=22, final=23 state)
      if (pythia.event[iPart].statusAbs()  == 22) {

	// Find Z_KK/gamma_KK or kk-gluon
	if( pythia.event[iPart].idAbs() == 5000023 
	 || pythia.event[iPart].idAbs() == 5100021
	 || pythia.event[iPart].idAbs() == 5100039){
	  hMass.fill( pythia.event[iPart].m() );
	}
	
      }	else if ( pythia.event[iPart].statusAbs()  == 23 ) {
	
	// Find graviton/unparticle
	if( pythia.event[iPart].idAbs() == 5000039){
	  hMass.fill( pythia.event[iPart].m() );
	}

	// Find mono-jets
	if (nTest == 1 || nTest == 2) {
	  if ( pythia.event[iPart].idAbs() <= 6 
            || pythia.event[iPart].idAbs() == 21 ){
	    if (tmp_monojet >= 0) {
	      std::cout << "More than one (hard process) mono-jet ! \n";
	    } else {
	      tmp_monojet  = iPart;
	    }
	  }
	}

      }
    }

    // Validation mono-jet wrt G.Giudice et al. paper [hep-ph/9811291v2]
    if (tmp_monojet >= 0) {
      double tmp_eta = pythia.event[tmp_monojet].eta();
      double tmp_et = pythia.event[tmp_monojet].eT();
      double tmp_et_cut = 1000;
      if ( tmp_et >=  tmp_et_cut && abs(tmp_eta) < 3 ) {
	hEtjet.fill( fabs(tmp_et) );
      }    
    }
    
  }  
 
  // Final statistics.
  pythia.stat(); 
  cout << hMass << hEtjet;

  return 0;
}
