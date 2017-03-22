// main07.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Illustration how to generate various two-body channels from 
// astroparticle processes, e.g. neutralino annihilation or decay.
// To this end a "blob" of energy is created with unit cross section,
// from the fictitious collision of two non-radiating incoming e+e-.
// In the accompanying main29.cmnd file the decay channels of this
// blob can be set up. Furthermore, only gamma, e+-, p/pbar and
// neutrinos are stable, everything else is set to decay.
// (The "single-particle gun" of main21.cc offers another possible 
// approach to the same problem.) 

#include "Pythia.h"

using namespace Pythia8; 
 
//==========================================================================

// A derived class for (e+ e- ->) GenericResonance -> various final states.

class Sigma1GenRes : public Sigma1Process {

public:

  // Constructor.
  Sigma1GenRes() {}

  // Evaluate sigmaHat(sHat): dummy unit cross section. 
  virtual double sigmaHat() {return 1.;}

  // Select flavour. No colour or anticolour.
  virtual void setIdColAcol() {setId( -11, 11, 999999);
    setColAcol( 0, 0, 0, 0, 0, 0);}

  // Info on the subprocess.
  virtual string name()    const {return "GenericResonance";}
  virtual int    code()    const {return 9001;}
  virtual string inFlux()  const {return "ffbarSame";}

};

//==========================================================================

int main() {

  // Pythia generator.
  Pythia pythia;

  // A class to generate the fictitious resonance initial state. 
  SigmaProcess* sigma1GenRes = new Sigma1GenRes();

  // Hand pointer to Pythia.
  pythia.setSigmaPtr( sigma1GenRes);

  // Read in the rest of the settings and data from a separate file.
  pythia.readFile("main07.cmnd");

  // Initialization.
  pythia.init();

  // Extract settings to be used in the main program.
  int nEvent  = pythia.mode("Main:numberOfEvents");
  int nAbort  = pythia.mode("Main:timesAllowErrors");

  // Histogram particle spectra.
  Hist eGamma("energy spectrum of photons",        100, 0., 250.);
  Hist eE(    "energy spectrum of e+ and e-",      100, 0., 250.);
  Hist eP(    "energy spectrum of p and pbar",     100, 0., 250.);
  Hist eNu(   "energy spectrum of neutrinos",      100, 0., 250.);
  Hist eRest( "energy spectrum of rest particles", 100, 0., 250.);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    // Loop over all particles and analyze the final-state ones.
    for (int i = 0; i < pythia.event.size(); ++i) 
    if (pythia.event[i].isFinal()) {
      int idAbs = pythia.event[i].idAbs();
      double eI = pythia.event[i].e();
      if (idAbs == 22) eGamma.fill(eI);
      else if (idAbs == 11) eE.fill(eI);
      else if (idAbs == 2212) eP.fill(eI);
      else if (idAbs == 12 || idAbs == 14 || idAbs == 16) eNu.fill(eI);
      else {
        eRest.fill(eI);
        cout << " Error: stable id = " << pythia.event[i].id() << endl;
      }
    }

  // End of event loop.
  }

  // Final statistics and histograms.
  pythia.stat();
  cout << eGamma << eE << eP << eNu << eRest;

  // Done.
  delete sigma1GenRes;
  return 0;
}
