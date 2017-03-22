// main22.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Simple illustration how to provide (a) your own resonance-width class,
// and (b) your own cross-section class, with instances handed in to Pythia.
// The hypothetical scenario is that top would have been so long-lived
// that a toponium resonance Theta could form. Then production could
// proceed via q qbar -> gamma*/Z* -> Theta, with decay either to
// a fermion pair or (dominantly) to three gluons.
// The implementation is not physically correct in any number of ways,
// but should exemplify the strategy needed for realistic cases.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// The ResonanceTheta class handles a toponium resonance.

class ResonanceTheta : public ResonanceWidths {

public:

  // Constructor.
  ResonanceTheta(int idResIn) {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double normTheta2qqbar, normTheta2llbar, normTheta2ggg;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  // Superfluous here, so skipped.
  //virtual void calcPreFac(bool = false);

  // Calculate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceTheta::initConstants() {

  // Dummy normalization of couplings to the allowed decay channels.
  normTheta2qqbar = 0.0001;
  normTheta2llbar = 0.0001;
  normTheta2ggg   = 0.001;
}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceTheta::calcWidth(bool) {

  // Expression for Theta -> q qbar (q up to b). Colour factor.
  if (id1Abs < 6) widNow = 3. * normTheta2qqbar * mHat;

  // Expression for Theta -> l lbar (l = e, mu, tau).
  else if (id1Abs == 11  || id1Abs == 13 || id1Abs == 15)
    widNow = normTheta2llbar * mHat;

  // Expression for Theta -> g g g. Colour factor.
  else if (id1Abs == 21) widNow = 8. * normTheta2ggg * mHat;

}

//==========================================================================

// A derived class for q qbar -> Theta (toponium bound state).

class Sigma1qqbar2Theta : public Sigma1Process {

public:

  // Constructor.
  Sigma1qqbar2Theta() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). Assumed flavour-independent so simple.
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return "q qbar -> Theta";}
  virtual int    code()       const {return 621;}
  virtual string inFlux()     const {return "qqbarSame";}
  virtual int    resonanceA() const {return 663;}

private:

  // Store flavour-specific process information and standard prefactor.
  int    idTheta;
  double mRes, GammaRes, m2Res, GamMRat, normTheta2qqbar, sigma;

  // Pointer to properties of Theta, to access decay width.
  ParticleDataEntry* particlePtr;

};

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1qqbar2Theta::initProc() {

  // Store Theta mass and width for propagator.
  idTheta  = 663;
  mRes     = particleDataPtr->m0(idTheta);
  GammaRes = particleDataPtr->mWidth(idTheta);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Same normlization as in ResonanceTheta for coupling strength.
  normTheta2qqbar = 0.0001;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(idTheta);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat); first step when inflavours unknown.

void Sigma1qqbar2Theta::sigmaKin() {

  // Incoming width with colour factor.
  double widthIn  = normTheta2qqbar * mH / 3.;

  // Breit-Wigner, including some (guessed) spin factors.
  double sigBW    = 12. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Outgoing width: only includes channels left open.
  double widthOut = particlePtr->resWidthOpen(663, mH);

  // Total answer.
  sigma = widthIn * sigBW * widthOut;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1qqbar2Theta::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idTheta);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 1, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for Theta -> g g g.

double Sigma1qqbar2Theta::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Should be Theta decay. (This is only option here, so overkill.)
  if (iResEnd != iResBeg || process[iResBeg].idAbs() != idTheta)
    return 1.;

  // Should be decay to three gluons.
  int i1 = process[iResBeg].daughter1();
  int i2 = i1 + 1;
  int i3 = i2 + 1;
  if (i3 != process[iResBeg].daughter2() || process[i1].id() != 21)
    return 1.;

  // Energy fractions x_i = 2 E_i/m_Theta of gluons in Theta rest frame.
  double x1 = 2. * process[i1].p() * process[iResBeg].p()
            / process[iResBeg].m2();
  double x2 = 2. * process[i2].p() * process[iResBeg].p()
            / process[iResBeg].m2();
  double x3 = 2. * process[i3].p() * process[iResBeg].p()
            / process[iResBeg].m2();

  // Matrix-element expression for Theta -> g g g.
  double wtME = pow2( (1. - x1) / (x2 * x3) )
    + pow2( (1. - x2) / (x1 * x3) ) + pow2( (1. - x3) / (x1 * x2) );
  double wtMEmax = 2.;
  return wtME / wtMEmax;

}

//==========================================================================

int main() {

  // Number of events to generate. Max number of errors.
  // Warning: generation of complete events is much slower than if you use
  // PartonLevel:all = off to only get cross sections, so adjust nEvent.
  int nEvent = 1000;
  int nAbort = 5;

  // Pythia generator.
  Pythia pythia;

  // Create the toponium resonance and a few production/decay channels.
  // Warning: many more exist, e.g. weak ones of one top quark.
  // Note: to obtain the correct width for the Breit-Wigner you must
  // include all channels, but you only need leave those on that you
  // want to study.
  pythia.readString("663:new = Theta void 3 0 0 342.0 0.2 300. 400. 0.");
  pythia.readString("663:addChannel = 1 0. 0 1 -1");
  pythia.readString("663:addChannel = 1 0. 0 2 -2");
  pythia.readString("663:addChannel = 1 0. 0 3 -3");
  pythia.readString("663:addChannel = 1 0. 0 4 -4");
  pythia.readString("663:addChannel = 1 0. 0 5 -5");
  pythia.readString("663:addChannel = 1 0. 0 11 -11");
  pythia.readString("663:addChannel = 1 0. 0 13 -13");
  pythia.readString("663:addChannel = 1 0. 0 15 -15");
  pythia.readString("663:addChannel = 1 0. 0 21 21 21");

  // Create instance of a class to calculate the width of Theta to the
  // above channels. Hand in pointer to Pythia.
  // Note: Pythia will automatically delete this pointer,
  // along with all other resonances.
  ResonanceWidths* resonanceTheta = new ResonanceTheta(663);
  pythia.setResonancePtr(resonanceTheta);

  // Create instance of a class to generate the q qbar -> Theta process
  // from an external matrix element. Hand in pointer to Pythia.
  SigmaProcess* sigma1Theta = new Sigma1qqbar2Theta();
  pythia.setSigmaPtr(sigma1Theta);

  // Optionally only compare cross sections.
  //pythia.readString("PartonLevel:all = off");
  pythia.readString("Check:nErrList = 2");

  // Initialization for LHC.
  pythia.init();

  // Book histogram.
  Hist mTheta("Theta mass", 100, 300., 400.);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Fill Theta mass. End of event loop.
    mTheta.fill( pythia.process[5].m() );
  }

  // Final statistics. Print histogram.
  pythia.stat();
  cout << mTheta;

  // Done.
  delete sigma1Theta;
  return 0;
}
