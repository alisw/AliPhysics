// SigmaExtraDim.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Stefan Ask for the *LED* routines.
// Header file for extra-dimensional-process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma(1/2)Process.

#ifndef Pythia8_SigmaExtraDim_H
#define Pythia8_SigmaExtraDim_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A derived class for g g -> G^* (excited graviton state).

class Sigma1gg2GravitonStar : public Sigma1Process {

public:

  // Constructor.
  Sigma1gg2GravitonStar() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for G* decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return "g g -> G*";}
  virtual int    code()       const {return 5001;}
  virtual string inFlux()     const {return "gg";}
  virtual int    resonanceA() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics.
  bool   eDsmbulk, eDvlvl;
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, sigma;

  // Couplings between graviton and SM (indexed by particle id).
  double eDcoupling[27];

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* gStarPtr;

};

//==========================================================================

// A derived class for f fbar -> G^* (excited graviton state).

class Sigma1ffbar2GravitonStar : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2GravitonStar() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for G* decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return "f fbar -> G*";}
  virtual int    code()       const {return 5002;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    resonanceA() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics.
  bool   eDsmbulk, eDvlvl;
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, sigma0;

  // Couplings between graviton and SM (indexed by particle id).
  double eDcoupling[27];

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* gStarPtr;

};

//==========================================================================

// A derived class for q qbar -> g^*/KK-gluon^* (excited kk-gluon state).

class Sigma1qqbar2KKgluonStar : public Sigma1Process {

public:

  // Constructor.
  Sigma1qqbar2KKgluonStar() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for g* decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return "q qbar -> g*/KK-gluon*";}
  virtual int    code()       const {return 5006;}
  virtual string inFlux()     const {return "qqbarSame";}
  virtual int    resonanceA() const {return idKKgluon;}

private:

  // Parameters set at initialization or for current kinematics.
  int    idKKgluon;
  double mRes, GammaRes, m2Res, GamMRat;
  double sumSM, sumInt, sumKK, sigSM, sigInt, sigKK;

  // Couplings between kk gluon and SM (indexed by particle id).
  // Helicity dependent couplings. Use vector/axial-vector
  // couplings internally, gv/ga = 0.5 * (gL +/- gR).
  double eDgv[10], eDga[10];

  // Interference parameter.
  int interfMode;

  // Pointer to properties of the particle species, to access decay
  // channels.
  ParticleDataEntry* gStarPtr;

};

//==========================================================================

// A derived class for g g -> G^* g (excited graviton state).

class Sigma2gg2GravitonStarg : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2GravitonStarg() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight: currently isotropic (except secondary top decay)..
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return "g g -> G* g";}
  virtual int    code()    const {return 5003;}
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics.
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, openFrac, sigma;

};

//==========================================================================

// A derived class for q g -> G^* q (excited graviton state).

class Sigma2qg2GravitonStarq : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2GravitonStarq() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight: currently isotropic (except secondary top decay).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return "q g -> G* q";}
  virtual int    code()    const {return 5004;}
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics.
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, openFrac, sigma;

};

//==========================================================================

// A derived class for q qbar -> G^* g (excited graviton state).

class Sigma2qqbar2GravitonStarg : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2GravitonStarg() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight: currently isotropic (except secondary top decay).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return "q qbar -> G* g";}
  virtual int    code()    const {return 5005;}
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics.
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, openFrac, sigma;

};

//==========================================================================

// NOAM: A derived class for, f fbar -> (gamma/Z)_KKTower -> F Fbar,
// for one heavy F.
// Process provided by N. Hod et al. and is described in arXiv:XXXX.YYYY

class Sigma2ffbar2TEVffbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2TEVffbar(int idIn, int codeIn) : idNew(idIn), codeSave(codeIn) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W decay angles in top decay (else inactive).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return nameSave;}
  virtual int    code()       const {return codeSave;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual bool   isSChannel() const {return true;}
  virtual int    idSChannel() const {return 5000023;}
  virtual int    resonanceA() const {return 23;}
  virtual int    resonanceB() const {return 5000023;}
  virtual int    id3Mass()    const {return idNew;}
  virtual int    id4Mass()    const {return idNew;}
  // Add phase-space sampling also around the Z_KK resonance.
  virtual int    resonanceA();
  virtual int    resonanceB();

private:

  // Values stored for process type.
  string  nameSave;
  int     idNew, gmZmode, codeSave, nexcitationmax;
  bool    isPhysical;
  double  gPlusf, gMinusf, gPlusF, gMinusF, gPlusTop, gMinusTop, gf, gF;
  double  mRes, m2Res, mStar, mTop, m2Top, mZKKn, m2ZKKn, m2gmKKn, mgmKKn,
          alphaemfixed;
  double  helicityME2, coefTot, coefAngular;
  double  mr, betaf, cosThe, openFracPair;
  double  wgmKKFactor, wgmKKn, wZKKn,
          wZ0, ttbarwZKKn, ttbarwgmKKn,
          ttbarwFactorA, ttbarwFactorB;
  double  phaseSpacemHatMin, phaseSpacemHatMax;
  complex gammaProp, resProp, gmPropKK, ZPropKK, totalProp;
  complex mI;
};

//==========================================================================

// A derived class for g g -> U/G g (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2gg2LEDUnparticleg : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2gg2LEDUnparticleg( bool Graviton ) : eDgraviton(Graviton) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return
    (eDgraviton ? "g g -> G g" : "g g -> U g") ;}
  virtual int    code()       const {return (eDgraviton ? 5021 : 5045);}
  virtual string inFlux()     const {return "gg";}
  virtual int    id3Mass()    const {return 5000039;}
  virtual int    id4Mass()    const {return 21;}

private:

  bool   eDgraviton;
  int    eDspin, eDnGrav, eDidG, eDcutoff;
  double mG, mGS, eDsigma0, eDdU, eDLambdaU, eDlambda, eDconstantTerm,
         eDtff, eDcf;

};

//==========================================================================

// A derived class for q g -> U/G q (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2qg2LEDUnparticleq : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2qg2LEDUnparticleq( bool Graviton) : eDgraviton(Graviton) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return
    (eDgraviton ? "q g -> G q" : "q g -> U q") ;}
  virtual int    code()       const {return (eDgraviton ? 5022 : 5046);}
  virtual string inFlux()     const {return "qg";}
  virtual int    id3Mass()    const {return 5000039;}

private:

  bool   eDgraviton;
  int    eDspin, eDnGrav, eDidG, eDcutoff;
  double mG, mGS, eDsigma0, eDdU, eDLambdaU, eDlambda, eDconstantTerm,
         eDtff, eDgf, eDcf;

};

//==========================================================================

// A derived class for q qbar -> U/G g (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2qqbar2LEDUnparticleg : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2qqbar2LEDUnparticleg( bool Graviton) : eDgraviton(Graviton) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return
    (eDgraviton ? "q qbar -> G g" : "q qbar -> U g") ;}
  virtual int    code()       const {return (eDgraviton ? 5023 : 5047);}
  virtual string inFlux()     const {return "qqbarSame";}
  virtual int    id3Mass()    const {return 5000039;}
  virtual int    id4Mass()    const {return 21;}

private:

  bool   eDgraviton;
  int    eDspin, eDnGrav, eDidG, eDcutoff;
  double mG, mGS, eDsigma0, eDdU, eDLambdaU, eDlambda, eDconstantTerm,
         eDtff, eDgf, eDcf;

};

//==========================================================================

// A derived class for f fbar -> U/G Z (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2ffbar2LEDUnparticleZ : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2ffbar2LEDUnparticleZ( bool Graviton) : eDgraviton(Graviton) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return
    (eDgraviton ? "f fbar -> G Z" : "f fbar -> U Z") ;}
  virtual int    code()       const {return (eDgraviton ? 5024 : 5041);}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    id3Mass()    const {return 5000039;}
  virtual int    id4Mass()    const {return 23;}
  virtual int    resonanceA() const {return 23;}
  virtual int    gmZmode()    const {return 2;}

private:

  // Constants: could only be changed in the code itself.
  static const double FIXRATIO;

  int    eDspin, eDnGrav, eDcutoff, eDidG;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDratio, eDlambdaPrime,
         eDtff, eDconstantTerm;
  double sHS, tHS, uHS, tHC, uHC, tHQ, uHQ, tHuH, mU, mUS, mZ, widZ,
         mZS, mwZS, eDsigma0;

};

//==========================================================================

// A derived class for f fbar -> U/G gamma (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2ffbar2LEDUnparticlegamma : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2ffbar2LEDUnparticlegamma( bool Graviton) : eDgraviton(Graviton) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return
    (eDgraviton ? "f fbar -> G gamma" : "f fbar -> U gamma") ;}
  virtual int    code()       const {return (eDgraviton ? 5025 : 5042);}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    id3Mass()    const {return 5000039;}
  virtual int    id4Mass()    const {return 22;}

private:

  // Constants: could only be changed in the code itself.
  static const double FIXRATIO;

  int    eDspin, eDnGrav, eDcutoff, eDidG;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDratio, eDlambdaPrime,
         eDtff, eDconstantTerm;
  double sHS, tHS, uHS, tHC, uHC, tHQ, uHQ, tHuH, mU, mUS, mZ,
         mZS, eDsigma0;

};

//==========================================================================

// A derived class for f fbar -> (LED G*/U*) -> gamma gamma
// (virtual graviton/unparticle exchange).

class Sigma2ffbar2LEDgammagamma : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2ffbar2LEDgammagamma( bool Graviton) : eDgraviton(Graviton) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return
    (eDgraviton ? "f fbar -> (LED G*) -> gamma gamma"
    : "f fbar -> (U*) -> gamma gamma") ;}
  virtual int    code()       const {return (eDgraviton ? 5026 : 5043);}
  virtual string inFlux()     const {return "ffbarSame";}

private:

  int    eDspin, eDcutoff, eDnGrav, eDnegInt;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDlambda2chi,
         eDterm1, eDterm2, eDterm3, eDtff;

};

//==========================================================================

// A derived class for g g -> (LED G*/U*) -> gamma gamma
// (virtual graviton/unparticle exchange).

class Sigma2gg2LEDgammagamma : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2gg2LEDgammagamma( bool Graviton) : eDgraviton(Graviton) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return (eDgraviton
    ? "g g -> (LED G*) -> gamma gamma" : "g g -> (U*) -> gamma gamma") ;}
  virtual int    code()       const {return (eDgraviton ? 5027 : 5044);}
  virtual string inFlux()     const {return "gg";}

private:

  int    eDspin, eDcutoff, eDnGrav;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDlambda2chi, eDsigma0, eDtff;

};

//==========================================================================

// A derived class for f fbar -> (LED G*/U*) -> l lbar
// (virtual graviton/unparticle exchange).
// Does not include t-channel contributions relevant for e^+e^- to e^+e^-

class Sigma2ffbar2LEDllbar : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2ffbar2LEDllbar( bool Graviton) : eDgraviton(Graviton) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return
    (eDgraviton ? "f fbar -> (LED G*) -> l l" : "f fbar -> (U*) -> l l") ;}
  virtual int    code()       const {return (eDgraviton ? 5028 : 5048);}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual bool   isSChannel() const {return true;}

private:

  int    eDspin, eDcutoff, eDnGrav,eDnxx, eDnxy, eDnegInt;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDlambda2chi, eDtff,
         eDmZ, eDmZS, eDGZ, eDGZS, eDabsMeU, eDdenomPropZ, eDrePropGamma,
         eDrePropZ, eDimPropZ, eDabsAS, eDreA, eDreABW, eDpoly1, eDpoly2,
         eDpoly3;

};

//==========================================================================

// A derived class for g g -> (LED G*/U*) -> l lbar
// (virtual graviton/unparticle exchange).

class Sigma2gg2LEDllbar : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2gg2LEDllbar( bool Graviton) : eDgraviton(Graviton) {}


  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat() {return eDsigma0;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return
    (eDgraviton ? "g g -> (LED G*) -> l l" : "g g -> (U*) -> l l") ;}
  virtual int    code()       const {return (eDgraviton ? 5029 : 5049);}
  virtual string inFlux()     const {return "gg";}

private:

  int    eDspin, eDcutoff, eDnGrav;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDlambda2chi, eDsigma0, eDtff;

};

//==========================================================================

// A derived class for g g -> (LED G*) -> g g.

class Sigma2gg2LEDgg : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2LEDgg() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "g g -> (LED G*) -> g g";}
  virtual int    code()   const {return 5030;}
  virtual string inFlux() const {return "gg";}

private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigTU, sigSum, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for g g -> (LED G*) -> q qbar.

class Sigma2gg2LEDqqbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2LEDqqbar() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "g g -> (LED G*) -> q qbar (uds)";}
  virtual int    code()   const {return 5031;}
  virtual string inFlux() const {return "gg";}

private:

  // Number of quarks to be considered in massless approximation.
  int    nQuarkNew;

  // Values stored for colour flow selection.
  int    idNew;
  double mNew, m2New, sigTS, sigUS, sigSum, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for q g -> (LED G*) -> q g.
// Use massless approximation also for Q since no alternative.

class Sigma2qg2LEDqg : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2LEDqg() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "q g -> (LED G*) -> q g";}
  virtual int    code()   const {return 5032;}
  virtual string inFlux() const {return "qg";}

private:

  // Values stored for colour flow selection.
  double sigTS, sigTU, sigSum, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for q q(bar)' -> (LED G*) -> q q(bar)'.

class Sigma2qq2LEDqq : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2LEDqq() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "q q(bar)' -> (LED G*) -> q q(bar)'";}
  virtual int    code()   const {return 5033;}
  virtual string inFlux() const {return "qq";}

 private:

  // Values stored for colour flow selection.
  double sigT, sigU, sigTU, sigST, sigSum;
  double sigGrT1, sigGrT2, sigGrU, sigGrTU, sigGrST;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for q qbar -> (LED G*) -> g g.

class Sigma2qqbar2LEDgg : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2LEDgg() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "q qbar -> (LED G*) -> g g";}
  virtual int    code()   const {return 5034;}
  virtual string inFlux() const {return "qqbarSame";}

 private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigSum, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for q qbar -> (LED G*) -> q' qbar'.

class Sigma2qqbar2LEDqqbarNew : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2LEDqqbarNew() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()  const {return "q qbar -> (LED G*) -> q' qbar' (uds)";}
  virtual int    code()   const {return 5035;}
  virtual string inFlux() const {return "qqbarSame";}

 private:

  // Number of quarks to be considered in massless approximation.
  int    nQuarkNew;

  // Values stored for colour flow selection.
  int    idNew;
  double mNew, m2New, sigS, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaExtraDim_H
