// SigmaLeptoquark.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for leptoquark-process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma(1/2)Process.
// Note: since leptoquark assumed scalar no need for decay-angles routines.

#ifndef Pythia8_SigmaLeptoquark_H
#define Pythia8_SigmaLeptoquark_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A derived class for q l -> LQ (leptoquark).

class Sigma1ql2LeptoQuark : public Sigma1Process {

public:

  // Constructor.
  Sigma1ql2LeptoQuark() : idQuark(), idLepton(), mRes(), GammaRes(), m2Res(),
    GamMRat(), kCoup(), widthIn(), sigBW(), LQPtr() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return "q l -> LQ (leptoquark)";}
  virtual int    code()       const {return 3201;}
  virtual string inFlux()     const {return "ff";}
  virtual int    resonanceA() const {return 42;}

private:

  // Parameters set at initialization or for current kinematics.
  int    idQuark, idLepton;
  double mRes, GammaRes, m2Res, GamMRat, kCoup, widthIn, sigBW;

  // Pointer to properties of the particle species, to access decay channel.
  ParticleDataEntry* LQPtr;

};

//==========================================================================

// A derived class for q g -> LQ l (leptoquark).

class Sigma2qg2LeptoQuarkl : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2LeptoQuarkl() : idQuark(), idLepton(), mRes(), GammaRes(), m2Res(),
    GamMRat(), kCoup(), openFracPos(), openFracNeg(), sigma0() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "q g -> LQ l (leptoquark)";}
  virtual int    code()    const {return 3202;}
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return 42;}

private:

  // Parameters set at initialization or for current kinematics.
  int    idQuark, idLepton;
  double mRes, GammaRes, m2Res, GamMRat, kCoup, openFracPos, openFracNeg,
         sigma0;

};

//==========================================================================

// A derived class for g g -> LQ LQbar (leptoquark).

class Sigma2gg2LQLQbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2LQLQbar() : mRes(), GammaRes(), m2Res(), GamMRat(), openFrac(),
    sigma() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "g g -> LQ LQbar (leptoquark)";}
  virtual int    code()    const {return 3203;}
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return 42;}
  virtual int    id4Mass() const {return 42;}

private:

  // Parameters set at initialization or for current kinematics.
  double mRes, GammaRes, m2Res, GamMRat, openFrac, sigma;

};

//==========================================================================

// A derived class for q qbar -> LQ LQbar (leptoquark).

class Sigma2qqbar2LQLQbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2LQLQbar() : idQuark(), mRes(), GammaRes(), m2Res(), GamMRat(),
    kCoup(), openFrac(), sigmaDiff(), sigmaSame() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {
    return (abs(id1) == idQuark) ? sigmaSame : sigmaDiff;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "q qbar -> LQ LQbar (leptoquark)";}
  virtual int    code()    const {return 3204;}
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return 42;}
  virtual int    id4Mass() const {return 42;}

private:

  // Parameters set at initialization or for current kinematics.
  int    idQuark;
  double mRes, GammaRes, m2Res, GamMRat, kCoup, openFrac, sigmaDiff,
         sigmaSame;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaLeptoquark_H
