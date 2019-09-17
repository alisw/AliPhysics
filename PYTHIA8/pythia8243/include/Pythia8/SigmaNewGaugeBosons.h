// SigmaNewGaugeBosons.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for new-gauge-boson-process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma1Process.

#ifndef Pythia8_SigmaNewGaugeBosons_H
#define Pythia8_SigmaNewGaugeBosons_H

#include "Pythia8/PythiaComplex.h"
#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// An intermediate class for f fbar -> Z'/W' -> WW/WZ -> 4 fermions.
// Copied from SigmaEW for gauge-boson-pair production.

class Sigma1ffbarZprimeWprime: public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbarZprimeWprime() {}

protected:

  // Internal products.
  Vec4    pRot[7];
  complex hA[7][7];
  complex hC[7][7];

  // Calculate and store internal products.
  void setupProd( Event& process, int i1, int i2, int i3, int i4,
    int i5, int i6);

  // Evaluate the F function of Gunion and Kunszt.
  complex fGK(int i1, int i2, int i3, int i4, int i5, int i6);

  // Evaluate the Xi function of Gunion and Kunszt.
  double xiGK( double tHnow, double uHnow, double s3now, double s4now);

  // Evaluate the Xj function of Gunion and Kunszt.
  double xjGK( double tHnow, double uHnow, double s3now, double s4now);

private:

};

//==========================================================================

// A derived class for f fbar -> gamma*/Z0/Z'0.

class Sigma1ffbar2gmZZprime : public Sigma1ffbarZprimeWprime {

public:

  // Constructor.
  Sigma1ffbar2gmZZprime() : gmZmode(), maxZpGen(), mRes(), GammaRes(),
    m2Res(), GamMRat(), sin2tW(), cos2tW(), thetaWRat(), mZ(), GammaZ(), m2Z(),
    GamMRatZ(), afZp(), vfZp(), coupZpWW(), anglesZpWW(), gamSum(), gamZSum(),
    ZSum(), gamZpSum(), ZZpSum(), ZpSum(), gamNorm(), gamZNorm(), ZNorm(),
    gamZpNorm(), ZZpNorm(), ZpNorm(), particlePtr() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for Z' decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return "f fbar -> gamma*/Z0/Zprime0";}
  virtual int    code()       const {return 3001;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    resonanceA() const {return 23;}
  virtual int    resonanceB() const {return 32;}

private:

  // Parameters set at initialization or for each new event.
  int    gmZmode, maxZpGen;
  double mRes, GammaRes, m2Res, GamMRat, sin2tW, cos2tW, thetaWRat,
         mZ, GammaZ, m2Z, GamMRatZ, afZp[20], vfZp[20], coupZpWW,
         anglesZpWW, gamSum, gamZSum, ZSum, gamZpSum, ZZpSum, ZpSum,
         gamNorm, gamZNorm, ZNorm, gamZpNorm, ZZpNorm, ZpNorm;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

// A derived class for f fbar' -> W'+-.

class Sigma1ffbar2Wprime : public Sigma1ffbarZprimeWprime {

public:

  // Constructor.
  Sigma1ffbar2Wprime() : mRes(), GammaRes(), m2Res(), GamMRat(), thetaWRat(),
    sigma0Pos(), sigma0Neg(), aqWp(), vqWp(), alWp(), vlWp(), coupWpWZ(),
    anglesWpWZ(), particlePtr() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return "f fbar' -> W'+-";}
  virtual int    code()       const {return 3021;}
  virtual string inFlux()     const {return "ffbarChg";}
  virtual int    resonanceA() const {return 34;}

private:

  // Parameters set at initialization.
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0Pos, sigma0Neg,
         aqWp, vqWp, alWp, vlWp, coupWpWZ, anglesWpWZ;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};
//==========================================================================

// A derived class for f fbar' -> R^0 (horizontal gauge boson).

class Sigma1ffbar2Rhorizontal : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2Rhorizontal() : mRes(), GammaRes(), m2Res(), GamMRat(),
    thetaWRat(), sigma0Pos(), sigma0Neg(), particlePtr() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return "f fbar' -> R^0";}
  virtual int    code()       const {return 3041;}
  virtual string inFlux()     const {return "ffbar";}
  virtual int    resonanceA() const {return 41;}

private:

  // Parameters set at initialization.
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0Pos, sigma0Neg;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia_SigmaNewGaugeBosons_H
