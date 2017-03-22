// SigmaTotal.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the class for cross section parametrizations.
// SigmaTotal: total and partial cross section in hadron-hadron collisions.

#ifndef Pythia8_SigmaTotal_H
#define Pythia8_SigmaTotal_H

#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The SigmaTotal class contains parametrizations of total, elastic and
// diffractive cross sections, and of the respective slope parameter.

class SigmaTotal {

public:

  // Constructor.
  SigmaTotal() : isCalc(false) {};

  // Store pointers and initialize data members.
  void init(Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn );

  // Calculate, or recalculate for new beams or new energy.
  bool calc(int idA, int idB, double eCM);

  // Confirm that initialization worked.
  bool   hasSigmaTot() const {return isCalc;}

  // Read out total and partial cross sections.
  double sigmaTot() const {return sigTot;}
  double sigmaEl()  const {return sigEl;}
  double sigmaXB()  const {return sigXB;}
  double sigmaAX()  const {return sigAX;}
  double sigmaXX()  const {return sigXX;}
  double sigmaAXB() const {return sigAXB;}
  double sigmaND()  const {return sigND;}

  // Calculate cross sections in MBR model.
  bool calcMBRxsecs(int idA, int idB, double eCM);

  // Get maximum of xi,dy distribution in MBR model (for event generation).
  double ddpMax()  const {return ddpmax;}
  double sdpMax()  const {return sdpmax;}
  double dpepMax() const {return dpepmax;}

  // Read out slope b in exp(b*t) dependence.
  double bSlopeEl()          const {return bEl;}
  double bSlopeXB(double sX) const { return 2.*bB + alP2 * log(s/sX) ;}
  double bSlopeAX(double sX) const { return 2.*bA + alP2 * log(s/sX) ;}
  double bSlopeXX(double sX1, double sX2) const {
    return alP2 * log( exp(4.) + s * s0 / (sX1 * sX2) ) ;}

  // Read out parameters of diffractive mass spectra.
  double mMinXB()  const {return mMinXBsave;}
  double mMinAX()  const {return mMinAXsave;}
  double mMinAXB() const {return mMinAXBsave;}
  double cRes()    const {return CRES;}
  double mResXB()  const {return mResXBsave;}
  double mResAX()  const {return mResAXsave;}
  double sProton() const {return SPROTON;}

  // Read out parameters of trial t spectra.
  double bMinSlopeXB() const { return max(2., 2. * bB);}
  double bMinSlopeAX() const { return max(2., 2. * bA);}
  double bMinSlopeXX() const { return alP2 * 4.;}

private:

  // Decide whether default or MBR diffractive cross sections.
  int    PomFlux;

  // Constants: could only be changed in the code itself.
  static const int    IHADATABLE[], IHADBTABLE[], ISDTABLE[], IDDTABLE[];
  static const double MMIN, EPSILON, ETA, X[], Y[], BETA0[], BHAD[],
                      ALPHAPRIME, CONVERTEL, CONVERTSD, CONVERTDD, MMIN0,
                      CRES, MRES0, CSD[10][8], CDD[10][9], SPROTON;

  // Integration of MBR cross sections and form factor approximation.
  static const int    NINTEG, NINTEG2;
  static const double HBARC2, FFA1, FFA2,FFB1, FFB2;

  // Initialization data, normally only set once.
  bool   isCalc, setTotal, zeroAXB, doDampen, setElastic;
  double sigAXB2TeV, sigTotOwn, sigElOwn, sigXBOwn, sigAXOwn, sigXXOwn,
         sigAXBOwn, maxXBOwn, maxAXOwn, maxXXOwn, maxAXBOwn, bSlope, rho,
         lambda, tAbsMin, alphaEM0, sigmaPomP, mPomP, pPomP;

  // Parameters of MBR model.
  double MBReps, MBRalpha, MBRbeta0, MBRsigma0, m2min, dyminSDflux,
         dyminDDflux, dyminCDflux, dyminSD, dyminDD, dyminCD,
         dyminSigSD, dyminSigDD, dyminSigCD, sdpmax, ddpmax, dpepmax;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Store values found by calc.
  double sigTot, sigEl, sigXB, sigAX, sigXX, sigAXB, sigND, bEl, s, bA, bB,
         alP2, s0, mMinXBsave, mMinAXsave, mMinAXBsave, mResXBsave,
         mResAXsave;

  // The error function erf(x) should normally be in your math library,
  // but if not uncomment this simple parametrization by Sergei Winitzki.
  //double erf(double x) { double x2 = x * x; double kx2 = 0.147 * x2;
  //  double tmp = sqrt(1. - exp(-x2 * (4./M_PI + kx2) / (1. + kx2)));
  //  return ((x >= 0.) ? tmp : -tmp); }

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaTotal_H
