// StandardModel.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file gives access to some Standard Model parameters.
// AlphaStrong: fix or first- or second-order running alpha_strong.

#ifndef Pythia8_StandardModel_H
#define Pythia8_StandardModel_H

#include "Pythia8/Basics.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The AlphaStrong class calculates the alpha_strong value at an arbitrary
// scale, given the value at m_Z, to zeroth, first or second order.

class AlphaStrong {

public:

  // Constructors.
  AlphaStrong() : isInit(false), order(0), nfmax(),
    Lambda3Save(0.), Lambda4Save(0.), Lambda5Save(0.), Lambda6Save(0.),
    Lambda3Save2(0.), Lambda4Save2(0.), Lambda5Save2(0.), Lambda6Save2(0.),
    scale2Min(0.), mc(0.), mb(0.), mt(0.), mc2(0.), mb2(0.), mt2(0.), useCMW(),
    lastCallToFull(false), valueRef(0.), valueNow(0.), scale2Now(0.) {}

  // Destructor.
  virtual ~AlphaStrong() {}

  // Initialization for given value at M_Z and given order.
  virtual void init(double valueIn = 0.12, int orderIn = 1, int nfmaxIn = 6,
    bool useCMWIn = false);

  // Set flavour threshold values: m_c, m_b, m_t.
  virtual void setThresholds(double mcIn, double mbIn, double mtIn) {
    mt=mtIn; mb=min(mt,mbIn); mc=min(mb,mcIn);}

  // alpha_S value and Lambda values.
  double alphaS(double scale2);
  double alphaS1Ord(double scale2);
  double alphaS2OrdCorr(double scale2);
  double Lambda3() const { return Lambda3Save; }
  double Lambda4() const { return Lambda4Save; }
  double Lambda5() const { return Lambda5Save; }
  double Lambda6() const { return (nfmax >= 6) ? Lambda6Save : Lambda5Save; }

  // Info: tell which scales we use for flavour thresholds.
  double muThres(int idQ);
  double muThres2(int idQ);

  // Return the CMW factor (for nF between 3 and 6).
  double facCMW( int nFin);

// Protected data members: accessible to derived classes.
protected:

  // Initialization data member.
  bool   isInit;

  // Running order and max number of flavours to use in running.
  int    order, nfmax;

  // Lambda values.
  double Lambda3Save, Lambda4Save, Lambda5Save, Lambda6Save;
  double Lambda3Save2, Lambda4Save2, Lambda5Save2, Lambda6Save2;

  // Smallest allowed renormalization scale.
  double scale2Min;

  // Flavour thresholds.
  static const double MZ;
  double mc, mb, mt;
  double mc2, mb2, mt2;

  // CMW rescaling factors.
  bool useCMW;
  static const double FACCMW3, FACCMW4, FACCMW5, FACCMW6;

  // Safety margins to avoid getting too close to LambdaQCD.
  static const double SAFETYMARGIN1, SAFETYMARGIN2;

// Private data members: not accessible to derived classes.
private:

  // Private constants: could only be changed in the code itself.
  static const int    NITER;

  // Private data members.
  bool   lastCallToFull;
  double valueRef, valueNow, scale2Now;

};

//==========================================================================

// The AlphaEM class calculates the alpha_electromagnetic value at an
// arbitrary scale, given the value at 0 and m_Z, to zeroth or first order.

class AlphaEM {

public:

  // Constructors.
  AlphaEM() : order(), alpEM0(), alpEMmZ(), mZ2(), bRun(), alpEMstep() {}

  // Initialization for a given order.
  void init(int orderIn, Settings* settingsPtr);

  // alpha_EM value.
  double alphaEM(double scale2);

private:

  // Constants: could only be changed in the code itself.
  static const double MZ, Q2STEP[5], BRUNDEF[5];

  // Data members.
  int    order;
  double alpEM0, alpEMmZ, mZ2, bRun[5], alpEMstep[5];

};

//==========================================================================

// The CoupSM class stores and returns electroweak couplings,
// including Cabibbo-Kobayashi-Maskawa mass mixing matrix elements.

class CoupSM {

public:

  // Constructor.
  CoupSM() : s2tW(), c2tW(), s2tWbar(), GFermi(), vfSave(), lfSave(), rfSave(),
    ef2Save(), vf2Save(), af2Save(), efvfSave(), vf2af2Save(), VCKMsave(),
    V2CKMsave(), V2CKMout(), rndmPtr() {}

  // Initialize, normally from Pythia::init().
  void init(Settings& settings, Rndm* rndmPtrIn);

  // alpha_S value and Lambda values.
  double alphaS(double scale2) {return alphaSlocal.alphaS(scale2);}
  double alphaS1Ord(double scale2) {return alphaSlocal.alphaS1Ord(scale2);}
  double alphaS2OrdCorr(double scale2) {
    return alphaSlocal.alphaS2OrdCorr(scale2);}
  double Lambda3() const {return alphaSlocal.Lambda3();}
  double Lambda4() const {return alphaSlocal.Lambda4();}
  double Lambda5() const {return alphaSlocal.Lambda5();}

  // Return alpha_EM value.
  double alphaEM(double scale2) {return alphaEMlocal.alphaEM(scale2);}

  // Return electroweak mixing angle and Fermi constant.
  double sin2thetaW() {return s2tW;}
  double cos2thetaW() {return c2tW;}
  double sin2thetaWbar() {return s2tWbar;}
  double GF() {return GFermi;}

  // Return electroweak couplings of quarks and leptons.
  double ef(int idAbs) {return efSave[idAbs];}
  double vf(int idAbs) {return vfSave[idAbs];}
  double af(int idAbs) {return afSave[idAbs];}
  double t3f(int idAbs) {return 0.5*afSave[idAbs];}
  double lf(int idAbs) {return lfSave[idAbs];}
  double rf(int idAbs) {return rfSave[idAbs];}

  // Return some squared couplings and other combinations.
  double ef2(int idAbs) {return ef2Save[idAbs];}
  double vf2(int idAbs) {return vf2Save[idAbs];}
  double af2(int idAbs) {return af2Save[idAbs];}
  double efvf(int idAbs) {return efvfSave[idAbs];}
  double vf2af2(int idAbs) {return vf2af2Save[idAbs];}

  // Return CKM value or square:
  // first index 1/2/3/4 = u/c/t/t', second 1/2/3/4 = d/s/b/b'.
  double VCKMgen(int genU, int genD) {return VCKMsave[genU][genD];}
  double V2CKMgen(int genU, int genD) {return V2CKMsave[genU][genD];}

  // Return CKM value or square for incoming flavours (sign irrelevant).
  double VCKMid(int id1, int id2);
  double V2CKMid(int id1, int id2);

  // Return CKM sum of squares for given inflavour, or random outflavour.
  double V2CKMsum(int id) {return V2CKMout[abs(id)];}
  int    V2CKMpick(int id);

protected:

  // Constants: could only be changed in the code itself.
  static const double efSave[20], afSave[20];

  // Couplings and VCKM matrix (index 0 not used).
  double s2tW, c2tW, s2tWbar, GFermi, vfSave[20], lfSave[20], rfSave[20],
         ef2Save[20], vf2Save[20], af2Save[20], efvfSave[20],
         vf2af2Save[20], VCKMsave[5][5], V2CKMsave[5][5], V2CKMout[20];

  // Pointer to the random number generator.
  Rndm*       rndmPtr;

  // An AlphaStrong instance for general use (but not MPI, ISR, FSR).
  AlphaStrong alphaSlocal;

  // An AlphaEM instance for general use (but not MPI, ISR, FSR).
  AlphaEM     alphaEMlocal;

};

//==========================================================================

// Generic couplings class

class Couplings : public CoupSM {

public:

  Couplings(bool isSUSYin = false) : isSUSY(isSUSYin) {}
  bool isSUSY;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_StandardModel_H
