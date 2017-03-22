// SigmaLeftRightSym.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for left-rights-symmetry differential cross sections.
// Contains classes derived from SigmaProcess via Sigma(1/2/3)Process.

#ifndef Pythia8_SigmaLeftRightSym_H
#define Pythia8_SigmaLeftRightSym_H

#include "SigmaProcess.h"

namespace Pythia8 {
 
//==========================================================================

// A derived class for f fbar -> Z_R^0 (righthanded gauge boson).

class Sigma1ffbar2ZRight : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2ZRight() {}

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
  virtual string name()       const {return "f fbar -> Z_R^0";}
  virtual int    code()       const {return 3101;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    resonanceA() const {return idZR;}

private:

  // Parameters set at initialization or for current kinematics. 
  int    idZR;
  double mRes, GammaRes, m2Res, GamMRat, sin2tW, sigma0;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* ZRPtr;

};
 
//==========================================================================

// A derived class for f fbar' -> W_R^+- (righthanded gauge boson).

class Sigma1ffbar2WRight : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2WRight() {}

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
  virtual string name()       const {return "f fbar' -> W_R^+-";}
  virtual int    code()       const {return 3102;}
  virtual string inFlux()     const {return "ffbarChg";}
  virtual int    resonanceA() const {return idWR;}

private:

  // Parameters set at initialization. 
  int    idWR;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0Pos, sigma0Neg;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};
 
//==========================================================================

// A derived class for l l -> H_L^++-- or H_R^++-- (doubly charged Higgs).

class Sigma1ll2Hchgchg : public Sigma1Process {

public:

  // Constructor.
  Sigma1ll2Hchgchg(int leftRightIn ) : leftRight(leftRightIn) {}

  // Initialize process. 
  virtual void initProc(); 

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()       const {return nameSave;}
  virtual int    code()       const {return codeSave;}
  virtual string inFlux()     const {return "ff";}
  virtual int    resonanceA() const {return idHLR;}

private:

  // Parameters set at initialization. 
  int    leftRight, idHLR, codeSave;
  string nameSave;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, yukawa[4][4];

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

}; 
 
//==========================================================================

// A derived class for l- gamma -> H_(L/R)^-- l+  (doubly charged Higgs).

class Sigma2lgm2Hchgchgl : public Sigma2Process {

public:

  // Constructor.
  Sigma2lgm2Hchgchgl(int leftRightIn, int idLepIn ) : leftRight(leftRightIn),
    idLep(idLepIn) {}

  // Initialize process. 
  virtual void initProc(); 

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()       const {return nameSave;}
  virtual int    code()       const {return codeSave;}
  virtual string inFlux()     const {return "fgm";}
  virtual int    resonanceA() const {return idHLR;}

private:

  // Parameters set at initialization. 
  int    leftRight, idHLR, idLep, codeSave;
  string nameSave;
  double yukawa[4], openFracPos, openFracNeg;

}; 
 
//==========================================================================

// A derived class for f_1 f_2 -> H_(L/R)^++-- f_3 f_4 (W+- W+- fusion).

class Sigma3ff2HchgchgfftWW : public Sigma3Process {

public:

  // Constructor.
  Sigma3ff2HchgchgfftWW(int leftRightIn) : leftRight(leftRightIn) {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "ff";}
  virtual int    id3Mass() const {return idHLR;}

  // Instructions for 3-body phase space with t-channel propagators.
  virtual int    idTchan1()        const {return 9900024;}
  virtual int    idTchan2()        const {return 9900024;}
  virtual double tChanFracPow1()   const {return 0.05;}
  virtual double tChanFracPow2()   const {return 0.9;}
  virtual bool   useMirrorWeight() const {return true;}

private:

  // Store standard prefactor.
  int    leftRight, idHLR, codeSave;
  string nameSave;
  double mWS, prefac, sigma0TU, sigma0T, openFracPos, openFracNeg;

};
 
//==========================================================================

// A derived class for f fbar -> H_(L/R)^++ H_(L/R)^--  (doubly charged Higgs).

class Sigma2ffbar2HchgchgHchgchg : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2HchgchgHchgchg(int leftRightIn) : leftRight(leftRightIn) {}

  // Initialize process. 
  virtual void initProc(); 

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()       const {return nameSave;}
  virtual int    code()       const {return codeSave;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    id3Mass()    const {return idHLR;}
  virtual int    id4Mass()    const {return idHLR;}
  virtual int    resonanceA() const {return 23;}

private:

  // Parameters set at initialization. 
  int    leftRight, idHLR, codeSave;
  string nameSave;
  double mRes, GammaRes, m2Res, GamMRat, sin2tW, preFac, yukawa[4][4], 
         openFrac;

}; 

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaLeftRightSym_H
