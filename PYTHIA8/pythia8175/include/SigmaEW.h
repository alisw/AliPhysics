// SigmaEW.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for electroweak process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma(1/2)Process.

#ifndef Pythia8_SigmaEW_H
#define Pythia8_SigmaEW_H

#include "PythiaComplex.h"
#include "SigmaProcess.h"

namespace Pythia8 {

 
//==========================================================================

// A derived class for q g -> q gamma (q = u, d, s, c, b).
// Use massless approximation also for Q since no alternative.

class Sigma2qg2qgamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2qgamma() {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "q g -> q gamma (udscb)";}
  virtual int    code()   const {return 201;}
  virtual string inFlux() const {return "qg";}

private:

  // Values stored for later use.
  double mNew, m2New, sigUS, sigma0;

};
 
//==========================================================================

// A derived class for q qbar -> g gamma.

class Sigma2qqbar2ggamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2ggamma() {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "q qbar -> g gamma";}
  virtual int    code()   const {return 202;}
  virtual string inFlux() const {return "qqbarSame";}

private:

  // Values stored for later use.
  double sigma0;

};
 
//==========================================================================

// A derived class for g g -> g gamma.

class Sigma2gg2ggamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2ggamma() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "g g -> g gamma";}
  virtual int    code()   const {return 203;}
  virtual string inFlux() const {return "gg";}

private:

  // Values stored for later use.  
  double chargeSum, sigma;

};
 
//==========================================================================

// A derived class for f fbar -> gamma gamma.

class Sigma2ffbar2gammagamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2gammagamma() {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "f fbar -> gamma gamma";}
  virtual int    code()   const {return 204;}
  virtual string inFlux() const {return "ffbarSame";}

private:

  // Values stored for later use.
  double sigTU, sigma0;

};
 
//==========================================================================

// A derived class for g g -> gamma gamma.

class Sigma2gg2gammagamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2gammagamma() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "g g -> gamma gamma";}
  virtual int    code()   const {return 205;}
  virtual string inFlux() const {return "gg";}

private:
  
  double charge2Sum, sigma;

};
 
//==========================================================================

// A derived class for f f' -> f f' via t-channel gamma*/Z0 exchange.

class Sigma2ff2fftgmZ : public Sigma2Process {

public:

  // Constructor.
  Sigma2ff2fftgmZ() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "f f' -> f f' (t-channel gamma*/Z0)";}
  virtual int    code()   const {return 211;}
  virtual string inFlux() const {return "ff";}

private:

  //  Z parameters for propagator.
  int    gmZmode;
  double mZ, mZS, thetaWRat, sigmagmgm, sigmagmZ, sigmaZZ;

};
 
//==========================================================================

// A derived class for f_1 f_2 -> f_3 f_4 via t-channel W+- exchange.

class Sigma2ff2fftW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ff2fftW() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()   const {return "f_1 f_2 -> f_3 f_4 (t-channel W+-)";}
  virtual int    code()   const {return 212;}
  virtual string inFlux() const {return "ff";}

private:

  //  W parameters for propagator.
  double mW, mWS, thetaWRat, sigma0;

};
 
//==========================================================================

// A derived class for q q' -> Q q" via t-channel W+- exchange.
// Related to Sigma2ff2fftW class, but with massive matrix elements.

class Sigma2qq2QqtW : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2QqtW(int idIn, int codeIn) : idNew(idIn), codeSave(codeIn) {}

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
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "ff";}
  virtual int    id3Mass() const {return idNew;}

private:

  // Values stored for process type. W parameters for propagator.
  int    idNew, codeSave;
  string nameSave;
  double mW, mWS, thetaWRat, sigma0, openFracPos, openFracNeg;

};
 
//==========================================================================

// A derived class for f fbar -> gamma*/Z0.

class Sigma1ffbar2gmZ : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2gmZ() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for Z decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()       const {return "f fbar -> gamma*/Z0";}
  virtual int    code()       const {return 221;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    resonanceA() const {return 23;}

private:

  // Parameters set at initialization or for each new event. 
  int    gmZmode;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, 
         gamSum, intSum, resSum, gamProp, intProp, resProp;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};
 
//==========================================================================

// A derived class for f fbar' -> W+-.

class Sigma1ffbar2W : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2W() {}

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
  virtual string name()       const {return "f fbar' -> W+-";}
  virtual int    code()       const {return 222;}
  virtual string inFlux()     const {return "ffbarChg";}
  virtual int    resonanceA() const {return 24;}

private:

  // Parameters set at initialization. 
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0Pos, sigma0Neg;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

// A derived class for f fbar -> gamma* -> f' fbar', summed over light f'.
// Allows pT-ordered evolution for multiparton interactions.

class Sigma2ffbar2ffbarsgm : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2ffbarsgm() {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {
    return "f fbar -> f' fbar' (s-channel gamma*)";}
  virtual int    code()       const {return 223;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual bool   isSChannel() const {return true;}

private:

  // Values stored for later use.
  int    idNew;
  double sigma0;

};

//==========================================================================

// A derived class for f fbar -> gamma*/Z0 -> F Fbar, for one heavy F.
// Allows pT cuts as for other 2 -> 2 processes.

class Sigma2ffbar2FFbarsgmZ : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2FFbarsgmZ(int idIn, int codeIn) : idNew(idIn), 
    codeSave(codeIn) {}

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
  virtual int    id3Mass()    const {return idNew;}
  virtual int    id4Mass()    const {return idNew;}
  virtual int    resonanceA() const {return 23;}

private:

  // Values stored for process type. Z parameters for propagator.
  int    idNew, codeSave, gmZmode;
  string nameSave;
  bool   isPhysical;
  double ef, vf, af, mRes, GammaRes, m2Res, GamMRat, thetaWRat,
         mr, betaf, cosThe, gamProp, intProp, resProp, openFracPair; 

};

//==========================================================================

// A derived class for f fbar' -> W+- -> F fbar", for one or two heavy F.
// Allows pT cuts as for other 2 -> 2 processes.

class Sigma2ffbar2FfbarsW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2FfbarsW(int idIn, int idIn2, int codeIn) : idNew(idIn), 
    idNew2(idIn2), codeSave(codeIn) {}

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
  virtual string inFlux()     const {return "ffbarChg";}
  virtual bool   isSChannel() const {return true;}
  virtual int    id3Mass()    const {return idNew;}
  virtual int    id4Mass()    const {return idPartner;}
  virtual int    resonanceA() const {return 24;}

private:

  // Values stored for process type. W parameters for propagator.
  int    idNew, idNew2, codeSave, idPartner;
  string nameSave;
  bool   isPhysical;
  double V2New, mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0,
         openFracPos, openFracNeg; 

};
 
//==========================================================================

// An intermediate class for f fbar -> gamma*/Z0/W+- gamma*/Z0/W-+.

class Sigma2ffbargmZWgmZW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbargmZWgmZW() {}

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
  double xiGK( double tHnow, double uHnow);

  // Evaluate the Xj function of Gunion and Kunszt.
  double xjGK( double tHnow, double uHnow);

private:

};
 
//==========================================================================

// A derived class for f fbar -> gamma*/Z0 gamma*/Z0.

class Sigma2ffbar2gmZgmZ : public Sigma2ffbargmZWgmZW {

public:

  // Constructor.
  Sigma2ffbar2gmZgmZ() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for simultaneous flavour choices.
  virtual double weightDecayFlav( Event& process); 

  // Evaluate weight for decay angles of the two gamma*/Z0.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "f fbar -> gamma*/Z0 gamma*/Z0";}
  virtual int    code()    const {return 231;}
  virtual string inFlux()  const {return "ffbarSame";}
  virtual int    id3Mass() const {return 23;}
  virtual int    id4Mass() const {return 23;}

private:

  // Parameters set at initialization or for each new event. 
  int    gmZmode, i1, i2, i3, i4, i5, i6;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0,  
         gamSum3, intSum3, resSum3, gamProp3, intProp3, resProp3,
         gamSum4, intSum4, resSum4, gamProp4, intProp4, resProp4,
         c3LL, c3LR, c3RL, c3RR, c4LL, c4LR, c4RL, c4RR, flavWt;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};
 
//==========================================================================

// A derived class for f fbar' -> Z0 W+-. (Here pure Z0, unfortunately.)

class Sigma2ffbar2ZW : public Sigma2ffbargmZWgmZW {

public:

  // Constructor.
  Sigma2ffbar2ZW() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for Z0 and W+- decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()       const {return "f fbar' -> Z0 W+- (no gamma*!)";}
  virtual int    code()       const {return 232;}
  virtual string inFlux()     const {return "ffbarChg";}
  virtual int    id3Mass()    const {return 23;}
  virtual int    id4Mass()    const {return 24;}
  virtual int    resonanceA() const {return 24;}

private:

  // Store W+- mass and width, and couplings.
  double mW, widW, mWS, mwWS, sin2thetaW, cos2thetaW, thetaWRat, cotT,
         thetaWpt, thetaWmm, lun, lde, sigma0, openFracPos, openFracNeg;

};
 
//==========================================================================

// A derived class for f fbar -> W+ W-.

class Sigma2ffbar2WW : public Sigma2ffbargmZWgmZW {

public:

  // Constructor.
  Sigma2ffbar2WW() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W+ and W- decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()       const {return "f fbar -> W+ W-";}
  virtual int    code()       const {return 233;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    id3Mass()    const {return 24;}
  virtual int    id4Mass()    const {return -24;}
  virtual int    resonanceA() const {return 23;}

private:

  // Store Z0 mass and width.
  double mZ, widZ, mZS, mwZS, thetaWRat, sigma0, cgg, cgZ, cZZ, cfg, 
    cfZ, cff, gSS, gTT, gST, gUU, gSU, openFracPair;

};
 
//==========================================================================

// An intermediate class for f fbar -> gamma*/Z0 g/gamma and permutations.

class Sigma2ffbargmZggm : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbargmZggm() {}

  // Initialize process. 
  virtual void initProc(); 

  // Evaluate weight for gamma&/Z0 decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

protected:

  // Parameters set at initialization or for each new event. 
  int    gmZmode;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, 
         gamSum, intSum, resSum, gamProp, intProp, resProp;

  // Evaluate current sum of flavour couplings times phase space. 
  void flavSum(); 

  // Evaluate current propagator terms of cross section. 
  void propTerm(); 

private:

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};
 
//==========================================================================

// A derived class for q qbar -> gamma*/Z0 g.

class Sigma2qqbar2gmZg : public Sigma2ffbargmZggm {

public:

  // Constructor.
  Sigma2qqbar2gmZg() {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "q qbar -> gamma*/Z0 g";}
  virtual int    code()    const {return 241;}
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return 23;}

private:

  // Values stored for later use.
  double sigma0;

};
 
//==========================================================================

// A derived class for q g -> gamma*/Z0 q.

class Sigma2qg2gmZq : public Sigma2ffbargmZggm {

public:

  // Constructor.
  Sigma2qg2gmZq() {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "q g-> gamma*/Z0 q";}
  virtual int    code()    const {return 242;}
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return 23;}

private:

  // Values stored for later use.
  double sigma0;

};
 
//==========================================================================

// A derived class for f fbar' -> gamma*/Z0 gamma.

class Sigma2ffbar2gmZgm : public Sigma2ffbargmZggm {

public:

  // Constructor.
  Sigma2ffbar2gmZgm() {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "f fbar -> gamma*/Z0 gamma";}
  virtual int    code()    const {return 243;}
  virtual string inFlux()  const {return "ffbarSame";}
  virtual int    id3Mass() const {return 23;}

private:

  // Values stored for later use.
  double sigma0;

};
 
//==========================================================================

// A derived class for f gamma -> gamma*/Z0 f.

class Sigma2fgm2gmZf : public Sigma2ffbargmZggm {

public:

  // Constructor.
  Sigma2fgm2gmZf() {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "f gamma -> gamma*/Z0 f";}
  virtual int    code()    const {return 244;}
  virtual string inFlux()  const {return "fgm";}
  virtual int    id3Mass() const {return 23;}

private:

  // Values stored for later use.
  double sigma0;

};
 
//==========================================================================

// An intermediate class for f fbar -> W+- g/gamma and permutations.

class Sigma2ffbarWggm : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbarWggm() {}

  // Evaluate weight for gamma&/Z0 decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

private:

};
 
//==========================================================================

// A derived class for q qbar' -> W+- g.

class Sigma2qqbar2Wg : public Sigma2ffbarWggm {

public:

  // Constructor.
  Sigma2qqbar2Wg() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "q qbar' -> W+- g";}
  virtual int    code()    const {return 251;}
  virtual string inFlux()  const {return "ffbarChg";}
  virtual int    id3Mass() const {return 24;}

private:

  // Values stored for later use.
  double sigma0, openFracPos, openFracNeg;

};
 
//==========================================================================

// A derived class for q g -> W+- q'.

class Sigma2qg2Wq : public Sigma2ffbarWggm {

public:

  // Constructor.
  Sigma2qg2Wq() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "q g-> W+- q'";}
  virtual int    code()    const {return 252;}
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return 24;}

private:

  // Values stored for later use.
  double sigma0, openFracPos, openFracNeg;

};
 
//==========================================================================

// A derived class for f fbar' -> W+- gamma.

class Sigma2ffbar2Wgm : public Sigma2ffbarWggm {

public:

  // Constructor.
  Sigma2ffbar2Wgm() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "f fbar' -> W+- gamma";}
  virtual int    code()    const {return 253;}
  virtual string inFlux()  const {return "ffbarChg";}
  virtual int    id3Mass() const {return 24;}

private:

  // Values stored for later use.
  double sigma0, openFracPos, openFracNeg;

};
 
//==========================================================================

// A derived class for f gamma -> W+- f'.

class Sigma2fgm2Wf : public Sigma2ffbarWggm {

public:

  // Constructor.
  Sigma2fgm2Wf() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "f gamma -> W+- f'";}
  virtual int    code()    const {return 254;}
  virtual string inFlux()  const {return "fgm";}
  virtual int    id3Mass() const {return 24;}

private:

  // Values stored for later use.
  double sigma0, openFracPos, openFracNeg;

};
//==========================================================================

// A derived class for gamma gamma -> f fbar.

class Sigma2gmgm2ffbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2gmgm2ffbar(int idIn, int codeIn) : idNew(idIn), codeSave(codeIn) {}

  // Initialize process. 
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "gmgm";}
  virtual int    id3Mass() const {return idMass;}
  virtual int    id4Mass() const {return idMass;}

private:

  // Member variables.
  int    idNew, codeSave, idMass, idNow;
  string nameSave;
  double ef4, s34Avg, sigTU, sigma, openFracPair;

};
 
//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaEW_H
