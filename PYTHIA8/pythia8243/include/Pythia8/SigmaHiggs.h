// SigmaHiggs.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// Part of code written by Marc Montull, CERN summer student 2007.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Higgs process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef Pythia8_SigmaHiggs_H
#define Pythia8_SigmaHiggs_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A derived class for f fbar -> H0 (SM), H1, H2 or A3 (BSM).

class Sigma1ffbar2H : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2H(int higgsTypeIn) : HResPtr(), mRes(), GammaRes(), m2Res(),
    GamMRat(), sigBW(), widthOut(), higgsType(higgsTypeIn), codeSave(),
    idRes() {}

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
  virtual string name()       const {return nameSave;}
  virtual int    code()       const {return codeSave;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    resonanceA() const {return idRes;}

private:

  // An H0, H1, H2 or A3 resonance object provides coupling
  // and propagator expressions.
  ParticleDataEntry* HResPtr;
  double mRes, GammaRes, m2Res, GamMRat, sigBW, widthOut;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for g g -> H0 (SM), H1, H2 or A3 (BSM).

class Sigma1gg2H : public Sigma1Process {

public:

  // Constructor.
  Sigma1gg2H(int higgsTypeIn) : HResPtr(), mRes(), GammaRes(), m2Res(),
    GamMRat(), sigma(), higgsType(higgsTypeIn), codeSave(), idRes() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return nameSave ;}
  virtual int    code()       const {return codeSave;}
  virtual string inFlux()     const {return "gg";}
  virtual int    resonanceA() const {return idRes;}

private:

  // A H0, H1, H2 or A3 resonance object provides coupling
  // and propagator expressions.
  ParticleDataEntry* HResPtr;
  double mRes, GammaRes, m2Res, GamMRat, sigma;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for gamma gamma -> H0 (SM Higgs), H1, H2 or A3 (BSM Higgs).

class Sigma1gmgm2H : public Sigma1Process {

public:

  // Constructor.
  Sigma1gmgm2H(int higgsTypeIn) : HResPtr(), mRes(), GammaRes(), m2Res(),
    GamMRat(), sigma(), higgsType(higgsTypeIn), codeSave(), idRes() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return nameSave;}
  virtual int    code()       const {return codeSave;}
  virtual string inFlux()     const {return "gmgm";}
  virtual int    resonanceA() const {return idRes;}

private:

  // A H0, H1, H2 or A3 resonance object provides coupling
  // and propagator expressions.
  ParticleDataEntry* HResPtr;
  double mRes, GammaRes, m2Res, GamMRat, sigma;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for f fbar -> H Z0.
// (H can be H0 SM or H1, H2, A3 from BSM).
class Sigma2ffbar2HZ : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2HZ(int higgsTypeIn) : mZ(), widZ(), mZS(), mwZS(), thetaWRat(),
    sigma0(), openFracPair(), coup2Z(), higgsType(higgsTypeIn), codeSave(),
    idRes() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return nameSave;}
  virtual int    code()       const {return codeSave;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual bool   isSChannel() const {return true;}
  virtual int    id3Mass()    const {return idRes;}
  virtual int    id4Mass()    const {return 23;}
  virtual int    resonanceA() const {return 23;}
  virtual int    gmZmode()    const {return 2;}

private:

  // Store Z0 mass and width.
  double mZ, widZ, mZS, mwZS, thetaWRat, sigma0, openFracPair, coup2Z;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for f fbar -> H W+- (Standard Model Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2ffbar2HW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2HW(int higgsTypeIn) : mW(), widW(), mWS(), mwWS(), thetaWRat(),
    sigma0(), openFracPairPos(), openFracPairNeg(), coup2W(),
    higgsType(higgsTypeIn), codeSave(), idRes() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return nameSave;}
  virtual int    code()       const {return codeSave;}
  virtual string inFlux()     const {return "ffbarChg";}
  virtual bool   isSChannel() const {return true;}
  virtual int    id3Mass()    const {return idRes;}
  virtual int    id4Mass()    const {return 24;}
  virtual int    resonanceA() const {return 24;}

private:

  // Store W+- mass and width, and couplings.
  double mW, widW, mWS, mwWS, thetaWRat, sigma0, openFracPairPos,
         openFracPairNeg, coup2W;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for f f' -> H f f' (Z0 Z0 fusion of SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma3ff2HfftZZ : public Sigma3Process {

public:

  // Constructor.
  Sigma3ff2HfftZZ(int higgsTypeIn) : mZS(), prefac(), sigma1(), sigma2(),
    openFrac(), coup2Z(), higgsType(higgsTypeIn), codeSave(), idRes() {}

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
  virtual int    id3Mass() const {return idRes;}

  // Instructions for 3-body phase space with t-channel propagators.
  virtual int    idTchan1()        const {return 23;}
  virtual int    idTchan2()        const {return 23;}
  virtual double tChanFracPow1()   const {return 0.05;}
  virtual double tChanFracPow2()   const {return 0.9;}
  virtual bool   useMirrorWeight() const {return true;}

private:

  // Store standard factors.
  double mZS, prefac, sigma1, sigma2, openFrac, coup2Z;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for f_1 f_2 -> H f_3 f_4 (W+ W- fusion of SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma3ff2HfftWW : public Sigma3Process {

public:

  // Constructor.
  Sigma3ff2HfftWW(int higgsTypeIn) : mWS(), prefac(), sigma0(), openFrac(),
    coup2W(), higgsType(higgsTypeIn), codeSave(), idRes() {}

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
  virtual int    id3Mass() const {return idRes;}

  // Instructions for 3-body phase space with t-channel propagators.
  virtual int    idTchan1()        const {return 24;}
  virtual int    idTchan2()        const {return 24;}
  virtual double tChanFracPow1()   const {return 0.05;}
  virtual double tChanFracPow2()   const {return 0.9;}
  virtual bool   useMirrorWeight() const {return true;}

private:

  // Store standard prefactor.
  double mWS, prefac, sigma0, openFrac, coup2W;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for g g -> H Q Qbar (Q Qbar fusion of SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma3gg2HQQbar : public Sigma3Process {

public:

  // Constructor.
  Sigma3gg2HQQbar(int idIn, int higgsTypeIn) : prefac(), sigma(),
    openFracTriplet(), coup2Q(), idNew(idIn), higgsType(higgsTypeIn),
    codeSave(), idRes() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return idRes;}
  virtual int    id4Mass() const {return idNew;}
  virtual int    id5Mass() const {return idNew;}

  // Instructions for 3-body phase space with t-channel propagators.
  virtual int    idTchan1()        const {return idNew;}
  virtual int    idTchan2()        const {return idNew;}
  virtual double tChanFracPow1()   const {return 0.4;}
  virtual double tChanFracPow2()   const {return 0.2;}
  virtual bool   useMirrorWeight() const {return false;}

private:

  // Store flavour-specific process information and standard prefactor.
  double prefac, sigma, openFracTriplet, coup2Q;
  int    idNew, higgsType, codeSave, idRes;
  string nameSave;

};

//==========================================================================

// A derived class for q qbar -> H Q Qbar (Q Qbar fusion of SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma3qqbar2HQQbar : public Sigma3Process {

public:

  // Constructor.
  Sigma3qqbar2HQQbar(int idIn, int higgsTypeIn) : prefac(), sigma(),
    openFracTriplet(), coup2Q(), idNew(idIn), higgsType(higgsTypeIn),
    codeSave(), idRes() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return idRes;}
  virtual int    id4Mass() const {return idNew;}
  virtual int    id5Mass() const {return idNew;}

  // Instructions for 3-body phase space with t-channel propagators.
  virtual int    idTchan1()        const {return idNew;}
  virtual int    idTchan2()        const {return idNew;}
  virtual double tChanFracPow1()   const {return 0.4;}
  virtual double tChanFracPow2()   const {return 0.2;}
  virtual bool   useMirrorWeight() const {return false;}

private:

  // Store flavour-specific process information and standard prefactor.
  double prefac, sigma, openFracTriplet, coup2Q;
  int    idNew, higgsType, codeSave, idRes;
  string nameSave;

};

//==========================================================================

// A derived class for q g -> H q (SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2qg2Hq : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2Hq(int idIn, int higgsTypeIn) : m2W(), thetaWRat(), sigma(),
    openFrac(), idNew(idIn), higgsType(higgsTypeIn), codeSave(), idRes() {}

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
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return idRes;}
  virtual int    id4Mass() const {return idNew;}

private:

  // Store flavour-specific process information and standard prefactor.
  double m2W, thetaWRat, sigma, openFrac;
  int    idNew, higgsType, codeSave, idRes;
  string nameSave;

};

//==========================================================================

// A derived class for g g -> H0 g (SM or BSM Higgs via heavy top loop).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2gg2Hglt : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2Hglt(int higgsTypeIn) : widHgg(), sigma(), openFrac(),
    higgsType(higgsTypeIn), codeSave(), idRes() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return idRes;}

private:

  // Store standard prefactor.
  double widHgg, sigma, openFrac;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for q g -> H q (SM or BSM Higgs via heavy top loop).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2qg2Hqlt : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2Hqlt(int higgsTypeIn) : widHgg(), sigma(), openFrac(),
    higgsType(higgsTypeIn), codeSave(), idRes() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return idRes;}

private:

  // Store standard prefactor.
  double widHgg, sigma, openFrac;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for q qbar -> H g (SM or BSM Higgs via heavy top loop).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2qqbar2Hglt : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2Hglt(int higgsTypeIn) : widHgg(), sigma(), openFrac(),
    higgsType(higgsTypeIn), codeSave(), idRes() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return idRes;}

private:

  // Store standard prefactor.
  double widHgg, sigma, openFrac;
  int    higgsType, codeSave, idRes;
  string nameSave;
};

//==========================================================================

// A derived class for f fbar' -> H+-.

class Sigma1ffbar2Hchg : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2Hchg() : HResPtr(), mRes(), GammaRes(), m2Res(), GamMRat(),
    m2W(), thetaWRat(), tan2Beta(), sigBW(), widthOutPos(), widthOutNeg() {}

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
  virtual string name()       const {return "f fbar' -> H+-";}
  virtual int    code()       const {return 1061;}
  virtual string inFlux()     const {return "ffbarChg";}
  virtual int    resonanceA() const {return 37;}

private:

  // A H0 resonance object provides coupling and propagator expressions.
  ParticleDataEntry* HResPtr;
  double mRes, GammaRes, m2Res, GamMRat, m2W, thetaWRat, tan2Beta, sigBW,
         widthOutPos, widthOutNeg;

};

//==========================================================================

// A derived class for q g -> H+- q'.

class Sigma2qg2Hchgq : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2Hchgq(int idIn, int codeIn, string nameIn) : idNew(idIn),
    codeSave(codeIn), idOld(), idUp(), idDn(), nameSave(nameIn), m2W(),
    thetaWRat(), tan2Beta(), sigma(), openFracPos(), openFracNeg() {}

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
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return 37;}
  virtual int    id4Mass() const {return idNew;}

private:

  // Store flavour-specific process information and standard prefactor.
  int    idNew, codeSave, idOld, idUp, idDn;
  string nameSave;
  double m2W, thetaWRat, tan2Beta, sigma, openFracPos, openFracNeg;

};

//==========================================================================

// A derived class for f fbar -> A0(H_3) h0(H_1) or A0(H_3) H0(H_2).

class Sigma2ffbar2A3H12 : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2A3H12(int higgsTypeIn) : higgsType(higgsTypeIn), higgs12(),
    codeSave(), coupZA3H12(), m2Z(), mGammaZ(), thetaWRat(), openFrac(),
    sigma0() {}

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
  virtual string inFlux()  const {return "ffbarSame";}
  virtual int    id3Mass() const {return 36;}
  virtual int    id4Mass() const {return higgs12;}

private:

  // Store flavour-specific process information and standard prefactor.
  int    higgsType, higgs12, codeSave;
  string nameSave;
  double coupZA3H12, m2Z, mGammaZ, thetaWRat, openFrac, sigma0;

};

//==========================================================================

// A derived class for f fbar -> H+- h0(H_1) or H+- H0(H_2).

class Sigma2ffbar2HchgH12 : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2HchgH12(int higgsTypeIn) : higgsType(higgsTypeIn), higgs12(),
    codeSave(), coupWHchgH12(), m2W(), mGammaW(), thetaWRat(), openFracPos(),
    openFracNeg(), sigma0() {}

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
  virtual string inFlux()  const {return "ffbarChg";}
  virtual int    id3Mass() const {return 37;}
  virtual int    id4Mass() const {return higgs12;}

private:

  // Store flavour-specific process information and standard prefactor.
  int    higgsType, higgs12, codeSave;
  string nameSave;
  double coupWHchgH12, m2W, mGammaW, thetaWRat, openFracPos, openFracNeg,
         sigma0;

};

//==========================================================================

// A derived class for f fbar -> H+ H-.

class Sigma2ffbar2HposHneg : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2HposHneg() : m2Z(), mGammaZ(), thetaWRat(), eH(), lH(),
    openFrac(), gamSig(), intSig(), resSig() {}

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
  virtual string name()    const {return "f fbar -> H+ H-";}
  virtual int    code()    const {return 1085;}
  virtual string inFlux()  const {return "ffbarSame";}
  virtual int    id3Mass() const {return 37;}
  virtual int    id4Mass() const {return 37;}

private:

  // Store flavour-specific process information and standard prefactor.
  double m2Z, mGammaZ, thetaWRat, eH, lH, openFrac, gamSig, intSig, resSig;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaHiggs_H
