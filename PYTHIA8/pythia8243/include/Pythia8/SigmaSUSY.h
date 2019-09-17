// SigmaSUSY.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Supersymmetric process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef Pythia8_SigmaSUSY_H
#define Pythia8_SigmaSUSY_H

#include "Pythia8/PhaseSpace.h"
#include "Pythia8/PythiaComplex.h"
#include "Pythia8/SigmaProcess.h"
#include "Pythia8/SusyCouplings.h"

namespace Pythia8 {

//==========================================================================

// An intermediate class for SUSY 2 -> 2 with nontrivial decay angles.

class Sigma2SUSY : public Sigma2Process {

public:

  // Constructor.
  Sigma2SUSY() { };

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

 };

//==========================================================================

// A derived class for q qbar -> neutralino_i neutralino_j.

class Sigma2qqbar2chi0chi0 : public Sigma2SUSY {

public:

  // Constructor.
  Sigma2qqbar2chi0chi0() : id3chi(), id4chi(), codeSave(), sigma0(),
    ui(), uj(), ti(), tj(), openFracPair(), coupSUSYPtr() {};

  // Constructor.
  Sigma2qqbar2chi0chi0(int id3chiIn, int id4chiIn, int codeIn) : sigma0(),
    ui(), uj(), ti(), tj(), openFracPair(), coupSUSYPtr() {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4chi   = id4chiIn;
    codeSave = codeIn;


    // Construct id codes from ordering indices.
    id3                  = 1000022;
    if (id3chi == 2) id3 = 1000023;
    if (id3chi == 3) id3 = 1000025;
    if (id3chi == 4) id3 = 1000035;
    if (id3chi == 5) id3 = 1000045;
    id4                  = 1000022;
    if (id4chi == 2) id4 = 1000023;
    if (id4chi == 3) id4 = 1000025;
    if (id4chi == 4) id4 = 1000035;
    if (id4chi == 5) id4 = 1000045;

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  //  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "ff";}
  virtual int    id3Mass() const {return abs(id3);}
  virtual int    id4Mass() const {return abs(id4);}
  virtual int    resonanceA() const {return 23;}
  virtual bool   isSUSY()  const {return true;}
  virtual double getSigma0() const {return sigma0;}

 protected:

  // Basic process information
  int     id3chi, id4chi, codeSave;
  string  nameSave;

  // Values stored for later use
  double  sigma0, ui, uj, ti, tj, openFracPair;
  complex propZ;

  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q qbar -> neutralino_i chargino_j.

class Sigma2qqbar2charchi0 : public Sigma2qqbar2chi0chi0 {

public:

  // Constructor.
  Sigma2qqbar2charchi0(int id3chiIn, int id4chiIn, int codeIn) {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4chi   = id4chiIn;
    codeSave = codeIn;

    // Construct id codes from ordering indices.
    id3 = (abs(id3chi) == 2) ? 1000037 : 1000024;
    if (id3chi < 0)  id3 = -id3;

    id4                  = 1000022;
    if (id4chi == 2) id4 = 1000023;
    if (id4chi == 3) id4 = 1000025;
    if (id4chi == 4) id4 = 1000035;
    if (id4chi == 5) id4 = 1000045;

  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  virtual int    resonanceA() const {return 24;}

protected :

  complex propW;

};

//==========================================================================

// A derived class for q qbar -> chargino+_i chargino-_j.

class Sigma2qqbar2charchar : public Sigma2qqbar2chi0chi0 {

public:

  // Constructor.
  Sigma2qqbar2charchar(int id3chiIn, int id4chiIn, int codeIn) {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4chi   = id4chiIn;
    codeSave = codeIn;

    // Construct id codes from ordering indices.
    id3 = (abs(id3chi) == 2) ?  1000037 :  1000024;
    id4 = (abs(id4chi) == 2) ? -1000037 : -1000024;

  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

};

//==========================================================================

// A derived class for q g -> neutralino_i squark_j (and cc)

class Sigma2qg2chi0squark : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2chi0squark() : id3chi(), id4sq(), codeSave(), sigma0(),
    ui(), uj(), ti(), tj(), openFracPair(), coupSUSYPtr() { };

  // Constructor.
  Sigma2qg2chi0squark(int id3chiIn, int id4sqIn, bool isUp, int codeIn) :
    sigma0(), ui(), uj(), ti(), tj(), openFracPair(), coupSUSYPtr() {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4sq    = id4sqIn;
    codeSave = codeIn;

    // Construct id codes from ordering indices.
    id3                  = 1000022;
    if (id3chi == 2) id3 = 1000023;
    if (id3chi == 3) id3 = 1000025;
    if (id3chi == 4) id3 = 1000035;
    if (id3chi == 5) id3 = 1000045;
    id4                  = 1000001 + (isUp ? 1 : 0);
    if (id4sq  == 2) id4 = 1000003 + (isUp ? 1 : 0);
    if (id4sq  == 3) id4 = 1000005 + (isUp ? 1 : 0);
    if (id4sq  == 4) id4 = 2000001 + (isUp ? 1 : 0);
    if (id4sq  == 5) id4 = 2000003 + (isUp ? 1 : 0);
    if (id4sq  == 6) id4 = 2000005 + (isUp ? 1 : 0);

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return abs(id3);}
  virtual int    id4Mass() const {return abs(id4);}
  virtual bool   isSUSY()  const {return true;}

 protected:

  // Basic process information
  int     id3chi, id4sq, codeSave;
  string  nameSave;

  // Values stored for later use
  double  sigma0, ui, uj, ti, tj, openFracPair;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q g -> chargino_i squark_j (incl cc)

class Sigma2qg2charsquark : public Sigma2qg2chi0squark {

public:

  // Constructor.
  Sigma2qg2charsquark(int id3chiIn, int id4sqIn, bool isUp, int codeIn) {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4sq    = id4sqIn;
    codeSave = codeIn;

    // Construct id codes from ordering indices.
    id3Sav                       = 1000024;
    if (abs(id3chi) == 2) id3Sav = 1000037;
    if (isUp)             id3Sav = -id3Sav;
    id4Sav                       = 1000001 + (isUp ? 1 : 0);
    if (id4sq  == 2) id4Sav      = 1000003 + (isUp ? 1 : 0);
    if (id4sq  == 3) id4Sav      = 1000005 + (isUp ? 1 : 0);
    if (id4sq  == 4) id4Sav      = 2000001 + (isUp ? 1 : 0);
    if (id4sq  == 5) id4Sav      = 2000003 + (isUp ? 1 : 0);
    if (id4sq  == 6) id4Sav      = 2000005 + (isUp ? 1 : 0);

    // Initial values, can be swapped to charge conjugates event by event.
    id3 = id3Sav;
    id4 = id4Sav;

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

 private:

  // Basic process information
  int id3Sav, id4Sav;

};

//==========================================================================

// A derived class for q q' -> ~q_i ~q_j

class Sigma2qq2squarksquark : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2squarksquark() : id3Sav(), id4Sav(), codeSave(), iGen3(), iGen4(),
    nNeut(), isUD(), onlyQCD(), m2Glu(), sigmaChar(), sigmaNeut(), sigmaGlu(),
    sigmaCharNeut(), sigmaCharGlu(), sigmaNeutGlu(), openFracPair(), tGlu(),
    uGlu(), sumCt(), sumCu(), sumNt(), sumNu(), sumGt(), sumGu(),
    sumInterference(), coupSUSYPtr() {}

  // Constructor.
  Sigma2qq2squarksquark(int id3In, int id4In, int codeIn) : iGen3(), iGen4(),
    nNeut(), isUD(), onlyQCD(), m2Glu(), sigmaChar(), sigmaNeut(), sigmaGlu(),
    sigmaCharNeut(), sigmaCharGlu(), sigmaNeutGlu(), openFracPair(), tGlu(),
    uGlu(), sumCt(), sumCu(), sumNt(), sumNu(), sumGt(), sumGu(),
    sumInterference(), coupSUSYPtr() {

    // Save ordering indices and process code
    id3Sav = id3In;
    id4Sav = id4In;
    codeSave = codeIn;
    // Initial values (flipped for c.c.)
    id3    = id3Sav;
    id4    = id4Sav;

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qq";}
  virtual int    id3Mass() const {return abs(id3Sav);}
  virtual int    id4Mass() const {return abs(id4Sav);}
  virtual bool   isSUSY()  const {return true;}

private:

  // Basic process information
  int     id3Sav, id4Sav, codeSave, iGen3, iGen4, nNeut;
  string  nameSave;
  bool    isUD, onlyQCD;

  // Storage of mass squares
  double m2Glu;
  vector<double> m2Neut, m2Char;

  // Flavor-independent prefactors.
  double sigmaChar, sigmaNeut, sigmaGlu;
  double sigmaCharNeut, sigmaCharGlu, sigmaNeutGlu;
  double openFracPair;

  // Point-by-point info
  double tGlu, uGlu;
  vector<double> tNeut, uNeut, tChar, uChar;
  double sumCt, sumCu, sumNt, sumNu, sumGt, sumGu, sumInterference;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;
};

//==========================================================================

// A derived class for q qbar' -> ~q_i ~q*_j

class Sigma2qqbar2squarkantisquark : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2squarkantisquark() : id3Sav(), id4Sav(), codeSave(), iGen3(),
    iGen4(), nNeut(), isUD(), isCC(), onlyQCD(), m2Glu(), xW(), openFracPair(),
    sigmaEW(), sigmaGlu(), sigmaEWG(), tGlu(), uGlu(), sumColS(), sumColT(),
    sumInterference(), coupSUSYPtr() {}

  // Constructor.
  Sigma2qqbar2squarkantisquark(int id3In, int id4In, int codeIn) : iGen3(),
    iGen4(), nNeut(), isUD(), isCC(), onlyQCD(), m2Glu(), xW(), openFracPair(),
    sigmaEW(), sigmaGlu(), sigmaEWG(), tGlu(), uGlu(), sumColS(), sumColT(),
    sumInterference(), coupSUSYPtr() {

    // Save ordering indices and process code
    // (always store squark first, antisquark second)
    id3Sav = abs(id3In);
    id4Sav = -abs(id4In);
    codeSave = codeIn;
    // Initial values
    id3    = id3Sav;
    id4    = id4Sav;

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qq";}
  virtual int    id3Mass() const {return abs(id3Sav);}
  virtual int    id4Mass() const {return abs(id4Sav);}
  virtual bool   isSUSY()  const {return true;}

private:

  // Basic process information
  int     id3Sav, id4Sav, codeSave, iGen3, iGen4, nNeut;
  string  nameSave;
  bool    isUD, isCC, onlyQCD;

  // Storage of mass squares
  double m2Glu;
  vector<double> m2Neut;

  // Flavor-independent prefactors: EW, strong, and interference
  double xW;
  double openFracPair;
  double sigmaEW, sigmaGlu, sigmaEWG;

  // Point-by-point info
  double tGlu, uGlu;
  vector<double> tNeut, uNeut;
  complex propZW;
  double sumColS, sumColT, sumInterference;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for g g -> ~q ~q*

class Sigma2gg2squarkantisquark : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2squarkantisquark() : id3Sav(), id4Sav(), codeSave(), sigma(),
    m2Sq(), openFracPair(), coupSUSYPtr() { }

  // Constructor.
  Sigma2gg2squarkantisquark(int id34In, int codeIn) : sigma(), m2Sq(),
    openFracPair(), coupSUSYPtr() {

    // Save ordering indices and process code
    // (always store squark first, antisquark second)
    id3Sav = abs(id34In);
    id4Sav = -abs(id34In);
    codeSave = codeIn;
    // Initial values
    id3    = id3Sav;
    id4    = id4Sav;

  }

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
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return abs(id3Sav);}
  virtual int    id4Mass() const {return abs(id4Sav);}
  virtual bool   isSUSY()  const {return true;}

private:

  // Basic process information
  int     id3Sav, id4Sav, codeSave;
  string  nameSave;
  double sigma, m2Sq, openFracPair;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q g -> ~q ~g

class Sigma2qg2squarkgluino : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2squarkgluino() : codeSave(), sigmaA(), sigmaB(), comFacHat(),
    m2Glu(), m2Sq(), openFracPair(), coupSUSYPtr() {}

  // Constructor.
  Sigma2qg2squarkgluino(int id3In, int codeIn) : sigmaA(), sigmaB(),
    comFacHat(), m2Glu(), m2Sq(), openFracPair(), coupSUSYPtr() {

    // Save ordering indices and process code
    codeSave = codeIn;
    // Initial values
    id3    = id3In;
    id4    = 1000021;

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return abs(id3);}
  virtual int    id4Mass() const {return 1000021;}
  virtual bool   isSUSY()  const {return true;}

private:

  // Basic process information
  int     codeSave;
  string  nameSave;
  double sigmaA, sigmaB, comFacHat, m2Glu, m2Sq, openFracPair;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for g g -> gluino gluino.

class Sigma2gg2gluinogluino : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2gluinogluino() : sigTS(), sigUS(), sigTU(), sigSum(), sigma(),
    openFracPair(), coupSUSYPtr() { }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "g g -> gluino gluino";}
  virtual int    code()    const {return 1201;}
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return 1000021;}
  virtual int    id4Mass() const {return 1000021;}
  virtual bool   isSUSY()  const {return true;}

private:

  // Values stored for process type and colour flow selection.
  double sigTS, sigUS, sigTU, sigSum, sigma, openFracPair;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q qbar -> gluino gluino.

class Sigma2qqbar2gluinogluino : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2gluinogluino() : openFracPair(), s34Avg(), sigS(), tHG(), uHG(),
    tHG2(), uHG2(), coupSUSYPtr() { }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "q qbar -> gluino gluino";}
  virtual int    code()    const {return 1202;}
  virtual string inFlux()  const {return "qq";}
  virtual int    id3Mass() const {return 1000021;}
  virtual int    id4Mass() const {return 1000021;}
  virtual bool   isSUSY()  const {return true;}

private:

  // Values stored for process type and colour flow selection.
  double openFracPair, s34Avg, sigS, tHG, uHG, tHG2, uHG2;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

class Sigma1qq2antisquark : public Sigma1Process {
public:

  // Constructor.
  Sigma1qq2antisquark() : mRes(), GammaRes(), m2Res(), sigBW(), widthOut(),
    codeSave(), idRes(), coupSUSYPtr() {}


  Sigma1qq2antisquark(int id3In) : mRes(), GammaRes(), m2Res(), sigBW(),
    widthOut(), codeSave(), coupSUSYPtr() {

    idRes = id3In;

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qq";}
  virtual bool   isSUSY()  const {return true;}
  virtual bool   isRPV()   const {return true;}
  virtual int    resonanceA() const {return idRes;}

private:

  // Values stored for process type and colour flow selection.
  double mRes, GammaRes, m2Res, sigBW, widthOut;
  int    codeSave, idRes;
  string nameSave;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};


//==========================================================================

// A derived class for q qbar -> neutralino_i gluino.

class Sigma2qqbar2chi0gluino : public Sigma2SUSY {

public:

  // Constructor.
  Sigma2qqbar2chi0gluino() : id3chi(), id4chi(), codeSave(), sigma0(), ui(),
    uj(), ti(), tj(), openFracPair(), coupSUSYPtr() {};

  // Constructor.
  Sigma2qqbar2chi0gluino(int id4chiIn, int codeIn) : id3chi(), sigma0(), ui(),
    uj(), ti(), tj(), openFracPair(), coupSUSYPtr() {

    // Save ordering indices and process code
    id3   = 1000021;
    id4chi   = id4chiIn;
    codeSave = codeIn;


    // Construct id codes from ordering indices.
    id4                  = 1000022;
    if (id4chi == 2) id4 = 1000023;
    if (id4chi == 3) id4 = 1000025;
    if (id4chi == 4) id4 = 1000035;
    if (id4chi == 5) id4 = 1000045;

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  //  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "ff";}
  virtual int    id3Mass() const {return abs(id3);}
  virtual int    id4Mass() const {return abs(id4);}
  virtual int    resonanceA() const {return 23;}
  virtual bool   isSUSY()  const {return true;}
  virtual double getSigma0() const {return sigma0;}

 protected:

  // Basic process information
  int     id3chi, id4chi, codeSave;
  string  nameSave;

  // Values stored for later use
  double  sigma0, ui, uj, ti, tj, openFracPair;

  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q qbar -> neutralino_i chargino_j.

class Sigma2qqbar2chargluino : public Sigma2qqbar2chi0gluino {

public:

  // Constructor.
  Sigma2qqbar2chargluino(int id4chiIn, int codeIn) {

    // Save ordering indices and process code
    id3   = 1000021;
    id4chi   = id4chiIn;
    codeSave = codeIn;

    // Construct id codes from ordering indices.
    id4 = (abs(id4chi) == 2) ? 1000037 : 1000024;
    if (id4chi < 0)  id4 = -id4;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  virtual int    resonanceA() const {return 24;}

protected :

  complex propW;

};

//==========================================================================

// A derived class for q qbar' -> ~q_i ~q*_j

class Sigma2qqbar2sleptonantislepton : public Sigma2qqbar2squarkantisquark {

public:

  // Constructor.
  Sigma2qqbar2sleptonantislepton() : id3Sav(), id4Sav(), codeSave(), iGen3(),
    iGen4(), nNeut(), isUD(), xW(), openFracPair(), sigmaEW(), sumColS(),
    sumColT(), sumInterference(), coupSUSYPtr() {}

  // Constructor.
  Sigma2qqbar2sleptonantislepton(int id3In, int id4In, int codeIn) : iGen3(),
    iGen4(), nNeut(), isUD(), xW(), openFracPair(), sigmaEW(), sumColS(),
    sumColT(), sumInterference(), coupSUSYPtr() {

    // Save ordering indices and process code
    // (always store squark first, antisquark second)
    id3Sav = abs(id3In);
    id4Sav = -abs(id4In);
    codeSave = codeIn;
    // Initial values
    id3    = id3Sav;
    id4    = id4Sav;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qq";}
  virtual int    id3Mass() const {return abs(id3Sav);}
  virtual int    id4Mass() const {return abs(id4Sav);}
  virtual bool   isSUSY()  const {return true;}

private:

  // Basic process information
  int     id3Sav, id4Sav, codeSave, iGen3, iGen4, nNeut;
  string  nameSave;
  bool    isUD;

  // Storage of mass squares
  vector<double> m2Neut;

  // Flavor-independent prefactors: EW, strong, and interference
  double xW;
  double openFracPair;
  double sigmaEW;

  // Point-by-point info
  vector<double> tNeut, uNeut;
  complex propZW;
  double sumColS, sumColT, sumInterference;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaSUSY_H
