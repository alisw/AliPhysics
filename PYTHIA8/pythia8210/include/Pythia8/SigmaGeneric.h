// SigmaGeneric.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Johan Bijnens,Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for various generic production processes, to be used as
// building blocks for some BSM processes.
// Currently represented by QCD pair production of colour triplet objects,
// with spin either 0, 1/2 or 1.

#ifndef Pythia8_SigmaGeneric_H
#define Pythia8_SigmaGeneric_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A derived class for g g -> qG qGbar (generic quark of spin 0, 1/2 or 1).

class Sigma2gg2qGqGbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2qGqGbar(int idIn, int codeIn, int spinIn,
    string nameIn = "g g -> qG qGbar") : idNew(idIn), codeSave(codeIn),
    spinSave(spinIn), nameSave(nameIn) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return idNew;}
  virtual int    id4Mass() const {return idNew;}

private:

  // Values stored for process type and colour flow selection.
  int    idNew, codeSave, spinSave, nCHV;
  string nameSave;
  bool   hasKappa;
  double openFracPair, sigma, sigTS, sigUS, sigSum, kappam1;

};

//==========================================================================

// A derived class for q qbar -> qG qGbar (generic quark of spin 0, 1/2 or 1).

class Sigma2qqbar2qGqGbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2qGqGbar(int idIn, int codeIn, int spinIn,
    string nameIn = "q qbar -> qG qGbar") : idNew(idIn), codeSave(codeIn),
    spinSave(spinIn), nameSave(nameIn) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return idNew;}
  virtual int    id4Mass() const {return idNew;}

private:

  // Values stored for process type and colour flow selection.
  int    idNew, codeSave, spinSave, nCHV;
  string nameSave;
  double openFracPair, sigma, sigSum, kappa;

};

//==========================================================================

// A derived class for f fbar -> fG fGbar (generic spin 0, 1/2 or 1 particle)
// via gamma^*/Z^* s-channel exchange. Still under development!! ??

class Sigma2ffbar2fGfGbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2fGfGbar(int idIn, int codeIn, int spinIn,
    string nameIn = "q qbar -> qG qGbar") : idNew(idIn), codeSave(codeIn),
    spinSave(spinIn), nameSave(nameIn) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "ffbarSame";}
  virtual int    id3Mass() const {return idNew;}
  virtual int    id4Mass() const {return idNew;}

private:

  // Values stored for process type and colour flow selection.
  int    idNew, codeSave, spinSave, nCHV;
  string nameSave;
  bool   hasColour;
  double eQHV2, openFracPair, sigma0, sigSum, kappa, colFac;

};

//==========================================================================

// A derived class for f fbar -> Zv, where Zv couples both to the SM and
// to a hidden sector. Primitive coupling structure.

class Sigma1ffbar2Zv : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2Zv() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat) for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return "f fbar -> Zv";}
  virtual int    code()       const {return 4941;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    resonanceA() const {return 4900023;}

private:

  // Store flavour-specific process information and standard prefactor.
  int    idZv;
  double mRes, GammaRes, m2Res, GamMRat, sigOut;

  // Pointer to properties of Zv, to access decay width.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaGeneric_H
