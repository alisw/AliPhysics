// SigmaOnia.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for charmonia/bottomonia process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef Pythia8_SigmaOnia_H
#define Pythia8_SigmaOnia_H

#include "SigmaProcess.h"

namespace Pythia8 {
 
//==========================================================================

// A derived class for g g -> QQbar[3S1(1)] g (Q = c or b).

class Sigma2gg2QQbar3S11g : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2QQbar3S11g(int idIn, int codeIn) : idNew(idIn), 
    codeSave(codeIn) {}

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
  virtual int    id3Mass() const {return idHad;}

 private:

  // Values stored for process type and colour flow selection.
  int    idNew, idHad, codeSave;
  string nameSave;
  double oniumME, sigma;

};
 
//==========================================================================

// A derived class for g g -> QQbar[3PJ(1)] g (Q = c or b, J = 0, 1 or 2).

class Sigma2gg2QQbar3PJ1g : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2QQbar3PJ1g(int idIn, int jIn, int codeIn) : idNew(idIn), 
    jSave(jIn), codeSave(codeIn) {}

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
  virtual int    id3Mass() const {return idHad;}

 private:

  // Values stored for process type and colour flow selection.
  int    idNew, idHad, jSave, codeSave;
  string nameSave;
  double oniumME, sigma;

};
 
//==========================================================================

// A derived class for q g -> QQbar[3PJ(1)] q (Q = c or b, J = 0, 1 or 2).

class Sigma2qg2QQbar3PJ1q : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2QQbar3PJ1q(int idIn, int jIn, int codeIn) : idNew(idIn), 
    jSave(jIn), codeSave(codeIn) {}

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
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return idHad;}

 private:

  // Values stored for process type and colour flow selection.
  int    idNew, idHad, jSave, codeSave;
  string nameSave;
  double oniumME, sigma;

};
 
//==========================================================================

// A derived class for q qbar -> QQbar[3PJ(1)] g (Q = c or b, J = 0, 1 or 2).

class Sigma2qqbar2QQbar3PJ1g : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2QQbar3PJ1g(int idIn, int jIn, int codeIn) : idNew(idIn), 
    jSave(jIn), codeSave(codeIn) {}

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
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return idHad;}

 private:

  // Values stored for process type and colour flow selection.
  int    idNew, idHad, jSave, codeSave;
  string nameSave;
  double oniumME, sigma;

};
 
//==========================================================================

// A derived class for g g -> QQbar[X(8)] g (Q = c or b, X = 3S1, 1S0 or 3PJ).

class Sigma2gg2QQbarX8g : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2QQbarX8g(int idIn, int stateIn, int codeIn) : idNew(idIn), 
    stateSave(stateIn), codeSave(codeIn) {}

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
  virtual int    id3Mass() const {return idHad;}

 private:

  // Values stored for process type and colour flow selection.
  int    idNew, idHad, stateSave, codeSave;
  string nameSave;
  double oniumME, sigma;

};
 
//==========================================================================

// A derived class for q g -> QQbar[X(8)] q (Q = c or b, X = 3S1, 1S0 or 3PJ).

class Sigma2qg2QQbarX8q : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2QQbarX8q(int idIn, int stateIn, int codeIn) : idNew(idIn), 
    stateSave(stateIn), codeSave(codeIn) {}

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
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return idHad;}

 private:

  // Values stored for process type and colour flow selection.
  int    idNew, idHad, stateSave, codeSave;
  string nameSave;
  double oniumME, sigma;

};
 
//==========================================================================

// A derived class for q qbar -> QQbar[X(8)] g (Q = c or b, 
//   X = 3S1, 1S0 or 3PJ).

class Sigma2qqbar2QQbarX8g : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2QQbarX8g(int idIn, int stateIn, int codeIn) : idNew(idIn), 
    stateSave(stateIn), codeSave(codeIn) {}

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
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return idHad;}

 private:

  // Values stored for process type and colour flow selection.
  int    idNew, idHad, stateSave, codeSave;
  string nameSave;
  double oniumME, sigma;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaOnia_H
