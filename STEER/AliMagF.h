#ifndef ALIMAGF_H
#define ALIMAGF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TNamed.h"

enum Field_t {kUndef=1, kConst=1, kConMesh=2, kDipoMap=3};

class AliMagF : public TNamed {

public:
  AliMagF(){}
  AliMagF(const char *name, const char *title, const Int_t integ, const Int_t map, 
	  const Float_t factor, const Float_t fmax);
  virtual ~AliMagF() {}
  virtual void Field(Float_t *x, Float_t *b);
  virtual Int_t Type() const {return fType;}
  virtual Float_t Max() const {return fMax;}
  virtual Int_t Map() const {return fMap;}
  virtual Int_t Integ() const {return fInteg;}
  virtual Float_t Factor() const {return fFactor;}
  virtual void ReadField() {}
  
protected:
  Int_t     fMap;    // Field Map identifier
  Int_t     fType;   // Mag Field type
  Int_t     fInteg;  // Integration method as indicated in Geant
  Float_t   fFactor; // Multiplicative factor
  Float_t   fMax;    // Max Field as indicated in Geant

  ClassDef(AliMagF,1)  //Base class for all Alice MagField
};

//ZDC part -------------------------------------------------------------------

  static const Float_t kG1=20.03;
  static const Float_t kFDIP=-37.34;
  static const Float_t kFDIMU=6.;
  static const Float_t kFCORN=11.72;
//
// ZBEG       Beginning of the inner triplet
// D1BEG      Beginning of separator dipole 1
// D2BEG      Beginning of separator dipole 2
// CORBEG     Corrector dipole beginning (because of dimuon arm)
//
  static const Float_t kCORBEG=1920,kCOREND=kCORBEG+190, kCORRA2=4.5*4.5;
//
  static const Float_t kZBEG=2300;
  static const Float_t kZ1BEG=kZBEG+   0,kZ1END=kZ1BEG+630,kZ1RA2=3.5*3.5;
  static const Float_t kZ2BEG=kZBEG+ 880,kZ2END=kZ2BEG+550,kZ2RA2=3.5*3.5;
  static const Float_t kZ3BEG=kZBEG+1530,kZ3END=kZ3BEG+550,kZ3RA2=3.5*3.5;
  static const Float_t kZ4BEG=kZBEG+2430,kZ4END=kZ4BEG+630,kZ4RA2=3.5*3.5;
  static const Float_t kD1BEG=5843.5    ,kD1END=kD1BEG+945,kD1RA2=4.5*4.5;
  static const Float_t kD2BEG=12113.2   ,kD2END=kD2BEG+945,kD2RA2=4.5*.5;

//ZDC part -------------------------------------------------------------------

#endif
