#ifndef ALIMAGF_H
#define ALIMAGF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//----------------------------------------------------------------------
// Basic magnetic field class
// Used in all the detectors, and also in the traking classes
// Author:
//----------------------------------------------------------------------

#include "TNamed.h"

enum Field_t {kUndef=1, kConst=1, kConMesh=2, kDipoMap=3};

class AliMagF : public TNamed {

public:
  AliMagF();
  AliMagF(const char *name, const char *title, Int_t integ, 
	  Float_t factor = 1., Float_t fmax = 10.);
  virtual ~AliMagF() {}
  virtual void Field(Float_t *x, Float_t *b) const;
  virtual Int_t Type() const {return fType;}
  virtual Float_t Max() const {return fMax;}
  virtual Int_t Map() const {return fMap;}
  virtual Int_t Integ() const {return fInteg;}
  virtual Float_t Factor() const {return fFactor;}
  virtual void ReadField() {}
  virtual void SetDebug(Int_t level=0) {fDebug=level;}
  virtual Float_t SolenoidField() const {return 2.;}
  virtual Int_t GetDebug() const {return fDebug;}
  
protected:
  Int_t     fMap;    // Field Map identifier
  Int_t     fType;   // Mag Field type
  Int_t     fInteg;  // Integration method as indicated in Geant
  Float_t   fFactor; // Multiplicative factor
  Float_t   fMax;    // Max Field as indicated in Geant
  Int_t     fDebug;  // Debug flag

  ClassDef(AliMagF,1)  //Base class for all Alice MagField
};

#endif
