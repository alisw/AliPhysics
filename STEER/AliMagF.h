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
  AliMagF(const AliMagF& maps);
  virtual ~AliMagF() {}
  AliMagF& operator=(const AliMagF& rhs);
  virtual void    Field(float *x, float *b)                  const;
  virtual void    Field(double *x, double *b)                const;
  virtual void    GetTPCInt(Float_t *xyz, Float_t *b)        const;
  virtual void    GetTPCIntCyl(Float_t *rphiz, Float_t *b)   const;
  virtual Int_t   Type() const {return fType;}
  virtual Float_t Max() const {return fMax;}
  virtual Int_t   Map() const {return fMap;}
  virtual Int_t   Integ() const {return fInteg;}
  virtual Int_t   PrecInteg() const {return fPrecInteg;}  
  virtual Float_t Factor() const {return fFactor;}
  virtual void    ReadField() {}
  virtual Float_t SolenoidField() const {return 2.;}
  virtual void    SetPrecInteg(Int_t integ);
  virtual void    SetReadField(Bool_t flag = kTRUE) {fReadField = flag;}
 protected:
  Int_t     fMap;       // Field Map identifier
  Int_t     fType;      // Mag Field type
  Int_t     fInteg;     // Default integration method as indicated in Geant
  Int_t     fPrecInteg; // Alternative integration method, e.g. for higher precision
  Float_t   fFactor;    // Multiplicative factor
  Float_t   fMax;       // Max Field as indicated in Geant
  Bool_t    fReadField; // Flag for reading the field from file (if available) 
  ClassDef(AliMagF,5)   //Base class for all Alice MagField
};

#endif
