#ifndef AliMagF_H
#define AliMagF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TNamed.h"
#include "TVector.h"

enum Field_t {Undef=1, Const=1, ConMesh=2};

class AliMagF : public TNamed {

protected:
  Int_t     fMap;    // Field Map identifier
  Int_t     fType;   // Mag Field type
  Int_t     fInteg;  // Integration method as indicated in Geant
  Float_t   fFactor; // Multiplicative factor
  Float_t   fMax;    // Max Field as indicated in Geant

public:
  AliMagF(){}
  AliMagF(const char *name, const char *title, const Int_t integ, const Int_t map, 
	  const Float_t factor, const Float_t fmax);
  virtual ~AliMagF() {}
  virtual void Field(Float_t *x, Float_t *b);
  virtual Int_t Type() {return fType;}
  virtual Float_t Max() const {return fMax;}
  virtual Int_t Map() const {return fMap;}
  virtual Int_t Integ() const {return fInteg;}
  virtual Float_t Factor() const {return fFactor;}
  virtual void ReadField() {}
  
  ClassDef(AliMagF,1)  //Base class for all Alice MagField
};

class AliMagFC  : public AliMagF
{
  //Alice Constant Magnetic Field
private:

public:
  AliMagFC(){}
  AliMagFC(const char *name, const char *title, const Int_t integ, const Int_t map, 
	   const Float_t factor, const Float_t fmax);
  virtual ~AliMagFC() {}
  virtual void Field(Float_t *x, Float_t *b);
  virtual void ReadField() {}
  
  ClassDef(AliMagFC,1)  //Class for all Alice Constant MagField 
};

class AliMagFCM : public AliMagF
{
  //Alice Magnetic Field with constan mesh
protected:

  Float_t    fXbeg;  // Start of mesh in x
  Float_t    fYbeg;  // Start of mesh in y
  Float_t    fZbeg;  // Start of mesh in z
  Float_t    fXdel;  // Mesh step in x
  Float_t    fYdel;  // Mesh step in y
  Float_t    fZdel;  // Mesh step in z
  Double_t   fXdeli; // Inverse of Mesh step in x
  Double_t   fYdeli; // Inverse of Mesh step in y
  Double_t   fZdeli; // Inverse of Mesh step in z
  Int_t      fXn;    // Number of mesh points in x
  Int_t      fYn;    // Number of mesh points in y
  Int_t      fZn;    // Number of mesh points in z
  TVector   *fB;     // Field map
public:
  AliMagFCM(){}
  AliMagFCM(const char *name, const char *title, const Int_t integ, const Int_t map, 
	   const Float_t factor, const Float_t fmax);
  virtual ~AliMagFCM() {delete fB;}
  virtual void Field(Float_t *x, Float_t *b);
  virtual void ReadField();

  inline Float_t Bx(const Int_t ix, const Int_t iy, const Int_t iz) {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix));
  }
  inline Float_t By(const Int_t ix, const Int_t iy, const Int_t iz) {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix)+1);
  }
  inline Float_t Bz(const Int_t ix, const Int_t iy, const Int_t iz) {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix)+2);
  }
  
  ClassDef(AliMagFCM,1)  //Class for all Alice MagField with Constant Mesh
};

#endif
